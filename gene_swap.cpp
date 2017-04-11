#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "stats.hpp" // to help with random number generation
#include "clone.hpp"
#include "dfe.hpp"

//#define DEBUG

class SNP{
    public:
        Mutation mutation;
        double freq;
};

typedef std::vector<SNP> SNPList;


class Measurement{
    public:
        int time;
        double total_size;
        double avg_fitness;
        std::vector<double> harvest_vector;
        std::vector<double> fitnesses; 
        std::vector<ExpressionVector> expression_vectors;
        std::vector<MutationList> clone_snps;
        SNPList snps;
};


int main(int argc, char * argv[]){

    if(argc<11){
        std::cerr << "Usage: ./program s NUb NUs pscramble Nm p cv k l\n";
        return 1;
    }

    // Create and seed random number generator
    Random random = create_random();

    const int num_replicates = 1;

    const int teq = 0;
    
    const double Ne = 1e07;
    const double clone_sample_size = 1000;
    
    double Un = 0; // neutral mutation rate
    //const double Ub = 10.0/Ne; // beneficial mutation rate
    //double s = 2e-02; // fitness effect of beneficial mutations
    //double Us = 10.0/Ne;
    
    const double s = atof(argv[1]);
    const double Ub = atof(argv[2])/Ne;    
    const double Us = atof(argv[3])/Ne;
    const double p_scramble = atof(argv[4]);
    const double m = atof(argv[5])/Ne;
    
    const int num_resources = atoi(argv[6]);
    const double resource_stddev = atof(argv[7]);
    const double alpha = 1.0/resource_stddev/resource_stddev;
    const double kavg = atof(argv[8]);
    const int lmax = atoi(argv[9]);
    const int barcode_dt = atoi(argv[10]);
    const bool energy_constraint = true;
    
    double p0 = 0.5; // probability that a gene is lost, vs modified

    int initial_num_migrants = 1;
   
    std::gamma_distribution<> draw_supply(alpha,1);
   
    std::vector<double> normalized_supply_vector;
    double total_supply;
    for(int i=0;i<num_resources;++i){
        double supply = draw_supply(random);
        normalized_supply_vector.push_back(supply);
        total_supply += supply;
    }
    for(auto & b : normalized_supply_vector){
        b /= total_supply;
    }
    
    std::vector<double> supply_vector;
    for(int i=0,imax=normalized_supply_vector.size();i<imax;++i){
        supply_vector.push_back( normalized_supply_vector[i]*Ne );
    }
    
    
    // set up times to sample things based on tmax and dt;
    std::vector<int> times;
    times.push_back(1);
    for(int i=1;i<=100;++i){
        times.push_back(i*20000);
    }
    
    const int num_timepoints = times.size();
    const bool sequence_population = false;
    const int sequencing_depth = 1000;

    
    // Create DFEs
    
    // Immigration "DFE" for seeding initial population
    //auto seed_immigration_dfe = UniformImmigrationDFE(0, supply_vector.size(), supply_vector.size(), 1, 0, true);
    //auto seed_immigration_dfe = UniformImmigrationDFE(0, supply_vector.size(), kavg, lmax, 0);
    auto seed_immigration_dfe = ExactUniformImmigrationDFE(0, supply_vector.size(), int(kavg), lmax, 0);
    
    // Immigration "DFE" for further immigrants over time
    auto immigration_dfe = UniformImmigrationDFE(m, supply_vector.size(), kavg, lmax, s);
    
    // DFE that modifies fitnes
    auto fitness_dfe = GaussianDFE(Ub,s,false);
    
    double gene_gain_mutation_rate = Us*p0/(p0+1-kavg/supply_vector.size());
    double gene_change_mutation_rate = Us*(1-kavg/supply_vector.size())*(1-p0)*(1-kavg/supply_vector.size())/kavg*(p0+1-kavg/supply_vector.size());
    double gene_loss_mutation_rate = Us*(1-kavg/supply_vector.size())*(p0)*(1-kavg/supply_vector.size())/kavg*(p0+1-kavg/supply_vector.size());
    
    auto gene_change_dfe = create_composite_dfe( LossOfFunctionDFE(gene_loss_mutation_rate), ChangeOfFunctionDFE(gene_change_mutation_rate, lmax) );
    //auto strategy_dfe = create_gain_loss_change_dfe(Us, supply_vector.size(), kavg, p0, lmax, p_scramble, s, true);
    auto strategy_dfe = create_nested_uncorrelated_dfe( SwapFunctionDFE(Us, supply_vector.size(), lmax, false), p_scramble, s);
    auto dfe = create_composite_dfe( fitness_dfe, strategy_dfe );
    
    
    ///////////////////////////////////////////////////////////////////////////
    // 
    // Evolve population!
    //
    ///////////////////////////////////////////////////////////////////////////
    for(int replicate_idx=0;replicate_idx<num_replicates;++replicate_idx){
    
        LabelGenerator label_generator; // for uniquely labeling mutations

        MutationList fixed_mutations; // for storing mutations that fixed (performance)
    
        // for storing mutations that are fixed in clades but not in entire population. 
        // (performance hack, not clear it does much for large # resources)
        std::unordered_map<Label, decltype(GenomePool().allocate())> clade_lines_of_descent;
        std::unordered_map<Label, decltype(PopulationPool().allocate())> clade_clones;
        std::unordered_map<Label, Mutation> clade_next_mutations;
        std::unordered_map<Label, bool> clade_mutation_mismatches; 
    
    
        Population population; // the population in the current generation
        Population new_population; // the population in the next generation

        // we use object pools to make allocation of things more efficient
        ClonePool clone_pool(1000); 
        PopulationPool population_pool(supply_vector.size()); 
        GenomePool genome_pool(supply_vector.size()); 
    
        ////////////////////////////////////////////////////////////////////   
        //
        // Initialize the population
        //
        ////////////////////////////////////////////////////////////////////
        double W0=1;
        double total_fitness = 0;
        double total_size = 0;
        double total_expected_size = 0;
        std::vector<double> harvest_vector(supply_vector.size(),0);
        GenePool gene_pool(supply_vector.size(),0); // not doing anything with this right now
    
        // Build population from initial migrant pool!
    
        for(int migrant_idx=0;migrant_idx<initial_num_migrants;++migrant_idx){
            auto new_clone_ptr = clone_pool.allocate();
            seed_immigration_dfe.form_migrant(random, label_generator, 0, *new_clone_ptr);
            new_clone_ptr->size = Ne/(initial_num_migrants);
            population.push_back(new_clone_ptr);
        }
    
    
        auto observation_time_ptr = times.begin();
        auto observation_time_end = times.end();

        int num_barcodings = 0;

        for(int t=1; observation_time_ptr < observation_time_end; ++t){
        
            /////////////////////////////////////////////////////////////////////////////
            //
            // Take measurement and report to stdout
            //
            /////////////////////////////////////////////////////////////////////////////
            if(t == *observation_time_ptr){
            
                Measurement measurement;
                
                // First make some popuation-wide measurements
                measurement.time = t;
                measurement.total_size = total_size;
                measurement.avg_fitness = W0;
                measurement.harvest_vector = harvest_vector;
            
            
                // Set up distribution to draw clones proportional to their size
                std::vector<double> weights;
                for(auto const & clone_ptr : population){
                    weights.push_back(clone_ptr->size);
                }
                std::discrete_distribution<int> draw_clone_idx(weights.begin(), weights.end());
    
                // Draw sample of clones to report on
                for(int i=0;i<clone_sample_size;++i){
    
                    auto clone_idx = draw_clone_idx(random);
                    auto clone_ptr = population[clone_idx];
                
                    // get list of polymorphic mutations in clone
                    MutationList clone_mutations;
                    // If any mutations are saved in a clade line of descent, 
                    // add these first
                    if( clade_lines_of_descent.find( clone_ptr->get_earliest_mutation().label ) != clade_lines_of_descent.end() ){
                        // if there are mutations stored in the clade list
                        if( clade_lines_of_descent[clone_ptr->get_earliest_mutation().label]->size()>1 ){
                    
                            clone_mutations.insert(clone_mutations.end(), clade_lines_of_descent[clone_ptr->get_earliest_mutation().label]->begin(), clade_lines_of_descent[clone_ptr->get_earliest_mutation().label]->end()-1);
                        }
                    
                    }
                
                    // then add the remaining mutations already in the clone mutation list
                    clone_mutations.insert(clone_mutations.end(), clone_ptr->mutations.begin(), clone_ptr->mutations.end());
                
                    // Save a copy of the clone's absolute (log) fitness
                    measurement.fitnesses.push_back(std::log(clone_ptr->fitness));
                    // Save a copy of the clone's expression vector
                    measurement.expression_vectors.push_back( ExpressionVector(clone_ptr->expression_levels) );
                    // Save mutation list
                    measurement.clone_snps.push_back(clone_mutations);
                
                } // done sampling clones!
                
                //
                // Measure SNP frequencies (slower)
                //
                if(sequence_population){
                
                    std::unordered_map<Label, int> snp_counts; 
                    std::unordered_map<Label, Mutation> snp_mutations;
        
                
                    // Draw sample of clones to report on
                    for(int i=0;i<sequencing_depth;++i){
    
                        auto clone_idx = draw_clone_idx(random);
                        auto clone_ptr = population[clone_idx];
                        
                        // Add mutations to snp_counts
                        for(auto & mutation : clone_ptr->mutations){
                            
                            // If mutation not in list of SNPs, add it
                            if( snp_counts.find( mutation.label ) == snp_counts.end() ){
                                snp_counts[mutation.label] = 0;
                                snp_mutations[mutation.label] = mutation;
                            }
                            snp_counts[mutation.label]+=1;
                        }
                    }
                    
                    // Add fixed mutations to list
                    for(auto & mutation : fixed_mutations){
                        measurement.snps.push_back(SNP{mutation, 1.0});
                    }
                    
                    // Done measuring SNP frequencies, now add them to list of SNPs
                    // now add these to the list
                    for(auto & label_count_pair : snp_counts){
                    
                        auto label = label_count_pair.first;
                        auto allele_count = label_count_pair.second;
                        
                        // Let's not include private mutations
                        if(allele_count < 2){
                            continue;
                        }
                        
                        double allele_frequency = allele_count*1.0/sequencing_depth;
                        
                        // If this mutation labels many clade fixations
                        // assign them the same frequency and append them
                        if( clade_lines_of_descent.find( label ) != clade_lines_of_descent.end() ){
                            
                            for(auto & mutation : *clade_lines_of_descent[label]){
                                measurement.snps.push_back(SNP{mutation, allele_frequency});
                            }
                        }
                        else{
                           // Mutation is not in a clade, so add it on its own
                           measurement.snps.push_back(SNP{snp_mutations[label], allele_frequency});
                        }
                    } // done adding SNPs
                } // done "sequencing" popluation
                
                //
                // Now report to STDOUT!
                //
                std::cout << "# Measurement\n";
                std::cout << "Time: " << measurement.time << std::endl;
                std::cout << "Total size: " << measurement.total_size << std::endl;
                std::cout << "Avg fitness: " << measurement.avg_fitness << std::endl;
    
                std::cout << "Scaled Harvest vector: ";
                for(int resource_idx=0;resource_idx<measurement.harvest_vector.size();++resource_idx){
                    std::cout << measurement.harvest_vector[resource_idx]/supply_vector[resource_idx] << " ";
                }
                std::cout << std::endl;
    
                std::cout << "Clone sample: " << measurement.clone_snps.size() << std::endl;
                for(int i=0,imax=measurement.clone_snps.size();i<imax;++i){
        
                    std::cout << i << " ; " << measurement.fitnesses[i] << " ; ";
        
                    for(auto & expression_level : measurement.expression_vectors[i]){
                        std::cout << expression_level.resource_idx << "," << expression_level.level << " ";
                    }
                    std::cout << " ; ";
                    for(auto & mutation : measurement.clone_snps[i]){
                        std::cout << mutation.label << "," << mutation.time << "," << mutation.fitness_effect << "," << mutation.resource_idx << "," << mutation.expression_level << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << "SNP sample: ";
                for(auto & snp : measurement.snps){
                    std::cout << snp.mutation.label << "," << snp.mutation.time << "," << snp.mutation.fitness_effect << "," << snp.mutation.resource_idx << "," << snp.mutation.expression_level << ";" << snp.freq << " ";
                }
                std::cout << std::endl;
                
                ++observation_time_ptr;
            } // done taking measurement!
        
            /////////////////////////////////////////////////////////////////////////////
            //
            // Birth and death cycle!
            //
            /////////////////////////////////////////////////////////////////////////////
        
            // First measure 
            // (1) total size
            // (2) mean fitness
            // (3) harvest vectors
            //std::cerr << "Measuring everything...\t";
            
            // First reset variables so that we can measure them for this generation!
            total_fitness = 0; // reset these variables so that we can 
            total_size = 0;    // remeasure them for the next generation
            total_expected_size = 0;
            for(auto & h : harvest_vector){
                h=0;
            }
            // Now do measurement
            for(auto & clone_ptr : population){
        
                double energy_penalty = 1;
                if(energy_constraint){
                    energy_penalty = (clone_ptr->total_expression*1.0/lmax/kavg)*std::exp(1-(clone_ptr->total_expression*1.0/lmax/kavg));
                }
        
                total_size += clone_ptr->size;
                total_fitness += clone_ptr->size*(clone_ptr->fitness*energy_penalty/W0);
                for(auto & expression_level : clone_ptr->expression_levels){
                    harvest_vector[expression_level.resource_idx] += (clone_ptr->size) * (clone_ptr->fitness*energy_penalty/W0) * expression_level.level*1.0/clone_ptr->total_expression;
                }
        
            }
            //std::cerr << "Done!\n";

            // Print status message
            if(t%10000==0){
                std::cerr << t << " " << total_size << " " << std::log(W0*total_fitness/total_size) << " " << fixed_mutations.size() << " " << population.size() << " " << clade_lines_of_descent.size() << std::endl;
            }
        
            double m = immigration_dfe.get_migration_rate();

            
            // Now generate offspring 
            //std::cerr << "Generating offspring for " << population.size() << "  lineages...\n";
            for(auto & clone_ptr : population){
            
                //std::cerr << "Calculating num offspring...\n";
                double expected_num_offspring = 0;
                for(auto & expression_level : clone_ptr->expression_levels){
                
                    double energy_penalty = 1;
                    if(energy_constraint){
                        energy_penalty =     (clone_ptr->total_expression*1.0/lmax/kavg)*std::exp(1-(clone_ptr->total_expression*1.0/lmax/kavg));
                    }
                
                    expected_num_offspring += supply_vector[expression_level.resource_idx]/harvest_vector[expression_level.resource_idx]*(clone_ptr->size) * (clone_ptr->fitness*energy_penalty/W0) * expression_level.level*1.0/clone_ptr->total_expression;
            
                }
                total_expected_size += expected_num_offspring;
                //std::cerr << "Done!\n";

            
                // Draw clonal and mutant offspring
                //std::cout << "Drawing offspring...\n";
                double U = dfe.get_mutation_rate(gene_pool, *clone_ptr);
            
                int num_clonal_offspring = sample_lineage_size(random, (1-U-m)*expected_num_offspring);
            
                int num_mutants = sample_lineage_size(random, U*expected_num_offspring);
            
            
                //std::cerr << "Propagating clonal offspring!" << std::endl;
                if(num_clonal_offspring > 0){
                    // progagate clones
                    clone_ptr->size = num_clonal_offspring;
                    new_population.push_back(clone_ptr);
                }
                
                //std::cerr << "Propagating mutants!" << std::endl;
                if(num_mutants > 0){
                    for(int i=0;i<num_mutants;++i){
                        auto new_clone_ptr = clone_pool.allocate();
                        new_clone_ptr->form_mutant(*clone_ptr, dfe.get_mutation(random, label_generator, t, gene_pool, *clone_ptr));
                        new_population.push_back(new_clone_ptr);
                    }
                }
                //std::cerr << "Done!" << std::endl;
            }
        
            // now do migrants!
            int num_migrants = sample_lineage_size(random, m*total_expected_size);
            if(num_migrants > 0){
                for(int i=0;i<num_migrants;++i){
                    auto new_clone_ptr = clone_pool.allocate();
                    immigration_dfe.form_migrant(random, label_generator, t, *new_clone_ptr);
                    new_population.push_back(new_clone_ptr);
                }
            }
            //std::cerr << "Done!\n";    
            // make the next generation the current generation
            population.swap(new_population);
            new_population.clear();
        
            // normalize fitness (why not)
            double Wavg = total_fitness/total_size;
            W0 *= Wavg;
        
            /////////////////////////////////////////////////////////////////////////////
            //
            // Barcode population (optional)
            //
            /////////////////////////////////////////////////////////////////////////////
            if(t % barcode_dt == 0){
            
                num_barcodings += 1;
            
                for(auto & clone_ptr : population){
                    clone_ptr->add_mutation( get_marker_mutation(label_generator, t) );
                }
            
            }
        
            /////////////////////////////////////////////////////////////////////////////
            //
            // Check for fixed mutations (performance) 
            //
            /////////////////////////////////////////////////////////////////////////////
             
            // First check for population-wide fixations
            //std::cerr << "Checking for population-wide fixations mutations!\n";
            auto mutation = population.front()->get_earliest_mutation();
            if(mutation != null_mutation){
                bool fixed = true;
                for(auto & clone_ptr : population){
                    if(clone_ptr->get_earliest_mutation() != mutation){
                        fixed = false;
                        break;
                    }
                }
                if(fixed){
            
                    // a mutation has fixed! 
                    //std::cout << "Fixed!\n";
                    if( clade_lines_of_descent.find( mutation.label ) != clade_lines_of_descent.end() ){
                        fixed_mutations.insert(fixed_mutations.end(), clade_lines_of_descent[mutation.label]->begin(), clade_lines_of_descent[mutation.label]->end());
                
                    }
                    else{
                        // otherwise, we have to push it ourselves
                        fixed_mutations.push_back(mutation);
                    }
            
                    clade_lines_of_descent.clear();
            
                    // now remove that mutation from everyone
                    for(auto & clone_ptr : population){
                        clone_ptr->remove_earliest_mutation();
                    }
                }
            } // done checking for population-wide fixations
        
            // Now check for fixations within top-level clades
            if(t%1000==0){ // slower, so don't do it as often...
            //if(false){ 
            
                bool something_fixed;
                
                do{
               
                    something_fixed = false;
               
                    //std::cout << "Checking for clade fixed mutations...\n";
            
                    clade_clones.clear();
                    clade_next_mutations.clear();
                    clade_mutation_mismatches.clear();
        
                    for (auto & pair : clade_lines_of_descent ) { 
                        clade_clones[pair.first] = population_pool.allocate();
                        clade_clones[pair.first]->clear();
                        clade_mutation_mismatches[pair.first] = false;
                    }
        
                    // Populate clade membership map
                    for(auto & clone_ptr : population){
            
                        auto mutation = clone_ptr->get_earliest_mutation(); 
                
                        if(mutation.label!=null_label){
                            // try to add it to the list
                
                            if( clade_lines_of_descent.find( mutation.label ) == clade_lines_of_descent.end() ){
                        
                                // not already present in the list! 
                                clade_lines_of_descent[mutation.label] = genome_pool.allocate();
                                clade_lines_of_descent[mutation.label]->clear();
                                clade_lines_of_descent[mutation.label]->push_back(mutation);
                        
                                clade_clones[mutation.label] = population_pool.allocate();
                                clade_clones[mutation.label]->clear();
                
                                clade_mutation_mismatches[mutation.label] = false;
                            }
                    
                            if( clade_next_mutations.find( mutation.label ) == clade_next_mutations.end() ){
                    
                                // no next mutation record yet!
                                clade_clones[mutation.label]->push_back(clone_ptr);
                                clade_next_mutations[mutation.label] = clone_ptr->get_second_earliest_mutation();
                        
                                if(clade_next_mutations[mutation.label]==null_mutation){
                                clade_mutation_mismatches[mutation.label] = true;
                                }   
            
                            }
                            else{
                                // already present in list! 
                                if(!clade_mutation_mismatches[mutation.label]){
                                    // not a mismatch yet
                                    if(clone_ptr->get_second_earliest_mutation() != clade_next_mutations[mutation.label]){
                                        // a mismatch!
                                        clade_mutation_mismatches[mutation.label] = true;
                                    }
                                    else{
                                        // a match, so add it to the list of clones
                                        clade_clones[mutation.label]->push_back(clone_ptr);
                                    }   
                                } 
                            }   
                        }
                    } // done checking for within-clade
            
                    //std::cerr << "Done checking!\n";
                    //std::cerr << "Seeing if anything fixed...\n";
                    // done looping over clones
                    // now see if anything fixed!   
                    for (auto & pair : clade_clones ) {
                
                        // if there are no clones corresponding to this clade, 
                        // delete it!
                        if(pair.second->empty()){
                            //std::cerr << "Deleting clade!\n";
                            clade_lines_of_descent.erase(pair.first);
                        }
                        else{
                            if(clade_mutation_mismatches[pair.first]==false){
                                // no mismatches, so clade fixation!
                                //std::cerr << "Found a clade fixation!\n";
                        
                                something_fixed=true;
                        
                                // add new mutation to line of descent
                                //std::cerr << "Switching labels!\n";
                                clade_lines_of_descent[pair.first]->push_back( clade_next_mutations[pair.first] );
                                // switch label to most recent mutation
                                clade_lines_of_descent[clade_next_mutations[pair.first].label] = genome_pool.allocate();  
                                clade_lines_of_descent[ clade_next_mutations[pair.first].label ]->swap( *clade_lines_of_descent[pair.first] );
                                clade_lines_of_descent.erase(pair.first);
                                //std::cerr << "Removing mutations from clones!";
                                // now remove the mutation from everyone
                                for(auto & clone_ptr : *pair.second){
                                    clone_ptr->remove_earliest_mutation();
                                }
                                //std::cerr << "Done with fixing!\n";
                            }   
                        }
                    }
            
                    clade_clones.clear();
                    clade_next_mutations.clear();
                    clade_mutation_mismatches.clear();
                }while(something_fixed);
        
                //std::cerr << "Done checking for clade fixations!\n";        
            } // done checking for clade fixations
        
                 
        } // done evolving popluation!
        //std::cout << "Done! " << total_size << std::endl;
        population.clear();
        new_population.clear();
        clade_lines_of_descent.clear();
    
    
    } // done looping over replicate populations
        
    return 0;
}


