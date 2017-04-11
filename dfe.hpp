#ifndef DFE_HPP
#define DFE_HPP

#include "stats.hpp"
#include "clone.hpp"
#include <cmath>
#include <memory>

inline Mutation get_marker_mutation(LabelGenerator & label_generator, int t){
    return Mutation{label_generator.get_next_label(),t,true,1,-1,0,-1,0};
}

// A dummy DFE that only returns neutral mutations
class NeutralDFE{
    public:
        const double U0;
        const bool track;
        NeutralDFE(double U, bool track=false): U0(U), track(track) {};
        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return U0; };
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, GenePool const & gene_pool, Clone const & individual){
            return Mutation{label_generator.get_next_label(),track,1,-1,0,-1,0};
        };
};

// A DFE that returns only fitness mutations
class DeltaDFE{
    public:
        const bool track;
        const double U0;
        const double W; 
        
        DeltaDFE(double U, double s, bool track=false): U0(U), W(std::exp(s)), track(track) {};
        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return U0; };
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
            return Mutation{label_generator.get_next_label(),t,track,W,-1,0,-1,0};
        };
};

// A DFE that returns only fitness mutations
class GaussianDFE{
    public:
        const bool track;
        const double U0;
        const double s; 
        
        GaussianDFE(double U, double s, bool track=false): U0(U), s(s), track(track) {};
        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return U0; };
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
            return Mutation{label_generator.get_next_label(),t,track,std::exp(s*sample_normal(random)),-1,0,-1,0};
        };
};

// A DFE that returns only fitness mutations
class UncorrelatedDFE{
    public:
        const bool track;
        const double U0;
        const double s; 
        
        UncorrelatedDFE(double U, double s, bool track=false): U0(U), s(s), track(track) {};
        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return U0; };
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
            
            double W = std::exp(s*sample_normal(random))/individual.fitness;
            
            return Mutation{label_generator.get_next_label(),t,track,W,-1,0,-1,0};
        };
};

// A DFE that returns only strategy mutations
class TwoResourceBetaDFE{
    public:
        const bool track;
        const double U0;
        const double alpha;
        const int lmax; // discretization
        
        TwoResourceBetaDFE(double U, double sigma, int lmax, bool track=false): U0(U), alpha(1.0/4/sigma/sigma-1), lmax(lmax), track(track) {
        
            for(int l=0;l<=lmax+1;++l){
                gamma_distributions.push_back( std::gamma_distribution<>(alpha*((l+1)*1.0/(lmax+2))) );
            }
        
        };
        
        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return U0; };
        
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
        
            int l = 0;
            for(auto & expression_level : individual.expression_levels){
                if(expression_level.resource_idx==0){
                    l = expression_level.level;
                    break;
                }
            }
            
            double a = gamma_distributions[l](random);
            double b = gamma_distributions[lmax-l](random);
            
            int new_l = 1+long(a/(a+b)*(lmax*(1-1e-09)));
            
            if(new_l<1 || new_l>lmax){
                std::cerr << new_l << std::endl;
            }
            
            return Mutation{label_generator.get_next_label(), t, track, 1, 0, new_l,1,lmax+1-new_l};
        };
        
    private:
        std::vector<std::gamma_distribution<>> gamma_distributions;
};


// A DFE that returns only strategy mutations
class GainOfFunctionDFE{
    public:
        const bool track;
        const double U0;
        const int num_resources;
        const int lmax; // max gene expression level
        
        GainOfFunctionDFE(double U, int num_resources, int lmax, bool track=false): U0(U), num_resources(num_resources), lmax(lmax), track(track), draw_resource_idx(0,num_resources-1), draw_expression_level(1,lmax) {};
        
        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return U0; };
        
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
        
            int resource_idx = draw_resource_idx(random);
            int expression_level = draw_expression_level(random);
        
            return Mutation{label_generator.get_next_label(), t, track, 1, resource_idx, expression_level,-1,0};
        };
        
    private:
        std::uniform_int_distribution<> draw_resource_idx;
        std::uniform_int_distribution<> draw_expression_level;
};

// A DFE that returns only strategy mutations
class LossOfFunctionDFE{
    public:
        const bool track;
        const double U0; // rate per gene
        
        LossOfFunctionDFE(double U, bool track=false): U0(U), track(track) {};
        
        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return U0*individual.expression_levels.size(); };
        
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
        
            auto draw_existing_resource_idx = std::uniform_int_distribution<>(0,individual.expression_levels.size()-1);
            
            int resource_idx = individual.expression_levels[draw_existing_resource_idx(random)].resource_idx;
            
            int expression_level = 0;
        
            return Mutation{label_generator.get_next_label(), t, track, 1, resource_idx, expression_level,-1,0};
        };
    
};

// A DFE that returns only strategy mutations
class ChangeOfFunctionDFE{
    public:
        const bool track;
        const double U0; // rate per gene
        const int lmax;
        ChangeOfFunctionDFE(double U, int lmax, bool track=false): U0(U), lmax(lmax), track(track), draw_expression_level(1,lmax) {};
        
        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return U0*individual.expression_levels.size(); };
        
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
        
            auto draw_existing_resource_idx = std::uniform_int_distribution<>(0,individual.expression_levels.size()-1);
            
            int resource_idx = individual.expression_levels[draw_existing_resource_idx(random)].resource_idx;
            
            int expression_level = draw_expression_level(random);
        
            return Mutation{label_generator.get_next_label(), t, track, 1, resource_idx, expression_level,-1,0};
        };
    
    private:
        std::uniform_int_distribution<> draw_expression_level;
};


// A DFE that returns only strategy mutations
class SwapFunctionDFE{
    public:
        const bool track;
        const double U0;
        const int num_resources;
        const int lmax; // max gene expression level
        
        SwapFunctionDFE(double U, int num_resources, int lmax, bool track=false): U0(U), num_resources(num_resources), lmax(lmax), track(track), draw_resource_idx(0,num_resources-1), draw_expression_level(1,lmax) {};
        
        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return U0; };
        
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
        
            
            int new_resource_idx = draw_resource_idx(random);
            int new_level = 0;
            for(auto & expression_level : individual.expression_levels){
                if(expression_level.resource_idx==new_resource_idx){
                    new_level = expression_level.level;
                    break;
                }
            }
            
            auto draw_existing_resource_idx = std::uniform_int_distribution<>(0,individual.expression_levels.size()-1);
            
            auto existing_expression_level = individual.expression_levels[draw_existing_resource_idx(random)];
            
            int existing_resource_idx = existing_expression_level.resource_idx;
            int existing_level = existing_expression_level.level;
            
            
            return Mutation{label_generator.get_next_label(), t, track, 1, new_resource_idx, existing_level, existing_resource_idx, new_level};
        
        };
        
    private:
        std::uniform_int_distribution<> draw_resource_idx;
        std::uniform_int_distribution<> draw_expression_level;
};



template <class DFE1, class DFE2>
class CompositeDFE{
    public:
        DFE1 dfe1;
        DFE2 dfe2;
        
        CompositeDFE(DFE1 dfe1, DFE2 dfe2): dfe1(dfe1), dfe2(dfe2) {};

        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return dfe1.get_mutation_rate(gene_pool, individual)+dfe2.get_mutation_rate(gene_pool, individual); };
        
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
            
            double U1 = dfe1.get_mutation_rate(gene_pool, individual);
            double U2 = dfe2.get_mutation_rate(gene_pool, individual);
            
            if( sample_uniform(random) < U1/(U1+U2) ){
                return dfe1.get_mutation(random, label_generator, t, gene_pool, individual);
            }
            else{
                return dfe2.get_mutation(random, label_generator, t, gene_pool, individual);
            }
        }
};

// A utility function for generating these DFEs with minimal syntactic overhead
template <class DFE1, class DFE2>
inline auto create_composite_dfe(DFE1 dfe1, DFE2 dfe2) -> CompositeDFE<DFE1, DFE2> {
    return CompositeDFE<DFE1, DFE2>(dfe1, dfe2);
}


template <class DFE>
class NestedUncorrelatedDFE{
    public:
        DFE dfe;
        double p; // probability to draw an uncorrelated fitness!
        double s;
        
        NestedUncorrelatedDFE(DFE dfe, double p, double s): dfe(dfe), p(p), s(s) {};

        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return dfe.get_mutation_rate(gene_pool, individual); };
        
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
            
            auto mutation = dfe.get_mutation(random, label_generator, t, gene_pool, individual);
            
            
            double W=1;
            if( sample_uniform(random) < p ){
                W = std::exp(s*sample_normal(random))/individual.fitness;
            }
            
            mutation.fitness_effect = W;
            return mutation;
        }
};

// A utility function for generating these DFEs with minimal syntactic overhead
template <class DFE>
inline auto create_nested_uncorrelated_dfe(DFE dfe, double p, double s) -> NestedUncorrelatedDFE<DFE> {
    return NestedUncorrelatedDFE<DFE>(dfe,p,s);
}

template <class DFE>
class NestedDeltaDFE{
    public:
        DFE dfe;
        double W;
        
        NestedDeltaDFE(DFE dfe, double s): dfe(dfe), W(std::exp(s)) {};

        double get_mutation_rate(GenePool const & gene_pool, Clone const & individual) const { return dfe.get_mutation_rate(gene_pool, individual); };
        
        Mutation get_mutation(Random & random, LabelGenerator & label_generator, int t, GenePool const & gene_pool, Clone const & individual){
            
            auto mutation = dfe.get_mutation(random, label_generator, t, gene_pool, individual);
            mutation.fitness_effect *= W;
            return mutation;
        }
};

// A utility function for generating these DFEs with minimal syntactic overhead
template <class DFE>
inline auto create_nested_delta_dfe(DFE dfe, double s) -> NestedDeltaDFE<DFE> {
    return NestedDeltaDFE<DFE>(dfe,s);
}




// A DFE for immigrants
class UniformImmigrationDFE{
    public:
        const bool tracked;
        const double m;
        const int num_total_resources;
        const double avg_used_resources;
        const int lmax; // max gene expression level
        const double s;
        
        UniformImmigrationDFE(double m, int num_total_resources, double avg_used_resources, int lmax, double s, bool tracked=false): tracked(tracked), m(m), s(s), num_total_resources(num_total_resources), avg_used_resources(avg_used_resources), lmax(lmax), draw_resource_use(avg_used_resources/num_total_resources), draw_expression_level(1,lmax) {};
        
        double get_migration_rate() const { return m; };
        
        void form_migrant(Random & random, LabelGenerator & label_generator, int t, Clone & clone, double W0=1){
            // size starts at one
            clone.size = 1;
            // fitness is random draw from some distribution
            clone.fitness = W0*std::exp( sample_normal(random)*s );
            clone.mutations.clear();
            // now draw random strategy
            clone.expression_levels.clear();
            clone.total_expression = 0;
            for(int resource_idx=0;resource_idx<num_total_resources;++resource_idx){
            
                if(draw_resource_use(random)){
                    int expression_level = draw_expression_level(random);
                    clone.expression_levels.push_back(ExpressionLevel{resource_idx, expression_level});
                    clone.total_expression += expression_level;
                    
                }
            
            }
            if(tracked){
                clone.mutations.push_back(get_marker_mutation(label_generator, t));
            }
        
        }
    private:
        std::bernoulli_distribution draw_resource_use;
        std::uniform_int_distribution<> draw_expression_level;
};

// A DFE for immigrants
class ExactUniformImmigrationDFE{
    public:
        const bool tracked;
        const double m;
        const int num_total_resources;
        const int k;
        const int lmax; // max gene expression level
        const double s;
        
        ExactUniformImmigrationDFE(double m, int num_total_resources, int k, int lmax, double s, bool tracked=false): tracked(tracked), m(m), s(s), num_total_resources(num_total_resources), k(k), lmax(lmax), draw_expression_level(1,lmax) {};
        
        double get_migration_rate() const { return m; };
        
        void form_migrant(Random & random, LabelGenerator & label_generator, int t, Clone & clone, double W0=1){
            // size starts at one
            clone.size = 1;
            // fitness is random draw from some distribution
            clone.fitness = W0*std::exp( sample_normal(random)*s );
            clone.mutations.clear();
            // now draw random strategy
            clone.expression_levels.clear();
            clone.total_expression = 0;
            
            resource_idxs.clear();
            for(int i=0;i<num_total_resources;++i){
                resource_idxs.push_back(i);
            }
            for(int i=0;i<k;++i){
            
                auto draw_existing_resource_idx = std::uniform_int_distribution<>(0,num_total_resources-i-1);
                
                auto resource_idx_idx = draw_existing_resource_idx(random);
                auto resource_idx = resource_idxs[resource_idx_idx];
                auto expression_level = draw_expression_level(random);
                clone.expression_levels.push_back(ExpressionLevel{resource_idx, expression_level});
                clone.total_expression+=expression_level;
                
                resource_idxs.erase(resource_idxs.begin()+resource_idx_idx);
            }
            
            if(tracked){
                clone.mutations.push_back(get_marker_mutation(label_generator, t));
            }
        
        }
    private:
        std::uniform_int_distribution<> draw_expression_level;
        std::vector<int> resource_idxs;
};

// A DFE for immigrants
class SpecialistImmigrationDFE{
    public:
        const bool tracked;
        const double m;
        const int num_total_resources;
        const int lmax; // max gene expression level
        const double s;
        
        SpecialistImmigrationDFE(double m, int num_total_resources, int lmax, double s, bool tracked=false): tracked(tracked), m(m), s(s), num_total_resources(num_total_resources), lmax(lmax), draw_resource_idx(0,num_total_resources-1) {};
        
        double get_migration_rate() const { return m; };
        
        void form_migrant(Random & random, LabelGenerator & label_generator, int t, Clone & clone, double W0=1){
            // size starts at one
            clone.size = 1;
            // fitness is random draw from some distribution
            clone.fitness = W0*std::exp( sample_normal(random)*s );
            clone.mutations.clear();
            // now draw random strategy
            clone.expression_levels.clear();
            clone.total_expression = 0;
            int special_resource_idx = draw_resource_idx(random);
            
            for(int resource_idx=0;resource_idx<num_total_resources;++resource_idx){
                
                int expression_level = 1;
                if(resource_idx==special_resource_idx){
                    expression_level=lmax;
                }
                
                clone.expression_levels.push_back(ExpressionLevel{resource_idx, expression_level});
                clone.total_expression += expression_level;
            }
            if(tracked){
                clone.mutations.push_back(get_marker_mutation(label_generator, t));
            }
        
        }
    private:
        std::uniform_int_distribution<> draw_resource_idx;
};

// A DFE for immigrants
class GeneralistImmigrationDFE{
    public:
        const bool tracked;
        const double m;
        const int num_total_resources;
        const int lmax; // max gene expression level
        const double s;
        
        GeneralistImmigrationDFE(double m, int num_total_resources, int lmax, double s, bool tracked=false): tracked(tracked), m(m), s(s), num_total_resources(num_total_resources), lmax(lmax), draw_resource_idx(0,num_total_resources-1) {};
        
        double get_migration_rate() const { return m; };
        
        void form_migrant(Random & random, LabelGenerator & label_generator, int t, Clone & clone, double W0=1){
            // size starts at one
            clone.size = 1;
            // fitness is random draw from some distribution
            clone.fitness = W0*std::exp( sample_normal(random)*s );
            clone.mutations.clear();
            // now draw random strategy
            clone.expression_levels.clear();
            clone.total_expression = 0;
            int special_resource_idx = draw_resource_idx(random);
            
            for(int resource_idx=0;resource_idx<num_total_resources;++resource_idx){
                
                int expression_level;
                
                if(resource_idx%2==0){
                    expression_level = ((lmax)/2);    
                }
                else{
                    expression_level = ((lmax)/2)+1;
                }
                
                clone.expression_levels.push_back(ExpressionLevel{resource_idx, expression_level});
                clone.total_expression += expression_level;
            }
            if(tracked){
                clone.mutations.push_back(get_marker_mutation(label_generator, t));
            }
        
        }
    private:
        std::uniform_int_distribution<> draw_resource_idx;
};

inline NestedUncorrelatedDFE<CompositeDFE<GainOfFunctionDFE, CompositeDFE<LossOfFunctionDFE, ChangeOfFunctionDFE>>> create_gain_loss_change_dfe(double Us, int p, double kavg, double p0, int lmax, double pscramble, double s, bool track=false){

    double gene_gain_mutation_rate = Us*p0/(p0+1-kavg/p);
    double gene_change_mutation_rate = Us*(1-kavg/p)*(1-p0)*(1-kavg/p)/kavg*(p0+1-kavg/p);
    double gene_loss_mutation_rate = Us*(1-kavg/p)*(p0)*(1-kavg/p)/kavg*(p0+1-kavg/p);

    auto gene_change_dfe = create_composite_dfe( LossOfFunctionDFE(gene_loss_mutation_rate, track), ChangeOfFunctionDFE(gene_change_mutation_rate, lmax, track) );
    auto strategy_dfe = create_nested_uncorrelated_dfe( create_composite_dfe( GainOfFunctionDFE(gene_gain_mutation_rate, p, lmax, track), gene_change_dfe ), pscramble, s);
    
    return strategy_dfe;
}


#endif
