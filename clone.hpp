#ifndef CLONE_HPP
#define CLONE_HPP

#include <vector>

#include "object_pool.hpp"

typedef int Label; // used for identifying unique mutations
Label null_label = -1;

typedef std::vector<double> GenePool;
typedef std::vector<double> SupplyVector;
typedef std::vector<double> HarvestVector;

 
class LabelGenerator{        
    public:
        LabelGenerator(): current_label(0) {};
        Label get_next_label() {
            if(current_label == null_label){
                ++current_label;
            }
            return current_label++;
        };
        Label current_label;
};


typedef std::vector<double> ResourceVector; // entry i gives percent of total input flux to resource i

class ExpressionLevel{
    public:
        int resource_idx;
        int level;
};

typedef std::vector<ExpressionLevel> ExpressionVector;


class Mutation{
    public:
        Label label;
        int time;
        bool track;
        double fitness_effect;
        int resource_idx;
        int expression_level;
        int resource_idx_2;
        int expression_level_2;
};

Mutation null_mutation{null_label,0,false,1,-1,0,-1,0};

inline bool operator<(Mutation const & m1, Mutation const & m2){ return m1.label < m2.label; }  
inline bool operator==(Mutation const & m1, Mutation const & m2){ return m1.label == m2.label; }
inline bool operator!=(Mutation const & m1, Mutation const & m2){ return !(m1==m2); }

typedef std::vector<Mutation> MutationList;


inline int sample_lineage_size(Random & random, double expected_size){
    return sample_poisson(random, expected_size);
}

class Clone{
    public:
       int size;
       double fitness; // W version
       
       // Strategy component 
       ExpressionVector expression_levels; // level of expression of used resources
       double total_expression;
       
       // Genome 
       MutationList mutations;
       
       Clone(){
           size=1;
           fitness=1;
           expression_levels.push_back(ExpressionLevel{0,1});
           total_expression=1;
       }

       void add_mutation(Mutation & mutation){
           
           // Update fitness
           fitness *= mutation.fitness_effect;
           
           // Update strategy
           if(mutation.resource_idx>=0){
           
               int old_expression_level = 0;
               
               for(auto expression_level_ptr=expression_levels.begin();expression_level_ptr!=expression_levels.end();++expression_level_ptr){
               
                   if(expression_level_ptr->resource_idx==mutation.resource_idx){
                    
                       old_expression_level = expression_level_ptr->level;
                        
                       if(mutation.expression_level>0){
                           // A modification
                           expression_level_ptr->level = mutation.expression_level;
                       }
                       else{
                           // A deletion
                           //std::cerr << "Deletion!" << std::endl;
                           expression_levels.erase(expression_level_ptr);
                       }
                    
                       break;        
                   }
               }
                
               if(old_expression_level==0){
                   // a gain of function mutation
                   expression_levels.push_back( ExpressionLevel{mutation.resource_idx, mutation.expression_level});
               }
                
               total_expression += mutation.expression_level-old_expression_level;      
           
           }
           
           // Update strategy (if double mutant)
           if(mutation.resource_idx_2>=0){
           
               int old_expression_level = 0;
               
               for(auto expression_level_ptr=expression_levels.begin();expression_level_ptr!=expression_levels.end();++expression_level_ptr){
               
                   if(expression_level_ptr->resource_idx==mutation.resource_idx_2){
                    
                       old_expression_level = expression_level_ptr->level;
                        
                       if(mutation.expression_level_2>0){
                           // A modification
                           expression_level_ptr->level = mutation.expression_level_2;
                       }
                       else{
                           // A deletion
                           //std::cerr << "Deletion!" << std::endl;
                           expression_levels.erase(expression_level_ptr);
                       }
                    
                       break;        
                   }
               }
                
               if(old_expression_level==0){
                   // a gain of function mutation
                   expression_levels.push_back( ExpressionLevel{mutation.resource_idx_2, mutation.expression_level_2});
               }
                
               total_expression += mutation.expression_level_2-old_expression_level;     
           
           }
           
           
           // Store mutation
           if(mutation.track){
               //std::cout << "Storing mutation!\n";
               mutations.push_back(mutation);
           }
       }

       void add_mutation(Mutation && mutation){
           
           // Update fitness
           fitness *= mutation.fitness_effect;
           
           // Update strategy
           if(mutation.resource_idx>=0){
           
               int old_expression_level = 0;
               
               for(auto expression_level_ptr=expression_levels.begin();expression_level_ptr!=expression_levels.end();++expression_level_ptr){
               
                   if(expression_level_ptr->resource_idx==mutation.resource_idx){
                    
                       old_expression_level = expression_level_ptr->level;
                        
                       if(mutation.expression_level>0){
                           // A modification
                           expression_level_ptr->level = mutation.expression_level;
                       }
                       else{
                           // A deletion
                           //std::cerr << "Deletion!" << std::endl;
                           expression_levels.erase(expression_level_ptr);
                       }
                    
                       break;        
                   }
               }
                
               if(old_expression_level==0){
                   // a gain of function mutation
                   expression_levels.push_back( ExpressionLevel{mutation.resource_idx, mutation.expression_level});
               }
                
               total_expression += mutation.expression_level-old_expression_level;
               
               // check total expression level
               /*int check_total_expression = 0;
               for(auto & expression_level : expression_levels){
                   check_total_expression += expression_level.level;
               }
           
               if(check_total_expression != total_expression){
                   std::cerr << "Total expression doesn't match up!" << std::endl;
                   std::cerr << mutation.expression_level << std::endl;
               }*/
                     
           }
           
           // Update strategy (if double mutant)
           if(mutation.resource_idx_2>=0){
           
               int old_expression_level = 0;
               
               for(auto expression_level_ptr=expression_levels.begin();expression_level_ptr!=expression_levels.end();++expression_level_ptr){
               
                   if(expression_level_ptr->resource_idx==mutation.resource_idx_2){
                    
                       old_expression_level = expression_level_ptr->level;
                        
                       if(mutation.expression_level_2>0){
                           // A modification
                           expression_level_ptr->level = mutation.expression_level_2;
                       }
                       else{
                           // A deletion
                           //std::cerr << "Deletion!" << std::endl;
                           expression_levels.erase(expression_level_ptr);
                       }
                    
                       break;        
                   }
               }
                
               if(old_expression_level==0){
                   // a gain of function mutation
                   expression_levels.push_back( ExpressionLevel{mutation.resource_idx_2, mutation.expression_level_2});
               }
                
               total_expression += mutation.expression_level_2-old_expression_level;     
           
           }
           
           // Store mutation
           if(mutation.track){
               //std::cout << "Storing mutation!\n";
               mutations.push_back(mutation);
           }
       }
       
       void form_mutant(Clone const & clone, Mutation & mutation){
           
           // First copy old clone to new one
           size = 1.0;
           fitness = clone.fitness;
           expression_levels = clone.expression_levels;
           total_expression = clone.total_expression;
           mutations = clone.mutations;
           // Then add mutation
           add_mutation(mutation);
       }

       void form_mutant(Clone const & clone, Mutation && mutation){
           // First copy old clone to new one
           size = 1.0;
           fitness = clone.fitness;
           expression_levels = clone.expression_levels;
           total_expression = clone.total_expression;
           mutations = clone.mutations;
           // Then add mutation
           add_mutation(mutation);
       }
       
       
       Mutation & get_earliest_mutation(){
           if(mutations.size() > 0){
               return mutations.front();
            }
            else{
               return null_mutation;
            }
       }
       
       Mutation & get_second_earliest_mutation(){
           if(mutations.size() > 1){
               return mutations[1];
            }
            else{
               return null_mutation;
            }
       }
       
       
       void remove_earliest_mutation(){
           mutations.erase(mutations.begin());
       }
};

typedef SharedObjectPool<Clone> ClonePool;
typedef std::vector<decltype(ClonePool().allocate())> Population;
typedef SharedObjectPool<Population> PopulationPool;
typedef SharedObjectPool<MutationList> GenomePool;

#endif
