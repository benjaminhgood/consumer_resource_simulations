import consumer_resource_sim
import numpy

s=1e-02
NUb=5
NUs=5
Nm=10
p=50
k=5
l=1
cv=0.05

immigration_measurements = consumer_resource_sim.run_simulation(s, 0,   0,   0, Nm, p, cv, k, l)
ecology_measurements =     consumer_resource_sim.run_simulation(s, 0,   Nm,  1, 0,  p, cv, k, l)
evolution_measurements =   consumer_resource_sim.run_simulation(s, NUb, NUs, 0, 0,  p, cv, k, l)

import pylab

pylab.figure()
    
pylab.plot(immigration_measurements[-1][3],'k-')
pylab.plot(ecology_measurements[-1][3],'r-')
pylab.plot(evolution_measurements[-1][3],'b-')

pylab.ylim([0,2])

pylab.figure()

immigration_ts = []
immigration_num_strains = []
immigration_homozygosities = []
for t,N,X,invh,clones,snps in immigration_measurements:
 
    strains = []
    for clone_idx, fitness, strategy_vector, mutations in clones:
        strains.append(tuple(strategy_vector))
    
    
    strains = set(strains)
    immigration_ts.append(t)
    immigration_num_strains.append(len(strains))
    immigration_homozygosities.append( consumer_resource_sim.calculate_strategy_homozygosity(clones))


immigration_homozygosities = numpy.array(immigration_homozygosities)


ecology_ts = []
ecology_num_strains = []
ecology_homozygosities = []
for t,N,X,invh,clones,snps in ecology_measurements:
 
    strains = []
    for clone_idx, fitness, strategy_vector, mutations in clones:
        strains.append(tuple(strategy_vector))
    
    
    strains = set(strains)
    ecology_ts.append(t)
    ecology_num_strains.append(len(strains))
    ecology_homozygosities.append( consumer_resource_sim.calculate_strategy_homozygosity(clones))

ecology_homozygosities = numpy.array(ecology_homozygosities)

evolution_ts = []
evolution_num_strains = []
evolution_homozygosities = []
for t,N,X,invh,clones,snps in evolution_measurements:
 
    strains = []
    for clone_idx, fitness, strategy_vector, mutations in clones:
        strains.append(tuple(strategy_vector))
    
        
    strains = set(strains)
    evolution_ts.append(t)
    evolution_num_strains.append(len(strains))
    evolution_homozygosities.append( consumer_resource_sim.calculate_strategy_homozygosity(clones))


evolution_homozygosities = numpy.array(evolution_homozygosities)


pylab.plot(immigration_ts, immigration_num_strains,'k.-') 
pylab.plot(immigration_ts, 1.0/immigration_homozygosities,'k.:') 
pylab.plot(ecology_ts, ecology_num_strains,'r.-')
pylab.plot(ecology_ts, 1.0/ecology_homozygosities,'r.:')
pylab.plot(evolution_ts, evolution_num_strains,'b.-')
pylab.plot(evolution_ts, 1.0/evolution_homozygosities,'b.:')

pylab.ylim([0,2*p])

pylab.figure()


for t,N,X,invh,clones,snps in immigration_measurements[-1:]:

    print "Calculating pairwise differences for t..."
    pi_matrix = consumer_resource_sim.calculate_pi_matrix(clones)
    flattened_pi_matrix = pi_matrix.flatten()*500
    flattened_pi_matrix.sort()
    
    print flattened_pi_matrix.min(), flattened_pi_matrix.max(), flattened_pi_matrix.mean()
    
    pylab.step(flattened_pi_matrix, 1.0-numpy.arange(0,len(flattened_pi_matrix))*1.0/len(flattened_pi_matrix),'k-',alpha=0.5)


for t,N,X,invh,clones,snps in ecology_measurements[-1:]:

    print "Calculating pairwise differences for t..."
    pi_matrix = consumer_resource_sim.calculate_pi_matrix(clones)
    flattened_pi_matrix = pi_matrix.flatten()*500
    flattened_pi_matrix.sort()
    
    print flattened_pi_matrix.min(), flattened_pi_matrix.max(), flattened_pi_matrix.mean()
    
    pylab.step(flattened_pi_matrix, 1.0-numpy.arange(0,len(flattened_pi_matrix))*1.0/len(flattened_pi_matrix),'r-',alpha=0.5)
    
for idx in [-30,-20,-10,-1]:
    
    t,N,X,invh,clones,snps = evolution_measurements[idx]

    print "Calculating pairwise differences for t..."
    pi_matrix = consumer_resource_sim.calculate_pi_matrix(clones)
    flattened_pi_matrix = pi_matrix.flatten()*500
    flattened_pi_matrix.sort()
    
    print flattened_pi_matrix.min(), flattened_pi_matrix.max(), flattened_pi_matrix.mean()
    
    pylab.step(flattened_pi_matrix, 1.0-numpy.arange(0,len(flattened_pi_matrix))*1.0/len(flattened_pi_matrix),'b-',alpha=0.5)
    

    
pylab.semilogx([1],[1],'k.')
pylab.xlim([1000,20000000])    
for clone_idx,fitness,strategy,mutations in evolution_measurements[-1][4][0:2]:
    print strategy

pylab.savefig('t2_distribution_1.pdf',bbox_inches='tight')

pylab.figure()

immigration_num_genes = []
for clone_idx,fitness,strategy_vector,mutations in immigration_measurements[-1][4]:
    immigration_num_genes.append((strategy_vector>0.5).sum())
immigration_num_genes.sort()

ecology_num_genes = []
for clone_idx,fitness,strategy_vector,mutations in ecology_measurements[-1][4]:
    ecology_num_genes.append((strategy_vector>0.5).sum())
ecology_num_genes.sort()

evolution_num_genes = []
for clone_idx,fitness,strategy_vector,mutations in evolution_measurements[-1][4]:
    evolution_num_genes.append((strategy_vector>0.5).sum())
evolution_num_genes.sort()

previous_evolution_num_genes = []
for clone_idx,fitness,strategy_vector,mutations in evolution_measurements[-25][4]:
    previous_evolution_num_genes.append((strategy_vector>0.5).sum())
previous_evolution_num_genes.sort()


pylab.step(immigration_num_genes, 1.0-numpy.arange(0,len(immigration_num_genes))*1.0/len(immigration_num_genes),color='k')
pylab.step(ecology_num_genes, 1.0-numpy.arange(0,len(ecology_num_genes))*1.0/len(ecology_num_genes),color='r')

pylab.step(evolution_num_genes, 1.0-numpy.arange(0,len(evolution_num_genes))*1.0/len(evolution_num_genes),color='b')
pylab.step(previous_evolution_num_genes, 1.0-numpy.arange(0,len(previous_evolution_num_genes))*1.0/len(previous_evolution_num_genes),color='b',alpha=0.5)

pylab.xlim([0,51])
    

#pylab.show()