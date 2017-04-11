import pylab
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy
import matplotlib.colors as colors
import matplotlib.cm as cmx
import sys

mpl.rcParams['font.size'] = 10.0
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


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

immigration_measurements = consumer_resource_sim.run_gene_swap_simulation(s, 0,   0,   0, Nm, p, cv, k, l, dt=10000)
immigration_t2s = consumer_resource_sim.calculate_pi_matrix(immigration_measurements[-1][4]).flatten()*10000/2
immigration_t2s.sort()
ts,immigration_num_strains = consumer_resource_sim.calculate_num_strain_timecourse(immigration_measurements)


ecology_measurements =     consumer_resource_sim.run_gene_swap_simulation(s, 0,   Nm,  1, 0,  p, cv, k, l, dt=10000)
ecology_t2s = consumer_resource_sim.calculate_pi_matrix(ecology_measurements[-1][4]).flatten()*10000/2
ecology_t2s.sort()
ts,ecology_num_strains = consumer_resource_sim.calculate_num_strain_timecourse(ecology_measurements)

evolution_measurements =   consumer_resource_sim.run_gene_swap_simulation(s, NUb, NUs, 0, 0,  p, cv, k, l, dt=100)
evolution_t2s = consumer_resource_sim.calculate_pi_matrix(evolution_measurements[-1][4]).flatten()*100/2
evolution_t2s.sort()
ts,evolution_num_strains = consumer_resource_sim.calculate_num_strain_timecourse(evolution_measurements)

replicate_evolution_measurements =  consumer_resource_sim.run_gene_swap_simulation(s, NUb, NUs, 0, 0,  p, cv, k, l, dt=100)
replicate_evolution_t2s = consumer_resource_sim.calculate_pi_matrix(replicate_evolution_measurements[-1][4]).flatten()*100/2
replicate_evolution_t2s.sort()
ts,replicate_evolution_num_strains = consumer_resource_sim.calculate_num_strain_timecourse(replicate_evolution_measurements)

pylab.figure(figsize=(5,1.7))

pylab.step(immigration_t2s, 1.0-numpy.arange(0,len(immigration_t2s))*1.0/len(immigration_t2s),'k-',label='immigration')
pylab.step(ecology_t2s, 1.0-numpy.arange(0,len(ecology_t2s))*1.0/len(ecology_t2s),'r-',label='strategy mutation')
pylab.step(evolution_t2s, 1.0-numpy.arange(0,len(evolution_t2s))*1.0/len(evolution_t2s),'b-',label='+ fitness mutation')
pylab.step(replicate_evolution_t2s, 1.0-numpy.arange(0,len(replicate_evolution_t2s))*1.0/len(replicate_evolution_t2s),'b-',alpha=0.5)

pylab.ylabel('Strain pairs >=$T$')
pylab.xlabel('Coalescence time, $T$')

pylab.xlim([1e03,1e07])
pylab.semilogx([1],[1],'k.')
pylab.legend(loc='lower left', frameon=False,fontsize=7)


pylab.savefig('gene_swap_t2_distribution.pdf',bbox_inches='tight')






pylab.figure(figsize=(2,1.7))
pylab.ylabel('Species number')

mean = immigration_num_strains[-20:].mean()
stddev = immigration_num_strains[-20:].std()
pylab.plot([1],[mean],'ko')
pylab.plot([1,1],[mean-stddev,mean+stddev],'k-')

mean = ecology_num_strains[-20:].mean()
stddev = ecology_num_strains[-20:].std()
pylab.plot([2],[mean],'ro')
pylab.plot([2,2],[mean-stddev,mean+stddev],'r-')

mean = evolution_num_strains[-20:].mean()
stddev = evolution_num_strains[-20:].std()
pylab.plot([3],[mean],'bo')
pylab.plot([3,3],[mean-stddev,mean+stddev],'b-')

mean = replicate_evolution_num_strains[-20:].mean()
stddev = replicate_evolution_num_strains[-20:].std()
pylab.plot([3.3],[mean],'bo',alpha=0.5)
pylab.plot([3.3,3.3],[mean-stddev,mean+stddev],'b-',alpha=0.5)
pylab.plot([0,3.5],[p,p],'k-',color='0.7',linewidth=0.25)
pylab.ylim([0,100])
pylab.xlim([0.5,3.5])
pylab.xticks([])
pylab.savefig('gene_swap_species_number.pdf',bbox_inches='tight')

pylab.figure(figsize=(5,1))

pylab.plot(consumer_resource_sim.calculate_last_harvest_vector(immigration_measurements),'k-')
pylab.plot(consumer_resource_sim.calculate_last_harvest_vector(ecology_measurements),'r-')

pylab.plot(consumer_resource_sim.calculate_last_harvest_vector(evolution_measurements),'b-')
pylab.plot(consumer_resource_sim.calculate_last_harvest_vector(replicate_evolution_measurements),'b-',alpha=0.5)

pylab.xlabel('Resource, $i$')
pylab.xticks([])
pylab.ylim([0,2])
pylab.savefig('gene_swap_resource_use.pdf',bbox_inches='tight')

