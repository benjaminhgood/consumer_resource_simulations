import pylab
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy
import matplotlib.colors as colors
import matplotlib.cm as cmx

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

upper_n0 = 1000
lower_n0 = 1

upper_immigration_measurements = consumer_resource_sim.run_evolve_generalist_simulation(s, 0,   0,   0, Nm, p, cv, k, l, upper_n0)
lower_immigration_measurements = consumer_resource_sim.run_evolve_generalist_simulation(s, 0,   0,   0, Nm, p, cv, k, l, lower_n0)

upper_ecology_measurements =     consumer_resource_sim.run_evolve_generalist_simulation(s, 0,   Nm,  1, 0,  p, cv, k, l, upper_n0)
lower_ecology_measurements =     consumer_resource_sim.run_evolve_generalist_simulation(s, 0,   Nm,  1, 0,  p, cv, k, l, lower_n0)

upper_evolution_measurements =         consumer_resource_sim.run_evolve_generalist_simulation(s, NUb, NUs, 0, 0,  p, cv, k, l, upper_n0)
lower_evolution_measurements =         consumer_resource_sim.run_evolve_generalist_simulation(s, NUb, NUs, 0, 0,  p, cv, k, l, lower_n0)


pylab.figure(figsize=(5,1))

pylab.plot(consumer_resource_sim.calculate_last_harvest_vector(upper_immigration_measurements),'k-')
pylab.plot(consumer_resource_sim.calculate_last_harvest_vector(lower_immigration_measurements),'k-')

pylab.plot(consumer_resource_sim.calculate_last_harvest_vector(upper_ecology_measurements),'r-')
pylab.plot(consumer_resource_sim.calculate_last_harvest_vector(lower_ecology_measurements),'r-')

pylab.plot(consumer_resource_sim.calculate_last_harvest_vector(upper_evolution_measurements),'b-')
pylab.plot(consumer_resource_sim.calculate_last_harvest_vector(lower_evolution_measurements),'b-')

pylab.xlabel('Resource, $i$')
pylab.xticks([])
pylab.ylim([0,2])
pylab.savefig('generalist_resource_use.pdf',bbox_inches='tight')

pylab.figure(figsize=(5,1.7))

pylab.xlabel('Generations, $t$')
pylab.ylabel('Species number')

pylab.loglog([1,1e06],[p,p],'k-',linewidth=0.25,color='0.7')


ts,num_strains = consumer_resource_sim.calculate_num_strain_timecourse(upper_immigration_measurements)

pylab.plot(ts,num_strains,'k-',label='Immigration')
ts,num_strains = consumer_resource_sim.calculate_num_strain_timecourse(lower_immigration_measurements)
pylab.plot(ts,num_strains,'k-')

ts,num_strains = consumer_resource_sim.calculate_num_strain_timecourse(upper_ecology_measurements)
pylab.plot(ts,num_strains,'r-',label='strategy mutations')
ts,num_strains = consumer_resource_sim.calculate_num_strain_timecourse(lower_ecology_measurements)
pylab.plot(ts,num_strains,'r-')

ts,num_strains = consumer_resource_sim.calculate_num_strain_timecourse(upper_evolution_measurements)
pylab.plot(ts,num_strains,'b-',label='+ fitness mutations')
ts,num_strains = consumer_resource_sim.calculate_num_strain_timecourse(lower_evolution_measurements)
pylab.plot(ts,num_strains,'b-')

pylab.legend(frameon=False,fontsize=7)

pylab.ylim([1,1e03])
pylab.xlim([1,1e06])


pylab.savefig('generalist_species_number.pdf',bbox_inches='tight')
#pylab.ylim([0,60])

#pylab.show()