import consumer_resource_sim
import numpy
import pylab
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib

##############################
#
# First set up the figure
#
##############################

mpl.rcParams['font.size'] = 4.0
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

generation_ticks = [10000*i for i in xrange(0,7)]
generation_ticklabels = ['%dk' % (10*i) for i in xrange(0,7)]
generation_ticklabels[0] = '0'

generation_ticks = [5000*i for i in xrange(0,13)]
generation_ticklabels = ['%dk' % (5*i) for i in xrange(0,13)]

frequency_ticks = [0,0.2,0.4,0.6,0.8,1.0]

pylab.figure(1,figsize=(5.5,3.5))
fig = pylab.gcf()

# make three panels panels
outer_grid  = gridspec.GridSpec(4, 1, height_ratios=[1,1,1,1], hspace=0.25)

axis1 = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(axis1)

#fixed_axis.set_title('%s basal clade' % parse_file.get_pretty_name(focal_population),loc='right',fontsize=4,y=0.91)
    

axis1.spines['top'].set_visible(False)
axis1.spines['right'].set_visible(False)
axis1.get_xaxis().tick_bottom()
axis1.get_yaxis().tick_left()
axis1.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
axis1.get_xaxis().set_tick_params(direction='out',length=3,pad=1)

axis1.set_xticks(generation_ticks)
axis1.set_xticklabels([])
axis1.set_yticks(frequency_ticks)

axis1.set_ylabel('Allele frequency, $f(t)$')


axis1.set_ylim([0,1.01])
axis1.set_xlim([0,60100])

axis2 = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(axis2)

#axis2.set_title('%s major clade' % parse_file.get_pretty_name(focal_population),loc='right',fontsize=4,y=0.91)


axis2.spines['top'].set_visible(False)
axis2.spines['right'].set_visible(False)
axis2.get_xaxis().tick_bottom()
axis2.get_yaxis().tick_left()
axis2.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
axis2.get_xaxis().set_tick_params(direction='out',length=3,pad=1)

axis2.set_xticks(generation_ticks)
axis2.set_xticklabels([])
axis2.set_yticks(frequency_ticks)
axis2.set_ylabel('Allele frequency, $f(t)$')

axis2.set_ylim([0,1.01])
axis2.set_xlim([0,60100])


axis3 = plt.Subplot(fig, outer_grid[2])
fig.add_subplot(axis3)

#minority_axis.set_title('%s minor clade' % parse_file.get_pretty_name(focal_population),loc='right',fontsize=4,y=0.91)


axis3.spines['top'].set_visible(False)
axis3.spines['right'].set_visible(False)
axis3.get_xaxis().tick_bottom()
axis3.get_yaxis().tick_left()
axis3.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
axis3.get_xaxis().set_tick_params(direction='out',length=3,pad=1)

axis3.set_xticks(generation_ticks)
axis3.set_xticklabels([])
axis3.set_yticks(frequency_ticks)

axis3.set_ylabel('Allele frequency, $f(t)$')


axis3.set_ylim([0,1.01])
axis3.set_xlim([0,60100])


axis4 = plt.Subplot(fig, outer_grid[3])
fig.add_subplot(axis4)

#axis4.set_title('%s extinct / unclassified' % parse_file.get_pretty_name(focal_population),loc='right',fontsize=4,y=0.91)


axis4.set_xlabel('Generation, $t$')

axis4.spines['top'].set_visible(False)
axis4.spines['right'].set_visible(False)
axis4.get_xaxis().tick_bottom()
axis4.get_yaxis().tick_left()
axis4.get_yaxis().set_tick_params(direction='out',length=3,pad=1)
axis4.get_xaxis().set_tick_params(direction='out',length=3,pad=1)

axis4.set_xticks(generation_ticks)
axis4.set_xticklabels(generation_ticklabels)
axis4.set_yticks(frequency_ticks)

axis4.set_ylabel('Allele frequency, $f(t)$')


axis4.set_ylim([0,1.01])
axis4.set_xlim([0,60100])

# Ecology parameters!
sstrategy = -0.02
NUs=0.01
cv=0.05
l=4

#################
#
# Pure ecology one!
#
##################
s=1e-02
NUb=0.005


measurements = consumer_resource_sim.run_tworesource_simulation(s, NUb, sstrategy, NUs, cv, l)

mutation_data = consumer_resource_sim.calculate_snp_timecourse(measurements)

for mutation,times,freqs in mutation_data:

    if (freqs<0.05).all():
        continue
        
    if consumer_resource_sim.is_strategy_mutation(mutation):
        color='r'
        zorder=2
    else:
        color='b'
        zorder=1
    axis1.plot(times,freqs,'-',color=color,zorder=zorder)

axis1.plot([1000000,1000001],'b-',label='General fitness mutation')
axis1.plot([1000000,1000001],'r-',label='Strategy mutation')
axis1.legend(frameon=False)

print measurements[0][4][0][2]
print measurements[-1][4][0][2]

#################
#
# Dynamic Evolution one
#
##################
s=5e-03
NUb=1

measurements = consumer_resource_sim.run_tworesource_simulation(s, NUb, sstrategy, NUs, cv, l)

mutation_data = consumer_resource_sim.calculate_snp_timecourse(measurements)

for mutation,times,freqs in mutation_data:

    if (freqs<0.05).all():
        continue
        
    if consumer_resource_sim.is_strategy_mutation(mutation):
        color='r'
        zorder=2
    else:
        color='b'
        zorder=1
    axis2.plot(times,freqs,'-',color=color,zorder=zorder)

print measurements[0][4][0][2]
print measurements[-1][4][0][2]

#################
#
# Evolution one
#
##################
s=5e-04
NUb=3000

measurements = consumer_resource_sim.run_tworesource_simulation(s, NUb, sstrategy, NUs, cv, l)

mutation_data = consumer_resource_sim.calculate_snp_timecourse(measurements)

for mutation,times,freqs in mutation_data:

    if (freqs<0.05).all():
        continue
        
    if consumer_resource_sim.is_strategy_mutation(mutation):
        color='r'
        zorder=2
    else:
        color='b'
        zorder=1
    axis3.plot(times,freqs,'-',color=color,zorder=zorder)

print measurements[0][4][0][2]
print measurements[-1][4][0][2]

#################
#
# Evolution >> ecology one
#
##################
s=2e-02
NUb=0.3

measurements = consumer_resource_sim.run_tworesource_simulation(s, NUb, sstrategy, NUs, cv, l)

mutation_data = consumer_resource_sim.calculate_snp_timecourse(measurements)

for mutation,times,freqs in mutation_data:

    if (freqs<0.05).all():
        continue
        
    if consumer_resource_sim.is_strategy_mutation(mutation):
        color='r'
        zorder=2
    else:
        color='b'
        zorder=1
    axis4.plot(times,freqs,'-',color=color,zorder=zorder)

print measurements[0][4][0][2]
print measurements[-1][4][0][2]


fig.savefig('tworesource_trafic_plots.pdf',bbox_inches='tight')
    