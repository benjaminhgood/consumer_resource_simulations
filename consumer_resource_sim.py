import subprocess
import numpy
from math import log,log10,exp
import sys
from numpy.random import randint

def parse_results(lines):
    
    measurements = []
    
    line_idx = 0
    while line_idx < len(lines) and lines[line_idx].strip()!="":
    
        line_idx+=1 # Header line
        time = float(lines[line_idx].split(":")[1])
        line_idx+=1
        N = float(lines[line_idx].split(":")[1])
        line_idx+=1
        avg_fitness = float(lines[line_idx].split(":")[1])
        line_idx+=1
        scaled_harvest_vector = numpy.array([float(item) for item in (lines[line_idx].split(":")[1].strip()).split()])
        line_idx+=1 
        num_clones = long(lines[line_idx].split(":")[1])
        line_idx+=1 
        
        clones = []
        for virtual_clone_idx in xrange(0,num_clones):
            
            items = lines[line_idx].split(";")
            clone_idx = long(items[0])
            fitness = float(items[1])
            
            subitems = (items[2].strip()).split()
            expression_vector = [0 for h in scaled_harvest_vector]
            for subitem in subitems:
                subsubitems = subitem.split(",")
                resource_idx = long(subsubitems[0])
                level = long(subsubitems[1])
                expression_vector[resource_idx] = level
            
            expression_vector = numpy.array(expression_vector)*1.0
            
            mutations = items[3].strip().split()
            
            clones.append((clone_idx, fitness, expression_vector, mutations))
              
            line_idx+=1
        
        if not lines[line_idx].startswith('SNP'):
            print "Big problem!", lines[line_idx][0:20]
            
        # Now get SNPs
        snps = {}
        for item in (lines[line_idx].split(":")[1].strip()).split():
            subitems = item.split(";")
            mutation = subitems[0].strip()
            freq = float(subitems[1])
            snps[mutation] = freq
        line_idx+=1        
        measurements.append((time,N,avg_fitness,scaled_harvest_vector,clones, snps))    
    
    return measurements


def run_simulation(s, NUb, NUs, pscramble, Nm, p, cv, k, l):

    exec_name = "./consumer_resource_sim"
    sys.stderr.write("Simulating s=%g NUb=%g NUs=%g ps=%g Nm=%g p=%g cv=%g k=%g l=%g...\n" % (s, NUb, NUs, pscramble, Nm, p, cv, k, l))
    
    args = [exec_name, str(s), str(NUb), str(NUs), str(pscramble), str(Nm), str(p), str(cv), str(k), str(l)]
    output = subprocess.Popen(args,shell=False,stdout=subprocess.PIPE).communicate()[0]
    lines = output.strip().split("\n")
    
    return parse_results(lines)

def run_evolve_generalist_simulation(s, NUb, NUs, pscramble, Nm, p, cv, k, l, n0):

    exec_name = "./evolve_generalist_sim"
    sys.stderr.write("Simulating s=%g NUb=%g NUs=%g ps=%g Nm=%g p=%g cv=%g k=%g l=%g n0=%g...\n" % (s, NUb, NUs, pscramble, Nm, p, cv, k, l, n0))
    
    args = [exec_name, str(s), str(NUb), str(NUs), str(pscramble), str(Nm), str(p), str(cv), str(k), str(l), str(n0)]
    output = subprocess.Popen(args,shell=False,stdout=subprocess.PIPE).communicate()[0]
    lines = output.strip().split("\n")
    
    return parse_results(lines)

def run_gene_swap_simulation(s, NUb, NUs, pscramble, Nm, p, cv, k, l, dt):

    exec_name = "./gene_swap_sim"
    sys.stderr.write("Simulating s=%g NUb=%g NUs=%g ps=%g Nm=%g p=%g cv=%g k=%g l=%g dt=%g...\n" % (s, NUb, NUs, pscramble, Nm, p, cv, k, l, dt))
    
    args = [exec_name, str(s), str(NUb), str(NUs), str(pscramble), str(Nm), str(p), str(cv), str(k), str(l), str(dt)]
    output = subprocess.Popen(args,shell=False,stdout=subprocess.PIPE).communicate()[0]
    lines = output.strip().split("\n")
    
    return parse_results(lines)


def run_tworesource_simulation(s,NUb,sstrategy,NUs,cv,l):

    exec_name = "./two_resource_sim"
    sys.stderr.write("Simulating sb=%g NUb=%g ss=%g NUs=%g cv=%g l=%g...\n" % (s, NUb, sstrategy, NUs, cv, l))
    
    args = [exec_name, str(s), str(NUb), str(sstrategy), str(NUs), str(cv), str(l)]
    output = subprocess.Popen(args,shell=False,stdout=subprocess.PIPE).communicate()[0]
    lines = output.strip().split("\n")
    
    return parse_results(lines)
    
    



def calculate_pi_matrix(clones):

    pi_matrix = numpy.zeros((len(clones),len(clones)))*1.0
    for i in xrange(0,len(clones)):
        for j in xrange(0, len(clones)):
            
            mutations1 = set(clones[i][3])
            mutations2 = set(clones[j][3])
            
            pi_matrix[i,j] = len(mutations1)+len(mutations2)-2*len(mutations1 & mutations2)
            pi_matrix[j,i] = pi_matrix[i,j]

    return pi_matrix
    

def calculate_strategy_matrix(clones):
    # probability that a random two clones will have same strategy
    strategy_matrix = numpy.zeros((len(clones),len(clones)))*1.0
    for i in xrange(0,len(clones)):
        for j in xrange(0, len(clones)):
            
            same_trategy = numpy.all(clones[i][2]==clones[j][2])*1.0
            
            strategy_matrix[i,j] = same_strategy
            strategy_matrix[j,i] = same_strategy

    return strategy_matrix
    
def calculate_strategy_homozygosity(clones):

    num_samples = 1000
    n = len(clones)
    num_same_strategy = 0
    for sample_ix in xrange(0,num_samples):
        
        i = randint(0,n)
        j = randint(0,n-1)
        if j>=i:
            j+=1
            
        if numpy.all(clones[i][2]==clones[j][2]):
            num_same_strategy+=1
            
    return num_same_strategy*1.0/num_samples
    
def calculate_snp_timecourse(measurements):

    ts = set([])
    snp_trajectories = {}
    for t,N,X,invh,clones,snps in measurements:
       
       ts.add(t)
       for mutation in snps.keys():
           if mutation not in snp_trajectories:
               snp_trajectories[mutation] = {}
           
           snp_trajectories[mutation][t]=snps[mutation]
    
    ts = [t for t in sorted(ts)]  
    
    final_snp_trajectories = []     
    for mutation in snp_trajectories.keys():
        
        frequencies = []
        for t in ts:
            if t in snp_trajectories[mutation]:
                frequencies.append(snp_trajectories[mutation][t])
            else:
                frequencies.append(0)
        times = numpy.array(ts)
        frequencies = numpy.array(frequencies)
        
        final_snp_trajectories.append((mutation, times, frequencies))
        
    return final_snp_trajectories
    
def is_strategy_mutation(mutation_str):
    items = mutation_str.split(",")
    resource_idx = long(items[3])
    if resource_idx >= 0:
        return True
    else:
        return False
        
 
def calculate_num_strain_timecourse(measurements):
    
    ts = []
    num_strains = []
       
    for t,N,X,invh,clones,snps in measurements:
 
        strains = []
        for clone_idx, fitness, strategy_vector, mutations in clones:
            strains.append(tuple(strategy_vector))
    
    
        strains = set(strains)
        ts.append(t)
        num_strains.append(len(strains))
        
    return numpy.array(ts), numpy.array(num_strains)
    
def calculate_last_harvest_vector(measurements):
    return measurements[-1][3]


