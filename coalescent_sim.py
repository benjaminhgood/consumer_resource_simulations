import subprocess
import numpy
from math import log,log10,exp
import sys

disabled=False
zeng_disabled=False


if sys.platform.startswith("linux") or sys.platform.startswith("darwin"):
    LINUX = True
else:
    LINUX = False

def calculate_time_spectrum(Ns,NUd,n,num_samples=1000):
    if disabled:
        return numpy.zeros(n-1)
    if LINUX:
        exec_name = "./simulate_time_spectrum"
    else:
        exec_name = "./simulate_time_spectrum.exe"
    print "Calculating", Ns, NUd, n
    args = [exec_name,str(Ns),str(NUd),str(n),str(num_samples)]
    output = subprocess.Popen(args,shell=False,stdout=subprocess.PIPE).communicate()[0]
    fs = numpy.array([float(item) for item in output.split()])
    return fs

def calculate_2d_time_spectrum(Ns1,NUd1,Ns2,NUd2,n,num_samples=1000):
    if disabled:
        return numpy.zeros(n-1)
    if LINUX:
        exec_name = "./simulate_2d_time_spectrum"
    else:
        exec_name = "./simulate_2d_time_spectrum.exe"
    
    print "Calculating", Ns1, NUd1, Ns2, NUd2, n
    args = [exec_name,str(Ns1),str(NUd1),str(Ns2),str(NUd2),str(n),str(num_samples)]
    output = subprocess.Popen(args,shell=False,stdout=subprocess.PIPE).communicate()[0]
    fs = numpy.array([float(item) for item in output.split()])
    return fs
    
def calculate_ttotn(Ns,NUd,n,num_samples=1000):
    return calculate_time_spectrum(Ns,NUd,n,num_samples).sum()

def calculate_2d_ttotn(Ns1,NUd1,Ns2,NUd2,n,num_samples=1000):
    return calculate_2d_time_spectrum(Ns1,NUd1,Ns2,NUd2,n,num_samples).sum()

def calculate_ttotn_ratio(Ns,NUd,n,num_samples=1000):
    if disabled:
        return numpy.zeros(n-1)
    if LINUX:
        exec_name = "./simulate_ratio"
    else:
        exec_name = "./simulate_ratio.exe"
    print "Calculating", Ns, NUd, n
    args = [exec_name,str(Ns),str(NUd),str(n),str(num_samples)]
    output = subprocess.Popen(args,shell=False,stdout=subprocess.PIPE).communicate()[0]
    return float(output) 

def calculate_t2(Ns,NUd,num_samples=1000):
    return calculate_ttotn(Ns,NUd,2,num_samples)/2

def calculate_2d_t2(Ns1,NUd1,Ns2,NUd2,num_samples=1000):
    return calculate_2d_ttotn(Ns1,NUd1,Ns2,NUd2,2,num_samples)/2


def calculate_recombination_time_spectrum(Ns,NUd,NR,n,num_samples=1000):
    if zeng_disabled:
        return numpy.zeros(n-1)
    if LINUX:
        exec_name = "./recombination_time_spectrum"
    else:
        exec_name = "./simulate_time_spectrum.exe"
    print "Calculating", Ns, NUd, NR, n
    args = [exec_name,str(Ns),str(NUd),str(NR),str(n),str(num_samples)]
    output = subprocess.Popen(args,shell=False,stdout=subprocess.PIPE).communicate()[0]
    fs = numpy.array([float(item) for item in output.split()])
    return fs

def calculate_recombination_ttotn(Ns,NUd,NR,n,num_samples=1000):
    fs = calculate_recombination_time_spectrum(Ns,NUd,NR,n,num_samples)
    return fs.sum()

def calculate_recombination_t2(Ns,NUd,NR,num_samples=1000):
    return calculate_recombination_ttotn(Ns,NUd,NR,2,num_samples)/2


    

