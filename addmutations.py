""" 
add mutations and compute sfs using msprime from 
slim tree sequence data
"""
import msprime
import pyslim
import numpy as np
import subprocess
import glob
import random

n=136

# first collect the tree file names :  ls treefiles* > treeseqfiles
# read in the tree files
# file = glob.glob('wfposa_*.trees')
# for l in file
### works but needs  treeseqfiles :
#with open('./treeseqfiles') as file:
#    f = [l.rstrip() for l in file]
afs = np.zeros(n-1, dtype=float)
# change the prefix  to read in the correct tree files
file = glob.glob('./wfnegm005shape2posh1s0_5Ne4me-12_*.trees')
f = [l.rstrip() for l in file]

s = 0
for i in f:
    # ts = pyslim.load("./prufatreeseq.trees")
    ts = pyslim.load(i)
    ts = ts.simplify()

# check if mrca reached at all trees
    for t in ts.trees():
        assert t.num_roots == 1, ("not coals on {} to {} ".format(t.interval[0], t.interval[1]))


    mutated = msprime.mutate(ts, rate=1e-6, random_seed=random.randint(4940,3930302902), keep=True)
    sample = np.arange(10000)
    np.random.shuffle(sample)
    print(sample[:136])
    m = mutated.allele_frequency_spectrum([sample[:n]], polarised=True,
                                                  span_normalise=False)
    print(m.size)
    s = s + np.sum(m[1:n])
    afs = afs + (m[1:n]/float(np.sum(m[1:n])))


print(afs/float(len(f)))
print(s/float(len(f)))
# record the sfs into a text file 
np.savetxt('./wfnegm005shape2posh1s0_5Ne4me-12_sfs_resout', afs/float(len(f)))
# np.savetxt('./skewedposC_sfs_resout', afs/float(len(f)) )
# mutated.dump("./mutationsover.trees")
