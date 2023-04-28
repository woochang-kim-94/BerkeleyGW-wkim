#!/usr/bin/env python

# Derek Vigil-Fowler, 02/23/17
# Usage: EXX.py x.dat kpt_list

# Script takes sigma_hp.log and kpt_list (just the output of kgrid.x for the sigma wavefunction) 
# files as input, computes the EXX total energy and writes this number out.  

import sys, numpy as np

try:
    sighpfilename=sys.argv[1]
    kptfilename=sys.argv[2]
except:
    sys.exit('Usage: EXX.py sigma_hp.log kpt_list')

print('What is the Fermi level in eV?')
EF=float(input())

print('How many spins in the calculation (1 or 2)?')
nspin=int(input())

# Grab the weights from the output file of kgrid.x, kgrid_list. 
kptfile=open(kptfilename,'r')

kptfile.readline()  # Top line is not useful
line=kptfile.readline()
nkpt=int(line.split()[0])
weightlist=[]
for ikpt in range(int(nkpt)):
   line=kptfile.readline()
   weightlist.append(line.split()[3])
weightlist=np.array(weightlist,dtype="float")
weightlist=weightlist/np.sum(weightlist)

kptfile.close()

# Grab the values of the bare exchange from sigma_hp.log for occupied states

sighpfile=open(sighpfilename,'r')
for iline in range(4):
   line=sighpfile.readline()
nband=int(line.split()[2])
sigx=np.zeros((nkpt,nband,nspin))
sigxweight=np.zeros((nkpt,nband,nspin))

for iline in range(5):
   line=sighpfile.readline()
line=sighpfile.readline()
ntran=int(line.split()[2])
for iline in range(ntran):
   line=sighpfile.readline()

for ikpt in range(nkpt):
   for ispin in range(nspin):
      for iline in range(4):
         line=sighpfile.readline()
      for iband in range(nband):
          line=sighpfile.readline()
          if float(line.split()[1])<EF:
             sigx[ikpt,iband,ispin]=float(line.split()[3])          
             sigxweight[ikpt,iband,ispin]=sigx[ikpt,iband,ispin]*weightlist[ikpt]


EXX=0.5*2.0/(nspin)*sigxweight.sum()

print 'EXX (eV): %13.9f' % (EXX)
print 'EXX (Ry): %13.9f' % (EXX/13.6057)
print 'EXX (Ha): %13.9f' % (EXX/(13.6057*2.0))
