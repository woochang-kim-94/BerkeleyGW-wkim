#!/usr/bin/env python

# 3/23/12 - B. D. Malone

# The purpose of this script is to break up the q-points needed in an epsilon
# calculation into smaller groups so that they can be computed in parallel.
# The script also sets up the appropriate input files and soft-links to 
# the epsmat files into a directory called 'put_together', from which you
# can easily merge together the resulting epsmat files using the epsmat_merge
# utility provided in BGW.
# An additional benefit of reducing the size of individual jobs is that the
# jobs will go through the queue much faster generally, and it may also become
# feasible for them to run in lower priority queues, saving CPU hours in your
# allocation.

# In the directory where you run this script you need 5 files:
# 1). WFN 
# 2). WFNq
# 3). your list of q-points (as given from kgrid.x). For the silicon example
# in the BGW distribution this would be the WFN_co.out file. This file can be
# called whatever you like, and you will specify it on the commandline 
# Note: the point at $\Gamma$ should be given as (0,0,0). This script will 
# take care of the small qshift. 
# 4). batch.temp - which is the batch file you want submitted to the 
# scheduler. The script will automatically detect the '#PBS -N' line and 
# change the job name to something useful for identification.
# 5). epsilon.inp.temp - This is a normal epsilon.inp input file that you want
# to run, but without the q-point information. For the silicon example in the 
# BGW distribution, the information you would want in this file would simply be
# the following:
# ----------------------------
# epsilon_cutoff 10.0
# number_bands 29
# band_occupation 4*1 25*0
# ----------------------------

# How to run:
# You simple need to do 'python split_up_epsilon_pts.py [your WFN_co.out file 
# from kgrid] [the # of points you want in each directory]'. For example, if 
# I have my points in a file called 'wfn0.out' and I want 1 q-point per 
# calculation, I would simply do:
# python split_up_epsilon_pts.py wfn0.out 1 

# line below makes compatible with Python versions less than 3.0
from __future__ import print_function

 


import sys,shutil,os,glob

# Note: The filerootname is simply the root name of what all your jobs
# will be called when sent to the scheduler. For example, with a filerootname
# of 'epsilon', your jobs will be called epsilon_1, epsilon_2,etc.
filerootname='epsilon'  


try:
    infilename=sys.argv[1]
    ptbreaks=int(sys.argv[2])
except:
    sys.exit('Error: please specify kpoint file and how many points you want per directory/calculation on commandline')

qshift = input("Enter your q-shift vector, assuming it is a small shift to the b3 vector: ")

print('Doing',ptbreaks,'points in each packet')            

if os.path.exists('batch.temp')==False:
    sys.exit('Error: You need a batch.temp file in this directory')
if os.path.exists('epsilon.inp.temp')==False:
    sys.exit('Error: You need an epsilon.inp.temp file in this directory')
if os.path.exists('WFN')==False and os.path.exists('WFN.h5')==False:
    sys.exit('Error: You need a WFN or WFN.h5 file in this directory')
if os.path.exists('WFNq')==False and os.path.exists('WFN.h5')==False:
    sys.exit('Error: You need a WFNq or WFNq.h5 file in this directory')
if os.path.exists(infilename)==False:
    sys.exit('Error: The file', infilename,'isn\'t really in this directory, is it?')


infile=open(infilename,'r')


pcounter=1
outfile=open('kpts.p'+str(pcounter)+'.dat','w')
infile.readline() # just the K_POINTS crystal header
numkpts=int(infile.readline().split()[0])
if numkpts>=ptbreaks:
    print('K_POINTS crystal', file=outfile)
    print(str(ptbreaks), file=outfile)
else:
    print('K_POINTS crystal', file=outfile)
    print(str(numkpts), file=outfile)
kptcounter=0
initvalue=0
while True:
    line=infile.readline()
    if line=='':  # EOF
        break
    kptcounter=kptcounter+1
    if numkpts-initvalue>=ptbreaks:
        print(line, end=' ', file=outfile)
    else:
        print(line, end=' ', file=outfile)
    if kptcounter%ptbreaks==0:
        outfile.close()
        initvalue=kptcounter
        if numkpts==initvalue:
            break
        pcounter=pcounter+1
        outfile=open('kpts.p'+str(pcounter)+'.dat','w')
        if numkpts-initvalue>=ptbreaks:
            print('K_POINTS crystal', file=outfile)
            print(ptbreaks, file=outfile)
        else:
            print('K_POINTS crystal', file=outfile)
            print(str(numkpts-initvalue), file=outfile)
outfile.close()


# now that we've created all the pcounter files, we just add these points 
# to all the chi0 files and create directories. 
totalp=pcounter

def set_up_chi0_file(iter):
    infile=open('epsilon.inp.temp','r')
    outfile=open('epsilon.inp','w')
    for line in infile:
        if line.find('epsilon_cutoff')!=-1:
            epscutoff=line.split()[1]
        print(line, end='', file=outfile)
    infile.close()
    infile2=open('kpts.p'+str(iter)+'.dat','r')
    infile2.readline()
    ptsinfile=infile2.readline().split()[0]
    print('begin qpoints', file=outfile)
    while True:
        line=infile2.readline().strip()
        if line=='':  # EOF
            break
        theline=line.split()
        thept=[theline[0],theline[1],theline[2]]
        if abs(float(thept[0]))<1E-10 and abs(float(thept[1]))<1E-10 and abs(float(thept[2]))<1E-10:
            # Modify small qshift in line below
            print('  {0} {1} {2:.12f} 1.0 1'.format(thept[0],thept[1],qshift), file=outfile)
        else:
            print('  {0} {1} {2} 1.0 0'.format(thept[0],thept[1],thept[2]), file=outfile)
    infile2.close()
    print('end', file=outfile)
    outfile.close()
    return epscutoff

def setup_epsmat_merge(epscutoff):
    print('Setting up epsmat_merge stuff')
    try:
        os.mkdir('put_together')
    except:
        pass # probably already there
    infile=open(infilename,'r')
    outfile=open(os.path.join('put_together','epsmat_merge.inp'),'w')
    infile.readline() # just header
    numkpts=int(infile.readline().split()[0])
    infile.readline()  # get rid of the Gamma-like point
    print(epscutoff,numkpts-1, file=outfile)
    while True:
        line=infile.readline()
        if line=='': # EOF
            break
        theline=line.split()
        print(theline[0],theline[1],theline[2],'1.0', file=outfile)
    infile.close()
    if ptbreaks==1:
        istart=2
        print(totalp-1, file=outfile)
    elif ptbreaks>1:
        istart=1
        print(totalp, file=outfile)
    else:
        sys.exit('Wait, what?')
    for i in range(istart,totalp+1):
        print('epsmat_p'+str(i)+'.h5', file=outfile)
        os.symlink('../p'+str(i)+'/epsmat.h5',os.path.join('put_together','epsmat_p'+str(i)+'.h5'))
    outfile.close()
        
        
def modify_batch(iter):
    infile=open('batch.temp','r')
    outfile=open('batch','w')
    for line in infile:
        if line.find('#PBS -N')!=-1:
            print('#PBS -N '+filerootname+'_'+str(iter), file=outfile)
        elif line.find('#SBATCH -J')!=-1:
            print('#SBATCH -J '+filerootname+'_'+str(iter), file=outfile)
        #elif line.find('cd $PBS_O_WORKDIR')!=-1:
        #    curdir=os.getcwd()
        #    path=os.path.join(curdir,'p'+str(iter))
        #    print >> outfile, 'cd '+path
        else:
            print(line, end='', file=outfile)
    infile.close()
    outfile.close()


for i in range(1,totalp+1):
    directname='p'+str(i)
    print('Setting up directory',directname)
    try:
        os.mkdir(directname)
    except:
        pass # probably already there
    modify_batch(i)
    # get new epsilon.inp file
    epscutoff=set_up_chi0_file(i)
    shutil.move('batch',directname)
    shutil.move('epsilon.inp',directname)
    if os.path.exists('WFN')==True and os.path.exists('WFN.h5')==False:
        os.symlink('../WFN',os.path.join(directname,'WFN'))
    else:
        os.symlink('../WFN.h5',os.path.join(directname,'WFN.h5'))
    if os.path.exists('WFNq')==True and os.path.exists('WFNq.h5')==False:
        os.symlink('../WFNq',os.path.join(directname,'WFNq'))
    else:
        os.symlink('../WFNq.h5',os.path.join(directname,'WFNq.h5'))
    

setup_epsmat_merge(epscutoff)

# delete unneeded kpt files

for filename in glob.glob('kpts.p*.dat'):
    os.remove(filename)


