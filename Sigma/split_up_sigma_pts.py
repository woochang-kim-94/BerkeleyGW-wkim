#!/usr/bin/env python

# 3/23/12 - B. D. Malone

# The point of this script is to break up the calculation of sigma into 
# a smaller number of points per calculation. To use it you simply do:

# python split_up_sigma_pts.py wfn0.out [number of kpoints per calculation]

# where wfn0.out is the list of kpoints you want to calculate sigma on
# (as generated from kgrid.x). Additionally, in your directory you will need
# the following files:
# 1). sigma.inp.temp - which is your sigma.inp file you want to use, with all
# the parameters necessary to run the calculation, with the exception of 
# the number_kpoints line and the begin_kpoints block (which this script
# fills in for you depending on your points and division of labor). 
# 2). batch.temp - which is the batch file you want to submit in each directory.
# 3). eps0mat, epsmat, WFN_inner, RHO, and vxc.dat - you should know what 
# these are.

# After a calculation, you will have a bunch of directories with your resulting
# data in them. If you run the script again, but giving the 'collect' option
# at startup, the script will go through each directory and collect the 
# information from the sigma_hp.log files, append them in the appropriate 
# manner, and write the resulting file to sigma_hop.log.combined in this
# directory. This file can be used with other utilities such as sig2wan.x, etc.

# FAQ:
# Q1: Why would I want to use this script?
# A1: For some systems sigma can be the most expensive part of the calculation,
# and breaking up the points can save time, get you into the queue faster, 
# qualify you for cheaper queues, etc.
# Q2: Why would I NOT want to use this script?
# A2: If you are doing a system with lots of G-vectors (say, in non-bulk 
# calculations or with semicore electrons) splitting up the calculation 
# can waste some time because currently sigma only reads the WFN_inner file
# on one processor. Thus you are forcing it to read these large wavefunctions
# over and over again, which can in some cases be a significant part of the 
# total wallclock time. Parallel IO (in development) would lessen this concern.
# Q3: I want my jobs to submit automatically, how do I do this?
# A3: Remove the commented-out line at the bottom where you see
# os.system(command) 

# An additional note is that when this script is run in 'collect' mode, it will
# compare eqp0 and eqp1 for different kpoints of interest and will write out
# a message if these things are greater than 10 meV. Currently you can see
# that this message only cares if you are at the Gamma, X, or L points in 
# zincblende materials. Is this a problem? That's for you to decide. In
# many materials the linearization of Sigma is pretty accurate, but if 
# eqp0 and eqp1 are off, you may want to recalculate sigma at a more 
# accurate quasiparticle energy using the eqp.py script so that you can 
# have highly accurate values, at least for the kpoints you care most about.

# line below makes compatible with Python versions less than 3.0
from __future__ import print_function


import sys,shutil,os,glob

# Note: The filerootname is simply the root name of what all your jobs
# will be called when sent to the scheduler. For example, with a filerootname
# of 'epsilon', your jobs will be called epsilon_1, epsilon_2,etc.
filerootname='sigma'  


try:
    infilename=sys.argv[1]
    ptbreaks=int(sys.argv[2])
except:
    print('Please specify kpoint file and how many points you want per directory/calculation on commandline');sys.exit(1)

print('Doing',ptbreaks,'points in each packet')            

if os.path.exists('batch.temp')==False:
    print('You need a batch.temp file in this directory');sys.exit(1)
if os.path.exists('sigma.inp.temp')==False:
    print('You need a sigma.inp.temp file in this directory');sys.exit(1)
if os.path.exists('WFN_inner')==False:
    print('You need a WFN_inner file in this directory');sys.exit(1)
if os.path.exists('RHO')==False:
    print('You need a RHO file in this directory');sys.exit(1)
if os.path.exists('vxc.dat')==False and os.path.exists('VXC')==False:
    print('You need a vxc.dat or VXC file in this directory');sys.exit(1)
if os.path.exists('epsmat')==False and os.path.exists('epsmat.h5')==False:
    print('You need an epsmat or epsmat.h5 file in this directory');sys.exit(1)
if os.path.exists('eps0mat')==False and os.path.exists('eps0mat.h5')==False:
    print('You need an eps0mat or eps0mat.h5 file in this directory');sys.exit(1)
if os.path.exists(infilename)==False:
    print('The file', infilename,'isn\'t really in this directory, is it?')

print('Do you want to run (1) or collect (2)?')
mode=int(input())

infile=open(infilename,'r')


pcounter=1
outfile=open('kpts.p'+str(pcounter)+'.dat','w')
infile.readline() # just the K_POINTS crystal header
numkpts=int(infile.readline().split()[0])
ptsfilelist=[]
if numkpts>=ptbreaks:
    print('K_POINTS crystal', file=outfile)
    print(str(ptbreaks), file=outfile)
    ptsfilelist.append(ptbreaks)
else:
    print('K_POINTS crystal', file=outfile)
    print(str(numkpts), file=outfile)
    ptsfilelist.append(numkpts)
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
            ptsfilelist.append(ptbreaks)
        else:
            print('K_POINTS crystal', file=outfile)
            print(str(numkpts-initvalue), file=outfile)
            ptsfilelist.append(numkpts-initvalue)
outfile.close()


# now that we've created all the pcounter files, we just add these points 
# to all the sigma files and create directories. 
totalp=pcounter

def set_up_sigma_file(iter):
    infile=open('sigma.inp.temp','r')
    outfile=open('sigma.inp','w')
    for line in infile:
        print(line, end='', file=outfile)
    infile.close()
    infile2=open('kpts.p'+str(iter)+'.dat','r')
    infile2.readline()
    ptsinfile=infile2.readline().split()[0]
    print('begin kpoints', file=outfile)
    while True:
        line=infile2.readline().strip()
        if line=='':  # EOF
            break
        theline=line.split()
        thept=[theline[0],theline[1],theline[2]]
        print('  {0} {1} {2} 1.0'.format(thept[0],thept[1],thept[2]), file=outfile)
    infile2.close()
    print('end', file=outfile)
    outfile.close()

def modify_batch(iter):
    infile=open('batch.temp','r')
    outfile=open('batch','w')
    for line in infile:
        if line.find('#PBS -N')!=-1:
            print('#PBS -N '+filerootname+'_'+str(iter), file=outfile)
        elif line.find('cd $PBS_O_WORKDIR')!=-1:
            curdir=os.getcwd()
            path=os.path.join(curdir,'p'+str(iter))
            print('cd '+path, file=outfile)
        else:
            print(line, end='', file=outfile)
    infile.close()
    outfile.close()

if mode==1:
    for i in range(1,totalp+1):
        directname='p'+str(i)
        print('Setting up directory',directname)
        try:
            os.mkdir(directname)
        except:
            pass # probably already there
        modify_batch(i)
        # get new sigma.inp file
        set_up_sigma_file(i)
        shutil.move('batch',directname)
        shutil.move('sigma.inp',directname)
        os.symlink('../WFN_inner',os.path.join(directname,'WFN_inner'))
        os.symlink('../RHO',os.path.join(directname,'RHO'))
        if os.path.exists('vxc.dat')==True:
            os.symlink('../vxc.dat',os.path.join(directname,'vxc.dat'))
        elif os.path.exists('VXC')==True:
            os.symlink('../VXC',os.path.join(directname,'VXC'))
        if os.path.exists('epsmat')==True:
            os.symlink('../epsmat',os.path.join(directname,'epsmat'))
        elif os.path.exists('epsmat.h5')==True:
            os.symlink('../epsmat.h5',os.path.join(directname,'epsmat.h5'))
        if os.path.exists('eps0mat')==True:
            os.symlink('../eps0mat',os.path.join(directname,'eps0mat'))
        elif os.path.exists('eps0mat.h5')==True:
            os.symlink('../eps0mat.h5',os.path.join(directname,'eps0mat.h5'))
        os.symlink('../eqp.py',os.path.join(directname,'eqp.py'))
        batchpath=os.path.join(directname,'batch')
        command='qsub '+batchpath
#        os.system(command)

elif mode==2:
    outfile=open('sigma_hp.log.combined','w')
    for i in range(1,totalp+1):
        directname='p'+str(i)
        print('Looking in directory',directname)
        filename=os.path.join(directname,'sigma_hp.log')
        numpts=ptsfilelist[i-1]
        infile=open(filename,'r')
        # we want the symmetries and header only for i=1
        while True:
            line=infile.readline()
            if line.find('band_index')!=-1:
                minband=int(line.split()[1])
                maxband=int(line.split()[2])
                numbands=maxband-minband+1
            if line.find('ntran')!=-1:
                numsym=int(line.split()[2])
                if i==1:
                    print(line, end=' ', file=outfile)
                break # go manual since we know what we want now
            if i==1:
                print(line, end=' ', file=outfile)
        for j in range(0,numsym+1):
            line=infile.readline()
            if i==1:
                print(line, end=' ', file=outfile)
        # now we grab k-points
        for k in range(0,numpts):
            for l in range(0,3):
                line=infile.readline()
                if l==0:
                    theline=line.split()
                    kptpos=[float(theline[2]),float(theline[3]),float(theline[4])]
                print(line, end=' ', file=outfile)
            for l in range(0,numbands+1):  # this includes a blank line
                line=infile.readline()     # after the kpoint
                if l==numbands:
                    pass
                else:  # we want to compare how close we are to 
                       # self-consistency in the linearization (only if it is
                       # the Gamma, L, or X point
                    ofinterest=[[0,0,0],[0,0.5,0.5],[0,0,0.5]]
                    for ptitem in ofinterest:
                        if abs(float(kptpos[0]-ptitem[0]))<1E-10 and abs(float(kptpos[1]-ptitem[1]))<1E-10 and abs(float(kptpos[2]-ptitem[2]))<1E-10:
                            theline=line.split()
                            eqp1=float(theline[9])
                            eqp0=float(theline[8])
                            if abs(eqp1-eqp0)>0.010:
                                print('Difference in eigenvalues of',abs(eqp1-eqp0),'for band', l+1,'kpt=',kptpos)
                print(line, end=' ', file=outfile)
        infile.close()
    outfile.close()


else:
    print('Don\'t take it so hard, not everyone is cut out for computational science')

# delete unneeded kpt files

for filename in glob.glob('kpts.p*.dat'):
    os.remove(filename)
