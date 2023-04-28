#!/usr/bin/env python

# 3/23/12 - B. D. Malone

# The purpose of this is to break up the kpoint files produced by kgrid.x so
# that wavefunctions can be calculated for smaller subsets of kpoints than all
# of the ones needed in one shot. This script is intended currently only for 
# Quantum-ESPRESSO input files.
# What the script does is the following:
# 1). Takes your kpoint list and appends these kpoints (split up however you 
#        like) onto your base QE input file and puts this in its own directory
# 2). It also puts in your *UPF pseudo files, a modified batch script, and your
#     pp.in file into that directory. The script will also set up your 
#     [prefix].save directory, link to your charge-density.dat and data-xml.dat
#     files. So everything will be set up for you to 'cd' to that directory
#     and submit your job (this script is easily modified to automatically
#     submit your jobs if you are brave). 
# 3). The script also sets up a 'put_together' directory, sets up your 
#     wfnmerge.x input file, and creates soft-links to your soon-to-exist 
#     WFN files. All you will have to do at this point is run wfnmerge.x.

# Q: So, what do I need in my directory to run this amazing script?
# A: Excellent question. You will need the following:
# 1) Your kgrid.x output file, let's call this kgrid.out
# 2). Your pseudopotential files, ending with 'UPF'
# 3). Your normal pp.in file for pw2bgw
# 4). Your base QE input file, which is the file you would want run in 
#     each directory but with no K-POINT information (i.e., just delete
#     everything beginning with the 'K_POINTS crystal' keyword). Call this file
#     [prefix].in.temp
# 5). Your batch input file, called batch.temp
# 6). Your charge-density.dat and data-xml.dat files that will be put into your
#        [prefix].save directories in each subdirectory.
# 

# HOW TO RUN:
# If you would like 2 k-points in each directory (i.e., 2 kpoints per batch 
# script submission), and you are running a graphene calculation with a prefix
# of 'gf' you would run the following command:
# python split_up_dft_kpts.py kgrid.out 2 gf

# And that's it. You just need to go submit your jobs or modify this script to 
# submit them for you.

# FAQ:
# Q: If -npools works, then why do I need this script?
# A: Let's say it does. Still, breaking up jobs explicitly into smaller k-point
# batches requires smaller jobs submitted to the queue for the same degree of
# parallelization, so you may go through the queue much faster than if you 
# simply submit a huge job to the queue using -npools.

# line below makes compatible with Python versions less than 3.0
from __future__ import print_function
import sys,glob,shutil,os


try:
    infilename=sys.argv[1]
    ptbreaks=int(sys.argv[2])
    prefix=sys.argv[3]
except:
    sys.exit('Error: Please specify kpoint file and how many points you want on commandline + the prefix')

print('Doing',ptbreaks,'points in each packet')            

filerootname='SIwfn06c'

if os.path.exists('charge-density.dat')==False:
    sys.exit('Error: You must have charge density in this directory')
if os.path.exists('data-file.xml')==False:
    sys.exit('Error: You must have data-file.xml in this directory')




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

totalp=pcounter

def set_up_pw_file(iter):
    infile=open(prefix+'.in.temp','r')
    outfile=open(prefix+'.in','w')
    for line in infile:
        if line.find('K_POINTS')!=-1:
            sys.exit('Error: The {0} file should not already have kpoints in it'.format(prefix+'.in.temp'))
        print(line, end=' ', file=outfile)
    infile.close()
    infile2=open('kpts.p'+str(iter)+'.dat','r')
    for line in infile2:
        print(line, end=' ', file=outfile)
    infile2.close()
    outfile.close()

def modify_batch(iter):
    infile=open('batch.temp','r')
    outfile=open('batch','w')
    for line in infile:
        if line.find('#PBS -N')!=-1:
            print('#PBS -N '+filerootname+'_'+str(iter), file=outfile)
        else:
            print(line, end=' ', file=outfile)
    infile.close()
    outfile.close()

def setup_wfn_merge():
    print('Setting up wfn merge')
    try:
        os.mkdir('put_together')
    except:
        pass
    outfile=open(os.path.join('put_together','wfnmerge.inp'),'w')
    infile=open('pp.in','r')
    for line in infile:
        if line.find('wfng_nk1')!=-1:
            grid1=line.split()[2]
        if line.find('wfng_nk2')!=-1:
            grid2=line.split()[2]
        if line.find('wfng_nk3')!=-1:
            grid3=line.split()[2]
        if line.find('wfng_dk1')!=-1:
            shift1=line.split()[2]
        if line.find('wfng_dk2')!=-1:
            shift2=line.split()[2]
        if line.find('wfng_dk3')!=-1:
            shift3=line.split()[2]
        if line.find('wfng_file')!=-1:
            wfnfilename=line.split()[2].replace("'","").replace('"','')
    infile.close()
    print('wfn.combined', file=outfile)
    print(grid1,grid2,grid3, file=outfile)
    print(shift1,shift2,shift3, file=outfile)
    print(totalp, file=outfile)
    print(numkpts, file=outfile)
    for i in range(1,totalp+1):
        print('wfg.'+str(i), file=outfile)
        os.symlink('../p'+str(i)+'/'+wfnfilename,os.path.join('put_together','wfg.'+str(i)))
    outfile.close()


for i in range(1,totalp+1):
    directname='p'+str(i)
    try:
        os.mkdir(directname)
    except:
        pass
    savedir=os.path.join(directname,prefix+'.save')
    try:
        os.mkdir(savedir)
    except:
        pass
    modify_batch(i)
    # get new pw file
    set_up_pw_file(i)
    shutil.move('batch',directname)
    shutil.move(prefix+'.in',directname)
    upffiles=glob.glob('*UPF')
    for thefile in upffiles:
        shutil.copy(thefile,directname)
    shutil.copy('pp.in',directname)
    os.symlink('../../charge-density.dat',os.path.join(savedir,'charge-density.dat'))
    os.symlink('../../data-file.xml',os.path.join(savedir,'data-file.xml'))

setup_wfn_merge()

# delete unneeded kpt files
for filename in glob.glob('kpts.p*.dat'):
    os.remove(filename)
