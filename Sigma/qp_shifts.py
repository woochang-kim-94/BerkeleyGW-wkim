#!/usr/bin/env python

# B. Malone, 10/18/09
# Usage: qp_shifts.py sigma_hp.log

# Script takes a sigma_hp.log as input and produces two files, 'cond' and 
# 'val'. These files contain the quasiparticle shifts for each band. 'cond' 
# contains all conduction bands and 'val' valence bands.

# Output format:
# E_lda E_qp-E_lda
#       .
#       .
#       . 

# You can use the files to determine the QP scissor shift parameters 
# ecs, evs, ecdel, evdel, ev0, and ec0. ec0, ev0, ecs, and evs can simply be read 
# off from the files as the lowest conduction band energy, highest valence 
# band energy, the shift of the lowest conduction band, and the shift of the
# highest valence band. To determine ecdel and evdel (the slopes) you will
# need to fit a line to the points with a program like gnuplot.

# For example:
# gnuplot > f(x)=a*x+b;a=1;b=1;     # a,b initial guesses, doesn't matter 
# gnuplot > fit f(x) 'cond' via a,b
# gnuplot > [output giving you the fitted parameters] 
# gnuplot > plot f(x), 'cond'       # let you see the fitted line along
                                    # with the data points


import sys


try:
    infilename=sys.argv[1]
except:
    sys.exit('Error: The sigma_hp.log file should be the first argument')

print('Which band is the top of the valence band?')
valtop=int(input())

infile=open(infilename,'r')

condlist=[]
vallist=[]

flag=0
kpcounter=0


while True:
    line=infile.readline()
    if line=='': # EOF
        break
    if line.find('k = ')!=-1 and line.find('spin =')!=-1: # found kpoint
        kpcounter=kpcounter+1
        bandcounter=0
        for i in range(1,3):   # read in two lines
            infile.readline()
        while True:        # reading in bands
            line=infile.readline()
            if line.find('.')!=-1: # there are energy values (not blank)
                bandcounter=bandcounter+1
                theline=line.split()
                evalue=float(theline[1])
                GWvalue=float(theline[9])
                diff=GWvalue-evalue
                addition=(evalue,diff)
                if int(theline[0])<=valtop:
                    vallist.append(addition)
                else:
                    condlist.append(addition)
            else:          # done reading all bands for this kpoint
                break

            
infile.close()

print 'There are',bandcounter,'bands and', kpcounter, 'kpoints in the file'


vallist.sort()
condlist.sort()

# now write to outfiles

outfile1=open('cond','w')
outfile2=open('val','w')

for i in condlist:
    print >>outfile1, i[0],i[1]

for i in vallist:
    print >>outfile2, i[0],i[1]

outfile1.close()
outfile2.close()

print('Done, info written to \'cond\' and \'val\'')
