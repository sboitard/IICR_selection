## simulates pairwise coalescence times (T2) with msms for several independent loci, each of them evolving under a single selective sweep model, and stores these times in a dataframe.

import subprocess
import os
import sys

# general parameters
sim_name=sys.argv[1] # file ID
nb_rep=int(sys.argv[2]) # number of simulated loci 
L=int(sys.argv[3]) # length of each locus (in bp)
theta=float(sys.argv[4]) # per site scaled mutation rate, 4Nmu
rho=float(sys.argv[5]) # per site scaled recombination rate, 4Nr
alpha=int(sys.argv[6]) # scaled selection intensity, 2Ns (fitness 1+s for homozygote mutants)
t0=float(sys.argv[7]) # end of the sweep (fixation time), time before present in 2N units

# call msms
cmd='msms/bin/msms 2 '+str(nb_rep)+' -t '+str(L*theta)+' -r '+str(L*rho)+' '+str(L)+' -T -Sp 0.5 -SAA '+str(alpha)+' -SaA '+str(alpha/2)+' -SF '+str(t0/2)+' -N 10000'
print cmd
p = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
(output, err) = p.communicate()
if not err == None:
    print 'error while running ms'
else:
    # parse msms output and store coalescence times
    lines=output.splitlines()
    coal=[]
    s=0
    for line in lines:
	if len(line)>0:
	    if line=='//':
	    	s+=1
	    elif line[0]=='[':
		buf=line.split(']')
		nb_sites=buf[0][1:]
		buf2=buf[1].split(',')
		t=buf2[0][3:]
		coal.append([s,nb_sites,t])
		# s: locus number (1<=s<=nb_rep)
		# nb_sites: number of sites of the current non-recombining segment
		# t: coalescence time for the current non-recombining segment
outfile=open(sim_name+'.times','w')
outfile.write('sample length time\n')
for t in coal:
    outfile.write(str(t[0])+' '+t[1]+' '+t[2]+'\n')

