# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:46:50 2013

@author: jeroen
"""

###############################################################################
### This script will check the progress of the job and write
### all finished events to a properly formatted data file.
### The output datafile can then be read in by a postprocessing script
###############################################################################


import os
import subprocess


###############################################################################
### BEGIN USER INPUT
###############################################################################

bindir = "/home/jmeidam/opt_tiger_ringdown/bin"
rundir = "/home/jmeidam/tiger_runs/TigerRingdownET/RingdownTDinj/GR"

injection = "GR"

#for now hardcoded:
dt = 1000.0

Nhypotheses = 8

lowcutSNR = 8.0
highcutSNR = 30.0

#Do you only want sources in the outfile that make the above SNR cut?
DoSNRCut = True


###############################################################################
### END USER INPUT
###############################################################################




XML="injections_%s_%s.xml" #formatted with %(injection,seed)

seeds = []
for i in os.listdir(rundir):
    #print i
    #print os.path.isdir(i)
    #if os.path.isdir(i):
    seeds.append(i)
    print "Found seed",i

Nevents=0

binfiles=os.listdir(bindir)
if 'ligolw_print' in binfiles:
    xmlprint = os.path.join(bindir,"ligolw_print")
else:
    xmlprint = os.path.join(bindir,"lwt_print")

print "using "+xmlprint

#function to find all occurrences of sub in a_str
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)

def NumsFromHyps(hypotheses):
    out = {}
    numlist = []
    nameint=0
    for hyp in hypotheses:
        dfreqlocs = list(find_all(hyp,'dfreq'))
        dtaulocs = list(find_all(hyp,'dtau'))
        print hyp

        name = ''
        if(dfreqlocs):
            for i in dfreqlocs:
                l = hyp[i+5]
                m = hyp[i+5+1]
                name+=(l+m)
        if(dtaulocs):
            for i in dtaulocs:
                l = hyp[i+4]
                m = hyp[i+4+1]
                name+=('9'+l+m)

        try:
            nameint = int(name)
            numlist.append(nameint)
            out[name]=hyp
        except ValueError:
            out['0']='GR'

    numlist.sort()
    numlist.append(0) #add the GR hypothesis last

    return out,numlist

def readxml(xml):
    tmpfile = "processparameters.txt"
    os.system(xmlprint+" "+xml+" -t process_params >> "+tmpfile)
    with open(tmpfile,'r') as f:
        processparameters = {}
        for line in f:
            row = line.split(',')
            key = row[2].strip() #strip() to remove any whitespace
            key = key.replace("\"","") #remove quotes
            value = row[4].strip() #strip() to remove any whitespace and '\n'
            value = value.replace("\"","") #remove quotes
            processparameters[key] = value

    os.system("rm "+tmpfile)
    return processparameters

def get_snr(filepath):
    with open(filepath,'r') as f:
        content = f.readlines()
        E1 = content[0].split()[1]
        E2 = content[1].split()[1]
        E3 = content[2].split()[1]
        NETWORK = content[3].split()[1]
    return ( E1, E2, E3, NETWORK )

out = []

for seed in seeds:
    print "Working on seed",seed
    seeddir = os.path.join(rundir,seed)
    xmlfile = os.path.join(seeddir,XML%(injection,seed))

    procparams = readxml(xmlfile)

    gpsend = float(procparams['--gps-end-time'])
    gpsstart = float(procparams['--gps-start-time'])

    Nevents = int( (gpsend - gpsstart)/dt )

    #get hypotheses from current directory
    hypotheses = [ d for d in os.listdir(seeddir) if os.path.isdir(os.path.join(seeddir,d)) ]
    hypdict,hypnums = NumsFromHyps(hypotheses)
    #print hypnums
    #print hypdict[str(hypnums[2])]

    #dictionary containing a list of all finished hypotheses for each event
    hyp_for_event = {}

    #essentially what will be written to the outfile
    out_thisseed = []

    #dictionary to check which hypothesis is complete (i.e. 1000 events)
    hypcount = {}
    for hyp in hypnums:
        hypcount[str(hyp)] = 0

    for event in range(Nevents):
        #create one set of hypotheses and their bayes factors for each event.

        gpstime = "%d"%int(gpsstart+event*dt)

        #output for current event
        out_thisevent = []

        #keep track of progress
        print event,gpstime

        #list of finished hypotheses for this event
        presenthypotheses = []

        for hyp in hypnums:
            #read all the available bayes files for this hypothesis

            hypname = hypdict[str(hyp)]
            bayesdir = os.path.join(seeddir,hypname,"rundir","engine")

            Bayesfile = "lalinferencenest_"+hypname+"_event"+str(event)+"-E1E2E3-"+gpstime+".0.0-"+str(event)+".dat_B.txt"

            if os.path.isfile(os.path.join(bayesdir,Bayesfile)):
                with open(os.path.join(bayesdir,Bayesfile), "r") as f:
                    line = f.readline()
                    LogBayes = line.split()[0]

                    #what are these again?
                    no1 = line.split()[1]
                    no2 = line.split()[2]
                    no3 = line.split()[3]

                snrfile = "snr_E1E2E3_"+gpstime+".0.dat"
                SNR = get_snr( os.path.join(seeddir,hypname,"rundir","SNR",snrfile) )

                string = "%s %s %s %d %s %s %s %s\n"%(seed,gpstime,SNR[3],hyp,LogBayes,no1,no2,no3)

                out_thisevent.append(string)
                presenthypotheses.append(hyp)

                count = hypcount[str(hyp)]
                count += 1
                hypcount[str(hyp)] = count

        hyp_for_event[str(event)] = presenthypotheses
        out_thisseed.append(out_thisevent)

    print "For seed "+seed+", "+str(len(out_thisseed))+" events were read in"

    #let us know how far along the runs are
    for hyp in hypcount:
        count = hypcount[hyp]
        percent = (float(count)/float(Nevents))*100.0
        if (count == Nevents):
            print "%s hypothesis %s is finished"%(seed,hyp)
        else: print "%s hypothesis %s is %.2f %% complete (%d/%d)"%(seed,hyp,percent,count,Nevents)

    #let us know which events are finished
    finished = 0
    for event in range(Nevents):
        line = hyp_for_event[str(event)]
        if (len(line)==len(hypnums)):
            finished+=1
    print "Currently %d events are completely finished"%finished


#make a list of events that are finished
finished = []
for event in range(Nevents):
    count = 0
    for i in hyp_for_event[str(event)]:
        count+=1
    if (count==Nhypotheses):
        finished.append(event)





#write out the dat file to read with the postprocessing script.
with open("./"+injection+".dat","w") as f:
    count = 0
    for event in finished:#range(Nevents):
        current_event = out_thisseed[event]
        SNR = float( ( current_event[0].split() )[2] )
        if (DoSNRCut):
            if (SNR >= lowcutSNR and SNR <= highcutSNR):
                count += 1
                for line in current_event:
                    f.write(line)
        else:
            for line in current_event:
                f.write(line)

        #end of event:
        f.write("\n")

if (DoSNRCut):
    print "%d sources made the SNR cut of [%.1f,%.1f]"%(count,lowcutSNR,highcutSNR)








