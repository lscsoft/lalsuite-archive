# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 10:25:53 2013

@author: jeroen
"""

import time
import os

from optparse import OptionParser

ligolw_print = 1

#######################################################################################
## Option Parser
#######################################################################################

usage=""" %prog [options]
Setup an html web page containing links to all posplots.html pages for all injection within one project folder.
"""

parser=OptionParser(usage)
parser.add_option("-p","--project-htmlpath",dest="projecthtmlpath",default=None, action="store", type="string",help="Public_html directory containing injections",metavar="PROJHTMLPATH")
parser.add_option("-r","--project-rundir",dest="projectrundir",default=None,action="store",type="string",help="Run directory containing injections",metavar="PROJRUNDIR")
parser.add_option("-o","--output-html",dest="outhtml",default=None,action="store",type="string",help="Location of the result page",metavar="OUTHTML")
parser.add_option("-u","--project-url",dest="projecturl",default=None,action="store",type="string",help="URL to project",metavar="PROJURL")
parser.add_option("-s","--gps-start-time",dest="starttime",action="store",type="string",default=None,help="Start time of analysis")
parser.add_option("-t","--time-step",dest="timestep",action="store",type="string",default=None,help="time step")
parser.add_option("-n","--no-postprocess",action="store_false", dest="nopost", default=True, help="postprocessing is not included in the dag")

#######################################################################################
## Create subfiles and dagfiles for each hypothesis
#######################################################################################

(opts,args)=parser.parse_args()
projectpath = opts.projecthtmlpath
outputfile = opts.outhtml
projectlink = opts.projecturl
rundir = opts.projectrundir
nopost = opts.nopost


#to be able to include SNR information I need the timestep and gps_start
#If these are incorrect, it will simply not find the snr files and omit them
gps_start = opts.starttime
dt = opts.timestep

injections = []
hypotheses = []

list_of_htmlfiles = {}
#html_files = {}

def get_snr(filepath):
    with open(filepath,'r') as f:
        content = f.readlines()
        E1 = content[0].split()[1]
        E2 = content[1].split()[1]
        E3 = content[2].split()[1]
        NETWORK = content[3].split()[1]
    return ( E1, E2, E3, NETWORK )

def dirname(directory):
    return str( os.path.basename(os.path.normpath(directory)) )

###############################################################################
## Write header of html file
###############################################################################
with open(outputfile,"w") as f:
    head=(
"<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n"
"<?xml version=\"1.0\" ?>\n"
"<html xmlns=\"http://www.w3.org/1999/xhtml\">\n"
"  <head>\n"
"    <title>\n"
"      "+dirname(projectpath)+"\n"
"    </title>\n"
"    <style type=\"text/css\">\n"
"\n"
"    p,h1,h2,h3,h4,h5\n"
"    { font-family:&quot;Trebuchet MS&quot;, Arial, Helvetica, sans-serif; }\n"
"\n"
"    p\n"
"    {font-size:14px;}\n"
"\n"
"    h1\n"
"    {font-size:20px;}\n"
"\n"
"    h2\n"
"    {font-size:18px;}\n"
"\n"
"    h3\n"
"    {font-size:16px;}\n"
"\n"
"    table\n"
"    {\n"
"    font-family:&quot;Trebuchet MS&quot;, Arial, Helvetica, sans-serif;\n"
#"    width:100%;\n"
"    border-collapse:collapse;\n"
"    }\n"
"    td,th\n"
"    {\n"
"    font-size:12px;\n"
"    border:1px solid #B5C1CF;\n"
"    padding:3px 7px 2px 7px;\n"
"    }\n"
"    th\n"
"    {\n"
"    font-size:14px;\n"
"    text-align:left;\n"
"    padding-top:5px;\n"
"    padding-bottom:4px;\n"
"    background-color:#B3CEEF;\n"
"    color:#ffffff;\n"
"    }\n"
"    #postable tr:hover\n"
"    {background: #DFF4FF;}\n"
"    #covtable tr:hover\n"
"    {background: #DFF4FF;}\n"
"    #statstable tr:hover\n"
"    {background: #DFF4FF;}\n"
"    \n"
"    .ppsection\n"
"    {border-bottom-style:double;}\n"
"    \n"
"    </style>\n"
"  </head>\n" )


    f.write(head)


def htmlappend(outname,text):
   with open(outname,'a') as f:
    f.write(text)

#with open(outputfile,'a') as f:
#    f.write(
#    "  <body> \n"
#    "    <h1> \n"
#    "    "+dirname(projectpath)+" \n"
#    "    </h1>\n"
#    )

def readlastline(fname):
    with open(fname, 'rb') as fh:
        first = next(fh)
        offs = -100
        while True:
            fh.seek(offs, 2)
            lines = fh.readlines()
            if len(lines)>1:
                last = lines[-1]
                break
            offs *= 2
    return last

def checkDAG(dagman_dot_out):

    if os.path.isfile(dagman_dot_out):
        daginfo = readlastline(dagman_dot_out)
        print daginfo
        return 1
        daginfo.find("****")

    else:
        return -2



def readxml(xml):
    tmpfile = "processparameters.txt"
    if (ligolw_print):
        os.system("ligolw_print "+xml+" -t process_params >> "+tmpfile)
    else:
        os.system("lwtprint "+xml+" -t process_params >> "+tmpfile)
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

def checklumdist(params):

    if( params.has_key('--min-distance') ):
        return 1
    else:
        return 0


def create_injectionparams_table(params,f,indent=""):

    LumDist = checklumdist(params)

    if (LumDist):
        dmin = "%.2f"%(float(params['--min-distance']))
        dmax = "%.2f"%(float(params['--max-distance']))
    else:
        dmin = "%.2f"%(float(params['--min-z']))
        dmax = "%.2f"%(float(params['--max-z']))

    qmin = "%.2f"%(float(params['--min-mratio']))
    qmax = "%.2f"%(float(params['--max-mratio']))
    f.write(indent+"<table id=\"injtable\" border=\"0\">")
    f.write(indent+"  <tr> <th colspan=\"2\" align=\"center\"> Injected parameters </th> </tr>")
    if (LumDist):
        f.write(indent+"  <tr> <td>Dmin</td>  <td>"+dmin+"</td> Mpc</tr>")
        f.write(indent+"  <tr> <td>Dmax</td>  <td>"+dmax+"</td> Mpc</tr>")
    else:
        f.write(indent+"  <tr> <td>zmin</td>  <td>"+dmin+"</td> </tr>")
        f.write(indent+"  <tr> <td>zmax</td>  <td>"+dmax+"</td> </tr>")
    f.write(indent+"  <tr> <td>qmin</td>  <td>"+qmin+"</td> </tr>")
    f.write(indent+"  <tr> <td>qmax</td>  <td>"+qmax+"</td> </tr>")
    if (params.has_key('--disable-spin')):
        f.write(indent+"  <tr> <td colspan=\"2\" align=\"center\">spin disabled</td> </tr>")
    else:
        #I want to include the spin settings here
        f.write(indent+"  <tr> <td colspan=\"2\" align=\"center\">spin enabled</td> </tr>")
    f.write(indent+"</table>")
    f.write(indent+"<p></p>")



text = ("  <body> \n"
        "    <h1> \n"
        "    "+dirname(projectpath)+" \n"
        "    </h1>\n" )

htmlappend(outputfile, text)

##This function adds a section to the page complete with table of htmls
def runproject(projectpath,outputfile):
    for injection in os.listdir(projectpath):
        print injection
        for seed in os.listdir(os.path.join(projectpath,injection)):
            print "  "+seed
            percent = (injection.split("_")[1]).replace("pc","",1)
            nonGRparam = injection.split("_")[0]

            jobstatus = checkDAG( os.path.join(rundir,injection,seed,"common_dag.dag.dagman.out") )
            if(jobstatus == 1):
                jobstring = "<span style=\"color: #008000\">Finished</span>"
            elif(jobstatus == -1):
                jobstring = "<span style=\"color: #800000\">Aborted</span>"
            elif(jobstatus == -2):
                jobstring = "<span style=\"color: #800080\">No dagman.out found</span>"
            else:
                jobstring = "<span style=\"color: #008080\">Ongoing</span>"
            xml = os.path.join(rundir,injection,seed)+"/injections_"+injection+"_"+seed+".xml"
            if (os.path.isfile(xml)):
                procparams = readxml(xml)
            else:
                print "Found no xml file for "+injection+"/"+seed+":"
                print xml
                break
            with open(outputfile,'a') as f:
                f.write("    <div class=\"ppsection\"> \n"
                    "      <h2> \n"
                    "        Injection ("+jobstring+")\n"
                    "      </h2> \n"
                    "      <p> shift in "+nonGRparam+" of "+percent+"% </p>\n"
                    "       \n")
                create_injectionparams_table(procparams,f,indent="      ")
                f.write("    </div> \n" )

            hypotheses = []
            tabledict = {}
            Highestevent = 0
            for hypothesis in os.listdir(os.path.join(projectpath,injection,seed)):
                print "    "+hypothesis
                hypotheses.append(hypothesis)
                events=[]
                for event in os.listdir(os.path.join(projectpath,injection,seed,hypothesis)):
                    print "      "+event
                    thisevent = event.split("_")[1]
                    n = int(thisevent.replace("event","",1))
                    events.append(n)
                    #Table should be as long as highest event
                    if (n>Highestevent):
                        Highestevent = n

                tabledict[hypothesis] = events

            #write the table
            #first the header
            f = open(outputfile,'a')
            f.write("    <div class=\"ppsection\">\n"
                    "    <h2>\n"
                    "      Table of web pages for seed "+seed+"\n"
                    "    </h2>\n")
            f.write("      <table border=\"1\" id=\"statstable\" width=\"100%\">\n")
            f.write("        <tr> <th>event #</th>")
            for hyp in hypotheses:
                f.write("<th>"+hyp+"</th>")
            f.write("</tr>\n")

            #then the rows
            for event in range(Highestevent+1):
                f.write("        <tr>")
                f.write("<td>"+str(event)+"</td>")
                gpstime = "%.1f"%float(float(procparams['--gps-start-time'])+event*float(dt))
                snrfilename = "snr_E1E2E3_%s.dat"%gpstime
                for hyp in hypotheses:
                    Bfile = "lalinferencenest_"+hyp+"_event"+str(event)+"-E1E2E3-"+gpstime+".0-"+str(event)+".dat_B.txt"
                    snrfile = os.path.join(rundir,injection,seed,hyp,"rundir","SNR",snrfilename)
                    #print snrfile
                    eventname = hyp+"_event"+str(event)
                    if (os.path.isfile(snrfile)):
                        SNR = get_snr(snrfile)[3]
                        #print 'SNR=',SNR
                        linkname = hyp+"_"+str(event)+" | "+SNR
                    else:
                        linkname = eventname
                    html = os.path.join(injection,seed,hyp,eventname,"posplots.html")
                    address = projectlink+html
                    link = "<a href=\""+address+"\"> "+linkname+" </a>"

                    if not ( os.path.isfile(os.path.join(projectpath,html)) ):
                        if (os.path.isfile(os.path.join(rundir,injection,seed,hyp,"rundir","engine",Bfile))):
                            f.write("<td> <span style=\"color: #008000\">Finished</span> </td> ")
                        else:
                            f.write("<td>"+"not finished"+"</td> ")
                    else:
                        f.write("<td>"+link+"</td> ")

                f.write("</tr>\n")
            f.write("      </table>\n")
            f.write("      <p></p>\n")
            f.write("    </div>\n")

            f.close()

runproject(projectpath, outputfile)

from datetime import datetime, timedelta

with open(outputfile,"a") as f:
    f.write("    <div>\n")
    current_time_utc = datetime.utcnow()
    time_cet = current_time_utc + timedelta(hours=2)
    datestring = str(time_cet.strftime("%d/%m/%Y"))
    timestring = str(time_cet.strftime("%H:%M:%S"))
    f.write("      <p>Document created on %s at %s CET</p>\n"%(datestring,timestring))
    f.write("    </div>\n")
    f.write("  </body>\n")
    f.write("</html>")

