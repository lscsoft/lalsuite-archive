#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       AppendToDatabase.py
#
#       Copyright 2012
#      Salvatore Vitale <salvatore.vitale@ligo.org>
#
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

#===============================================================================
# Preamble
#===============================================================================

#standard library imports
import sys
import os
import getpass

import numpy as np

#local application/library specific imports
#from pylal import SimInspiralUtils
#from pylal import bayespputils as bppu
from pylal import git_version

__author__="Salvatore  Vitale <salvatore.vitale@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date


def vararg_callback(option, opt_str, value, parser):
    assert value is None
    value = []
    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False
    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(arg)
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)        

def my_sort(a,b):
    if a=='GR':
        return -1
    elif b=='GR':
        return 1
    elif len(a)<len(b):
        return -1
    elif len(a)>len(b):
        return 1
    elif str(a)<str(b):
        return -1
    elif str(a)>str(b):
        return 1
    else:
        print "Could not order %s and %s"%(a,b)
    

def AppendToDatabase(   path,
                        subhyp,
                        remote_script,
                        remote_database,
                        time,
                        testP,
                        seed,
                        local_pc
                        #misc. optional
                        #injfile=None,eventnum=None,skyres=None,
                       
                    ):
    
    """
    """
    #clusters={"pcdev1.phys.uwm.edu":"Nemo","marlin.phys.uwm.edu":"Marlin","hydra.phys.uwm.edu":"Hydra","atlas1.atlas.aei.uni-hannover.de":"Atlas1","atlas2.atlas.aei.uni-hannover.de":"Atlas2","atlas3.atlas.aei.uni-hannover.de":"Atlas3","atlas4.atlas.aei.uni-hannover.de":"Atlas4","titan1.atlas.aei.uni-hannover.de":"Titan1","titan2.atlas.aei.uni-hannover.de":"Titan2","titan3.atlas.aei.uni-hannover.de":"Titan3","ldas-grid.ligo.caltech.edu":"Cit","sugar.phy.sys.edu":"Sugar","ldas-grid.ligo-la.caltech.edu":"LLO","ldas-grid.ligo-wa.caltech.edu":"LHO1","ldas-pcdev1.ligo-wa.caltech.edu":"LHO2"}
    #for un in os.uname():
    #    if un in clusters.keys():
    #        local_pc=clusters[un]
    #        break
    #    else:
    #        local_pc=os.uname()[1]
            
    header=" "
    string_to_write=""
    string_to_write+=str(seed)+" "
    header+="InspinjSeed "
    header+="User "
    string_to_write+=str(getpass.getuser())+" "
    string_to_write+=str(local_pc)+" "
    header+=str("cluster  time" )+ " "
    string_to_write+=str(time)+" "
    found=0
    failed={}

    for hyp in sorted(subhyp,my_sort):
        path_to_snr=os.path.join(path,hyp,'SNR',"snr_H1L1V1_"+str(time)+".0.dat")
        path_to_Bfile=os.path.join(path,hyp,'nest',"outfile_"+str(time)+".000000_H1L1V1.dat_B.txt")
        if os.path.isfile(path_to_Bfile):
            bfile=np.loadtxt(path_to_Bfile)
            logB=bfile[0]
            string_to_write+="%.3f"%logB+" "
            header+=str(hyp)+" "
            found+=1
        else:
            print "WARNING: The Bfile %s was not found"%str(path_to_Bfile)
            failed[hyp]=str(path_to_Bfile)
            string_to_write+=str(path_to_Bfile)+" "
               
    if found < len(subhyp):
        print "ERROR, some of the B files were not found. This is probably due to some failed inspnest run(s)"
        remote_database=remote_database+"_failed"
 
    if os.path.isfile(path_to_snr):
        snrfile=np.loadtxt(path_to_snr,skiprows=3,usecols=(1,1))
        snr=snrfile[0]
        string_to_write+="%.2f"%snr+" "
        header+=str("NetSNR ")

    string_to_write+=str(testP)+" "
    header+=str("testParameter_ShiftPc ")
    remote_server=remote_script[0:remote_script.index(":")]
    path_remote_script=remote_script[remote_script.index(":")+1:]

    #print "ssh "+remote_server + ' "'+ path_remote_script + " '" +string_to_write +"' "+ remote_database+'"'
    #os.system("ssh "+remote_server + ' "'+ path_remote_script + " '" +string_to_write +"' "+ remote_database+" '" +header +"' "+'"')    
    os.system("ssh "+remote_server + ' "'+ path_remote_script + " '" +string_to_write +"' "+ remote_database+'"')  

if __name__=='__main__':

    from optparse import OptionParser
    parser=OptionParser()
    parser.add_option("-p","--path", dest="path",type="string",help="This dir contain all the sub-hypotheses dirs", metavar="DIR")
    parser.add_option("-s","--subhyp",dest="subhyp",action="callback",callback=vararg_callback,help="A space-separed list of the subhypothesis (comprised the GR)",metavar="GR dphi1 dphi2 dphi1dphi2")
    parser.add_option("-R","--remote-script",type="string",action="store",help="The path to the remote script which writes on the datafile", metavar="svitale@login.nikhef.nl:/project/gravwav/safe_append.sh")
    parser.add_option("-r","--remote-database",type="string",action="store",help="The remote database on which to append", metavar="datafile.txt")
    parser.add_option("-T","--testParam-ShiftPc",type="string",action="store",help="The shifted parameter underscore the percent value of the shift, or GR", metavar="dphi6_1pc")
    parser.add_option("-Q","--inspinj-seed",default=None,action="store",type="string",help="The unique value of the inspinj seed", metavar="700000")
    parser.add_option("-t","--event-time",type="string",action="store",help="The time of injection", metavar="939936910")
    parser.add_option("-c","--cluster",type="string",action="store",help="The cluster from which the pipeline is being run", metavar="Atlas1")

    (opts,args)=parser.parse_args()


    AppendToDatabase( 
                        opts.path,
                        opts.subhyp,
                        opts.remote_script,
                        opts.remote_database,
                        opts.event_time,
                        opts.testParam_ShiftPc,
                        opts.inspinj_seed,
                        opts.cluster
                    )
#
