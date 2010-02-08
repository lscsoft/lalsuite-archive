#!/usr/bin/python

##############################################################################
#
#          options
#
##############################################################################

usage = """
usage: %prog [options] 

Estimate the sky position from a coincident trigger.

"""


def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser( usage=usage, version="%prog CVS $Id$ " )

  # options related to input and output
  parser.add_option("-g","--glob",action="store",type="string",\
      default=None, metavar=" GLOB",help="GLOB of thinca files to read" )
  parser.add_option("-c","--coarse-resolution",action="store",type="float",\
      default=4.0, metavar=" COARSE_RESOLUTION",help="number of sky points to throw in round one" )
  parser.add_option("-i","--fine-resolution",action="store",type="float",\
      default=0.5, metavar=" FINE_RESOLUTION",help="number of sky points to throw in round two" )
  parser.add_option("-d","--plotpoints",action="store_true",\
      default=False, help="make a color coded plot of the sky" )
  parser.add_option("-u","--usecatalog",action="store",type="string",\
      default=None, metavar=" CATALOG_NAME", help="galaxy catalog to use; must be specified if --listgalaxies option is used")
  parser.add_option("-D","--Deffcut",action="store_true",\
      default=True,help="only consider galaxies from here to the minimum effective distance measured. this is a very strong cut!") 
  parser.add_option("-f","--reference-frequency",action="store",type="float", default=0.0, metavar=" REFERENCE_FREQUENCY", \
    help="reference frequency for signal timing" )
  parser.add_option("--output-path", help="root of the HTML output")
  parser.add_option("-n","--narrow-threshold",action="store",type="float",\
      default=4, metavar=" L_NARROW_THRESH",help="threshold on L for 68% area" )
  parser.add_option("-w","--wide-threshold",action="store",type="float",\
      default=15, metavar=" L_WIDE_THRESH",help="threshold on L for 95% area" )
  parser.add_option("-o","--output",action="store",type="string",default=None,\
                    help="name of the xml file to store output in")
  parser.add_option("-z","--input_type",action="store",default="coinctable",\
                    help="specify the type of input in the glob.  valid options are coinctable (DEFAULT) and coire")
  parser.add_option("-y","--timing-only",action="store_true",default=False,\
                    help="only use timing information for sky localization")
  parser.add_option("-H", "--snr-threshold", action="store_true",default=True,\
                    help="use SNR dependent threshold")
  (options,args) = parser.parse_args()

  return options, sys.argv[1:]


##############################################################################
#
#          convenience functions
#
##############################################################################


##############################################################################
#
#          i/o setup
#
##############################################################################


##############################################################################
#
#          main program
#
##############################################################################



