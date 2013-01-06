import pylal.bayespputils as bppu
from optparse import OptionParser
from math import pi

parser=OptionParser()
parser.add_option("-i", "--input",dest="i")
parser.add_option("-o", "--output",dest="o")
parser.add_option("-j","--injFile",dest="j")
parser.add_option("-e","--event",type="int",dest="e")
parser.add_option("-u","--url",action="store_true",dest="u",default=False)
parser.add_option("-n","--name",dest="n")
parser.add_option("-p","--password",dest="p")

(options, args) = parser.parse_args()

if (options.i == None or options.o == None or options.j == None or options.e == None):
    parser.error("not enough arguments")
input=options.i
output=options.o
injFile=options.j
event=options.e
urlUsed=options.u
userName = options.n
userPassword = options.p

###############                                                                                                                            
def open_url_wget(url,folder,un=userName,pw=userPassword,eventNum=0, args=[]):
    import subprocess
    import urlparse
    name = folder+"/posteriors/posterior_"+str(eventNum)+".dat"
    if un is not None and pw is not None:
        args+=["-O",name,"--user",un,"--password",pw,"--no-check-certificate"]
    retcode=subprocess.call(['wget']+[url]+args)

    return retcode
###############load posterior#########                                                                                                     

data=bppu.PEOutputParser('common')
if(urlUsed):
    open_url_wget(input,output,eventNum=event)
    inputFileObj = open(output+'/posteriors/posterior_'+str(event)+'.dat')
else:
    inputFileObj = open(input)
dataObj0 = data.parse(inputFileObj)

                                                                                                                          
from pylal import SimInspiralUtils
injections = SimInspiralUtils.ReadSimInspiralFromFiles([injFile])
if(len(injections)<event):
    print "Error: You asked for event %d, but %s contains only %d injections" %(event,injfile,len(injections))
    sys.exit(1)
else:
    injection = injections[event]

posterior = bppu.Posterior(dataObj0,injection)
outFile = open(output+'/kdresult' + str(event), 'w')
outFile.write('label injection_cl injection_area\n')
confidenceLevels = [0.68,0.9]

def mapping(ra,dec):
    return (ra,dec)

#set up evrything for running kd algorithm
if 'ra' and 'dec' not in posterior.names:
    if 'rightascension' and 'declination' in posterior.names:
        new_param_names = ['ra','dec']
        post_names = ['rightascension','declination']
        posterior.append_multiD_mapping(new_param_names, mapping , post_names)
    else:
        print 'param not found'
print 'out test --------'
areas, nodeList, injInfo = bppu.kdtree_bin_sky_area(posterior,confidenceLevels,samples_per_bin = 50)
print 'details'
print injInfo
print areas

for cl in confidenceLevels:
    temp_area = areas[cl]*(180/pi)**2.
    outFile.write('kd_areaCL' +str(cl) + ' ' + str(temp_area) + '\n')

outFile.write('kd_injCL ' + str(injInfo[3])+' \n')                                                                                                      
temp_area = injInfo[4]*(180/pi)**2.
outFile.write('kd_injCLarea ' + str(temp_area) + '\n')
