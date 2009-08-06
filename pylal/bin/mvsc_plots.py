#!/usr/bin/python
"""Contains functions to create various plots from lalapps_mvsc_player"""

__author__ = 'Tristan Miller <tmiller@caltech.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'
__prog__ = 'mvsc_plots'

##############################################################################
# import modules
import math,re
from matplotlib import pyplot
from matplotlib import cbook

##############################################################################

def patread(filename,station=None,flag=False,headeronly=False):
    """Reads in a file from SprOutputWriterApp.

    If flag is false, reads a dat file.  Otherwise, reads a pat file."""
    
    try:
        f = open(filename)
    except IOError:
        print '***Error!*** Trouble opening file', filename
        return

    p = re.compile(r'\S+')

    #the flag is part of a workaround dealing with rounding issues in
    #SprOutputWriterApp
    # if flag is true, reads a pat file (rather than a dat file)
    # and outputs strings, without rounding
    if flag:
        h = f.readline()
        
    h = f.readline()
    header = p.findall(h)
    
    if headeronly:
        f.close()
        return header
        
    data = []
    cols = {}
    for i in range(len(header)):
        cols[header[i]] = i

    #read the rest of the file
    n = f.readline()
    m = p.findall(n)
    if m:
        for i in range(len(m)):
            data.append([float(m[i])])

    if flag:
        while m:
            n = f.readline()
            m = p.findall(n)
            if m:
                for i in range(len(m)):
                    data[i].append(m[i])
    else:
        while m:
            n = f.readline()
            m = p.findall(n)
            if m:
                for i in range(len(m)):
                    data[i].append(float(m[i]))
    
    f.close()

    if not station:
        print 'Finished reading file', filename
        return data, cols
    
    #modify cols dictionary to be more uniform (ie get rid of H1s and L1s)
    iter = cols.iterkeys()
    cols2 = {}
    
    while True:
        try:
            key = iter.next()
            origcopy = key
            
            for j in range(0,6,2):
                for i in range(0,len(station),2):
                    if key[j:j+2] == station[i:i+2]:
                        key = key[:j]+'t'+str(i/2+1)+key[j+2:]
                        break

            cols2[key] = cols[origcopy]

        except StopIteration:
            break
        
    print 'Finished reading file', filename
    return data, cols2, cols

##############################################################################

def inforead(filename):
    """Reads the printout of SprBaggerDecisionTreeApp and returns info.

    Info returned: List of FOM values
                   Table of number of tree splits
                   Number of timeslides in training set
                   Number of injections in training set
                   Number of timeslides in validation set
                   Number of injections in validation set"""

    try:
        f = open(filename)
    except IOError:
        print '***Error!*** Trouble opening file', filename
        return

    output = []    
    p1 = re.compile(r'Points in class [0(1):]{5}   (\S+)')
    counter = 0
    
    while True:
        n = f.readline()
        m = p1.match(n)
        if m:
            output.append(int(m.group(1)))
            counter += 1
        elif counter >= 4:
            break
        elif not n:
            print '***Error!*** Missing info in ', filename
            return
    
    p2 = re.compile(r'Validation FOM=(\S+)')
    fom = []
    flag = False
    while True:
        n = f.readline()
        m = p2.match(n)
        if m:
            fom.append(float(m.group(1)))
            flag = True
        elif flag:
            break
        elif not n:
            print '***Error!*** No FOM in ', filename
            return

    treesplits=[]
    p3 = re.compile(r'Variable\s+(\S+)\s+Splits\s+(\S+)\s+Delta FOM\s+(\S+)')
    while True:
        n = f.readline()
        m = p3.match(n)
        if m:
            treesplits.append([])
            treesplits[-1].append(m.group(1))
            treesplits[-1].append(int(m.group(2)))
            treesplits[-1].append(float(m.group(3)))
        elif not n:
            break
    
    f.close()
    return fom,treesplits,output[0],output[1],output[2],output[3]

##############################################################################

def ROC(timeslides,injections):
    """Computes true and false positive rates.

    May change the order of timeslides and injections."""

    timeslides.sort()
    injections.sort()
    n = len(injections)
    m = len(timeslides)

    truepos = []
    falsepos = []

    for i in range(n):
        if injections[i] != injections[i-1]:
            for y in range(n):
                if injections[y]>injections[i]:
                    truepos.append(1-(y-1)/float(n))
                    break
            else:
                truepos.append(0.0)

            #Note that the minimum false positive rate is set to 1/(m+1)
            #which is the best upper bound of the wilson confidence interval
            
            if timeslides[-1]<=injections[i]:
                falsepos.append(1/float(m+1))
            else:
                for x in range(2,m+1):
                    if timeslides[-x]<=injections[i]:
                        falsepos.append((x-1)/float(m))
                        break
                else:
                    falsepos.append(1.0)

    return truepos, falsepos

##############################################################################

def wilson(p,n):
    """Calculates the Wilson interval, the confidence interval for a binomial
    distribution.

    Returns the appropriate upper and lower error bars.
    The confidence level used is always 68%."""

    n = float(n)
    diff = math.sqrt(max(p*(1-p)/n + 0.25/n**2,0)) / (1+1/n)
    mid = (-p/n + 0.5/n) / (1+1/n)
    upper = mid+diff
    lower = diff-mid

    return upper,lower

##############################################################################

def stairs(x,y):
    """Transforms lists (x,y) into (x1,y1) such that plotting them will
    create a stairs graph."""

    x1 = []
    y1 = [y[0]]
    for i in range(len(x)-1):
        x1.append(x[i])
        x1.append(x[i])
        y1.append(y[i+1])
        y1.append(y[i+1])
    x1.append(x[-1])

    return x1,y1

##############################################################################

def top_events(data,cols,n,stations,patfile):
    """Finds the n events with the highest MVSC values.

    If multiple events tie for last, all will be included."""

    #must transpose data in order to sort it
    data2 = []
    for i in range(len(data[0])):
        data2.append([])
        for j in range(len(data)):
            data2[-1].append( data[j][i] )
        
    sorter = cbook.Sorter()
    sorter(data2,cols['Bagger'])
    mvsc_cutoff = data2[-min(len(data2),n)][cols['Bagger']]

    events = []
    index = []
    
    for i in range(1,len(data2)+1):
        if data2[-i][cols['Bagger']] >= mvsc_cutoff:
            index.append(i)
            #events.append([])
            #for j in range(len(data2[0])):
            #    events[-1].append(data2[-i][j])
        else:
            break

    #problem: SprOutputWriterApp rounds all numbers to six digits
    #workaround solution: below
    data3,cols,cols2 = patread(patfile,stations,True)

    for i in index:
        events.append([])
        
        for j in range(3):
            events[-1].append(data2[-i][j])
        for j in range(len(data3)-1):
            events[-1].append(data3[j][-i])
        events[-1].append(data2[-i][-1])

    return events

##############################################################################

def ROCplot(data,cols,op_point = 1,ts_trig_ratio = 25):
    """Creates an ROC plot from one file.

    Returns the confidence interval of the resulting efficiency."""
    
    timeslides = []
    injections = []
    timeslides_snr = []
    injections_snr = []

    for i in range(len(data[0])):
        if data[1][i] == 0:
            timeslides.append(data[cols['Bagger']][i])
            timeslides_snr.append(data[cols['t1get_effective_snr()']][i]**2+ \
                                  data[cols['t2get_effective_snr()']][i]**2)
        else:
            injections.append(data[cols['Bagger']][i])
            injections_snr.append(data[cols['t1get_effective_snr()']][i]**2+ \
                                  data[cols['t2get_effective_snr()']][i]**2)

    truepos,falsepos = ROC(timeslides,injections)
    truepos_snr,falsepos_snr = ROC(timeslides_snr,injections_snr)

    binomerru = []
    binomerrl = []
    for i in range(len(truepos)):
        upper,lower = wilson( truepos[i], len(injections) )
        binomerru.append(upper)
        binomerrl.append(lower)
        binomerru.append(0)
        binomerrl.append(0)
    t = binomerru.pop()
    t = binomerrl.pop()
    
    binomerru_snr = []
    binomerrl_snr = []
    for i in range(len(truepos_snr)):
        upper,lower = wilson( truepos_snr[i], len(injections) )
        binomerru_snr.append(upper)
        binomerrl_snr.append(lower)
        binomerru_snr.append(0)
        binomerrl_snr.append(0)
    t = binomerru_snr.pop()
    t = binomerrl_snr.pop()

    falsepos2 = []
    falsepos2_snr = []
    for i in range(len(falsepos)):
        falsepos2.append( falsepos[i]*float(len(timeslides))/ts_trig_ratio )
    for i in range(len(falsepos_snr)):
        falsepos2_snr.append( falsepos_snr[i]*float(len(timeslides))/ \
                              ts_trig_ratio )
        
    falsepos1,truepos1 = stairs(falsepos2,truepos)
    falsepos1_snr,truepos1_snr = stairs(falsepos2_snr,truepos_snr)

    pyplot.figure()
    pyplot.errorbar(falsepos1,truepos1,yerr=[binomerrl,binomerru] \
                    ,label='MVSC',color='green',linewidth=2)
    pyplot.errorbar(falsepos1_snr,truepos1_snr,yerr=[binomerrl_snr, \
        binomerru_snr], label='Effective SNR', \
                    color='blue',linewidth=2)
    pyplot.semilogx()
    pyplot.xlabel('False alarm number')
    pyplot.ylabel('True positive rate')
    pyplot.title('ROC curve')
    
    op_point = max(op_point,1/(float(len(timeslides))))
    pyplot.plot([op_point,op_point],[.85,1.01],'r',label='Operating point')
    
    for i in range(len(falsepos)):
        if falsepos2[i] <= op_point:
            break
            
    xmin,xmax = pyplot.xlim()
    mid = truepos[i]
    upper = mid + binomerru[2*i]
    lower = mid - binomerrl[2*i]
    pyplot.plot([xmin,xmax],[mid,mid],'r')
    pyplot.plot([xmin,xmax],[lower,lower],'r:')
    pyplot.plot([xmin,xmax],[upper,upper],'r:')
    
    pyplot.xlim(xmin,xmax)

    pyplot.ylim(.85,1.01)
    pyplot.legend(loc='lower right')

    return lower,upper

##############################################################################

def ROCcompare(array,cols,op_point=None):
    """Overlays several ROC plots to compare results."""

    sum_truepos = 0
    total_inj = 0
    
    pyplot.figure()
    for k in range(len(array)):
        data = array[k]
        
        timeslides = []
        injections = []
    
        for i in range(len(data[0])):
            if data[1][i] == 0:
                timeslides.append(data[cols['Bagger']][i])
            else:
                injections.append(data[cols['Bagger']][i])

        truepos,falsepos = ROC(timeslides,injections)

        falsepos1,truepos1 = stairs(falsepos,truepos)
        
        pyplot.semilogx(falsepos1,truepos1)
        
        if op_point:
        
            for i in range(len(falsepos)):
                if falsepos[i] <= op_point:
                    break

            sum_truepos += truepos[i]
            total_inj += len(injections)

    if op_point:
        pyplot.plot([op_point,op_point],[.85,1],'r',label='Operating point')
        
        xmin,xmax = pyplot.xlim()
        mid = sum_truepos / len(array)
        low,high = wilson(mid,total_inj)
        
        upper = mid+high
        lower = mid-low
        pyplot.plot([xmin,xmax],[mid,mid],'r')
        pyplot.plot([xmin,xmax],[lower,lower],'r:')
        pyplot.plot([xmin,xmax],[upper,upper],'r:')
        
        pyplot.xlim(xmin,xmax)

    pyplot.ylim(.85,1)
    pyplot.xlabel('False positive rate')
    pyplot.ylabel('True positive rate')
    pyplot.title('ROC curves of multiple identical runs')
    pyplot.legend(loc='lower right')

    if op_point:
        return mid,upper,lower

##############################################################################

def FOMplot(fom):
    """Plots FOM vs bagger cycle."""

    pyplot.figure()
    pyplot.plot(fom)
    pyplot.xlabel('Bagger cycle')
    pyplot.ylabel('Figure of Merit')
    pyplot.title('Figure of Merit vs number of Bagger cycles')
    pyplot.xlim(-len(fom)*0.05,len(fom))

##############################################################################

def FOMmean(fom):
    """Plots average FOM over several runs."""

    mean = []
    std = []
    
    for i in range(len(fom[0])):
        m = []
        s = []
        for j in range(len(fom)):
            m.append(fom[j][i])
            s.append(fom[j][i]**2)

        sum1 = sum(m)
        sumsq = sum(s)
        mean.append( sum1/len(m) )
        std.append( math.sqrt( max( \
            (sumsq/len(m) - (sum1/len(m))**2)/(len(m)-1),0 ) ) )

    pyplot.figure()
    pyplot.errorbar(range(1,len(mean)+1),mean,yerr=std)
    pyplot.xlabel('Bagger cycle')
    pyplot.ylabel('Figure of Merit')
    pyplot.title('Figure of Merit vs number of Bagger cycles')
    pyplot.xlim(-len(fom[0])*0.05,len(fom[0]))

##############################################################################

def mvsc_vs_effsnr(data,cols, mvsc_cutoff=None,effsnr_cutoff=None,\
                   zerodata=None):
    """Plots mvsc values vs the sum of the squares of effective snr."""

    timeslides = []
    injections = []
    timeslides_snr = []
    injections_snr = []

    for i in range(len(data[0])):
        if data[1][i] == 0:
            timeslides.append(data[cols['Bagger']][i])
            timeslides_snr.append(data[cols['t1get_effective_snr()']][i]**2+ \
                                  data[cols['t2get_effective_snr()']][i]**2)
        else:
            injections.append(data[cols['Bagger']][i])
            injections_snr.append(data[cols['t1get_effective_snr()']][i]**2+ \
                                  data[cols['t2get_effective_snr()']][i]**2)

    pyplot.figure()
    pyplot.semilogy(timeslides,timeslides_snr,'kx',label='Timeslides')
    pyplot.plot(injections,injections_snr,'r+',mec="red",label='Injections')
    pyplot.xlabel('MVSC')
    pyplot.ylabel('Combined effective SNR squared')
    pyplot.title('MVSC vs Effective SNR')

    if zerodata:
        zerolag_snr = []
        for i in range(len(zerodata[0])):
            zerolag_snr.append(zerodata[cols['t1get_effective_snr()']][i]**2+ \
            zerodata[cols['t2get_effective_snr()']][i]**2 )
            
        pyplot.plot(zerodata[cols['Bagger']],zerolag_snr,'*g',mec='g',\
                    label='Zero Lag' )
            
    if mvsc_cutoff:
        ymin,ymax = pyplot.ylim()
        pyplot.plot([mvsc_cutoff,mvsc_cutoff],[ymin,ymax],'b',\
                    label='Operating point cutoff')
    if effsnr_cutoff:
        pyplot.plot([-.1,1.1],[effsnr_cutoff,effsnr_cutoff],'b')

    pyplot.xlim(-.1,1.1)
    pyplot.legend(loc='upper left')

##############################################################################

def fraction_detected(data, cols, afar = None,mvsc_cutoff=None, \
                      distance=False ):
    """Graphs the fraction detected vs snr or distance.

    Returns the false alarm rate (afar), the mvsc cutoff and the effsnr
    cutoff.  Graphs vs snr by default, but if distance=True, will graph
    against distance instead.  Requires either the afar or mvsc_cutoff
    in order to choose operating point."""

    if distance:
        minbin = 1.
        stepfactor = 3
        maxbin = 1000000000
    else:
        #Make histograms based on following snr bins:
        minbin = 60.
        stepfactor = 1.2
        maxbin = 1000
    
    snrbins = [ minbin ]
    i = 1
    while snrbins[-1] <= maxbin:
        snrbins.append( stepfactor**(i) * snrbins[0] )
        i += 1
    
    inj = []
    ts = []

    if distance:
        for i in range(len(data[0])):
            if data[1][i] == 0:
                ts.append( ( data[cols['Bagger']][i], \
                             data[cols['t1get_effective_snr()']][i]**2+ \
                             data[cols['t2get_effective_snr()']][i]**2 ) )
            else:
                inj.append( ( data[cols['Bagger']][i], \
                       data[cols['t1get_effective_snr()']][i]**2+ \
                       data[cols['t2get_effective_snr()']][i]**2, \
                       (data[cols['t1eff_distance']][i]+ \
                       data[cols['t2eff_distance']][i])**3/8 ) )
    else:
        for i in range(len(data[0])):
            if data[1][i] == 0:
                ts.append( ( data[cols['Bagger']][i], \
                         data[cols['t1get_effective_snr()']][i]**2+ \
                         data[cols['t2get_effective_snr()']][i]**2 ) )
            else:
                inj.append( ( data[cols['Bagger']][i], \
                    data[cols['t1get_effective_snr()']][i]**2+ \
                    data[cols['t2get_effective_snr()']][i]**2, \
                    data[cols['t1snr']][i]**2+data[cols['t2snr']][i]**2 ) )

    sorter = cbook.Sorter()
    ts_mvsc = sorter(ts,0)
    ts_effsnr = sorter(ts,1)
    
    if afar:
    #Determine cutoff values which allow only the given number of false alarms
    
        cutoff = int(afar*len(ts)) + 1
        mvsc_cutoff = ts_mvsc[-cutoff][0]
    elif mvsc_cutoff:
        for i in range(1,len(ts_mvsc)+1):
            if ts_mvsc[-i][0] < mvsc_cutoff:
                break

        cutoff = i
        afar = float(cutoff - 1)/len(ts)
    else:
        print '***Error*** Must give fraction_detected either an afar value '+\
              'or a mvsc_cutoff value'

    effsnr_cutoff = ts_effsnr[-cutoff][1]

    #Determine which injections were detected
    det_mvsc = []
    notdet_mvsc = []
    det_effsnr = []
    notdet_effsnr = []
    for i in range(len(inj)):
        if inj[i][0] >= mvsc_cutoff:
            det_mvsc.append(inj[i][2])
        else:
            notdet_mvsc.append(inj[i][2])

        if inj[i][1] >= effsnr_cutoff:
            det_effsnr.append(inj[i][2])
        else:
            notdet_effsnr.append(inj[i][2])

    #Count the number in each bin
    mvsc_frac = []
    effsnr_frac = []
    mvsc_erru = []
    mvsc_errl = []
    effsnr_erru = []
    effsnr_errl = []
    for i in range(len(snrbins)-1):
        det_mvsc_total = 0
        mvsc_total = 0

        for j in range(len(det_mvsc)):
            if (det_mvsc[j] > snrbins[i]) & (det_mvsc[j] <= snrbins[i+1]):
                det_mvsc_total += 1
        for j in range(len(notdet_mvsc)):
            if (notdet_mvsc[j] > snrbins[i]) & \
                   (notdet_mvsc[j] <= snrbins[i+1]):
                mvsc_total += 1

        mvsc_total += det_mvsc_total

        det_effsnr_total = 0
        effsnr_total = 0
        
        for j in range(len(det_effsnr)):
            if (det_effsnr[j] > snrbins[i]) & (det_effsnr[j] <= snrbins[i+1]):
                det_effsnr_total += 1
        for j in range(len(notdet_effsnr)):
            if (notdet_effsnr[j] > snrbins[i]) & \
                   (notdet_effsnr[j] <= snrbins[i+1]):
                effsnr_total += 1

        effsnr_total += det_effsnr_total

        #Calculate the fraction detected in each bin:
        try:
            mvsc_frac.append(float(det_mvsc_total)/mvsc_total)
            upper,lower = wilson(float(det_mvsc_total)/mvsc_total, \
                mvsc_total )
            mvsc_erru.append(upper)
            mvsc_errl.append(lower)
        except ZeroDivisionError:
            mvsc_frac.append(0.5)
            mvsc_erru.append(0.5)
            mvsc_errl.append(0.5)
        try:
            effsnr_frac.append(float(det_effsnr_total)/effsnr_total)
            upper,lower = wilson( float(det_effsnr_total)/effsnr_total, \
                effsnr_total )
            effsnr_erru.append(upper)
            effsnr_errl.append(lower)
        except ZeroDivisionError:
            effsnr_frac.append(0.5)
            effsnr_erru.append(0.5)
            effsnr_errl.append(0.5)

    pyplot.figure()
    pyplot.semilogx()
    pyplot.errorbar(snrbins[:-1],effsnr_frac,marker='s',linestyle='-', \
        yerr=[effsnr_errl,effsnr_erru],label='Effective SNR',mec='blue')
    pyplot.errorbar(snrbins[:-1],mvsc_frac,marker='*',linestyle='-', \
        yerr=[mvsc_errl,mvsc_erru],label='MVSC',mec='green')
    pyplot.xlim(snrbins[0]/stepfactor**2,snrbins[-1]*stepfactor**2)
    pyplot.ylim(-.1,1.1)

    if distance:
        pyplot.xlabel('Effective Volume')
        pyplot.legend(loc='lower left')
    else:
        pyplot.xlabel('Combined SNR squared')
        pyplot.legend(loc='lower right')
        
    pyplot.ylabel('Fraction detected')
    pyplot.title('Fraction of injections detected using different measures')

    return afar,mvsc_cutoff,effsnr_cutoff

##############################################################################

def snr_vs_chisqr(data,cols,afar = 1.0/2000,zerodata=None ):
    """Plots SNR vs Chi Squared, indicating which triggers were correctly
    classified, and which were incorrectly classified.

    afar is the allowed false alarm rate in the classification"""
    
    inj = []
    ts = []

    #separate into timeslides and injections
    for i in range(len(data[0])):
        if data[1][i] == 0:
            ts.append( ( data[cols['Bagger']][i], \
                        data[cols['t1chisq']][i]**2+ \
                        data[cols['t2chisq']][i]**2, \
                        data[cols['t1snr']][i]**2+ \
                        data[cols['t2snr']][i]**2 ) )
        else:
            inj.append( ( data[cols['Bagger']][i], \
                        data[cols['t1chisq']][i]**2+ \
                        data[cols['t2chisq']][i]**2, \
                        data[cols['t1snr']][i]**2+ \
                        data[cols['t2snr']][i]**2 ) )

    cutoff = int(afar*len(ts)) + 1
    
    #Determine cutoff value which allows only the given number of false alarms
    sorter = cbook.Sorter()
    ts_mvsc = sorter(ts,0)
    mvsc_cutoff = ts_mvsc[-cutoff][0]

    #Classify all triggers
    falsepos_chi = []
    falsepos_snr = []
    trueneg_chi = []
    trueneg_snr = []
    truepos_chi = []
    truepos_snr = []
    falseneg_chi = []
    falseneg_snr = []
    for i in range(len(ts)):
        if ts[i][0] >= mvsc_cutoff:
            falsepos_chi.append(ts[i][1])
            falsepos_snr.append(ts[i][2])
        else:
            trueneg_chi.append(ts[i][1])
            trueneg_snr.append(ts[i][2])
    for i in range(len(inj)):
        if inj[i][0] >= mvsc_cutoff:
            truepos_chi.append(inj[i][1])
            truepos_snr.append(inj[i][2])
        else:
            falseneg_chi.append(inj[i][1])
            falseneg_snr.append(inj[i][2])

    pyplot.figure()
    pyplot.loglog(trueneg_snr,trueneg_chi,'xk',label='Filtered timeslides')
    pyplot.plot(truepos_snr,truepos_chi,'+r',mec='r', \
                label='Detected injections')
    pyplot.plot(falseneg_snr,falseneg_chi,'oc',mec='c', \
                label='Missed injections')
    pyplot.plot(falsepos_snr,falsepos_chi,'sm',mec='m',label='False alarms')

    if zerodata:
        zero = []
        for i in range(len(zerodata[0])):
            zero.append( ( zerodata[cols['Bagger']][i],
                        zerodata[cols['t1chisq']][i]**2+ \
                        zerodata[cols['t2chisq']][i]**2, \
                        zerodata[cols['t1snr']][i]**2+ \
                        zerodata[cols['t2snr']][i]**2 ) )

        zeropos_chi = []
        zeropos_snr = []
        zeroneg_chi = []
        zeroneg_snr = []
        for i in range(len(zero)):
            if zero[i][0] >= mvsc_cutoff:
                zeropos_chi.append(zero[i][1])
                zeropos_snr.append(zero[i][2])
            else:
                zeroneg_chi.append(zero[i][1])
                zeroneg_snr.append(zero[i][2])

        pyplot.plot(zeroneg_snr,zeroneg_chi,'xb', \
            label='Zero lag (noise)')
        pyplot.plot(zeropos_snr,zeropos_chi,'*g',mec='g', \
            label='Zero lag (signal)')
        
    pyplot.legend(loc='upper left')
    pyplot.ylabel('Combined Chisq squared')
    pyplot.xlabel('Combined SNR squared')
    pyplot.title('SNR vs Chi Squared')
 
##############################################################################

def FARplot(data,cols,zerodata=None,ts_trig_ratio = 25):
    """Graphs the cumulative number of detections vs MVSC threshold."""
    
    ts_mvsc = []

    for i in range(len(data[0])):
        if data[1][i] == 0:
            ts_mvsc.append(data[cols['Bagger']][i])

    ts_mvsc.sort()
    n = len(ts_mvsc)
    
    far = []
    mvsc = []
    mvsc_cutoff = 1
    flag = True
    
    for i in range(1,n):
        if ts_mvsc[i] != ts_mvsc[i-1]:
            far.append(float(n-i)/ts_trig_ratio)
            mvsc.append(ts_mvsc[i])
            
            if flag & (far[-1] < 1):
                mvsc_cutoff = mvsc[-1]
                flag = False

    far2,mvsc2 = stairs(far,mvsc)

    binomerru = []
    binomerrl = []
    for i in range(len(far)):
        upper,lower = wilson( (far[i]/n)*ts_trig_ratio, n )
        binomerru.append((upper*n)/ts_trig_ratio)
        binomerrl.append((lower*n)/ts_trig_ratio)
        binomerru.append(0)
        binomerrl.append(0)
    t = binomerru.pop()
    t = binomerrl.pop()
    
    pyplot.figure()
    pyplot.errorbar(mvsc2,far2,yerr=[binomerrl,binomerru], \
                    label='Expected FAN (based on timeslides)')
    pyplot.loglog()

    if zerodata:
        zero_mvsc = zerodata[cols['Bagger']][:]
        zero_mvsc.sort()
        m = len(zero_mvsc)
        
        zerorate = []
        zero_mvsc1 = []
        
        for i in range(1,m):
            if zero_mvsc[i] != zero_mvsc[i-1]:
                zerorate.append(float(m-i))
                zero_mvsc1.append(zero_mvsc[i])
            
        zerorate2,zero_mvsc2 = stairs(zerorate,zero_mvsc1)

        zbinomerru = []
        zbinomerrl = []
        for i in range(len(zerorate)):
            upper,lower = wilson( zerorate[i]/m, m )
            zbinomerru.append(upper*m)
            zbinomerrl.append(lower*m)
            zbinomerru.append(0)
            zbinomerrl.append(0)
        t = zbinomerru.pop()
        t = zbinomerrl.pop()

        pyplot.errorbar(zero_mvsc2,zerorate2,yerr=[zbinomerrl,zbinomerru],\
                        label='cumulative number of zero lags')

    pyplot.legend(loc='lower left')
    
    xmin,xmax = pyplot.xlim()
    ymin,ymax = pyplot.ylim()
    pyplot.plot([xmin,xmax],[1,1],'r')
    pyplot.plot([mvsc_cutoff,mvsc_cutoff],[ymin,ymax],'r')
    pyplot.xlim(xmin,xmax)
    pyplot.ylim(ymin,ymax)
    
    pyplot.xlabel('MVSC value cutoff')
    pyplot.title('FAR plot')

    return mvsc_cutoff

##############################################################################

def IFANplot(data,cols,zerodata,ts_trig_ratio=25):
    """Plots the inverse false alarm number vs cumulative number of detections.
    """
    
    ts_mvsc = []

    for i in range(len(data[0])):
        if data[1][i] == 0:
            ts_mvsc.append(data[cols['Bagger']][i])

    ts_mvsc.sort()
    n = len(ts_mvsc)

    zero_mvsc = zerodata[cols['Bagger']][:]
    zero_mvsc.sort()
    m = len(zero_mvsc)

    cumnumber = []
    cumnumber_erru = []
    cumnumber_errl = []
    ifan = []
    ifan_erru = []
    ifan_errl = []
    for i in range(1,m):
        if zero_mvsc[i] != zero_mvsc[i-1]:
            cumnumber.append(float(m-i))
            upper,lower = wilson(cumnumber[-1]/m,m)
            cumnumber_erru.append(upper*m)
            cumnumber_errl.append(lower*m)

            for j in range(n):
                if ts_mvsc[j] >= zero_mvsc[i]:
                    break
                
            ifan.append(ts_trig_ratio/float(n-j))
            upper1,lower1 = wilson(float(n-j)/n,n)
            upper2 = lower1*n/(float(n-j)*(float(n-j)-lower1*n))
            lower2 = upper1*n/(float(n-j)*(float(n-j)+upper1*n))
            ifan_erru.append(upper2*ts_trig_ratio)
            ifan_errl.append(lower2*ts_trig_ratio)

    pyplot.figure()
    pyplot.errorbar(ifan,cumnumber,xerr=[ifan_errl,ifan_erru], \
                    yerr=[cumnumber_errl,cumnumber_erru])
    pyplot.loglog()

    ymin,ymax = pyplot.ylim()
    xmin,xmax = pyplot.xlim()
    pyplot.plot([xmin,xmax],[-xmin,-xmax],'r',alpha=0.5)
    pyplot.ylim(ymin,ymax)
    pyplot.xlim(xmin,xmax)
    
    pyplot.xlabel('Inverse False Alarm Number')
    pyplot.ylabel('Cumulative Number')
    pyplot.title('IFAN plot')
