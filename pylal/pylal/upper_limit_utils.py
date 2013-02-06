import numpy
from scipy import random
from scipy import interpolate
import bisect
import sys

from glue.ligolw import lsctables
from glue.ligolw import dbtables
from pylal.xlal import constants
from pylal.imr_utils import make_sim_inspiral_row_from_columns_in_db, time_within_segments
from pylal import rate

import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot


def margLikelihoodMonteCarlo(VTs, lambs, mu, mcerrs=None):
    '''
    This function marginalizes the loudest event likelihood over unknown
    Monte Carlo errors, assumed to be independent between each experiment.
    '''
    if mcerrs is None:
        mcerrs = [0]*len(VTs)

    # combine experiments, propagating statistical uncertainties
    # in the measured efficiency
    likely = 1
    for vA,dvA,mc in zip(VTs,lambs,mcerrs):
        if mc == 0:
            # we have perfectly measured our efficiency in this mass bin
            # so the posterior is given by eqn (11) in Biswas et al. [arXiv:0710.0465]
            likely *= (1+mu*vA*dvA)*numpy.exp(-mu*vA)
        else:
            # we have uncertainty in our efficiency in this mass bin and
            # want to marginalize it out using eqn (24) of Biswas et al.
            k = (vA/mc)**2 # k is 1./fractional_error**2
            likely *= (1+mu*vA*(1/k+dvA))*(1+mu*vA/k)**(-(k+1))

    return likely


def margLikelihood(VTs, lambs, mu, calerr=0, mcerrs=None):
    '''
    This function marginalizes the loudest event likelihood over unknown
    Monte Carlo and calibration errors. The vector VTs is the sensitive
    volumes for independent searches and lambs is the vector of loudest
    event likelihood. The statistical errors are assumed to be independent
    between each experiment while the calibration errors are applied
    the same in each experiment.
    '''
    if calerr == 0:
        return margLikelihoodMonteCarlo(VTs,lambs,mu,mcerrs)

    std = calerr
    mean = 0 # median volume = 1

    fracerrs = numpy.linspace(0.33,3,1e2) # assume we got the volume to a factor of three or better
    errdist = numpy.exp(-(numpy.log(fracerrs)-mean)**2/(2*std**2))/(fracerrs*std) # log-normal pdf
    errdist /= errdist.sum() #normalize

    likely = sum([ pd*margLikelihoodMonteCarlo(delta*VTs,lambs,mu,mcerrs) for delta, pd in zip(fracerrs,errdist)]) #marginalize over errors

    return likely


def normalize_pdf(mu, pofmu):
    """
    Takes a function pofmu defined at rate sample values mu and 
    normalizes it to be a suitable pdf. Both mu and pofmu must be
    arrays or lists of the same length. 
    """
    if min(pofmu) < 0:
        raise ValueError, "Probabilities cannot be negative, please don't \
        ask me to normalize a function with negative values!"
    if min(mu) < 0:
        raise ValueError, "Rates cannot be negative, please don't \
        ask me to normalize a function over a negative domain!"
    # simple trapezoidal integral
    dmu = mu[1:] - mu[:-1]
    mean_in_bin = (pofmu[1:] + pofmu[:-1]) /2
    dp = dmu*mean_in_bin

    return mu, pofmu/sum(dp)


def compute_upper_limit(mu_in, post, alpha = 0.9):
    """
    Returns the upper limit mu_high of confidence level alpha for a
    posterior distribution post on the given parameter mu.
    """
    # simple trapezoidal integral
    dmu = mu_in[1:] - mu_in[:-1]
    mean_post = (post[1:] + post[:-1]) /2
    dp = dmu*mean_post

    if 0 < alpha < 1:
        high_idx = bisect.bisect_left( dp.cumsum()/dp.sum(), alpha )
        # if alpha is in (0,1] and post is non-negative, bisect_left
        # will always return an index in the range of mu since
        # post.cumsum()/post.sum() will always begin at 0 and end at 1
        mu_high = mu_in[high_idx]
    elif alpha == 1:
        mu_high = numpy.max(mu_in[post>0])
    else:
        raise ValueError, "Confidence level must be in (0,1]."

    return mu_high


def compute_lower_limit(mu_in, post, alpha = 0.9):
    """
    Returns the lower limit mu_low of confidence level alpha for a
    posterior distribution post on the given parameter mu.
    """
    # simple trapezoidal integral
    dmu = mu_in[1:] - mu_in[:-1]
    mean_post = (post[1:] + post[:-1]) /2
    dp = dmu*mean_post

    if 0 < alpha < 1:
        low_idx = bisect.bisect_right( dp.cumsum()/dp.sum(), 1-alpha )
        # if alpha is in [0,1) and post is non-negative, bisect_right
        # will always return an index in the range of mu since
        # post.cumsum()/post.sum() will always begin at 0 and end at 1
        mu_low = mu_in[low_idx]
    elif alpha == 1:
        mu_low = numpy.min(mu_in[post>0])
    else:
        raise ValueError, "Confidence level must be in (0,1]."

    return mu_low


def confidence_interval( mu, post, alpha = 0.9 ):
    '''
    Returns the minimal-width confidence interval [mu_low,mu_high] of
    confidence level alpha for a distribution post on the parameter mu.
    '''
    if not 0 < alpha < 1:
        raise ValueError, "Confidence level must be in (0,1)."

    # choose a step size for the sliding confidence window
    alpha_step = 0.01

    # initialize the lower and upper limits
    mu_low = numpy.min(mu)
    mu_high = numpy.max(mu)

    # find the smallest window (by delta-mu) stepping by dalpha
    for ai in numpy.arange( 0, 1-alpha, alpha_step ):
        ml = compute_lower_limit( mu, post, 1 - ai )
        mh = compute_upper_limit( mu, post, alpha + ai)
        if mh - ml < mu_high - mu_low:
            mu_low = ml
            mu_high = mh

    return mu_low, mu_high


def integrate_efficiency(dbins, eff, err=0, logbins=False):

    if logbins:
        logd = numpy.log(dbins)
        dlogd = logd[1:]-logd[:-1]
        dreps = numpy.exp( (numpy.log(dbins[1:])+numpy.log(dbins[:-1]))/2) # log midpoint
        vol = numpy.sum( 4*numpy.pi *dreps**3 *eff *dlogd )
        verr = numpy.sqrt(numpy.sum( (4*numpy.pi *dreps**3 *err *dlogd)**2 )) #propagate errors in eff to errors in v
    else:
        dd = dbins[1:]-dbins[:-1]
        dreps = (dbins[1:]+dbins[:-1])/2 #midpoint
        vol = numpy.sum( 4*numpy.pi *dreps**2 *eff *dd )
        verr = numpy.sqrt(numpy.sum( (4*numpy.pi *dreps**2 *err *dd)**2 )) #propagate errors in eff to errors in v

    return vol, verr


def compute_efficiency(f_dist,m_dist,dbins):
    '''
    Compute the efficiency as a function of distance for the given sets of found
    and missed injection distances.
    Note that injections that do not fit into any dbin get lost :(.
    '''
    efficiency = numpy.zeros( len(dbins)-1 )
    error = numpy.zeros( len(dbins)-1 )
    for j, dlow in enumerate(dbins[:-1]):
        dhigh = dbins[j+1]
        found = numpy.sum( (dlow <= f_dist)*(f_dist < dhigh) )
        missed = numpy.sum( (dlow <= m_dist)*(m_dist < dhigh) )
        if found+missed == 0: missed = 1.0 #avoid divide by 0 in empty bins
        efficiency[j] = 1.0*found /(found + missed)
        error[j] = numpy.sqrt(efficiency[j]*(1-efficiency[j])/(found+missed))

    return efficiency, error


def mean_efficiency_volume(found, missed, dbins):

    if len(found) == 0: # no efficiency here
        return numpy.zeros(len(dbins)-1),numpy.zeros(len(dbins)-1), 0, 0

    # only need distances
    f_dist = numpy.array([l.distance for l in found])
    m_dist = numpy.array([l.distance for l in missed])

    # compute the efficiency and its variance
    eff, err = compute_efficiency(f_dist,m_dist,dbins)
    vol, verr = integrate_efficiency(dbins, eff, err)

    return eff, err, vol, verr


def filter_injections_by_mass(injs, mbins, bin_num , bin_type, bin_num2=None):
    '''
    For a given set of injections (sim_inspiral rows), return the subset
    of injections that fall within the given mass range.
    '''

    if bin_type == "Mass1_Mass2":
        m1bins = numpy.concatenate((mbins.lower()[0],numpy.array([mbins.upper()[0][-1]])))
        m1lo = m1bins[bin_num]
        m1hi = m1bins[bin_num+1]
        m2bins = numpy.concatenate((mbins.lower()[1],numpy.array([mbins.upper()[1][-1]])))
        m2lo = m2bins[bin_num2]
        m2hi = m2bins[bin_num2+1]
        newinjs = [l for l in injs if ( (m1lo<= l.mass1 <m1hi and m2lo<= l.mass2 <m2hi) or (m1lo<= l.mass2 <m1hi and m2lo<= l.mass1 <m2hi))]
        return newinjs

    mbins = numpy.concatenate((mbins.lower()[0],numpy.array([mbins.upper()[0][-1]])))
    mlow = mbins[bin_num]
    mhigh = mbins[bin_num+1]
    if bin_type == "Chirp_Mass":
        newinjs = [l for l in injs if (mlow <= l.mchirp < mhigh)]
    elif bin_type == "Total_Mass":
        newinjs = [l for l in injs if (mlow <= l.mass1+l.mass2 < mhigh)]
    elif bin_type == "Component_Mass": #it is assumed that m2 is fixed
        newinjs = [l for l in injs if (mlow <= l.mass1 < mhigh)]
    elif bin_type == "BNS_BBH":
        if bin_num == 0 or bin_num == 2: #BNS/BBH case
            newinjs = [l for l in injs if (mlow <= l.mass1 < mhigh and mlow <= l.mass2 < mhigh)]
        else:
            newinjs = [l for l in injs if (mbins[0] <= l.mass1 < mbins[1] and mbins[2] <= l.mass2 < mbins[3])] #NSBH
            newinjs += [l for l in injs if (mbins[0] <= l.mass2 < mbins[1] and mbins[2] <= l.mass1 < mbins[3])] #BHNS

    return newinjs


def compute_volume_vs_mass(found, missed, mass_bins, bin_type, dbins=None, ploteff=False):
    """
    Compute the average luminosity an experiment was sensitive to given the sets
    of found and missed injections and assuming luminosity is unformly distributed
    in space.
    """
    # mean and std estimate for luminosity (in L10s)
    volArray = rate.BinnedArray(mass_bins)
    vol2Array = rate.BinnedArray(mass_bins)

    # found/missed stats
    foundArray = rate.BinnedArray(mass_bins)
    missedArray = rate.BinnedArray(mass_bins)

    #
    # compute the mean luminosity in each mass bin
    #
    effvmass = []
    errvmass = []

    if bin_type == "Mass1_Mass2": # two-d case first
        for j,mc1 in enumerate(mass_bins.centres()[0]):
            for k,mc2 in enumerate(mass_bins.centres()[1]):
                newfound = filter_injections_by_mass( found, mass_bins, j, bin_type, k)
                newmissed = filter_injections_by_mass( missed, mass_bins, j, bin_type, k)

                foundArray[(mc1,mc2)] = len(newfound)
                missedArray[(mc1,mc2)] = len(newmissed)

                # compute the volume using this injection set
                meaneff, efferr, meanvol, volerr = mean_efficiency_volume(newfound, newmissed, dbins)
                effvmass.append(meaneff)
                errvmass.append(efferr)
                volArray[(mc1,mc2)] = meanvol
                vol2Array[(mc1,mc2)] = volerr

        return volArray, vol2Array, foundArray, missedArray, effvmass, errvmass


    for j,mc in enumerate(mass_bins.centres()[0]):

        # filter out injections not in this mass bin
        newfound = filter_injections_by_mass( found, mass_bins, j, bin_type)
        newmissed = filter_injections_by_mass( missed, mass_bins, j, bin_type)

        foundArray[(mc,)] = len(newfound)
        missedArray[(mc,)] = len(newmissed)

        # compute the volume using this injection set
        meaneff, efferr, meanvol, volerr = mean_efficiency_volume(newfound, newmissed, dbins)
        effvmass.append(meaneff)
        errvmass.append(efferr)
        volArray[(mc,)] = meanvol
        vol2Array[(mc,)] = volerr

    return volArray, vol2Array, foundArray, missedArray, effvmass, errvmass


def log_volume_derivative_fit(x, vols):
    '''
    Performs a linear least squares to log(vols) as a function of x.
    '''
    if numpy.min(vols) == 0:
        print >> sys.stderr, "Warning: cannot fit log volume derivative, one or more volumes are zero!"
        print >> sys.stderr, vols
        return (0,0)

    coeffs, resids, rank, svs, rcond = numpy.polyfit(x, numpy.log(vols), 1, full=True)
    if coeffs[0] < 0: #negative derivatives may arise from rounding error
        print >> sys.stderr, "Warning: Derivative fit resulted in Lambda =", coeffs[0]
        coeffs[0] = 0
        print >> sys.stderr, "The value Lambda = 0 was substituted"

    return coeffs


#
# FIXME: this function is only a slight variant on a function that
# already exists in imr_utils.py. The difference here make the present
# function integrate more smoothly into svim. In the future, however,
# we would like to switch over to using the function in imr_utils.py
# so that there is only one such function.
#
def get_all_min_far_inspiral_injections(connection, segments = None, table_name = "coinc_inspiral"):
	"""
	This function returns the found injections from a database and the
	minimum far associated with them as tuple of the form (far, sim). It also tells
	you all of the injections that should have been injected.  Subtracting the two
	outputs	should tell you the missed injections
	"""

	if table_name == dbtables.lsctables.CoincInspiralTable.tableName:
		found_query = 'SELECT sim_inspiral.*, coinc_inspiral.combined_far FROM sim_inspiral JOIN coinc_event_map AS mapA ON mapA.event_id == sim_inspiral.simulation_id JOIN coinc_event_map AS mapB ON mapB.coinc_event_id == mapA.coinc_event_id JOIN coinc_inspiral ON coinc_inspiral.coinc_event_id == mapB.event_id JOIN coinc_event on coinc_event.coinc_event_id == coinc_inspiral.coinc_event_id WHERE mapA.table_name = "sim_inspiral" AND mapB.table_name = "coinc_event" AND injection_in_segments(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)'

	elif table_name == dbtables.lsctables.CoincRingdownTable.tableName:
		found_query = 'SELECT sim_inspiral.*, coinc_ringdown.false_alarm_rate FROM sim_inspiral JOIN coinc_event_map AS mapA ON mapA.event_id == sim_inspiral.simulation_id JOIN coinc_event_map AS mapB ON mapB.coinc_event_id == mapA.coinc_event_id JOIN coinc_ringdown ON coinc_ringdown.coinc_event_id == mapB.event_id JOIN coinc_event on coinc_event.coinc_event_id == coinc_ringdown.coinc_event_id WHERE mapA.table_name = "sim_inspiral" AND mapB.table_name = "coinc_event" AND injection_in_segments(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)'

	elif table_name == dbtables.lsctables.MultiBurstTable.tableName:
		found_query = 'SELECT sim_inspiral.*, multi_burst.false_alarm_rate FROM sim_inspiral JOIN coinc_event_map AS mapA ON mapA.event_id == sim_inspiral.simulation_id JOIN coinc_event_map AS mapB ON mapB.coinc_event_id == mapA.coinc_event_id JOIN multi_burst ON multi_burst.coinc_event_id == mapB.event_id JOIN coinc_event on coinc_event.coinc_event_id == multi_burst.coinc_event_id WHERE mapA.table_name = "sim_inspiral" AND mapB.table_name = "coinc_event" AND injection_in_segments(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)'

	else:
		raise ValueError("table must be in " + " ".join(allowed_analysis_table_names()))

	def injection_was_made(end_time, end_time_ns, segments = segments):
		return time_within_segments(end_time, end_time_ns, segments)

	# restrict the found injections to only be within certain segments
	connection.create_function("injection_in_segments", 2, injection_was_made)

	# get the mapping of a record returned by the database to a sim
	# inspiral row. Note that this is DB dependent potentially, so always
	# do this!
	make_sim_inspiral = make_sim_inspiral_row_from_columns_in_db(connection)

	found_injections = {}

	total_query = 'SELECT * FROM sim_inspiral WHERE injection_in_segments(geocent_end_time, geocent_end_time_ns)'
	for values in connection.cursor().execute(total_query):
		sim = make_sim_inspiral(values)
		found_injections.setdefault(sim.simulation_id, (float('inf'),sim))

	for values in connection.cursor().execute(found_query):
		# all but the last column is used to build a sim inspiral object
		sim = make_sim_inspiral(values[:-1])
		far = values[-1]
		# update with the minimum far seen until now
		this_inj = found_injections[sim.simulation_id]
		if far < this_inj[0]:
			found_injections[sim.simulation_id] = (far, sim)

	return found_injections.values()
