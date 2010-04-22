MultiNest v 2.7
Farhan Feroz, Mike Hobson
f.feroz@mrao.cam.ac.uk
arXiv:0704.3704 & arXiv:0809.3437
Released June 2009

---------------------------------------------------------------------------

MultiNest library requires lapack library and the Makefile uses Intel compilers so you might have to update the
Makefiles accordingly.

The code is MPI compatible. In order to disable the MPI parallelization, remove -DMPI compilation flag.

You might need to change the stacksize in order for MultiNest to run interactively. 256MB should be adequate.
The command to do it on Unix systems with different shells is;

sh & ksh:
ulimit -s [stacksize in KB]

csh:
limit stacksize [stacksize in KB]

---------------------------------------------------------------------------

gfortran compiler:

You might need to use the following flag while compiling MultiNest with gfortran compiler to remove the
restriction imposed by gfortran on line length.

-ffree-line-length-none

---------------------------------------------------------------------------

The subtoutine to begin MultiNest are as follows:

subroutine nestRun(mmodal, ceff, nlive, tol, efr, ndims, nPar, nCdims, maxModes, updInt, nullZ, root, seed,
pWrap, feedback, resume, loglike, context)

logical mmodal 	 		!do mode separation?
integer nlive 	 		!number of live points
logical ceff			!run in constant efficiency mode
real*8 tol 		 	!evidence tolerance factor
real*8 efr 		 	!sampling efficiency
integer ndims	 		!number of dimensions
integer nPar	 		!total no. of parameters
integer nCdims			!no. of parameters on which clustering should be performed (read below)
integer maxModes		!maximum no. of modes (for memory allocation)
integer updInt			!iterations after which the output files should be written
real*8 nullZ			!null log-evidence (read below)
character(LEN=100) root  	!root for MultiNest output files
integer seed 		 	!random no. generator seed, -ve value for seed from the sys clock
integer pWrap[ndims]		!wraparound parameters?
logical feedback		!need update on sampling progress?
logical resume			!resume from a previous run?
loglike(Cube,ndims,nPar,lnew) 	!subroutine which gives lnew=loglike(Cube(ndims))
	integer ndims
      	real*8 Clube(nPar)
      	real*8 lnew
integer context			!dummy integer for the C interface
      
--------------------------------------------------------------------------- 

likelihood routine: slikelihood(Cube,ndims,nPar,lnew)

Cube(1:nPar) has nonphysical parameters
scale Cube(1:n_dim) & return the scaled parameters in Cube(1:n_dim) & additional parameters that you want to
returned by MultiNest along with the actual parameters in Cube(n_dim+1:nPar)
Return the log-likelihood in lnew

---------------------------------------------------------------------------

Checkpointing:

MultiNest is able to checkpoint. It creates [root]resume.dat file & stores information in it after every
updInt iterations to checkpoint, where updInt is set by the user. If you don't want to resume your program from
the last run run then make sure that you either delete [root]resume.dat file or set the parameter resume to F
before starting the sampling.

---------------------------------------------------------------------------

Constant Efficiency Mode:

If ceff is set to T, then the enlargement factor of the bounding ellipsoids are tuned so that the sampling
efficiency is as close to the target efficiency (set by efr) as possible. This does mean however, that the
evidence value may not be accurate.

---------------------------------------------------------------------------

Periodic Boundary Conditions:

In order to sample from parameters with periodic boundary conditions (or wraparound parameters), set pWrap[i],
where i is the index of the parameter to be wraparound, to a non-zero value. If pWrap[i] = 0, then the ith
parameter is not wraparound.

---------------------------------------------------------------------------

Constant Efficiency Mode:

If ceff is set to T, then the enlargement factor of the bounding ellipsoids are tuned so that the sampling
efficiency is as close to the target efficiency (set by efr) as possible. This does mean however, that the
evidence value may not be accurate.

---------------------------------------------------------------------------

Sampling Parameters:

Here I describe the recommended paramter values to be used with the MultiNest. For detailed description please
refer to the paper arXiv:0809.3437

nPar: 
Total no. of parameters, should be equal to ndims in most cases but if you need to store some additional
parameters with the actual parameters then you need to pass them through the likelihood routine.


efr:
defines the sampling efficiency. 0.8 and 0.3 are recommended for parameter estimation & evidence evalutaion
respectively.


tol:
A value of 0.5 should give good enough accuracy.

           
nCdims:
If mmodal is T, MultiNest will attempt to separate out the modes. Mode separation is done through a clustering
algorithm. Mode separation can be done on all the parameters (in which case nCdims should be set to ndims) & it
can also be done on a subset of parameters (in which case nCdims < ndims) which might be advantageous as
clustering is less accurate as the dimensionality increases. If nCdims < ndims then mode separation is done on
the first nCdims parameters.


nullZ:
If mmodal is T, MultiNest can find multiple modes & also specify which samples belong to which mode. It might be
desirable to have separate samples & mode statistics for modes with local log-evidence value greater than a
particular value in which case nullZ should be set to that value. If there isn't any particulrly interesting
nullZ value, then nullZ should be set to a very large negative number (e.g. -1.d90).

---------------------------------------------------------------------------

Progress Monitoring:

MultiNest produces [root]physlive.dat & [root]ev.dat files after every updInt iterations which can be used to
monitor the progress. The format & contents of  these two files are as follows:

[root]physlive.dat:
This file contains the current set of live points. It has nPar+2 columns. The first nPar columns are the ndim
parameter values along with the (nPar-ndim)  additional parameters that are being passed by the likelihood
routine for MultiNest to save along with the ndim parameters. The nPar+1 column is the log-likelihood value &
the last column is the node no. (used for clustering).

[root]ev.dat:
This file contains the set of rejected points. It has nPar+3 columns. The first nPar columns are the ndim
parameter values along with the (nPar-ndim)  additional parameters that are being passed by the likelihood
routine for MultiNest to save along with the ndim parameters. The nPar+1 column is the log-likelihood value,
nPar+2 column is the log(prior mass) & the last column  is the node no. (used for clustering).

---------------------------------------------------------------------------

Posterior Samples:

These files are created after every updInt*10 iterations of the algorithm & at the end of sampling.

MultiNest will produce five posterior sample files in the root, given by the user, as following


[root].txt
Compatable with getdist with 2+nPar columns. Columns have sample probability, -2*loglikehood, samples. Sample
probability is the sample prior mass multiplied by its likelihood & normalized by the evidence.


[root]post_separate.dat
This file is only created if mmodal is set to T. Posterior samples for modes with local log-evidence value
greater than nullZ, separated by 2 blank lines. Format is the same as [root].txt file.


[root]stats.dat
Contains the global log-evidence, its error & local log-evidence with error & parameter means & standard
deviations as well as the  best fit & MAP parameters of each of the mode found with local log-evidence > nullZ.


[root]post_equal_weights.dat
Contains the equally weighted posterior samples

---------------------------------------------------------------------------

C Interface:

A C interface to MultiNest is provided in example_eggbox directory.

---------------------------------------------------------------------------

Toy Problems

There are 3 toy programs included with MultiNest.

example_obj_detect: The object detection problem discussed in arXiv:0704.3704. The positions, amplitudes &
widths of the Gaussian objects can be modified through params.f90 file. Sampling parameters are also set in
params.f90.

example_gauss_shell: The Gaussian shells problem discussed in arXiv:0809.3437. The dimensionality, positions and
thickness of these shells can be modified through params.f90 file. Sampling parameters are also set in
params.f90.

example_eggbox: The C/C++ interface includes the egg box problem discussed in arXiv:0809.3437. The toy 
problem and sampling parameters are set in cnest.cc file.
