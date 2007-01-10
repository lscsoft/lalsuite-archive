Python frgetvect for Pylal
--------------------------
This module is intended to allow the reading and writing of gravitational
wave Frame files in Python.  It interfaces lalframe.

Written by Nickolas Fotopoulos (nvf@gravity.phys.uwm.edu)

Prerequisites
-------------
numpy (http://www.scipy.org/NumPy)
LAL

Usage
-----
In a Python script or interactive session, simply "import pylal.Fr" and call
pylal.Fr.frgetvect and pylal.Fr.putvect as indicated below.

Function docstrings
-----------------
frgetvect(filename, channel, start=-1, span=-1, verbose=False)
    
    Python adaptation of frgetvect (based on Matlab frgetvect).
    Reads a vector from a Fr frame to a numpy array.

    The input arguments are:
       1) filename - accepts name of frame file or FFL.
       2) channel - channel name (ADC, SIM or PROC channel)
       3) start - starting frequency or GPS time (default = -1)
                  A value <=0 will read from the first available time/frequency.       4) span - quantity of data in seconds (default = -1)
                 A value <=0 will return the whole vector of the first frame
                 (and in the first file, for the case of an FFL).
       5) verbose - Verbose (True) or silent (False) (default = False)
    
    Returned data (in a tuple):
       1) ADC, SIM, or PROC data
       2) start frequency/GPS time
       3) x-axis spacing (dt for time-series, df for frequency-series)
       4) Unit of x-axis as a string
       5) Unit of y-axis as a string

frputvect(filename, channellist, history='', verbose=False)
    
    The inverse of frgetvect -- write numpy arrays to a Fr frame file.
    
    The input arguments are:
        1) filename - name of file to write.
        2) channellist - list of dictionaries with the fields below:
            1) name - channel name - string
            2) data - list of one-dimensional vectors to write
            3) start - lower limit of the x-axis in GPS time or Hz
            4) dx - spacing of x-axis in seconds or Hz
            5) x_unit - unit of x-axis as a string (default = '')
            6) y_unit - unit of y-axis as a string (default = '')
            7) kind - 'PROC', 'ADC', or 'SIM' (default = 'PROC')
            8) type - type of data (default = 0):
                   0 - Unknown/undefined
                   1 - Time series
                   2 - Frequency series
                   3 - Other 1-D series
                   4 - Time-frequency
                   5 - Wavelets
                   6 - Multi-dimensional
            9) subType - sub-type of frequency series (default = 0):
                   0 - Unknown/undefined
                   1 - DFT
                   2 - Amplitude spectral density
                   3 - Power spectral density
                   4 - Cross spectral density
                   5 - Coherence
                   6 - Transfer function
       3) history - history string (default = '')
       4) verbose - Verbose (True) or silent (False) (default = False)
    
    Returns None

Tested configurations
---------------------
Tested on Debian GNU/Linux (x86), Red Hat Linux (x86,x86_64) and Mac OS X (PPC).
