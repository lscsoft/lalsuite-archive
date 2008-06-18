#!/bin/bash
source /etc/profile
source /usr1/grid/.bash_profile
source /usr1/grid/chrism/suggested_env.txt

/usr1/grid/chrism/GEOpublishstate.py --run A5 --version 4 --logfile /usr1/grid/chrism/GEOpublishstate --database seg_cit --server ldas-cit.ligo.caltech.edu --frametype RDS_C01_L3 --start 835920013 --finder /usr1/grid/chrism/LSCdataFind --fetcher /usr1/grid/chrism/FrStateFetcherSPARC32


