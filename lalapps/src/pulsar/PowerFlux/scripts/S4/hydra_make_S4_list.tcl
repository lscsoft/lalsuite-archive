#!/usr/bin/env tclsh

set skip_list {}

foreach instrument {L1 H1 H2} {
	set FILENAME "/scratch4/volodya/S4.$instrument.list.ht.txt"
	file delete $FILENAME
	for { set i 1 } { $i < 297 } { incr i } {
		if { [lsearch -exact $skip_list $i] >= 0 } { continue }
		set node [format "%03d" $i]
		catch {
			exec find /netdata/s$node/S4/ -name \*${instrument}_RDS_C01_LX\*.gwf >> $FILENAME 
			} err
		if { $err != "" } { puts "$err" }

                catch {
                        exec find /netdatc/s$node/S4/ -name \*${instrument}_RDS_C01_LX\*.gwf >> $FILENAME
                        } err
                if { $err != "" } { puts "$err" }
		}
	}
