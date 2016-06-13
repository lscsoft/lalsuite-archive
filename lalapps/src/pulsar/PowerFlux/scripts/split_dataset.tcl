#!/usr/bin/env tclsh

foreach {var value} {
	SFT_ROOT	"/archive/home/volodya/SFT-3/"
	DETECTOR	"L1"
	DATASET		"S5.${DETECTOR}.ht.3.geo"
	NODES_START	1
	NODES_STOP	300
	ACTION	"help"
	} {
	global $var
	set $var [subst -nocommands $value]
	}

foreach {var value} $argv {
	global $var
	set $var $value
	}

set NODE_LIST {}

proc assign_nodes {} {
global NODES_START NODES_STOP NODE NGROUPS ALLGROUPS
set k 0
set ALLGROUPS {}
for { set i $NODES_START } { $i < $NODES_STOP } { incr i } {
	catch {
		set dir "/data/node$i/volodya/"
		file mkdir $dir
		set NODE($k) $dir
		lappend ALLGROUPS $k
		incr k
		puts stderr $dir
		}
	}

set NGROUPS $k
}


proc output_dataset {} {
global DATASET NGROUPS NODE DETECTOR

set D [open "${DATASET}.dst" "w"]
puts $D "new_dataset \"$DATASET\""
switch $DETECTOR {
	H1  { puts $D "detector  \"LHO\"" }
	H2  { puts $D "detector  \"LHO\"" }
	L1  { puts $D "detector  \"LLO\"" }
	}
for {set i 0 } { $i < $NGROUPS } { incr i } {
	set dir $NODE($i)
	puts $D "lock_file \"$dir/read.lock\""
	puts $D "directory \"$dir/$DATASET\""
	}
close $D
}

proc assign_groups {} {
global NGROUPS GROUP FILES SFT_ROOT DATASET
set FILES [glob "$SFT_ROOT/$DATASET/*"]

set STEP [expr ceil([llength $FILES]/$NGROUPS)]

puts "STEP=$STEP"

set k 0
foreach file $FILES {
	set GROUP($file) $k
	puts "$file -> $k"

	incr k	
	if { $k >= $NGROUPS } { set k 0 }
	}
}

proc remove_groups { groups } {
global NODE DATASET
foreach group $groups {
	set dir "$NODE($group)/$DATASET/"
	puts stderr "Cleaning $dir"
	file delete -force -- $dir
	}
}

proc create_dirs { groups } {
global NODE DATASET
foreach group $groups {
	set dir "$NODE($group)/$DATASET/"
	puts stderr "Creating $dir"
	while { [ catch { file mkdir $dir } err ] } {
		puts stderr "ERROR creating $dir : $err, retrying."
		after 1000
		}
	}
}

proc copy_groups { groups } {
global NODE GROUP FILES DATASET

remove_groups $groups
create_dirs $groups

foreach file $FILES {
	set g $GROUP($file)
	if { [lsearch -integer $groups $g] < 0 } {
		continue
		}
	set dir "$NODE($GROUP($file))/$DATASET/"
	puts "$file -> $dir"
	while { [catch {file copy $file "$dir/[file tail $file]" } err ]  } {
		puts stderr "ERROR copying $file -> $dir: $err, retrying."
		after 1000
		}
	}
}


proc save_status { filename } {
set F [open $filename "w"]

foreach var {NGROUPS DATASET DETECTOR } {
	global $var
	puts $F "global $var"
	puts $F "set $var {[set $var]}"
	}

foreach arr {NODE GROUP} {
	global $arr
	puts $F "global $arr"
	foreach {key val} [array get $arr] {
		puts -nonewline $F "set ${arr}"
		puts $F "($key) {$val}"
		}
	}

foreach L {ALLGROUPS FILES} {
	global $L
	puts $F "global $L"
	puts $F "set $L {{[join [set $L] "} {"]}}"
	}
puts $F "set LOAD_COMPLETE 1"
close $F
}

proc prepare_plan {} {
global DATASET
assign_nodes
assign_groups
save_status "${DATASET}.plan.tcl"
}

proc execute_plan {} {
global NGROUPS ALLGROUPS DATASET
source "${DATASET}.plan.tcl"
puts "NGROUPS=$NGROUPS"
output_dataset
copy_groups $ALLGROUPS
}

proc help {} {
puts "usage: split_dataset.tcl VAR1 VALUE1 ..."
puts "\tACTION:  help, prepare_plan, execute_plan"
}

foreach {var value} $argv {
	global $var
	set $var $value
	}

$ACTION