############################################################
#
# This script converts LAMMPS trj into DCD
#
# usage:
# set trjFile SMTH.lammpstrj;
# set outName SMTH;
# trj2dcd $trjFile $outName
#

package provide wccnt 0.1


namespace eval ::wccnt:: {
    namespace export wccnt*

    proc wccnttrj2dcd { args } {
	
	# Info
	proc usage {} {
	    vmdcon -info {usage: wccnt trj2dcd [args...]
		-trjFile  : trajectory in LAMMPS format
		-outName : output name		
	    }
	    return
	}
	if { [llength $args] < 1 } then { usage; return }
	
	
	# Parse options
	for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {
	    set arg [ lindex $args $argnum ]
	    set val [ lindex $args [expr $argnum + 1]]
	    switch -- $arg {
		"-trjFile"  { set trjFile  $val; incr argnum; }
		"-outName"  { set outName  $val; incr argnum; }
		default { error "error: trj2dcd: unknown option: $arg" }
	    }
	}

	
	# check non-default variables    
	set checkTRJFILE   [ info exists trjFile ];
	set checkOUTNAME   [ info exists outName ];
	
	if { $checkTRJFILE < 1 } {
	    error "error: trj2dcd: need to define variable -trjFile"
	}
	if { $checkOUTNAME < 1 } {
	    error "error: trj2dcd: need to define variable -outName"
	}


	########### MAIN ############

	# load molecule
	mol new $trjFile type lammpstrj first 0 last -1 step 1 waitfor all;
	set molID  [ molinfo top ];
	
	# export dcd file
	set selAll [ atomselect $molID all ];
	animate write dcd $outName.dcd beg 0 end -1 sel $selAll waitfor all $molID;
	
	# clean
	$selAll delete;		
	mol delete $molID;    	
    }
}

