####################################################################
#
# This script creates a graphene sheet : PSF, PDB, TOP and PAR files
#
#

package provide wccnt 0.1

namespace eval ::wccnt:: {
    namespace export wccnt*

    proc wccntGRSbuild { args } {

	# Info
	proc usage {} {
	    vmdcon -info {usage: wccnt GRSbuild [args...]
		-lx        : length in x direction (nm)
		-ly        : length in y direction (nm)
		-type      : <armchair|zigzag> [ default : "armchair" ]
		-cc 	   : bond length (nm) [ default : "0.1418" ]
		-ma	   : <C-C|B-N> carbon/ boron-nitride pairs [ default: C-C ]
		-nlayers   : stacking layers of graphene [ default : 1 ]
		-perXY      : periodic bonds along x-y-axis. 0=no, 1=yes [ default : 0 ]
		-segName   : segname   [ default : GRS ]
		-resName   : resname   [ default : GRS ]
		-atomType  : atom type [ default : CA ]
		-shiftXYZ  : move graphene sheet by a 3D vector  (in nm) [ default : "0 0 0" ]
		-unitGRSnm : length of GRS resid unit (in nm) [ default : 1 ]
		-outName   : output name		
	    }
	    return
	}
	if { [llength $args] < 1 } then { usage; return }


	# Set the defaults
	set perXY      0;
	set type      "armchair";
	set ma        "C-C";
	set cc	      0.1418;
	set nlayers   1;
	set segName   GRS;
	set resName   GRS;
	set atomType  CA;
	set shiftXYZ  "0 0 0";
	set unitGRSnm 1;

	
	# Parse options
	for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {
	    set arg [ lindex $args $argnum ]
	    set val [ lindex $args [expr $argnum + 1]]
	    switch -- $arg {
		"-lx"   { set lx   $val; incr argnum; }
		"-ly"   { set ly   $val; incr argnum; }
		"-type" { set type $val; incr argnum; }
		"-cc"     { set cc   $val; incr argnum; }
		"-ma"     { set ma   $val; incr argnum; }
		"-nlayers"     { set nlayers   $val; incr argnum; }
		"-perXY"     { set perXY   $val; incr argnum; }
		"-segName"  { set segName  $val; incr argnum; }
		"-resName"  { set resName  $val; incr argnum; }
		"-atomType" { set atomType $val; incr argnum; }
		"-shiftXYZ" { set shiftXYZ $val; incr argnum; }
		"-unitGRSnm"    { set unitGRSnm    $val; incr argnum; }
		"-outName"  { set outName  $val; incr argnum; }
		default { error "error: GrapheneBuild: unknown option: $arg" }
	    }
	}
	# check non-default variables    
	set checkLX    [ info exists lx ];
	set checkLY    [ info exists ly ];
	set checkOUTNAME   [ info exists outName ];
	
	if { $checkLX < 1 } {
	    error "error: GRSbuild: need to define variable -lx"
	}    
	if { $checkLY < 1 } {
	    error "error: GRSbuild: need to define variable -ly"
	}
	if { $checkOUTNAME < 1 } {
	    error "error: GRSbuild: need to define variable -outName"
	}

	proc GRSper { lx ly type nlayers perXY shiftXYZ unitGRSnm segName resName atomType outName } {
	    
	    #
	    # This script creates a periodic graphene sheet
	    # the periodicity is in the x- and y-direction
	    #
	    # lx   : length - x for GRS - real
	    # ly   : length - y for GRS - real
	    # type : armchair or zigzag
	    # nlayers  :  number of stacked GRS layers
	    # perXY     :  0 for non periodic system; 1 for periodic system - integer
	    # shiftXYZ :  move GRS by a 3D vector 
	    # unitGRSnm    :  resid unit, lenght in NM - integer
	    # outName  : output name
	    #
	    # test:
	    # -----
	    # set lx 20.;
	    # set ly 20;
	    # set type "armchair";
	    # set nlayers 2;
	    # set perXY 1; 
	    # set shiftXYZ "0 0 0";
	    # set unitGRSnm 0.5;
	    # set outName testper
	    # GRSper  $lx $ly $type $perXY $shiftXYZ $unitGRSnm $outName
	    #
	    
	    
	    ################## MAIN #################

	    # 0.- procedures
	    # ---------------
	    
	    # procedure to produce string from A to ZZZ
	    proc chainName3 { numChain2 } {
		
		# 475253 = ZZZZ
		# 18277 = ZZZ
		# 701 = ZZ 
		# 26 = Z
		# more numbers require new design
		
		proc chainName { numChain } {
		    
		    # template names
		    set listTemplate "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z";
		    set lengthTemplate 26;
		    
		    set chainName [ lindex $listTemplate  [ expr $numChain%26 ] ];
		    
		    set prefixIndex [ expr ($numChain/26) - 1 ];
		    set prefix [ lindex $listTemplate $prefixIndex ];
		    
		    return "$prefix$chainName"	
		}
		
		if { $numChain2 <= 701 } {
		    # two letters
		    chainName $numChain2;	
		} else {
		    # more than two letters
		    set num1 [ expr $numChain2 - 702 ];
		    set perNum [ expr 702 - 26 ];
		    set num2 [ expr $num1/$perNum ];
		    set num3  [ expr $num1 - ( $num2 * $perNum ) + 26 ];
	    
		    set preLeft [ chainName $num2 ];
		    set preRight [ chainName $num3 ];
		    
		    return "$preLeft$preRight";	
		}    
	    }
	    
	    
	    
	    # 1.- create nanotube
	    # --------------------

	    # small nanotube for resid counting
	    package require nanotube
	    graphene -lx $unitGRSnm -ly $unitGRSnm -type $type -nlayers $nlayers ;
	    set molID0 [ molinfo top ];
	    set selAll [ atomselect $molID0 all ]
	    set atomsInRes [ $selAll num ];
	    $selAll delete;
	    mol delete $molID0;
	    
	    # create GRS
	    package require nanotube
	    graphene -lx $lx -ly $ly -type $type -nlayers $nlayers ;  
	    # move GRS
	    set molID1   [ molinfo top ];
	    set selAll   [ atomselect $molID1 all ]
	    set shiftXYZ [ vecscale 10 $shiftXYZ ]
	    $selAll moveby $shiftXYZ;

	    # output PSF/PDB	    
	    animate write psf $outName.NonPer.psf sel $selAll waitfor all $molID1;
	    animate write pdb $outName.NonPer.pdb sel $selAll waitfor all $molID1;
	    
	    # pbc information
	    set aPBC     [ molinfo $molID1 get a ];
	    set bPBC     [ molinfo $molID1 get b ];
	    set cPBC     [ molinfo $molID1 get c ];
	    set alphaPBC [ molinfo $molID1 get alpha ];
	    set betaPBC  [ molinfo $molID1 get beta ];
	    set gammaPBC [ molinfo $molID1 get gamma ];
	    
	    # clean
	    $selAll delete;
	    mol delete $molID1;
	    
	    # 2.- renaming	    
	    # -------------------
  	    
	    # load structure
	    mol new  $outName.NonPer.psf type psf waitfor all;
	    set molID2  [ molinfo top ];
	    mol addfile $outName.NonPer.pdb type pdb molid $molID2 waitfor all;	    
	    
	    # common features
	    set selAll [ atomselect $molID2 all ];
	    $selAll set segname $segName;
	    $selAll set segid   $segName;
	    $selAll set type    $atomType;
	    $selAll set beta 0;
	    $selAll set occupancy 0;
	    
	    
	    # atom names and resids    
	    set listIndex [ $selAll get index ];

	    set i 0;
	    set iName 0;
	    set iRes 1;

	    foreach index $listIndex {
		# atom name		
		if { $iName == $atomsInRes } {
		    set iName 0;
		}
		set newName [ chainName3 $iName ];

		# resid
		set iRes [ expr ($i/$atomsInRes) +1 ]

		# change
		set selOne [ atomselect $molID2 "index $i" ];
		$selOne set name $newName;
		$selOne set resid $iRes;
		$selOne delete;

		incr i;
		incr iName;
	    }


	    # resnames
	    set residList [ lsort -unique -increasing -integer [ $selAll get resid ] ];
	    
	    # number of atoms in resids
	    set numResidList "";
	    foreach resid $residList {
		set selOneRes [ atomselect $molID2 "resid $resid" ];
		lappend numResidList [ $selOneRes num ];
		$selOneRes delete;
	    }	    	    	    
	    set numFirstResid [ lindex $numResidList 0 ];

	    # add X if there are two different residues
	    set xChar "X";
	    set resNameEnd "$resName$xChar"

	    foreach numResid $numResidList resid $residList {
		set selOneRes [ atomselect $molID2 "resid $resid" ];
		if { $numResid == $numFirstResid } {		    
		    $selOneRes set resname $resName;
		} else {
		    $selOneRes set resname $resNameEnd;
		}
		$selOneRes delete;
	    }


	    # output molecule
	    animate write psf $outName.NonPer.psf sel $selAll waitfor all $molID2;
	    animate write pdb $outName.NonPer.pdb sel $selAll waitfor all $molID2;
	    $selAll delete;
	    mol delete $molID2;
	    unset listIndex;
	    	    
	    

	    # 3.- remove improper terms
	    # --------------------------
	    
	    # load molecule
	    mol new $outName.NonPer.psf type psf waitfor all;
	    set molID3   [ molinfo top ];
	    mol addfile  $outName.NonPer.pdb type pdb molid $molID3 waitfor all;
	    set selAll   [ atomselect $molID3 all ];
	    
	    # remove improper terms
	    package require topotools
	    topo -molid $molID3 -sel $selAll clearimpropers
    
	    # output PSF/PDB
	    animate write psf $outName.NonImpr.psf sel $selAll waitfor all $molID3;
	    animate write pdb $outName.NonImpr.pdb sel $selAll waitfor all $molID3;	    
	    
	    # clean
	    $selAll delete;
	    mol delete $molID3;
	    file delete -force $outName.NonPer.psf;
	    file delete -force $outName.NonPer.pdb;
	    	    


	    # 4.- periodic bonds
	    # -------------------
	    
	    if { $perXY > 0 } {
		
		# load molecule
		mol new $outName.NonImpr.psf type psf waitfor all;
		set molID4 [ molinfo top ];
		mol addfile  $outName.NonImpr.pdb type pdb molid $molID4 waitfor all;
		
		# up and down rings
		set selAll [ atomselect $molID4 all ];
		foreach { cenX cenY cenZ } [ measure center $selAll ] { break };
		$selAll delete;		
		set selRingY   [ atomselect $molID4 "(numbonds == 2) and (y < $cenY)" ];
		set selRingYUp   [ atomselect $molID4 "(numbonds == 2) and (y > $cenY)" ];
		set selRingX [ atomselect $molID4 "(numbonds == 2) and (x < $cenX)" ];
		set selRingXUp [ atomselect $molID4 "(numbonds == 2) and (x > $cenX)" ];

		# move down ring close to up ring
		set xPer     [ molinfo $molID4 get a ];
		set yPer     [ molinfo $molID4 get b ];
		set moveX   "$xPer 0 0";
		set moveXmin [ vecscale -1 $moveX ];
		$selRingX moveby $moveX;
		
		# list of bonds between up and down rings
		set indexRingX [ $selRingX get index ];
		set indexRingXUp [ $selRingXUp get index ];
		
		set cutoffC 1.418;
		set perBondsX "";
		
		foreach index $indexRingX {
		    set nearSel [ atomselect $molID4 "(all within $cutoffC of index $index ) and (index $indexRingXUp)" ];
		    set nearIndex [ $nearSel get index ];
		    
		    foreach index2 $nearIndex {
			lappend perBondsX "$index $index2";
		    }		    
		    $nearSel delete;
		}
		
		# move down ring back
		$selRingX moveby $moveXmin;
		
		# add bonds
		foreach eachPair $perBondsX {
		    foreach { indexLeft indexRight } $eachPair { break };
		    topo -molid $molID4 addbond $indexLeft $indexRight
		}

		set moveY   "0 $yPer 0";
		set moveYmin [ vecscale -1 $moveY ];
		$selRingY moveby $moveY;
		
		set indexRingY [ $selRingY get index ];
		set indexRingYUp [ $selRingYUp get index ];
		set perBondsY "";

		foreach index $indexRingY {
		    set nearSel [ atomselect $molID4 "(all within $cutoffC of index $index ) and (index $indexRingYUp)" ];
		    set nearIndex [ $nearSel get index ];
		    
		    foreach index2 $nearIndex {
			lappend perBondsY "$index $index2";
		    }		    
		    $nearSel delete;
		}

		# move down ring back
		$selRingY moveby $moveYmin;

		# add bonds
		foreach eachPair $perBondsY {
		    foreach { indexLeft indexRight } $eachPair { break };
		    topo -molid $molID4 addbond $indexLeft $indexRight
		}		

		# output PSF/PDB
		set selAll [ atomselect $molID4 all ];
		animate write psf $outName.Per.psf sel $selAll waitfor all $molID4;
		animate write pdb $outName.Per.pdb sel $selAll waitfor all $molID4;
		
		# clean
		$selAll delete;
		$selRingX delete;
		$selRingXUp delete;
		$selRingY delete;
		$selRingYUp delete;
		unset perBondsX;
		unset perBondsY;
		mol delete $molID4;
		file delete -force $outName.NonImpr.psf;
		file delete -force $outName.NonImpr.pdb;
		
	    } else {	
		file rename -force $outName.NonImpr.psf $outName.Per.psf;
		file rename -force $outName.NonImpr.pdb $outName.Per.pdb;
	    }
    
	    
	    
	    # 5.- regenerate angles/dihedrals
	    # --------------------------------
	    package require psfgen
	    resetpsf
	    psfcontext reset;
	    readpsf  $outName.Per.psf;
	    coordpdb $outName.Per.pdb;    
	    regenerate angles dihedrals    
	    writepsf $outName.psf
	    writepdb $outName.pdb
	    
	    # clean
	    resetpsf
	    psfcontext reset;
	    file delete -force $outName.Per.psf;
	    file delete -force $outName.Per.pdb;
	    
	    
    
	    # 6.- add PBC info back
	    # -----------------------
	    
	    # load molecule
	    mol new $outName.psf type psf waitfor all;
	    set molID4 [ molinfo top ];
	    mol addfile  $outName.pdb type pdb molid $molID4 waitfor all;
	    
	    # set values
	    molinfo $molID4 set a $aPBC;
	    molinfo $molID4 set b $bPBC;
	    molinfo $molID4 set c $cPBC;
	    molinfo $molID4 set alpha $alphaPBC;
	    molinfo $molID4 set beta  $betaPBC;
	    molinfo $molID4 set gamma $gammaPBC;
	    
	    # output molecule
	    set selAll [ atomselect $molID4 all ];
	    animate write pdb $outName.pdb sel $selAll waitfor all $molID4;

	    # clean
	    $selAll delete;	    
	    mol delete $molID4
	    
	}
	
	
	
	proc GRStoppar { psfFile pdbFile outName } {

	    #
	    # This script generates a GRS topology
	    # from a PSF/PDB structure to be used in CHARMM2LAMMPS
	    #
	    # psfFile : PSF file
	    # pdbFile : PDB file
	    # outName : output name
	    #
	    # test:
	    # -----
	    # set psfFile test04per.psf
	    # set pdbFile test04per.pdb
	    # set outName test10;
	    # GRStop $psfFile $pdbFile $outName;
	    #
	    
	    
	    ############## MAIN ################
	    
	    # 1.- molecule info
	    # --------------------
	    
	    # load structure
	    mol new $psfFile type psf waitfor all;
	    set molID1  [ molinfo top ];
	    mol addfile $pdbFile type pdb molid $molID1 waitfor all;	    

	    # PBC info
	    set aPBC     [ molinfo $molID1 get a ];
	    set bPBC     [ molinfo $molID1 get b ];
	    set cPBC     [ molinfo $molID1 get c ];
	    set alphaPBC [ molinfo $molID1 get alpha ];
	    set betaPBC  [ molinfo $molID1 get beta ];
	    set gammaPBC [ molinfo $molID1 get gamma ];
	    
	    # segname info
	    # currently, the script only works with one segment, one atom type
	    set selAll  [ atomselect $molID1 all ];
	    set segName [ $selAll get segname ];
	    set segName [ lsort -unique $segName ];
	    set segName [ lindex $segName 0 ];   
	    set atomType [ $selAll get type ];
	    set atomType [ lsort -unique $atomType ];
	    set atomType [ lindex $atomType 0 ];
	    $selAll delete;

	    # >>>>>>>>>>>>>>>>

	    # resname/resid info	    
	    set selAll       [ atomselect $molID1 all ];
	    set listResname  [ lsort -unique [ $selAll get resname ] ];
	    set lListResname [ llength $listResname ];
	    set listResid    [ lsort -unique [ $selAll get resid ] ];
	    set listResid    [ lsort -unique -integer -increasing $listResid ];
	    $selAll delete;

	    # first resid
	    set residFirst [ lindex $listResid 0 ];
	    set selFirst   [ atomselect $molID1 "resid $residFirst" ];
	    set indexFirst [ $selFirst get index ];
	    set resnFirst  [ $selFirst get resname ];
	    set resnFirst  [ lsort -unique $resnFirst ];
	    $selFirst delete;

	    set listTopIndex "";
	    lappend listTopIndex $indexFirst;
	    set listTopName "";
	    lappend listTopName $resnFirst;
	    
	    # last resid
	    if { $lListResname > 1 } { 		
		set residLast [ lindex $listResid end ];
		set selLast   [ atomselect $molID1 "resid $residLast" ];
		set indexLast [ $selLast get index ];
		set resnLast  [ $selLast get resname ];
		set resnLast  [ lsort -unique $resnLast ];
		$selLast delete;

		lappend listTopIndex $indexLast;
		lappend listTopName  $resnLast;
	    }
	    # >>>>>>>>>>>>>>>>


	    # 2.- feed topology
	    # -------------------
	    
	    # topology file
	    set outTOP [ open $outName.top w ];
	    
	    # mass
	    puts $outTOP "31  1";
	    puts $outTOP " ";
	    puts $outTOP "MASS     1 $atomType    12.01100 C";
	    puts $outTOP " ";
	    puts $outTOP "AUTO ANGLES DIHE";
	    puts $outTOP " ";

	    # residues
	    foreach currIndex $listTopIndex currResName $listTopName {
		
		# header
		puts $outTOP "RESI $currResName          0.00";
				
		# atoms
		foreach index $currIndex {
		    set  selOne  [ atomselect $molID1 "index $index" ];
		    set  nameOne [ $selOne get name ];
		    puts $outTOP "ATOM $nameOne $atomType      0.00";
		    $selOne delete;
		}
				
		# create bond list
		set bondList "";
		
		foreach index $currIndex {

		    # find neigh around each index
		    set selOne     [ atomselect $molID1 "index $index" ];
		    set neighIndex [ $selOne getbonds ];
		    set neighIndex [ lindex $neighIndex 0 ];
		    $selOne delete;

		    set selOne     [ atomselect $molID1 "(index $neighIndex) and (index $currIndex)" ];
		    set neighIndex [ $selOne get index ];
		    $selOne delete;
		
		    # create a list of all indexes
		    foreach indexB $neighIndex { 
			set selA  [ atomselect $molID1 "index $index" ];
			set selB  [ atomselect $molID1 "index $indexB" ];
			set nameA [ $selA get name ];
			set nameB [ $selB get name ];
			
			if { $index < $indexB } {
			    lappend bondList "$nameA $nameB";
			} elseif { $index > $indexB } {
			    lappend bondList "$nameB $nameA";
			} else {
			}
			
			$selA delete;
			$selB delete;
		    }
		}
		
		# sort bond list
		set bondList [ lsort -unique $bondList ];
		
		# output bonds
		foreach valPAIRunique $bondList {
		    set valA [  lindex $valPAIRunique 0 ];
		    set valB [  lindex $valPAIRunique 1 ];
		    puts $outTOP "BOND $valA $valB";
		}
		
		puts $outTOP " ";
		puts $outTOP " ";

		unset bondList;
	    }
	    	    
	    # close topology file;
	    puts  $outTOP " ";
	    close $outTOP;
	    	    	    
	    # clear
	    mol delete $molID1;
	    
	    
	    
	    # 3.- create molecule from topology file
	    # ---------------------------------------    
	    package require psfgen
	    resetpsf
	    psfcontext reset;
	    topology $outName.top;
	    segment  $segName { pdb $pdbFile } 
	    coordpdb $pdbFile $segName;
	    regenerate angles dihedrals;
	    writepsf $outName.TOP.psf;
	    writepdb $outName.TOP.pdb;
	    resetpsf;
	    psfcontext reset;
	    

	    
	    # 4.- add PBC info back
	    # -----------------------
	    
	    # load molecule
	    mol new $outName.TOP.psf type psf waitfor all;
	    set molID2 [ molinfo top ];
	    mol addfile  $outName.TOP.pdb type pdb molid $molID2 waitfor all;
	    
	    # set PBC values
	    molinfo $molID2 set a $aPBC;
	    molinfo $molID2 set b $bPBC;
	    molinfo $molID2 set c $cPBC;
	    molinfo $molID2 set alpha $alphaPBC;
	    molinfo $molID2 set beta  $betaPBC;
	    molinfo $molID2 set gamma $gammaPBC;
    
	    # output molecule
	    set selAll [ atomselect $molID2 all ];
	    animate write pdb $outName.TOP.pdb sel $selAll waitfor all $molID2;
	    
	    # clean
	    mol delete $molID2
	    $selAll delete;


	    
	    # 5.- create parameter file
	    # ---------------------------

	    set outPAR [ open $outName.par w ];
	    	    
	    puts $outPAR "* This is a reduced version of the CHARMM22 parameter file for GRS";
	    puts $outPAR "*";
	    puts $outPAR " ";
	    puts $outPAR "BONDS";
	    puts $outPAR "!V(bond) = Kb(b - b0)**2 ";
	    puts $outPAR "!";
	    puts $outPAR "!atom type Kb          b0";
	    puts $outPAR "$atomType $atomType 305.000  1.3750 ";
	    puts $outPAR " ";
	    puts $outPAR " ";
	    puts $outPAR "ANGLES";
	    puts $outPAR "!V(angle) = Ktheta(Theta - Theta0)**2";
	    puts $outPAR "!V(Urey-Bradley) = Kub(S - S0)**2 ";
	    puts $outPAR "! ";
	    puts $outPAR "!atom types     Ktheta    Theta0   Kub     S0 ";
	    puts $outPAR "$atomType $atomType $atomType  40.000    120.00   35.00   2.41620";
	    puts $outPAR " ";
	    puts $outPAR " ";
	    puts $outPAR "DIHEDRALS";
	    puts $outPAR "!V(dihedral) = Kchi(1 + cos(n(chi) - delta)) ";
	    puts $outPAR "! ";
	    puts $outPAR "!atom types             Kchi    n   delta ";
	    puts $outPAR "$atomType $atomType $atomType $atomType  3.1000  2   180.00 ";
	    puts $outPAR " ";
	    puts $outPAR " ";
	    puts $outPAR "NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch - ";
	    puts $outPAR "cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5  ";
	    puts $outPAR "!V(Lennard-Jones) = Eps,i,j\[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6\] ";
	    puts $outPAR "! ";
	    puts $outPAR "!atom  ignored    epsilon      Rmin/2 ";
	    puts $outPAR "$atomType   0.000000  -0.0565     1.8178 ";
	    puts $outPAR " ";
	    puts $outPAR " ";
	    puts $outPAR "END";
	    puts $outPAR " ";

	    close $outPAR;
	    
	}
	
	
	
	# ---------------
	# run procedures
	# ---------------
	
	# create GRS
	GRSper $lx $ly $type $nlayers $perXY $shiftXYZ $unitGRSnm $segName $resName $atomType $outName

	# create topology
	GRStoppar $outName.psf $outName.pdb $outName;

	# clean
	file delete -force $outName.TOP.psf;
	file delete -force $outName.TOP.pdb;
		        
    }   
}
