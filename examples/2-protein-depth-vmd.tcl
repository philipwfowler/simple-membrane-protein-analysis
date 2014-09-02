
# load the trajectory
mol load gro trajectory-files/peptso-1a.gro
mol addfile trajectory-files/peptso-1a-100ns-dt1ns.xtc type xtc first 0 last -1 waitfor all

# open a file for writing the data to
set OUTPUT [open "dat/2-protein-depth-vmd.dat" w]

# identify the protein
set protein [atomselect top "protein"]

# identify the bilayer
set bilayer [atomselect top "resname POPC"]


# iterate through the trajectory, frame by frame
for {set timestep 0} {$timestep < [molinfo top get numframes]} {incr timestep 1} {

	# update the atom selections
	$protein frame $timestep
	$bilayer frame $timestep
	
	set proteinZ [lindex [measure center $protein weight mass] 2]
	set bilayerZ [lindex [measure center $bilayer weight mass] 2]

	set proteinCoord [expr $proteinZ - $bilayerZ]
			
	# so we can see what is happening, write the output to the screen for the moment
	puts "$timestep $proteinCoord"

	# write the data for this timestep to the file (the bit in the brackets formats the data nicely)
	puts $OUTPUT [format "%7i %7.3f" $timestep $proteinCoord]

}
	
# close the file 	
close $OUTPUT
