
# load trajectory
mol load gro trajectory-files/peptso-1a.gro
mol addfile trajectory-files/peptso-1a-100ns-dt1ns.xtc type xtc first 0 last -1 waitfor all

# open a file for writing the data to
set OUTPUT [open "dat/2-protein-depth-vmd.dat" w]

# identify the protein
set protein [atomselect top "protein"]

# identify the bilayer
set bilayer [atomselect top "resname POPC"]

# loop over all frames
for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame 1} {

	# update the atom selections
	$protein frame $frame
	$bilayer frame $frame
	
	# find out the z coordinate of the centre of mass of the protein
	set proteinZ [lindex [measure center $protein weight mass] 2]
	
	# find out the z coordinate of the centre of mass of the bilayer
	set bilayerZ [lindex [measure center $bilayer weight mass] 2]

	# calculate the relative depth of the protein
	set proteinCoord [expr $proteinZ - $bilayerZ]
			
	# so we can see what is happening, write the output to the screen for the moment
	puts "$frame $proteinCoord"

	# write out the frame number and the relative depth of the protein to the file (the bit in the quotes formats the data nicely)
	puts $OUTPUT [format "%7i %7.3f" $frame $proteinCoord]

}
	
# close the file 	
close $OUTPUT
