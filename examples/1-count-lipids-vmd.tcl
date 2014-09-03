
# load the trajectory
mol load gro trajectory-files/peptso-1a.gro
mol addfile trajectory-files/peptso-1a-100ns-dt1ns.xtc type xtc first 0 last -1 waitfor all

# open a file for writing the data to
set OUTPUT [open "dat/1-count-lipids-vmd.dat" w]

# loop over all frames
for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame 1} {

	# identify all lipid atoms within 0.36nm of any protein atom
	set lipids [atomselect top "resname POPC and within 3.6 of protein" frame $frame]
	
    # count the lipid atoms
	set lipidsNumber [$lipids num]
		
	# so we can see what is happening, write the output to the screen for the moment
	puts "$frame $lipidsNumber"

	# write out the frame number and number of lipid atoms to file  (the bit in the quotes formats the data nicely)
	puts $OUTPUT [format "%7i %7i" $frame $lipidsNumber]

}
	
# close the file 	
close $OUTPUT
