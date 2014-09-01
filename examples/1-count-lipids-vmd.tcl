
# load the trajectory
mol load gro trajectory-files/peptso-1a.gro
mol addfile trajectory-files/peptso-1a-100ns-dt1ns.xtc type xtc first 0 last -1 waitfor all

# open a file for writing the data to
set OUTPUT [open "dat/1-vmd-lipid-contacts.dat" w]

# iterate through the trajectory, frame by frame
for {set timestep 0} {$timestep < [molinfo top get numframes]} {incr timestep 1} {

	# identify the lipid atoms within 0.36 nm of the protein
	set lipids [atomselect top "resname POPC and within 3.6 of protein"]
	
    # how many atoms is that?
	set lipidsNumber [$lipids num]
		
	# so we can see what is happening, write the output to the screen for the moment
	puts "$timestep $lipidsNumber"

	# write the data for this timestep to the file (the bit in the brackets formats the data nicely)
	puts $OUTPUT "$timestep $lipidsNumber"

}
	
# close the file 	
close $OUTPUT
