###############################
## COUNT SOLUTE NEAR PROTEIN ##
###############################

file delete -force Results/
file mkdir Results
file mkdir Results/datToPlotCSNP

set CSNP_all_filenamesoutfile [open Results/datToPlotCSNP/CSNP_all_filenames.dat w]

if {[info exists miMD]} {
  mol delete $miMD
}

# load MD with frames
mol new {../../GMXeval/CPGMX/npt_PR.gro} type {gro}
set miMD [lindex [molinfo list] end]
mol addfile {../../GMXeval/CenterSuperimposeProteinInBox/MS/123pmdcfMS.trr} type {trr} waitfor all $miMD 
animate delete  beg 0 end 0 skip 0 $miMD
set nfMD [molinfo $miMD get numframes]
set nfMDmo [expr $nfMD-1]
set eqt 9

for {set r 1} {$r < 27} {incr r 1} {
set CSNPoutfile [open "Results/datToPlotCSNP/CSNP_WofP$r.dat" w]
puts $CSNP_all_filenamesoutfile "[format %-50s CSNP_WofP$r]"

  for {set f 0} {$f < $nfMD} {incr f} {
    set selW_O [atomselect top "water and oxygen within $r of protein" frame $f]
    set selNA [atomselect top "name NA and within $r of protein" frame $f]
    set selCL [atomselect top "name CL and within $r of protein" frame $f]
    puts $CSNPoutfile " $f\t[llength [$selW_O get resid]]\t[llength [$selNA get resid]]\t[llength [$selCL get resid]]\t$r"
  }
close $CSNPoutfile 
}
close $CSNP_all_filenamesoutfile


############################
## GRAPHIC REPRESENTATION ##
############################
display projection Orthographic
display rendermode GLSL
color scale method BGR
display depthcue on 
display resize 1024 1024
scale by 0.833000
scale by 0.833000
axes location Off
color Display Background white
color Labels Bonds violet2
color Name C black
#set volmap resolution 
set vmr 0.25

mol addrep $miMD
set mnr 0
mol modselect $mnr $miMD "protein and resid 1 to 21"
mol modstyle $mnr $miMD NewCartoon 0.200000 10.000000 4.100000 0 
mol modcolor $mnr $miMD ColorID 32
mol addrep $miMD
incr mnr
mol modselect $mnr $miMD "protein and resid 22 to 51"
mol modstyle $mnr $miMD NewCartoon 0.200000 10.000000 4.100000 0 
mol modcolor $mnr $miMD ColorID 7

mol addrep $miMD
incr mnr
mol modselect $mnr $miMD "protein"
mol modstyle $mnr $miMD CPK 0.500000 0.200000 10.000000 10.000000
mol modmaterial $mnr $miMD Transparent
mol modcolor $mnr $miMD Name

mol addrep $miMD
incr mnr
mol modselect $mnr $miMD "water"
mol modstyle $mnr $miMD CPK 0.500000 0.200000 10.000000 10.000000
mol modmaterial $mnr $miMD Transparent
mol modcolor $mnr $miMD Name
mol addrep $miMD

mol addrep $miMD
incr mnr
mol modselect $mnr $miMD "(ions and name NA)"
mol modstyle $mnr $miMD CPK 1.300000 0.300000 12.000000 12.000000
mol modcolor $mnr $miMD ColorID 23
mol addrep $miMD

mol addrep $miMD
incr mnr
mol modselect $mnr $miMD "(ions and name CL)"
mol modstyle $mnr $miMD CPK 0.850000 0.300000 12.000000 12.000000
mol modcolor $mnr $miMD ColorID 18
mol addrep $miMD
mol delrep $mnr $miMD

##END if NA and/or Cl##

display resetview
mol top $miMD

###################
## END OF SCRIPT ##
###################
