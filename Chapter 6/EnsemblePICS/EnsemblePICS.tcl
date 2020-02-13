###################
## ENSEMBLE PICS ##
###################

if {[info exists miMD]} {
  mol delete $miMD
}

# load MD with frames
mol new {../GMXeval/CPGMX/npt_PR.gro} type {gro}
set miMD [lindex [molinfo list] end]
mol addfile {../GMXeval/CenterSuperimposeProteinInBox/MS/123pmdcfMS.trr} type {trr} waitfor all $miMD 
animate delete  beg 0 end 0 skip 0 $miMD

############################
## GRAPHIC REPRESENTATION ##
############################
display projection Orthographic
display depthcue off 
display resize 1024 1024
scale by 0.833000
scale by 0.833000
axes location Off
color Display Background white
color Labels Bonds violet2
#set volmap resolution 
set vmr 8

mol addrep $miMD
mol modselect 0 $miMD "(protein and resid 1 to 21)"
mol modstyle 0 $miMD NewCartoon 0.200000 10.000000 4.100000 0 
mol modcolor 0 $miMD ColorID 4
mol modmaterial 0 $miMD Translucent
mol addrep $miMD
mol modselect 1 $miMD "(protein and resid 22 to 51)"
mol modstyle 1 $miMD NewCartoon 0.200000 10.000000 4.100000 0 
mol modcolor 1 $miMD ColorID 7
mol modmaterial 1 $miMD Translucent

mol addrep $miMD
mol modselect 2 $miMD "protein and (not hydrogen) and (not sidechain)"
mol modstyle 2 $miMD NewCartoon 0.200000 10.000000 4.100000 0 
mol modstyle 2 $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modmaterial 2 $miMD Goodsell

mol addrep $miMD
mol modselect 3 $miMD "protein and (not name H1 H2 H3 HA HA1 HN N C O OT1 OT2) and (not hydrogen)"
mol modstyle 3 $miMD CPK 0.500000 0.200000 12.000000 12.000000
color Name C black

mol addrep $miMD
mol modselect 4 $miMD "water within 20 of (protein and resid 1 to 51)"
mol modstyle 4 $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modmaterial 4 $miMD Transparent
mol modcolor 4 $miMD Name
mol addrep $miMD
mol modselect 5 $miMD "water within 20 of (protein and resid 1 to 51)"
volmap occupancy [atomselect $miMD "water within 20 of (protein and resid 1 to 51)"] -res $vmr -allframes -combine avg -mol $miMD
display resetview
mol modstyle 5 $miMD Isosurface 0.450000 0 0 0 1 1
mol modmaterial 5 $miMD Transparent
mol modcolor 5 $miMD ColorID 9
#mol showrep $miMD 5 1

##if no LIG(acetic acid ...)
if {0} {
mol addrep $miMD
mol modselect 6 $miMD "(resname LIG) and within 20 of (protein and resid 1 to 51)"
mol modstyle 6 $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modmaterial 6 $miMD Transparent
mol modcolor 6 $miMD Name
mol addrep $miMD
mol modselect 7 $miMD "(resname LIG) and within 20 of (protein and resid 1 to 51)"
volmap occupancy [atomselect $miMD "(resname LIG) and within 20 of (protein and resid 1 to 51)"] -res $vmr -allframes -combine avg -mol $miMD
display resetview
mol modstyle 7 $miMD Isosurface 0.460000 0 0 0 1 1
mol modmaterial 7 $miMD Transparent
mol modcolor 7 $miMD ColorID 27
}
#endif0 no LIG

##if NA and/or Cl##

mol addrep $miMD
mol modselect 6 $miMD "(ions and name NA) and within 20 of (protein and resid 1 to 51)"
mol modstyle 6 $miMD CPK 1.300000 0.300000 12.000000 12.000000
mol modcolor 6 $miMD ColorID 23
mol addrep $miMD
mol modselect 7 $miMD "(ions and name NA) and within 20 of (protein and resid 1 to 51)"
volmap occupancy [atomselect $miMD "(ions and name NA) and within 20 of (protein and resid 1 to 51)"] -res $vmr -allframes -combine avg -mol $miMD
display resetview
mol modstyle 7 $miMD Isosurface 0.003 1 0 0 1 1
mol modmaterial 7 $miMD Transparent
mol modcolor 7 $miMD ColorID 23
#mol showrep $miMD 7 1

mol addrep $miMD
mol modselect 8 $miMD "(ions and name CL) and within 20 of (protein and resid 1 to 51)"
mol modstyle 8 $miMD CPK 0.850000 0.300000 12.000000 12.000000
mol modcolor 8 $miMD ColorID 18
mol addrep $miMD
mol modselect 9 $miMD "(ions and name CL) and within 20 of (protein and resid 1 to 51)"
volmap occupancy [atomselect $miMD "(ions and name CL) and within 20 of (protein and resid 1 to 51)"] -res $vmr -allframes -combine avg -mol $miMD
display resetview
mol modstyle 9 $miMD Isosurface 0.003 2 0 0 1 1
mol modmaterial 9 $miMD Transparent
mol modcolor 9 $miMD ColorID 18
#mol showrep $miMD 9 1

mol delrep 10 $miMD

##END if NA and/or Cl##

axes location off
mol top $miMD

###################
## END OF SCRIPT ##
###################
