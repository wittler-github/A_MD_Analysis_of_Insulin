####################################
## HB VISUALIZATION ENSEMBLE PICS ##
####################################

if {[info exists miMD]} {
  mol delete $miMD
}

label delete Bonds all

# load MD with frames
mol new {../../../../GMXeval/CPGMX/npt_PR.gro} type {gro}
set miMD [lindex [molinfo list] end]
mol addfile {../../../../GMXeval/CenterSuperimposeProteinInBox/MS/123pmdcfMS.trr} type {trr} waitfor all $miMD 
animate delete  beg 0 end 0 skip 0 $miMD
set nf [expr double([molinfo $miMD get numframes])]

set unsortedHBlist [open "../Results/SortCountHB/Sorted_HB_ht_25%_nrHBs_28.dat" r]
set HBlines [read $unsortedHBlist]
close $unsortedHBlist
set HBlines [split $HBlines "\n"]
set sHB {}
set hbcnr 0
foreach line $HBlines {
  if { [regexp {D\(\d+\s+[A-Z]{3}\s+[A-Z]+[0-9]*\)H\(\d+\s+[A-Z]{3}\s+[A-Z]+[0-9]*\)\s+A\(\d+\s+[A-Z]{3}\s+[A-Z]+[0-9]*\)} $line match]} { 
    lappend sHB $line   
    puts "$line"
    incr hbcnr
    puts "HB nr $hbcnr"
  }
}
for {set i 0} {$i < [llength $sHB]} {incr i} {
  regexp {D\(\d+\s+[A-Z]{3}\s+[A-Z]+[0-9]*\)H\((\d+)\s+[A-Z]{3}\s+([A-Z]+[0-9]*)\)\s+A\((\d+)\s+[A-Z]{3}\s+([A-Z]+[0-9]*)\)} [lindex $sHB $i] wm rpname pname raname aname
  set p_index [[atomselect $miMD "protein and resid $rpname and name $pname"] get index]
  set a_index [[atomselect $miMD "protein and resid $raname and name $aname"] get index]
  label add Bonds ${miMD}/$p_index ${miMD}/$a_index
}

animate goto 508

############################
## GRAPHIC REPRESENTATION ##
############################
display projection Orthographic
display rendermode GLSL
display depthcue on 
display resize 1024 1024
scale by 0.833000
scale by 0.833000
axes location Off
color Display Background white
color Labels Bonds purple
label textthickness 0.00001
label textsize 0.0001
color Name C black

set mnr 0
mol addrep $miMD
mol modselect $mnr $miMD "protein and resid 1 to 21"
mol modstyle $mnr $miMD NewCartoon 0.180000 10.000000 4.100000 0 
mol modcolor $mnr $miMD ColorID 4
mol modmaterial $mnr $miMD Transparent

incr mnr
mol addrep $miMD
mol modselect $mnr $miMD "protein and resid 22 to 51"
mol modstyle $mnr $miMD NewCartoon 0.180000 10.000000 4.100000 0 
mol modcolor $mnr $miMD ColorID 7
mol modmaterial $mnr $miMD Transparent

incr mnr
mol addrep $miMD
mol modselect $mnr $miMD "protein"
mol modstyle $mnr $miMD NewCartoon 0.200000 10.000000 4.100000 0 
mol modstyle $mnr $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modmaterial $mnr $miMD Opaque
mol showrep $miMD $mnr 0

incr mnr
mol addrep $miMD
mol modselect $mnr $miMD "protein and (not sidechain)"
mol modstyle $mnr $miMD NewCartoon 0.200000 10.000000 4.100000 0 
mol modstyle $mnr $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modmaterial $mnr $miMD Goodsell

incr mnr
mol addrep $miMD
#mol modselect $mnr $miMD "protein and (not name H1 H2 H3 HA HA1 HN N C O OT1 OT2) and (not hydrogen) and (resname PHE LEU ILE ALA VAL PRO GLY)"
mol modselect $mnr $miMD "protein and (not name H1 H2 H3 HA HA1 HN N C O OT1 OT2) and (resname PHE LEU ILE ALA VAL PRO GLY)"
mol modstyle $mnr $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modcolor $mnr $miMD ColorID 5
mol modmaterial $mnr $miMD Transparent

incr mnr
mol addrep $miMD
#mol modselect $mnr $miMD "protein and (not name H1 H2 H3 HA HA1 C O OT1 OT2) and (not hydrogen) and (resname PRO)"
mol modselect $mnr $miMD "protein and (not name H1 H2 H3 HA HA1 C O OT1 OT2) and (resname PRO)"
mol modstyle $mnr $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modcolor $mnr $miMD ColorID 5
mol modmaterial $mnr $miMD Transparent

incr mnr
mol addrep $miMD
#mol modselect $mnr $miMD "protein and (not name H1 H2 H3 HA HA1 HN N C O OT1 OT2) and (not hydrogen) and (resname GLU ASP)"
mol modselect $mnr $miMD "protein and (not name H1 H2 H3 HA HA1 HN N C O OT1 OT2) and (resname GLU ASP)"
mol modstyle $mnr $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modcolor $mnr $miMD ColorID 1
mol modmaterial $mnr $miMD Transparent

incr mnr
mol addrep $miMD
#mol modselect $mnr $miMD "protein and (not name H1 H2 H3 HA HA1 HN N C O OT1 OT2) and (not hydrogen) and (resname LYS ARG HIS)"
mol modselect $mnr $miMD "protein and (not name H1 H2 H3 HA HA1 HN N C O OT1 OT2) and (resname LYS ARG HIS)"
mol modstyle $mnr $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modcolor $mnr $miMD ColorID 0
mol modmaterial $mnr $miMD Transparent

incr mnr
mol addrep $miMD
#mol modselect $mnr $miMD "protein and (not name H1 H2 H3 HA HA1 HN N C O OT1 OT2) and (not hydrogen) and (resname TYR THR GLN ASN SER CYS)"
mol modselect $mnr $miMD "protein and (not name H1 H2 H3 HA HA1 HN N C O OT1 OT2) and (resname TYR THR GLN ASN SER CYS)"
mol modstyle $mnr $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modcolor $mnr $miMD ColorID 3
mol modmaterial $mnr $miMD Transparent

incr mnr
mol addrep $miMD
mol modselect $mnr $miMD "water"
mol modstyle $mnr $miMD CPK 0.500000 0.200000 12.000000 12.000000
mol modmaterial $mnr $miMD Transparent
mol modcolor $mnr $miMD Name
mol showrep $miMD $mnr 0

incr mnr
mol delrep $mnr $miMD

mol top $miMD

puts "END OF SCRIPT"
###################
## END OF SCRIPT ##
###################
