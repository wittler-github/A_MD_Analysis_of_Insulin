###################################################
## SALTBRIDGES, HYDROGENBONDS, RESIDUE DISTANCES ##
###################################################

file delete -force Results
file mkdir Results

if {[info exists miMD]} {
  mol delete $miMD
}
label delete Bonds 0
label delete Bonds all

# equlibration time in frames(1ns)
set eqt 9.0
set eqtf [expr int($eqt)]

# load MD with frames
mol new {../../../GMXeval/CPGMX/npt_PR.gro} type {gro}
set miMD [lindex [molinfo list] end]
mol addfile {../../../GMXeval/CenterSuperimposeProteinInBox/MS/123pmdcfMS.trr} type {trr} waitfor all $miMD 
animate delete  beg 0 end 0 skip 0 $miMD
set nf [expr double([molinfo $miMD get numframes])]
set nfmo [expr $nf-1.0]

set selCA [atomselect $miMD "protein and name CA"]
set listselCA [$selCA get index]   

#################
## SALTBRIDGES ##
#################
saltbr -sel [atomselect $miMD "protein"] -frames $eqtf:last -ondist 3.2 -outdir Results

##############################
## HYDROGEN BONDS DISTANCES ##
##############################
label delete Atoms all
label delete Bonds all

file mkdir Results/datToPlotHB
file mkdir Results/PlotHB

set critHBdis 3.5
set critHBang 30.0

hbonds -sel1 [atomselect $miMD "protein"] -dist $critHBdis -ang $critHBang -polar yes -type unique -frames $eqtf:1:end -upsel yes -writefile yes -outfile Results/datToPlotHB/HBnumber_vs_frame.dat -detailout Results/datToPlotHB/detailedHBs.dat -log LOG.dat -plot no

set unsortedHBlist [open "Results/datToPlotHB/detailedHBs.dat" r]
set HBlines [read $unsortedHBlist]
close $unsortedHBlist

set HBlines [split $HBlines "\n"]
set sHB {}

foreach line $HBlines {
  if { [regexp {\s*[A-Z]{3}\d+[-]+(Main|Side)+[-]+[A-Z]+\d*\s+[A-Z]{3}\d+[-]+(Main|Side)+[-]+[A-Z]+\d*\s+\d+\.\d+} $line match]} { 
    lappend sHB $line   
  }
}

set cldpa {}

for {set i 0} {$i < [llength $sHB]} {incr i} {
  set HB_index {}
  regexp {\d+\.\d+} [lindex $sHB $i 2] pc
  regexp {\s*[A-Z]{3}(\d+)[-]+(Main|Side)+[-]+([A-Z]+\d*)} [lindex $sHB $i 0] d d1 d2 d3
  set d_index [[atomselect $miMD "protein and resid $d1 and name $d3"] get index]
  set p_index [lindex [[atomselect $miMD "protein and resid $d1 and hydrogen within 1.2 of (resid $d1 and name $d3)"] get index] 0]
  regexp {\s*[A-Z]{3}(\d+)[-]+(Main|Side)+[-]+([A-Z]+\d*)} [lindex $sHB $i 1] a a1 a2 a3
  set a_index [[atomselect $miMD "protein and resid $a1 and name $a3"] get index]

  lappend HB_index $pc
  lappend HB_index $d_index
  lappend HB_index $p_index
  lappend HB_index $a_index 
  lappend cldpa $HB_index
}

#Sort HB in order from 0-100% presence
set Scldpa [lsort -real -index 0 $cldpa] 

set HBoutfile [open Results/datToPlotHB/HB_all_filenames.dat w]
set HBpercoutfile [open Results/datToPlotHB/HBperc_all_filenames.dat w]

set ProteinselHB [atomselect $miMD "protein"]
set lProteinselHB [$ProteinselHB get index]
set listHBinoduplicates {}

for {set i [expr [llength $Scldpa]-1]} {$i >= 0} {incr i -1} { 

  set listP {}
  set listA {}

  #set Dli [lindex [lindex $Scldpa $i] 1]
  set Pli [lindex [lindex $Scldpa $i] 2]
  set selHBPli [atomselect $miMD "index $Pli"]
  set Ali [lindex [lindex $Scldpa $i] 3]
  set selHBAli [atomselect $miMD "index $Ali"]

  #If proton P includes no numbers eg. HN
  if {[string is alpha "[$selHBPli get {name}]"] == 1} {
    lappend listP $Pli
  }
  #If proton P includes 2 numbers eg. HH21
  if {([regexp {[[:digit:]]{2}} "[$selHBPli get {name}]" match] == 1)} { 
    if {[regexp {11} "[$selHBPli get {name}]" match] == 1} {
      set Pwc [string trimright "resid [$selHBPli get {resid}] and name [$selHBPli get {name}]" "11"]
    }
    if {[regexp {12} "[$selHBPli get {name}]" match] == 1} {
      set Pwc [string trimright "resid [$selHBPli get {resid}] and name [$selHBPli get {name}]" "12"]
    }
    if {[regexp {13} "[$selHBPli get {name}]" match] == 1} {
      set Pwc [string trimright "resid [$selHBPli get {resid}] and name [$selHBPli get {name}]" "13"]
    }
    if {[regexp {21} "[$selHBPli get {name}]" match] == 1} {
      set Pwc [string trimright "resid [$selHBPli get {resid}] and name [$selHBPli get {name}]" "21"]
    }
    if {[regexp {22} "[$selHBPli get {name}]" match] == 1} {
      set Pwc [string trimright "resid [$selHBPli get {resid}] and name [$selHBPli get {name}]" "22"]
    }
    if {[regexp {23} "[$selHBPli get {name}]" match] == 1} {
      set Pwc [string trimright "resid [$selHBPli get {resid}] and name [$selHBPli get {name}]" "23"]
    }
    append Pwc "11"
    set selPwc [atomselect $miMD "$Pwc"]
    if {-1 != [ lsearch $lProteinselHB [$selPwc get index] ]} {
      lappend listP [$selPwc get index]
    }
    set Pwc [string trimright "$Pwc" "11"]
    append Pwc "12"
    set selPwc [atomselect $miMD "$Pwc"]
    if {-1 != [ lsearch $lProteinselHB [$selPwc get index] ]} {
      lappend listP [$selPwc get index]
    }
    set Pwc [string trimright "$Pwc" "12"]
    append Pwc "13"
    set selPwc [atomselect $miMD "$Pwc"]
    if {-1 != [ lsearch $lProteinselHB [$selPwc get index] ]} {
      lappend listP [$selPwc get index]
    }
    set Pwc [string trimright "$Pwc" "13"]
    append Pwc "21"
    set selPwc [atomselect $miMD "$Pwc"]
    if {-1 != [ lsearch $lProteinselHB [$selPwc get index] ]} {
      lappend listP [$selPwc get index]
    }
    set Pwc [string trimright "$Pwc" "21"]
    append Pwc "22"
    set selPwc [atomselect $miMD "$Pwc"]
    if {-1 != [ lsearch $lProteinselHB [$selPwc get index] ]} {
      lappend listP [$selPwc get index]
    }
    set Pwc [string trimright "$Pwc" "22"]
    append Pwc "23"
    set selPwc [atomselect $miMD "$Pwc"]
    if {-1 != [ lsearch $lProteinselHB [$selPwc get index] ]} {
      lappend listP [$selPwc get index]
    }
  }

  #If proton P includes 1 number eg. H1 (and not 2 numbers)
  if {([regexp {[[:digit:]]{1}} "[$selHBPli get {name}]" match] == 1) && ([regexp {[[:digit:]]{2}} "[$selHBPli get {name}]" match] == 0)} {
    if {[regexp {1} "[$selHBPli get {name}]" match] == 1} {
      set Pwc [string trimright "resid [$selHBPli get {resid}] and name [$selHBPli get {name}]" "1"]        
    }
    if {[regexp {2} "[$selHBPli get {name}]" match] == 1} {
      set Pwc [string trimright "resid [$selHBPli get {resid}] and name [$selHBPli get {name}]" "2"]
    }
    if {[regexp {3} "[$selHBPli get {name}]" match] == 1} {
      set Pwc [string trimright "resid [$selHBPli get {resid}] and name [$selHBPli get {name}]" "3"]
    }
    append Pwc "1"
    set selPwc [atomselect $miMD "$Pwc"]
    if {-1 != [ lsearch $lProteinselHB [$selPwc get index] ]} {
      lappend listP [$selPwc get index]
    }
    set Pwc [string trimright "$Pwc" "1"]
    append Pwc "2"
    set selPwc [atomselect $miMD "$Pwc"]
    if {-1 != [ lsearch $lProteinselHB [$selPwc get index] ]} {
      lappend listP [$selPwc get index]
    }
    set Pwc [string trimright "$Pwc" "2"]
    append Pwc "3"
    set selPwc [atomselect $miMD "$Pwc"]
    if {-1 != [ lsearch $lProteinselHB [$selPwc get index] ]} {
      lappend listP [$selPwc get index]
    }
  }
  #If acceptor atom A has no number e.g. O or N
  if {[string is alpha "[$selHBAli get {name}]"] == 1} {
    lappend listA $Ali
  }
  #If acceptor A has 1 number, eg. OT1
  if {[regexp {[[:digit:]]{1}} "[$selHBAli get {name}]" match] == 1} {
    if {[regexp {1} "[$selHBAli get {name}]" match] == 1} {
      set Awc [string trimright "resid [$selHBAli get {resid}] and name [$selHBAli get {name}]" "1"]
    }
    if {[regexp {2} "[$selHBAli get {name}]" match] == 1} {
      set Awc [string trimright "resid [$selHBAli get {resid}] and name [$selHBAli get {name}]" "2"]
    }
    append Awc "1"
    set selAwc [atomselect $miMD "$Awc"]
    if {-1 != [ lsearch $lProteinselHB [$selAwc get index] ]} {
      lappend listA [$selAwc get index]
    }
    set Awc [string trimright "$Awc" "1"]
    append Awc "2"
    set selAwc [atomselect $miMD "$Awc"]
    if {-1 != [ lsearch $lProteinselHB [$selAwc get index] ]} {
      lappend listA [$selAwc get index]
    }
    set Awc [string trimright "$Awc" "2"]
    append Awc "3"
    set selAwc [atomselect $miMD "$Awc"]
    if {-1 != [ lsearch $lProteinselHB [$selAwc get index] ]} {
      lappend listA [$selAwc get index]
    }
  }

  if {(0 != [llength $listP]) && (0 != [llength $listA])} {
    for {set p 0} {$p < [llength $listP]} {incr p} {
      for {set a 0} {$a < [llength $listA]} {incr a} {
        set Plip [lindex $listP $p]
        set selHBPlip [atomselect $miMD "index $Plip"]
        set Alia [lindex $listA $a]
        set selHBAlia [atomselect $miMD "index $Alia"]
        set Dlip [[atomselect $miMD "index $Plip"] getbonds] 
        set selHBDlip [atomselect $miMD "index $Dlip"]
        set HBi "D([lindex [$selHBDlip get {resid resname name}] 0])H([lindex [$selHBPlip get {resid resname name}] 0]) A([lindex [$selHBAlia get {resid resname name}] 0])"
        if {-1 != [lsearch $listHBinoduplicates "$HBi"]} {
          continue
        }
        lappend listHBinoduplicates "$HBi" 
        label add Bonds ${miMD}/$Plip ${miMD}/$Alia
        puts $HBoutfile "$HBi"
        set HBioutfile [open Results/datToPlotHB/$HBi.dat w]
    	set PercHBDHBAC 0.0
    	set PercHBDC 0.0 
    	set PercHBAC 0.0 
    	for {set f 0} {$f < $nf} {incr f} { 
          $selHBDlip frame $f
          $selHBPlip frame $f
          $selHBAlia frame $f
      	  set cD [lindex [$selHBDlip get {x y z}] 0]  
          set cP [lindex [$selHBPlip get {x y z}] 0]
          set cA [lindex [$selHBAlia get {x y z}] 0]
      	  set HBdis($f.r) [veclength [vecsub $cA $cD]]
      	  set varpsiHNHA($f.r) [expr 180-(180/$M_PI)*acos(([vecdot [vecsub $cA $cP] [vecsub $cD $cP]])/([veclength [vecsub $cA $cP]]*[veclength [vecsub $cD $cP]]))] 
      	  if {(($f>=$eqtf) && ($f<$nf)) && (($HBdis($f.r)<$critHBdis) && ($varpsiHNHA($f.r)<$critHBang))} {
            set PercHBDHBAC [expr $PercHBDHBAC+100.0*(1.0/($nf-$eqtf))]
          }
      	  if {(($f>=$eqtf) && ($f<$nf)) && ($HBdis($f.r)<$critHBdis)} {
            set PercHBDC [expr $PercHBDC+100.0*(1.0/($nf-$eqtf))] 
          }
      	  if {(($f>=$eqtf) && ($f<$nf)) && ($varpsiHNHA($f.r)<$critHBang)} {
            set PercHBAC [expr $PercHBAC+100.0*(1.0/($nf-$eqtf))]
      	  }
          puts $HBioutfile " [format "%-f" $f]\t[format "%-.16f" $HBdis($f.r)]\t[format "%-.16f" $varpsiHNHA($f.r)]"
        } 

        #check format  
        puts $HBpercoutfile "[format %50s $HBi]\t[format "%- f" [lindex [lindex $Scldpa $i] 0]]\t[format "%- f" $PercHBDHBAC]\t[format "%- f" $PercHBDC]\t[format "%- f" $PercHBAC]"
        close $HBioutfile
      }
    }
  }  
}

close $HBoutfile
close $HBpercoutfile

###############################################
## CARBON ALPHA & GEOMETRIC CENTER DISTANCES ##
###############################################
file mkdir Results/datToPlotCAGC
file mkdir Results/PlotCAGC
set CAGCoutfile [open Results/datToPlotCAGC/CAGC_all_filenames.dat w]   
 
for {set i 0} {$i < [llength $listselCA]} {incr i} { 
  set selCAi [atomselect $miMD "protein and (resid [expr $i+1]) and alpha"]
  set selMCi [atomselect $miMD "protein and (resid [expr $i+1]) and not sidechain"]
  set selSCi [atomselect $miMD "protein and (resid [expr $i+1]) and sidechain"]
  set selWCi [atomselect $miMD "protein and (resid [expr $i+1])"]
  for {set j [expr $i+1]} {$j < [llength $listselCA]} {incr j} { 
    set selCAj [atomselect $miMD "protein and (resid [expr $j+1]) and alpha"]
    set selMCj [atomselect $miMD "protein and (resid [expr $j+1]) and not sidechain"]
    set selSCj [atomselect $miMD "protein and (resid [expr $j+1]) and sidechain"]
    set selWCj [atomselect $miMD "protein and (resid [expr $j+1])"]

    set ijCAGC "[lindex [$selCAi get {resid resname}] 0] to [lindex [$selCAj get {resid resname}] 0]"
    puts $CAGCoutfile "$ijCAGC"
    set ijCAGCoutfile [open "Results/datToPlotCAGC/$ijCAGC.dat" w] 

    for {set f 0} {$f < $nf} {incr f} { 
      $selCAi frame $f
      $selCAj frame $f  
      set cCAi [lindex [$selCAi get {x y z}] 0]
      set cCAj [lindex [$selCAj get {x y z}] 0]
      set CAdis($f.r) [veclength [vecsub $cCAi $cCAj]]

      $selMCi frame $f
      $selMCj frame $f   
      set mci {0 0 0}
      set mcj {0 0 0}
      foreach coord [$selMCi get {x y z}] {
        set mci [vecadd $mci $coord]
      }
      set mci [vecscale [expr 1.0/[$selMCi num]] $mci]
      foreach coord [$selMCj get {x y z}] {
        set mcj [vecadd $mcj $coord]
      }
      set mcj [vecscale [expr 1.0/[$selMCj num]] $mcj]
      set MCdis($f.r) [veclength [vecsub $mci $mcj]] 

      $selSCi frame $f
      $selSCj frame $f   
      set sci {0 0 0}
      set scj {0 0 0}
      foreach coord [$selSCi get {x y z}] {
        set sci [vecadd $sci $coord]
      }
      set sci [vecscale [expr 1.0/[$selSCi num]] $sci]
      foreach coord [$selSCj get {x y z}] {
        set scj [vecadd $scj $coord]
      }
      set scj [vecscale [expr 1.0/[$selSCj num]] $scj]
      set SCdis($f.r) [veclength [vecsub $sci $scj]]

      $selWCi frame $f
      $selWCj frame $f   
      set wci {0 0 0}
      set wcj {0 0 0}
      foreach coord [$selWCi get {x y z}] {
        set wci [vecadd $wci $coord]
      }
      set wci [vecscale [expr 1.0/[$selWCi num]] $wci]
      foreach coord [$selWCj get {x y z}] {
        set wcj [vecadd $wcj $coord]
      }
      set wcj [vecscale [expr 1.0/[$selWCj num]] $wcj]
      set WCdis($f.r) [veclength [vecsub $wci $wcj]]

      $selSCi frame $f
      $selMCj frame $f   
      set sci {0 0 0}
      set mcj {0 0 0}
      foreach coord [$selSCi get {x y z}] {
        set sci [vecadd $sci $coord]
      }
      set sci [vecscale [expr 1.0/[$selSCi num]] $sci]
      foreach coord [$selMCj get {x y z}] {
        set mcj [vecadd $mcj $coord]
      }
      set mcj [vecscale [expr 1.0/[$selMCj num]] $mcj]
      set SCMCdis($f.r) [veclength [vecsub $sci $mcj]]

      $selSCi frame $f
      $selCAj frame $f   
      set sci {0 0 0}
      set caj {0 0 0}
      foreach coord [$selSCi get {x y z}] {
        set sci [vecadd $sci $coord]
      }
      set sci [vecscale [expr 1.0/[$selSCi num]] $sci]
      foreach coord [$selCAj get {x y z}] {
        set caj [vecadd $caj $coord]
      }
      set caj [vecscale [expr 1.0/[$selCAj num]] $caj]
      set SCCAdis($f.r) [veclength [vecsub $sci $caj]]

      puts $ijCAGCoutfile " [format "%-f" $f]\t[format "%-.16f" $CAdis($f.r)]\t[format "%-.16f" $MCdis($f.r)]\t[format "%-.16f" $SCdis($f.r)]\t[format "%-.16f" $WCdis($f.r)]\t[format "%-.16f" $SCMCdis($f.r)]\t[format "%-.16f" $SCCAdis($f.r)]"
    }
    close $ijCAGCoutfile    
    $selCAj delete
    $selMCj delete
    $selSCj delete
    $selWCj delete
  }
  $selCAi delete
  $selMCi delete
  $selSCi delete
  $selWCi delete
}
close $CAGCoutfile

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
light 0 on
light 1 off
light 2 on
light 3 off

set mnr 0
mol addrep $miMD
mol modselect $mnr $miMD "protein"
mol modstyle $mnr $miMD NewCartoon 0.180000 10.000000 4.100000 0 
mol modcolor $mnr $miMD ColorID 4
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
mol modselect $mnr $miMD "not protein"
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
