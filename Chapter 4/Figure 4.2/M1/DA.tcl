##########################################
## DIHEDRAL ANGLES PHI/PSI/CHI1,2,3,4,5 ##
##########################################

file delete -force Results
file mkdir Results

file mkdir Results/datToPlotDA
file mkdir Results/PlotDA

if {[info exists miMD]} {
  mol delete $miMD
}

# load MD with frames
mol new {renumberedBDvmdwritten4ins7_M1chainCD.gro} type {gro}
set miMD [lindex [molinfo list] end]
set nf [molinfo $miMD get numframes]
set nfmo [expr $nf -1]
set MSf 0

set FNDAoutfile [open Results/datToPlotDA/DA_all_filenames.dat w] 
set AllresidDAoutfile [open Results/datToPlotDA/AllresidDA.dat w] 

set selCA [atomselect $miMD "protein and alpha"]
set listselCA [$selCA get index] 
for {set i 0} {$i < [llength $listselCA]} {incr i} {

  set ir [expr $i+1]
  set selCAi [atomselect $miMD "index [lindex $listselCA $i]"]     
  set DAi "[lindex [$selCAi get {resid resname}] 0]"
  puts $DAi
  puts $FNDAoutfile "$DAi"
  set DAoutfile [open Results/datToPlotDA/$DAi.dat w]

  if {[string match GLY [$selCAi get {resname}]]} {
    puts "GLY"
    for {set f 0} {$f < $nf} {incr f} {
      $selCAi frame $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\tNaN\tNaN\tNaN\tNaN\tNaN"
      }
    }
  }
  if {[string match ALA [$selCAi get {resname}]]} {
    puts "ALA"
    for {set f 0} {$f < $nf} {incr f} {
      $selCAi frame $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\tNaN\tNaN\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match VAL [$selCAi get {resname}]]} {
    puts "VAL"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG1"]
    set indexCHI1 [$selCHI1 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\tNaN\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match LEU [$selCAi get {resname}]]} {
    puts "LEU"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG CD1"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match ILE [$selCAi get {resname}]]} {
    puts "ILE"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG1"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG1 CD"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match PRO [$selCAi get {resname}]]} {
    puts "PRO"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG CD"]
    set selCHI3 [atomselect $miMD "(protein and resid $ir) and name CB CG CD N"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    set indexCHI3 [$selCHI3 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      set DAchi3($f.r) [measure dihed $indexCHI3]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\t[format "%+ f" $DAchi3($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\t[format "%+ f" $DAchi3($f.r)]\tNaN\tNaN"
      }
    }
  }

  if {[string match PHE [$selCAi get {resname}]]} {
    puts "PHE"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG CD1"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match CYS [$selCAi get {resname}]]} {
    puts "CYS"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB SG"]
    set indexCHI1 [$selCHI1 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\tNaN\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match SER [$selCAi get {resname}]]} {
    puts "SER"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB OG"]
    set indexCHI1 [$selCHI1 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\tNaN\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match THR [$selCAi get {resname}]]} {
    puts "THR"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB OG1"]
    set indexCHI1 [$selCHI1 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\tNaN\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match TYR [$selCAi get {resname}]]} {
    puts "TYR"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG CD1"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match ASN [$selCAi get {resname}]]} {
    puts "ASN"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG OD1"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match GLN [$selCAi get {resname}]]} {
    puts "GLN"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG CD"]
    set selCHI3 [atomselect $miMD "(protein and resid $ir) and name CB CG CD OE1"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    set indexCHI3 [$selCHI3 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      set DAchi3($f.r) [measure dihed $indexCHI3]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\t[format "%+ f" $DAchi3($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\t[format "%+ f" $DAchi3($f.r)]\tNaN\tNaN"
      }
    }
  }

  if {[string match ASP [$selCAi get {resname}]]} {
    puts "ASP"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG OD1"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match GLU [$selCAi get {resname}]]} {
    puts "GLU"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG CD"]
    set selCHI3 [atomselect $miMD "(protein and resid $ir) and name CB CG CD OE1"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    set indexCHI3 [$selCHI3 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      set DAchi3($f.r) [measure dihed $indexCHI3]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\t[format "%+ f" $DAchi3($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\t[format "%+ f" $DAchi3($f.r)]\tNaN\tNaN"
      }
    }
  }

  if {[string match HIS [$selCAi get {resname}]]} {
    puts "HIS"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG ND1"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\tNaN\tNaN\tNaN"
      }
    }
  }

  if {[string match LYS [$selCAi get {resname}]]} {
    puts "LYS"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG CD"]
    set selCHI3 [atomselect $miMD "(protein and resid $ir) and name CB CG CD CE"]
    set selCHI4 [atomselect $miMD "(protein and resid $ir) and name CG CD CE NZ"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    set indexCHI3 [$selCHI3 get index]
    set indexCHI4 [$selCHI4 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      set DAchi3($f.r) [measure dihed $indexCHI3]
      set DAchi4($f.r) [measure dihed $indexCHI4]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\t[format "%+ f" $DAchi3($f.r)]\t[format "%+ f" $DAchi4($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\t[format "%+ f" $DAchi3($f.r)]\t[format "%+ f" $DAchi4($f.r)]\tNaN"
      }
    }
  }

  if {[string match ARG [$selCAi get {resname}]]} {
    puts "ARG"
    set selCHI1 [atomselect $miMD "(protein and resid $ir) and name N CA CB CG"]
    set selCHI2 [atomselect $miMD "(protein and resid $ir) and name CA CB CG CD"]
    set selCHI3 [atomselect $miMD "(protein and resid $ir) and name CB CG CD NE"]
    set selCHI4 [atomselect $miMD "(protein and resid $ir) and name CG CD NE CZ"]
    set selCHI5 [atomselect $miMD "(protein and resid $ir) and name CD NE CZ NH1"]
    set indexCHI1 [$selCHI1 get index]
    set indexCHI2 [$selCHI2 get index]
    set indexCHI3 [$selCHI3 get index]
    set indexCHI4 [$selCHI4 get index]
    set indexCHI5 [$selCHI5 get index]
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) [$selCAi get {phi}]
      set DApsi($f.r) [$selCAi get {psi}]
      set DAchi1($f.r) [measure dihed $indexCHI1]
      set DAchi2($f.r) [measure dihed $indexCHI2]
      set DAchi3($f.r) [measure dihed $indexCHI3]
      set DAchi4($f.r) [measure dihed $indexCHI4]
      set DAchi5($f.r) [measure dihed $indexCHI5]
      puts $DAoutfile "$f\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\t[format "%+ f" $DAchi3($f.r)]\t[format "%+ f" $DAchi4($f.r)]\t[format "%+ f" $DAchi5($f.r)]"
      if {$f==$MSf} {
        puts $AllresidDAoutfile "$ir\t[format "%+ f" $DAphi($f.r)]\t[format "%+ f" $DApsi($f.r)]\t[format "%+ f" $DAchi1($f.r)]\t[format "%+ f" $DAchi2($f.r)]\t[format "%+ f" $DAchi3($f.r)]\t[format "%+ f" $DAchi4($f.r)]\t[format "%+ f" $DAchi5($f.r)]"
      }
    }
  }
  close $DAoutfile
}

close $FNDAoutfile
close $AllresidDAoutfile
puts "Finished Generating Dihedral Angle Data"
###################
## END OF SCRIPT ##
###################
