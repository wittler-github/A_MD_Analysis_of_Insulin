######################################################
## DA PHI/PSI/CHI1,2,3,4,5 COMPARISON TO EXP BOUNDS ##
######################################################

file delete -force Results
file mkdir Results

file mkdir Results/datToPlotDA
file mkdir Results/PlotDA

if {[info exists miMD]} {
  mol delete $miMD
}

# load MD with frames
mol new {../renumberedBDvmdwritten4ins7_M1chainCD.gro} type {gro}
set miMD [lindex [molinfo list] end]

set nf [molinfo $miMD get numframes]
set nfmo [expr $nf -1]
set MSf 0

set eqt 0.0
set eqtf [expr int($eqt)]

#Open (PDB entry).mr file with DA restraints sort them into lists
set nmrdotmr [open 2kjj.mr]
set read_nmrdotmr [read $nmrdotmr]
close $nmrdotmr
set mExpDA {}; 
foreach {rep hej} [regexp -all -inline {(\s*assign\s*\(\s*resid\s*\d*\s*and\s*name\s*[A-Z]+\d*\s*\)\s*\n\s*\(\s*resid\s*\d*\s*and\s*name\s*[A-Z]+\d*\s*\)\s*\n\s*\(\s*resid\s*\d*\s*and\s*name\s*[A-Z]+\d*\s*\)\s*\n\s*\(\s*resid\s*\d*\s*and\s*name\s*[A-Z]+\d*\s*\)\s*\d+.\d+\s*[-+]*\d+\s*\d+.\d+\s*\d+)} $read_nmrdotmr]  {
  lappend mExpDA $rep
}

set list_mExpDAid {}
set expDAnr 0

foreach DAid $mExpDA {
  regexp {\s*assign\s*\(\s*(resid\s*\d*\s*and\s*name\s*[A-Z]+\d*)\s*\)\s*\n\s*\(\s*(resid\s*\d*\s*and\s*name\s*[A-Z]+\d*)\s*\)\s*\n\s*\(\s*(resid\s*\d*\s*and\s*name\s*[A-Z]+\d*)\s*\)\s*\n\s*\(\s*(resid\s*\d*\s*and\s*name\s*[A-Z]+\d*)\s*\)\s*\d+.\d+\s*([-+]*\d+)\s*(\d+.\d+)\s*\d+} $DAid Aall_id A1_id A2_id A3_id A4_id angle bound
  set tempExpDAid [concat "{$A1_id}" "{$A2_id}" "{$A3_id}" "{$A4_id}" "{$angle}" "{$bound}"]	
  #puts "$tempExpDAid"
  lappend list_mExpDAid $tempExpDAid
  incr expDAnr
}

for {set i_eDA 0} {$i_eDA < [llength $list_mExpDAid]} {incr i_eDA} {
  #puts "[lindex $list_mExpDAid $i_eDA 0] [lindex $list_mExpDAid $i_eDA 1] [lindex $list_mExpDAid $i_eDA 2] [lindex $list_mExpDAid $i_eDA 3] [lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5]"
}

##Calculate DA and compare with exp DA

set FNDAoutfile [open Results/datToPlotDA/E_DA_all_filenames.dat w] 
set AllresidDAoutfile [open Results/datToPlotDA/E_AllresidDA.dat w] 
set Cid_expDAoutfile [open Results/datToPlotDA/Cid_expDA.dat w]

set selCA [atomselect $miMD "protein and alpha"]
set listselCA [$selCA get index] 
#set Vda 0
for {set i 0} {$i < [llength $listselCA]} {incr i} {
  set ir [expr $i+1]
  set irpo [expr $ir+1]
  set selCAi [atomselect $miMD "index [lindex $listselCA $i]"]     
  set DAi "[lindex [$selCAi get {resid resname}] 0]"
  #puts $DAi
  puts $FNDAoutfile $DAi
  set DAoutfile [open Results/datToPlotDA/$DAi.dat w]

  set phi_LB "NOTexp"  
  set phi_UB "NOTexp"
  set psi_LB "NOTexp"  
  set psi_UB "NOTexp"
  set chi1_LB "NOTexp"  
  set chi1_UB "NOTexp"
  set chi2_LB "NOTexp"  
  set chi2_UB "NOTexp"
  set chi3_LB "NOTexp"  
  set chi3_UB "NOTexp"
  set chi4_LB "NOTexp"  
  set chi4_UB "NOTexp"
  set chi5_LB "NOTexp"  
  set chi5_UB "NOTexp"

  set phi_lAi "NOTEXP"
  set psi_lAi "NOTEXP"
  set chi1_lAi "NOTEXP"
  set chi2_lAi "NOTEXP"
  set chi3_lAi "NOTEXP"
  set chi4_lAi "NOTEXP"
  set chi5_lAi "NOTEXP"

  set t_phi 'NAN'
  set t_psi 'NAN'
  set t_chi1 'NAN'
  set t_chi2 'NAN'
  set t_chi3 'NAN'
  set t_chi4 'NAN'
  set t_chi5 'NAN'

  #set V_loub 0
  puts $Cid_expDAoutfile "\n\n## rID $ir LOOP EXP DA ##\n"  

  for {set i_eDA 0} {$i_eDA < [llength $list_mExpDAid]} {incr i_eDA} {
    if {("resid $i and name C"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name C"==[lindex $list_mExpDAid $i_eDA 3])} {
      set A1_c [atomselect $miMD "resid $i and name C"]
      set A2_c [atomselect $miMD "resid $ir and name N"]
      set A3_c [atomselect $miMD "resid $ir and name CA"]
      set A4_c [atomselect $miMD "resid $ir and name C"]
      set phi_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
      set t_phi [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
      if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
        set phi_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
        set phi_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
      } else {
        set phi_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
        set phi_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
      }
      puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $phi_LB $phi_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_phi\n"
    }
    if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name C"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $irpo and name N"==[lindex $list_mExpDAid $i_eDA 3])} {
      set A1_c [atomselect $miMD "resid $ir and name N"]
      set A2_c [atomselect $miMD "resid $ir and name CA"]
      set A3_c [atomselect $miMD "resid $ir and name C"]
      set A4_c [atomselect $miMD "resid $irpo and name N"]
      set psi_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
      set t_psi [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
      if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
        set psi_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
        set psi_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
      } else {
        set psi_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
        set psi_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
      }
      puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $psi_LB $psi_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_psi\n"
    }
    if {[string match VAL [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG1"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG1"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
    }
    if {[string match ILE [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG1"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG1"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG1"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG1"]
        set A4_c [atomselect $miMD "resid $ir and name CD"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
    }
    if {[string match LEU [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CD1"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name CD1"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
    }
    if {[string match PRO [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name CD"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
      if {("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CB"]
        set A2_c [atomselect $miMD "resid $ir and name CG"]
        set A3_c [atomselect $miMD "resid $ir and name CD"]
        set A4_c [atomselect $miMD "resid $ir and name N"]
        set chi3_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi3 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi3_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi3_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi3_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi3_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi3_LB $chi3_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi3\n"
      }
    }
    if {[string match PHE [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CD1"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name CD1"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
    }
    if {[string match CYS [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name SG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name SG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
    }
    if {[string match SER [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name OG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name OG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
    }
    if {[string match THR [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name OG1"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name OG1"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
    }
    if {[string match TYR [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CD1"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name CD1"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
    }
    if {[string match ASN [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name OD1"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name OD1"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
    }
    if {[string match GLN [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name CD"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
      if {("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name OE1"==[lindex $list_mExpDAid $i_eDA 3])} { 
        set A1_c [atomselect $miMD "resid $ir and name CB"]
        set A2_c [atomselect $miMD "resid $ir and name CG"]
        set A3_c [atomselect $miMD "resid $ir and name CD"]
        set A4_c [atomselect $miMD "resid $ir and name OE1"]
        set chi3_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi3 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi3_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi3_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi3_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi3_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi3_LB $chi3_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi3\n"
      }
    }
    if {[string match ASP [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name OD1"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name OD1"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
    }
    if {[string match GLU [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name CD"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
      if {("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name OE1"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CB"]
        set A2_c [atomselect $miMD "resid $ir and name CG"]
        set A3_c [atomselect $miMD "resid $ir and name CD"]
        set A4_c [atomselect $miMD "resid $ir and name OE1"]
        set chi3_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi3 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi3_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi3_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi3_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi3_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi3_LB $chi3_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi3\n"
      }
    }
    if {[string match HIS [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} { 
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name ND1"==[lindex $list_mExpDAid $i_eDA 3])} { 
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name ND1"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
    }
    if {[string match LYS [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 3])} { 
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name CD"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
      if {("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CE"==[lindex $list_mExpDAid $i_eDA 3])} { 
        set A1_c [atomselect $miMD "resid $ir and name CB"]
        set A2_c [atomselect $miMD "resid $ir and name CG"]
        set A3_c [atomselect $miMD "resid $ir and name CD"]
        set A4_c [atomselect $miMD "resid $ir and name CE"]
        set chi3_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi3 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi3_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi3_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi3_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi3_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi3_LB $chi3_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi3\n"
      }
      if {("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CE"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name NZ"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CG"]
        set A2_c [atomselect $miMD "resid $ir and name CD"]
        set A3_c [atomselect $miMD "resid $ir and name CE"]
        set A4_c [atomselect $miMD "resid $ir and name NZ"]
        set chi4_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi4 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi4_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi4_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi4_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi4_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi4_LB $chi4_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi4\n"
      }
    }
    if {[string match ARG [$selCAi get {resname}]]} {
      if {("resid $ir and name N"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name N"]
        set A2_c [atomselect $miMD "resid $ir and name CA"]
        set A3_c [atomselect $miMD "resid $ir and name CB"]
        set A4_c [atomselect $miMD "resid $ir and name CG"]
        set chi1_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi1 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi1_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi1_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi1_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi1_LB $chi1_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi1\n"
      }
      if {("resid $ir and name CA"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CA"]
        set A2_c [atomselect $miMD "resid $ir and name CB"]
        set A3_c [atomselect $miMD "resid $ir and name CG"]
        set A4_c [atomselect $miMD "resid $ir and name CD"]
        set chi2_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi2 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi2_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi2_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi2_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi2_LB $chi2_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi2\n"
      }
      if {("resid $ir and name CB"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name NE"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CB"]
        set A2_c [atomselect $miMD "resid $ir and name CG"]
        set A3_c [atomselect $miMD "resid $ir and name CD"]
        set A4_c [atomselect $miMD "resid $ir and name NE"]
        set chi3_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi3 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi3_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi3_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi3_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi3_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi3_LB $chi3_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi3\n"
      }
      if {("resid $ir and name CG"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name NE"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name CZ"==[lindex $list_mExpDAid $i_eDA 3])} {
        set A1_c [atomselect $miMD "resid $ir and name CG"]
        set A2_c [atomselect $miMD "resid $ir and name CD"]
        set A3_c [atomselect $miMD "resid $ir and name NE"]
        set A4_c [atomselect $miMD "resid $ir and name CZ"]
        set chi4_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi4 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi4_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi4_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi4_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi4_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi4_LB $chi4_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi4\n"
      }
      if {("resid $ir and name CD"==[lindex $list_mExpDAid $i_eDA 0])&&("resid $ir and name NE"==[lindex $list_mExpDAid $i_eDA 1])&&("resid $ir and name CZ"==[lindex $list_mExpDAid $i_eDA 2])&&("resid $ir and name NH1"==[lindex $list_mExpDAid $i_eDA 3])} {  
        set A1_c [atomselect $miMD "resid $ir and name CD"]
        set A2_c [atomselect $miMD "resid $ir and name NE"]
        set A3_c [atomselect $miMD "resid $ir and name CZ"]
        set A4_c [atomselect $miMD "resid $ir and name NH1"]
        set chi5_lAi [list [concat [$A1_c get index] [$A2_c get index] [$A3_c get index] [$A4_c get index]]]
        set t_chi5 [expr abs([lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5])]
        if {([expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]>=180)} {
          set chi5_UB [expr -[lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]   
          set chi5_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]     
        } else {
          set chi5_UB [expr [lindex $list_mExpDAid $i_eDA 4]+[lindex $list_mExpDAid $i_eDA 5]]
          set chi5_LB [expr [lindex $list_mExpDAid $i_eDA 4]-[lindex $list_mExpDAid $i_eDA 5]]
        }
        puts $Cid_expDAoutfile "[lindex $list_mExpDAid $i_eDA 0]\n[lindex $list_mExpDAid $i_eDA 1]\n[lindex $list_mExpDAid $i_eDA 2]\n[lindex $list_mExpDAid $i_eDA 3] $chi5_LB $chi5_UB ([lindex $list_mExpDAid $i_eDA 4] [lindex $list_mExpDAid $i_eDA 5])\t$t_chi5\n"
      }
    }
  }

  if {$phi_lAi!="NOTEXP"} {
    #puts "phi_lAi exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
	 if {$t_phi>=180} {
        if {($phi_LB <= [measure dihed [lindex $phi_lAi 0]])||([measure dihed [lindex $phi_lAi 0]]<=$phi_UB)} {
          set DAphi($f.r) 1.0
        } else {
          set DAphi($f.r) 0.0
        }	
      }
	 if {$t_phi<=180} {
        if {($phi_LB <= [measure dihed [lindex $phi_lAi 0]])&&([measure dihed [lindex $phi_lAi 0]]<=$phi_UB)} {
          set DAphi($f.r) 1.0
        } else {
          set DAphi($f.r) 0.0
        }
      }
      puts $Cid_expDAoutfile "$f\t[format "%- .2f" [measure dihed [lindex $phi_lAi 0]]]\t\t$phi_LB $phi_UB\t\t\t$t_phi"
    }
    puts $Cid_expDAoutfile "\n"
  }
  if {$phi_lAi=="NOTEXP"} {
   #puts "phi_lAi not exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAphi($f.r) "NaN"
    }
  }
  if {$psi_lAi!="NOTEXP"} {
    #puts "psi_lAi exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
	 if {$t_psi>=180} {
        if {($psi_LB <= [measure dihed [lindex $psi_lAi 0]])||([measure dihed [lindex $psi_lAi 0]]<=$psi_UB)} {
          set DApsi($f.r) 1.0
        } else {
          set DApsi($f.r) 0.0
        }	
      }
	 if {$t_psi<=180} {
        if {($psi_LB <= [measure dihed [lindex $psi_lAi 0]])&&([measure dihed [lindex $psi_lAi 0]]<=$psi_UB)} {
          set DApsi($f.r) 1.0
        } else {
          set DApsi($f.r) 0.0
        }
      }
      puts $Cid_expDAoutfile "$f\t[format "%- .2f" [measure dihed [lindex $psi_lAi 0]]]\t\t$psi_LB $psi_UB\t\t\t$t_psi"
    }
    puts $Cid_expDAoutfile "\n"
  }
  if {$psi_lAi=="NOTEXP"} {
    #puts "psi_lAi not exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DApsi($f.r) "NaN"
    }
  }  
  if {$chi1_lAi!="NOTEXP"} {
    #puts "chi1_lAi exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
	 if {$t_chi1>=180} {
        if {($chi1_LB <= [measure dihed [lindex $chi1_lAi 0]])||([measure dihed [lindex $chi1_lAi 0]]<=$chi1_UB)} {
          set DAchi1($f.r) 1.0
        } else {
          set DAchi1($f.r) 0.0
        }	
      }
	 if {$t_chi1<=180} {
        if {($chi1_LB <= [measure dihed [lindex $chi1_lAi 0]])&&([measure dihed [lindex $chi1_lAi 0]]<=$chi1_UB)} {
          set DAchi1($f.r) 1.0
        } else {
          set DAchi1($f.r) 0.0
        }
      }
      puts $Cid_expDAoutfile "$f\t[format "%- .2f" [measure dihed [lindex $chi1_lAi 0]]]\t\t$chi1_LB $chi1_UB\t\t\t$t_chi1"
    }
    puts $Cid_expDAoutfile "\n"
  }
  if {$chi1_lAi=="NOTEXP"} {
    #puts "chi1_lAi not exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAchi1($f.r) "NaN"
    }
  } 
  if {$chi2_lAi!="NOTEXP"} {
    #puts "chi2_lAi exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
	 if {$t_chi2>=180} {
        if {($chi2_LB <= [measure dihed [lindex $chi2_lAi 0]])||([measure dihed [lindex $chi2_lAi 0]]<=$chi2_UB)} {
          set DAchi2($f.r) 1.0
        } else {
          set DAchi2($f.r) 0.0
        }	
      }
	 if {$t_chi2<=180} {
        if {($chi2_LB <= [measure dihed [lindex $chi2_lAi 0]])&&([measure dihed [lindex $chi2_lAi 0]]<=$chi2_UB)} {
          set DAchi2($f.r) 1.0
        } else {
          set DAchi2($f.r) 0.0
        }
      }
      puts $Cid_expDAoutfile "$f\t[format "%- .2f" [measure dihed [lindex $chi2_lAi 0]]]\t\t$chi2_LB $chi2_UB\t\t\t$t_chi2"
    }
    puts $Cid_expDAoutfile "\n"
  }
  if {$chi2_lAi=="NOTEXP"} {
    #puts "chi2_lAi not exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAchi2($f.r) "NaN"
    }
  } 
  if {$chi3_lAi!="NOTEXP"} {
    #puts "chi3_lAi exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
	 if {$t_chi3>=180} {
        if {($chi3_LB <= [measure dihed [lindex $chi3_lAi 0]])||([measure dihed [lindex $chi3_lAi 0]]<=$chi3_UB)} {
          set DAchi3($f.r) 1.0
        } else {
          set DAchi3($f.r) 0.0
        }	
      }
	 if {$t_chi3<=180} {
        if {($chi3_LB <= [measure dihed [lindex $chi3_lAi 0]])&&([measure dihed [lindex $chi3_lAi 0]]<=$chi3_UB)} {
          set DAchi3($f.r) 1.0
        } else {
          set DAchi3($f.r) 0.0
        }
      }
      puts $Cid_expDAoutfile "$f\t[format "%- .2f" [measure dihed [lindex $chi3_lAi 0]]]\t\t$chi3_LB $chi3_UB\t\t\t$t_chi3"
    }
    puts $Cid_expDAoutfile "\n"
  }
  if {$chi3_lAi=="NOTEXP"} {
    #puts "chi3_lAi not exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAchi3($f.r) "NaN"
    }
  } 
  if {$chi4_lAi!="NOTEXP"} {
    #puts "chi4_lAi exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
	 if {$t_chi4>=180} {
        if {($chi4_LB <= [measure dihed [lindex $chi4_lAi 0]])||([measure dihed [lindex $chi4_lAi 0]]<=$chi4_UB)} {
          set DAchi4($f.r) 1.0
        } else {
          set DAchi4($f.r) 0.0
        }	
      }
	 if {$t_chi4<=180} {
        if {($chi4_LB <= [measure dihed [lindex $chi4_lAi 0]])&&([measure dihed [lindex $chi4_lAi 0]]<=$chi4_UB)} {
          set DAchi4($f.r) 1.0
        } else {
          set DAchi4($f.r) 0.0
        }
      }
      puts $Cid_expDAoutfile "$f\t[format "%- .2f" [measure dihed [lindex $chi4_lAi 0]]]\t\t$chi4_LB $chi4_UB\t\t\t$t_chi4"
    }
    puts $Cid_expDAoutfile "\n"
  }
  if {$chi4_lAi=="NOTEXP"} {
    #puts "chi4_lAi not exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAchi4($f.r) "NaN"
    }
  } 
  if {$chi5_lAi!="NOTEXP"} {
    #puts "chi5_lAi exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
	 if {$t_chi5>=180} {
        if {($chi5_LB <= [measure dihed [lindex $chi5_lAi 0]])||([measure dihed [lindex $chi5_lAi 0]]<=$chi5_UB)} {
          set DAchi5($f.r) 1.0
        } else {
          set DAchi5($f.r) 0.0
        }	
      }
	 if {$t_chi5<=180} {
        if {($chi5_LB <= [measure dihed [lindex $chi5_lAi 0]])&&([measure dihed [lindex $chi5_lAi 0]]<=$chi5_UB)} {
          set DAchi5($f.r) 1.0
        } else {
          set DAchi5($f.r) 0.0
        }
      }
      puts $Cid_expDAoutfile "$f\t[format "%- .2f" [measure dihed [lindex $chi5_lAi 0]]]\t\t$chi5_LB $chi5_UB\t\t\t$t_chi5"
    }
    puts $Cid_expDAoutfile "\n"
  }
  if {$chi5_lAi=="NOTEXP"} {
    #puts "chi5_lAi not exist"
    for {set f 0} {$f < $nf} {incr f} {
      animate goto $f 
      set DAchi5($f.r) "NaN"
    }
  } 

  for {set f 0} {$f < $nf} {incr f} {
    puts $DAoutfile "$f\t$DAphi($f.r)\t$DApsi($f.r)\t$DAchi1($f.r)\t$DAchi2($f.r)\t$DAchi3($f.r)\t$DAchi4($f.r)\t$DAchi5($f.r)"
    if {$f==$MSf} {
      puts $AllresidDAoutfile "$ir\t$DAphi($f.r)\t$DApsi($f.r)\t$DAchi1($f.r)\t$DAchi2($f.r)\t$DAchi3($f.r)\t$DAchi4($f.r)\t$DAchi5($f.r)"
    }
  }  
  close $DAoutfile
}

close $FNDAoutfile
close $AllresidDAoutfile
close $Cid_expDAoutfile

puts "Finished Generating Dihedral Angle Data"
###################
## END OF SCRIPT ##
###################
