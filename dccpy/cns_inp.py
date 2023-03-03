#===========================================================
# this module is for all the script-input file for CNS
#===========================================================


##########################################################
def xplo2d():
    '''the input to generate script input
    '''
    inp='''auto
N
N
N
N
Y
Y
Y
Y
INSERT pdb file here
N
N
QUIT
'''
    
    out='xplo2d.inp'
    fw=open(out, "w")
    fw.write ("%s" % inp)
    fw.close()

    return out
    

##########################################################


def generate():
    '''the input to generate the generate.inp
    '''

    
    inp='''{+ file: generate_easy.inp +}
{+ directory: general +}
{+ description: Generate coordinate and structure file for simple models +}
{+ comment:
           This is designed to be a means of generating a coordinate
           and structure file for commonly encountered models: protein
           and/or DNA/RNA plus waters and maybe ligands. The coordinates
           are provided by the user in a single input PDB file.
           Disulphide bonds will be automatically determined by distance.
           If required generate hydrogens. Any atoms with unknown
           coordinates can be automatically generated +}
{+ authors: Paul Adams and Axel Brunger +}
{+ copyright: Yale University +}

{- Guidelines for using this file:
   - all strings must be quoted by double-quotes
   - logical variables (true/false) are not quoted
   - do not remove any evaluate statements from the file -}

{- Special patches will have to be entered manually at the relevant points
   in the file - see comments throughout the file -}

{- begin block parameter definition -} define(

{============================== important =================================}

{* Different chains in the structure must have either unique segid or
   chainid records. If this is no the case, the end of a chain must
   be delimited by a TER card. *}

{* A break in a chain can be detected automatically or should be delimited
   by a BREAK card. In this case no patch (head, tail or link) will be 
   applied between the residues that bound the chain break. *}

{* NB. The input PDB file must finish with an END statement *}

{=========================== coordinate files =============================}

{* coordinate file *}
{===>} coordinate_infile= "pdb100d.ent_NA";

{* convert chainid to segid if chainid is non-blank *}
{+ choice: true false +}
{===>} convert=true;

{* separate chains by segid - a new segid starts a new chain *}
{+ choice: true false +}
{===>} separate=true;

{============================ renaming atoms ===============================}

{* some atoms may need to be renamed in the topology database to conform
   to what is present in the coordinate file *}

{* delta carbon in isoleucine is named CD in CNS
   what is it currently called in the coordinate file? *}
{* this will not be changed if left blank *}
{===>} ile_CD_becomes="CD1";

{* terminal oxygens are named OT1 and OT2 in CNS
   what are they currently called in the coordinate file? *}
{* these will not be changed if left blank *}
{===>} OT1_becomes="O";
{===>} OT2_becomes="OXT";

{======================= automatic mainchain breaks ========================}

{* automatically detect mainchain breaks in proteins based on distance *}
{* the peptide link at break points will be removed *}
{+ choice: true false +}
{===>} auto_break=true;

{* cutoff distance in Angstroms for identification of breaks *}
{* the default of 2.5A should be reasonable for most cases. If the input
   structure has bad geometry it may be necessary to increase this distance *}
{===>} break_cutoff=2.5;

{======================= automatic disulphide bonds ========================}

{* cutoff distance in Angstroms for identification of disulphides *}
{* the default of 3.0A should be reasonable for most cases. If the input
   structure has bad geometry it may be necessary to increase this distance *}
{===>} disulphide_dist=3.0;

{========================= RNA to DNA conversion  ==========================}

{* All nucleic acid residues initially have ribose sugars (rather than
   deoxyribose). A patch must be applied to convert the ribose to deoxyribose
   for DNA residues. Select those residues which need to have the patch
   applied to make them DNA. *}
{* Make sure that the atom selection is specific for the nucleic acid
   residues *}
{===>} dna_sele=(none);
	    
{========================= generate parameters =============================}

{* hydrogen flag - determines whether hydrogens will be output *}
{* must be true for NMR, atomic resolution X-ray crystallography 
   or modelling.  Set to false for most X-ray crystallographic 
   applications at resolution > 1A *}
{+ choice: true false +}
{===>} hydrogen_flag=false;

{* which hydrogens to build *}
{+ choice: "all" "unknown" +}
{===>} hydrogen_build="all";

{* selection of atoms other than hydrogens for which coordinates
   will be generated *}
{* to generate coordinates for all unknown atoms use: (not(known)) *}
{===>} atom_build=(none);

{* selection of atoms to be deleted *}
{* to delete no atoms use: (none) *}
{===>} atom_delete=(none);

{* set bfactor flag *}
{+ choice: true false +}
{===>} set_bfactor=false;

{* set bfactor value *}
{===>} bfactor=15.0;

{* set occupancy flag *}
{+ choice: true false +}
{===>} set_occupancy=false;

{* set occupancy value *}
{===>} occupancy=1.0;

{============================= output files ================================}

{* output structure file *}
{===>} structure_outfile="generate_15902.mtf"; 

{* output coordinate file *}
{===>} coordinate_outfile="generate_15902.pdb"; 

{* format output coordinates for use in o *}
{* if false then the default CNS output coordinate format will be used *}
{+ choice: true false +}
{===>} pdb_o_format=true;

{================== protein topology and parameter files ===================}

{* protein topology file *}
{===>} prot_topology_infile="CNS_TOPPAR:protein.top";

{* protein linkage file *}
{===>} prot_link_infile="CNS_TOPPAR:protein.link";

{* protein parameter file *}
{===>} prot_parameter_infile="CNS_TOPPAR:protein_rep.param";

{================ nucleic acid topology and parameter files =================}

{* nucleic acid topology file *}
{===>} nucl_topology_infile="CNS_TOPPAR:dna-rna.top";

{* nucleic acid linkage file *}
{* use CNS_TOPPAR:dna-rna-pho.link for 5'-phosphate *}
{===>} nucl_link_infile="CNS_TOPPAR:dna-rna.link";

{* nucleic acid parameter file *}
{===>} nucl_parameter_infile="CNS_TOPPAR:dna-rna_rep.param";

{=================== water topology and parameter files ====================}

{* water topology file *}
{===>} water_topology_infile="CNS_TOPPAR:water.top";

{* water parameter file *}
{===>} water_parameter_infile="CNS_TOPPAR:water_rep.param";

{================= carbohydrate topology and parameter files ===============}

{* carbohydrate topology file *}
{===>} carbo_topology_infile="CNS_TOPPAR:carbohydrate.top";

{* carbohydrate parameter file *}
{===>} carbo_parameter_infile="CNS_TOPPAR:carbohydrate.param";

{============= prosthetic group topology and parameter files ===============}

{* prosthetic group topology file *}
{===>} prost_topology_infile="";

{* prosthetic group parameter file *}
{===>} prost_parameter_infile="";

{=================== ligand topology and parameter files ===================}

{* ligand topology file *}
{===>} ligand_topology_infile="LIGAND.top"; 

{* ligand parameter file *}
{===>} ligand_parameter_infile="LIGAND.param"; 

{===================== ion topology and parameter files ====================}

{* ion topology file *}
{===>} ion_topology_infile="CNS_TOPPAR:ion.top";

{* ion parameter file *}
{===>} ion_parameter_infile="CNS_TOPPAR:ion.param";

{===========================================================================}
{         things below this line do not need to be changed                  }
{===========================================================================}

 ) {- end block parameter definition -}

 checkversion 1.3

 evaluate ($log_level=quiet)

 topology
   if ( &BLANK%prot_topology_infile = false ) then
     @@&prot_topology_infile
   end if
   if ( &BLANK%nucl_topology_infile = false ) then
     @@&nucl_topology_infile
   end if
   if ( &BLANK%water_topology_infile = false ) then
     @@&water_topology_infile
   end if
   if ( &BLANK%carbo_topology_infile = false ) then
     @@&carbo_topology_infile
   end if
   if ( &BLANK%prost_topology_infile = false ) then
     @@&prost_topology_infile
   end if
   if ( &BLANK%ligand_topology_infile = false ) then
     @@&ligand_topology_infile
   end if
   if ( &BLANK%ion_topology_infile = false ) then
     @@&ion_topology_infile
   end if
 end

 topology
   if ( &BLANK%prot_break_infile = false ) then
     @@&prot_break_infile
   end if
 end

 parameter
   if ( &BLANK%prot_parameter_infile = false ) then
     @@&prot_parameter_infile
   end if
   if ( &BLANK%nucl_parameter_infile = false ) then
     @@&nucl_parameter_infile
   end if
   if ( &BLANK%water_parameter_infile = false ) then
     @@&water_parameter_infile
   end if
   if ( &BLANK%carbo_parameter_infile = false ) then
     @@&carbo_parameter_infile
   end if
   if ( &BLANK%prost_parameter_infile = false ) then
     @@&prost_parameter_infile
   end if
   if ( &BLANK%ligand_parameter_infile = false ) then
     @@&ligand_parameter_infile
   end if
   if ( &BLANK%ion_parameter_infile = false ) then
     @@&ion_parameter_infile
   end if
 end

 segment
   chain
     if ( &convert = true ) then
       convert=true
     end if
     if ( &separate = true ) then
       separate=true
     end if
     @@&prot_link_infile
     @@&nucl_link_infile
     coordinates @@&coordinate_infile
   end
 end

 if ( &BLANK%ile_CD_becomes = false ) then
   do (name=&ile_CD_becomes) (resname ILE and name CD)
 end if
 if ( &BLANK%OT1_becomes = false ) then
   do (name=&OT1_becomes) (name OT1)
 end if
 if ( &BLANK%OT2_becomes = false ) then
   do (name=&OT2_becomes) (name OT2)
 end if

 coordinates 
   if ( &convert = true ) then
     convert=true
   end if
   @@&coordinate_infile

 set echo=off end
 show sum(1) ( not(hydrogen) and not(known) )
 if ( $select = 0 ) then
   display  %INFO: There are no coordinates missing for non-hydrogen atoms
 end if
 set echo=on end

 if ( $log_level = verbose ) then
   set message=normal echo=on end
 else
   set message=off echo=off end
 end if

 if ( &auto_break = true ) then

   evaluate ($break=0)

   for $id1 in id ( name C and bondedto(name CA) and bondedto(name O) ) loop break

     show (segid) (id $id1)
     evaluate ($segid1=$result)
     show (resid) (id $id1)
     evaluate ($resid1=$result)
     show (resname) (id $id1)
     evaluate ($resname1=$result)

     show sum(1) (id $id1 and known)
     if ( $result = 0 ) then
       display unknown coordinates for segid $segid1 resname $resname1 resid $resid1 name C
       display this coordinate must be known for automatic chain break detection
       abort
     end if

     identity (store1) ( name N and bondedto( segid $segid1 and resid $resid1 and name c ) )

     if ( $select = 1 ) then
       show element (store1) (attribute store1 > 0)
       evaluate ($id2=$result)
       show (segid) (id $id2)
       evaluate ($segid2=$result)
       show (resid) (id $id2)
       evaluate ($resid2=$result)
       show (resname) (id $id2)
       evaluate ($resname2=$result)

       show sum(1) (id $id2 and known)
       if ( $result = 0 ) then
         display unknown coordinates for segid $segid2 resname $resname2 resid $resid2 name N
         display this coordinate must be known for automatic chain break detection
         abort
       end if

       pick bond
         (name c and segid $segid1 and resid $resid1)
         (name n and segid $segid2 and resid $resid2)
         geometry

       if ( $result > &break_cutoff ) then
         evaluate ($break=$break+1)
         evaluate ($seg1.$break=$segid1)
         evaluate ($res1.$break=$resid1)
         evaluate ($seg2.$break=$segid2)
         evaluate ($res2.$break=$resid2)
         if ( $resname2 = PRO ) then
           evaluate ($patch.$break=DPPP)
         elseif ( $resname2 = CPR ) then
           evaluate ($patch.$break=DPPP)
         else
           evaluate ($patch.$break=DPEP)
         end if
       end if
     end if

   end loop break

   evaluate ($counter=1)

   while ($counter <= $break) loop delete
     patch $patch.$counter
       reference=-=(segid $seg1.$counter and resid $res1.$counter)
       reference=+=(segid $seg2.$counter and resid $res2.$counter)
     end
     buffer message
       display peptide link removed (applied $patch.$counter): from \
$seg1.$counter[a4] $res1.$counter[a4] to $seg2.$counter[a4] $res2.$counter[a4]
     end
     evaluate ($counter=$counter+1)
   end loop delete

 end if

 evaluate ($disu=0)

 for $id1 in id ( resname CYS and name SG ) loop dis1

   show (segid) (id $id1)
   evaluate ($segid1=$result)
   show (resid) (id $id1)
   evaluate ($resid1=$result)

   identity (store1) (all)

   for $id2 in id ( resname CYS and name SG and 
                  ( attr store1 > $id1 ) ) loop dis2

     show (segid) (id $id2)
     evaluate ($segid2=$result)
     show (resid) (id $id2)
     evaluate ($resid2=$result)

     pick bond (id $id1) (id $id2) geometry

     if ( $result <= &disulphide_dist ) then
       evaluate ($disu=$disu+1)
       evaluate ($seg1.$disu=$segid1)
       evaluate ($seg2.$disu=$segid2)
       evaluate ($res1.$disu=$resid1)
       evaluate ($res2.$disu=$resid2)
     end if

   end loop dis2

 end loop dis1

 evaluate ($counter=1)
 while ( $counter <= $disu ) loop disu
   patch disu
     reference=1=(segid $seg1.$counter and resid $res1.$counter)
     reference=2=(segid $seg2.$counter and resid $res2.$counter)
   end
   buffer message
     display disulphide added: from \
$seg1.$counter[a4] $res1.$counter[a4] to $seg2.$counter[a4] $res2.$counter[a4]
   end
   evaluate ($counter=$counter+1)
 end loop disu

 {- patching of RNA to DNA -}
 evaluate ($counter=0)
 for $id in id ( tag and (&dna_sele) ) loop dna
   evaluate ($counter=$counter+1)
   show (segid) (id $id)
   evaluate ($dna.segid.$counter=$result)
   show (resid) (id $id)
   evaluate ($dna.resid.$counter=$result)
 end loop dna
 evaluate ($dna.num=$counter)

 evaluate ($counter=0)
 while ($counter < $dna.num) loop dnap
   evaluate ($counter=$counter+1)
   patch deox reference=nil=(segid $dna.segid.$counter and
                             resid $dna.resid.$counter) end
 end loop dnap

 set messages=normal end
 set echo=on end

 if (&hydrogen_flag=false) then
   delete selection=( hydrogen ) end
 end if

 delete selection=( &atom_delete ) end

 identity (store1) (none)

 identity (store1) (&atom_build)
 if ( &hydrogen_build = "all" ) then
   identity (store1) (store1 or hydrogen)
 elseif ( &hydrogen_build = "unknown" ) then
   identity (store1) (store1 or (not(known) and hydrogen))
 end if

 show sum(1) (store1)
 evaluate ($tobuild=$result)

 if ( $tobuild > 0 ) then

   fix selection=(not(store1)) end

   show sum(1) (store1)
   evaluate ($moving=$result)

   if ( $moving > 0 ) then
     for $id in id (tag and byres(store1)) loop avco

       show ave(x) (byres(id $id) and known)
       evaluate ($ave_x=$result)
       show ave(y) (byres(id $id) and known)
       evaluate ($ave_y=$result)
       show ave(z) (byres(id $id) and known)
       evaluate ($ave_z=$result)

       do (x=$ave_x) (byres(id $id) and store1)
       do (y=$ave_y) (byres(id $id) and store1)
       do (z=$ave_z) (byres(id $id) and store1)
 
     end loop avco 

     do (x=x+random(2.0)) (store1)
     do (y=y+random(2.0)) (store1)
     do (z=z+random(2.0)) (store1)

     {- start parameter for the side chain building -}
     parameter
       nbonds
         rcon=20. nbxmod=-2 repel=0.9  wmin=0.1 tolerance=1.
         rexp=2 irexp=2 inhibit=0.25
       end
     end

     {- Friction coefficient, in 1/ps. -}
     do (fbeta=100) (store1)

     evaluate ($bath=300.0)
     evaluate ($nstep=500)
     evaluate ($timestep=0.0005)

     do (refy=mass) (store1)

     do (mass=20) (store1)

     igroup interaction 
       (store1) (store1 or known)
     end

     {- turn on initial energy terms -}
     flags exclude * include bond angle vdw end
 
     minimize powell nstep=50  nprint=10 end

     do (vx=maxwell($bath)) (store1)
     do (vy=maxwell($bath)) (store1)
     do (vz=maxwell($bath)) (store1)

     flags exclude vdw include impr end

     dynamics cartesian
       nstep=50
       timestep=$timestep
       tcoupling=true temperature=$bath
       nprint=$nstep
       cmremove=false
     end

     flags include vdw end

     minimize powell nstep=50 nprint=10 end

     do (vx=maxwell($bath)) (store1)
     do (vy=maxwell($bath)) (store1)
     do (vz=maxwell($bath)) (store1)

     dynamics cartesian
       nstep=50
       timestep=$timestep
       tcoupling=true temperature=$bath
       nprint=$nstep
       cmremove=false
     end

     parameter
       nbonds
         rcon=2. nbxmod=-3 repel=0.75
       end
     end

     minimize powell nstep=100 nprint=25 end

     do (vx=maxwell($bath)) (store1)
     do (vy=maxwell($bath)) (store1)
     do (vz=maxwell($bath)) (store1)

     dynamics cartesian
       nstep=$nstep
       timestep=$timestep
       tcoupling=true temperature=$bath
       nprint=$nstep
       cmremove=false
     end

     {- turn on all energy terms -}
     flags include dihe ? end

     {- set repel to ~vdw radii -}
     parameter
       nbonds
         repel=0.89
       end
     end

     minimize powell nstep=500 nprint=50 end

     flags exclude * include bond angl impr dihe vdw end

     {- return masses to something sensible -}
     do (mass=refy) (store1)

     do (vx=maxwell($bath)) (store1)
     do (vy=maxwell($bath)) (store1)
     do (vz=maxwell($bath)) (store1)

     dynamics cartesian
       nstep=$nstep
       timestep=$timestep
       tcoupling=true temperature=$bath
       nprint=$nstep
       cmremove=false
     end

     {- some final minimisation -}
     minimize powell
       nstep=500
       drop=40.0
       nprint=50
     end

     print thres=0.02 bonds
     print thres=5. angles

   end if
  
   fix selection=( none ) end

 end if

 set echo=false end
 show sum(1) (not(known))
 if ( $result < 100 ) then
   for $id in id (not(known)) loop print
     show (segid) (id $id)
     evaluate ($segid=$result)
     show (resname) (id $id)
     evaluate ($resname=$result)
     show (resid) (id $id)
     evaluate ($resid=$result)
     show (name) (id $id)
     evaluate ($name=$result)
     buffer message
       display unknown coordinates for atom: $segid[a4] $resname[a4] $resid[a4] $name[a4]
     end
   end loop print
 else
   buffer message
     display unknown coordinates for more than 100 atoms
   end
 end if
 set echo=true end

 if (&set_bfactor=true) then
   do (b=&bfactor) ( all )
 else
   show ave(b) (known and not(store1))
   do (b=$result) (store1 and (attr b < 0.01))
 end if

 if (&set_occupancy=true) then
   do (q=&occupancy) ( all )
 end if

 set echo=false end
 show sum(1) (store1)
 if ( $result < 100 ) then
   for $id in id (store1) loop print
     show (segid) (id $id)
     evaluate ($segid=$result)
     show (resname) (id $id)
     evaluate ($resname=$result)
     show (resid) (id $id)
     evaluate ($resid=$result)
     show (name) (id $id)
     evaluate ($name=$result)
     buffer message
       display coordinates built for atom: $segid[a4] $resname[a4] $resid[a4] $name[a4]
     end 
   end loop print
 else
   buffer message
     display coordinates built for more than 100 hundred atoms
   end
 end if
 set echo=true end

 set remarks=reset end

 buffer message
   to=remarks
   dump
 end

 write structure output=&structure_outfile end

 if ( &pdb_o_format = true ) then
   write coordinates format=PDBO output=&coordinate_outfile end
 else
   write coordinates output=&coordinate_outfile end
 end if

 stop
'''

    out='generate.inp'
    fw=open(out, 'w')
    fw.write ("%s" %inp )
    fw.close()
    return out
    
#########################################################
    
def model_stat():
    
    inp='''{+ file: model_stats.inp +}
{+ directory: xtal_refine +}
{+ description: Crystallographic model statistics +}
{+ authors: Axel T. Brunger, and Paul D. Adams +}
{+ copyright: Yale University +}

{- Guidelines for using this file:
   - all strings must be quoted by double-quotes
   - logical variables (true/false) are not quoted
   - do not remove any evaluate statements from the file 
   - the selections store3 through store9 are available for general use -}

{- begin block parameter definition -} define(

{============================ coordinates ============================}

{* coordinate file *}
{===>} coordinate_infile= "generate_15902.pdb";

{==================== molecular information ==========================}

{* topology files *}
{===>} topology_infile_1="CNS_TOPPAR:protein.top";
{===>} topology_infile_2="CNS_TOPPAR:dna-rna.top";
{===>} topology_infile_3="CNS_TOPPAR:water.top";
{===>} topology_infile_4="CNS_TOPPAR:ion.top";
{===>} topology_infile_5="CNS_TOPPAR:carbohydrate.top";
{===>} topology_infile_6="";
{===>} topology_infile_7="";
{===>} topology_infile_8="";

{* linkage files for linear, continuous polymers (protein, DNA, RNA) *}
{===>} link_infile_1="CNS_TOPPAR:protein.link";
{===>} link_infile_2="CNS_TOPPAR:dna-rna-pho.link";
{===>} link_infile_3="";

{* parameter files *}
{===>} parameter_infile_1= "CNS_TOPPAR:protein_rep.param";
{===>} parameter_infile_2= "CNS_TOPPAR:water_rep.param";
{===>} parameter_infile_3= "CNS_TOPPAR:dna-rna_rep.param";
{===>} parameter_infile_4= "CNS_TOPPAR:ion.param";
{===>} parameter_infile_5= "LIGAND.param";
{===>} parameter_infile_6="";
{===>} parameter_infile_7="";
{===>} parameter_infile_8="";

{* molecular topology file: optional (leave blank for auto generation) *}
{* 
   Auto generation of the molecular topology from the coordinates should only 
   be used if:
   (1) Each distinct protein, DNA, or RNA chain must have a separate segid 
       (or chainid if the chainid is non-blank). 
   (2) Each contiguous protein, RNA, or RNA chain must not be disrupted by 
       other types of residues or ligands.  Rather, these other residues 
       should be listed after protein, RNA/DNA chains. 
   (3) Disulphides are automatically detected based on distances between the sulfur atoms
      (must be less than 3 A apart).
   (4) Broken protein/RNA/DNA chains without terminii must be more than 2.5 A apart to be recognized as such.
   (5) N-linked glycan links are automatically recognized if the bonded atoms are less than 2.5 A apart.
   (6) Automatic generation cannot be used with alternate conformations. 
   For ligands, the user must make suitable topology and parameter files.
   For non-standard covalent linkages, the custom patch file should be used.
   Alternatively, the generate.inp or generate_easy.inp task files
   can be used to generated the mtf prior to running this task file.
    *}
{===>} structure_infile= "generate_15902.mtf";

{* for auto generation: extra linkages and modifications by custom patches *}
{===>} patch_infile="";

{* force field settings file *}
{===>} force_field_infile="";

{====================== crystallographic data ========================}

{* space group *}
{* use International Table conventions with subscripts substituted
   by parenthesis *}
{===>} sg="P2(1)2(1)2(1)";

{* unit cell parameters in Angstroms and degrees *}
{+ table: rows=1 "cell" cols=6 "a" "b" "c" "alpha" "beta" "gamma" +}
{===>} a=23.98;
{===>} b=40.77;
{===>} c=44.84;
{===>} alpha=90.00;
{===>} beta=90.00;
{===>} gamma=90.00;

{* anomalous f' f'' library file *}
{* If a file is not specified, no anomalous contribution will be included *}
{+ choice: "CNS_XRAYLIB:anom_cu.lib" "CNS_XRAYLIB:anom_mo.lib" "" user_file +}
{===>} anom_library="";

{* reflection files *}
{* specify non-anomalous reflection files before anomalous reflection files. *}
{* files must contain unique array names otherwise errors will occur *}
{===>} reflection_infile_1="r100dsf.ent.cns";
{===>} reflection_infile_2="";
{===>} reflection_infile_3="";
{===>} reflection_infile_4="";

{* reciprocal space array containing observed amplitudes: required *}
{===>} obs_f="fobs";

{* reciprocal space array containing sigma values for amplitudes: required *}
{===>} obs_sigf="sigma";

{* reciprocal space array containing test set for cross-validation: optional *}
{* required for the calculation of cross-validated sigmaA values *}
{* required for the maximum likelihood targets *}
{===>} test_set="test";

{* number for selection of test reflections: required for cross-validation *}
{* ie. reflections with the test set array equal to this number will be
       used for cross-validation, all other reflections form the working set *}
{===>} test_flag=1;

{* reciprocal space array containing weighting scheme for observed 
   amplitudes: optional *}
{* only used for the "residual" and "vector" targets - this will
   default to a constant value of 1 if array is not present *}
{===>} obs_w="";

{* reciprocal space array containing observed intensities: optional *}
{* required for the "mli" target *}
{===>} obs_i="";

{* reciprocal space array containing sigma values for intensities: optional *}
{* required for the "mli" target *}
{===>} obs_sigi="";

{* experimental phase probability distribution: optional *}
{* required for the "mlhl" target *}
{* reciprocal space arrays with Hendrickson-Lattman coefficients A,B,C,D *}
{+ table: rows=1 "HL coefficients" cols=4 "A" "B" "C" "D" +}
{===>} obs_pa="";
{===>} obs_pb="";
{===>} obs_pc="";
{===>} obs_pd="";

{* complex reciprocal space array containing experimental phases: optional *}
{* required for the "mixed" and "vector" targets *}
{===>} obs_phase="";

{* reciprocal space array containing experimental figures of merit: optional *}
{* required for the "mixed" target *}
{===>} obs_fom="";

{* resolution limits *}
{* the full resolution range of observed data should be used in refinement.
   A bulk solvent correction should be applied to allow the use of low
   resolution terms. If no bulk solvent correction is applied, data must
   be truncated at a lower resolution limit of between 8 and 6 Angstrom. *}
{+ table: rows=1 "resolution" cols=2 "lowest" "highest" +}
{===>} low_res=8.00;
{===>} high_res=1.90;

{* apply rejection criteria to amplitudes or intensities *}
{+ choice: "amplitude" "intensity" +}
{===>} obs_type="amplitude";

{* Observed data cutoff criteria: applied to amplitudes or intensities *}
{* reflections with magnitude(Obs)/sigma < cutoff are rejected. *}
{===>} sigma_cut=0.0;

{* rms outlier cutoff: applied to amplitudes or intensities *}
{* reflections with magnitude(Obs) > cutoff*rms(Obs) will be rejected *}
{===>} obs_rms=10000;

{=================== non-crystallographic symmetry ===================}

{* NCS-restraints/constraints file *}
{* see auxiliary/ncs.def *}
{===>} ncs_infile="";

{============ overall B-factor and bulk solvent corrections ==========}

{* overall B-factor correction *}
{+ choice: "no" "isotropic" "anisotropic" +}
{===>} bscale="anisotropic";

{* bulk solvent correction *}
{* a mask is required around the molecule(s). The region
   outside this mask is the solvent region *}
{+ choice: true false +}
{===>} bulk_sol=true;

{* bulk solvent mask file *}
{* mask will be read from O type mask file if a name is given
   otherwise calculated from coordinates of selected atoms *}

{===>} bulk_mask_infile="";

{* automatic bulk solvent parameter optimization for e-density level sol_k (e/A^3) and B-factor sol_b (A^2) *}
{+ choice: true false +}
{===>} sol_auto=true;

{* fixed solvent parameters (used if the automatic option is turned off) *}
{+ table: rows=1 "bulk solvent" cols=2 "e-density level sol_k (e/A^3)" "B-factor sol_b (A^2) " +}
{===>} sol_k=0.3;
{===>} sol_b=50.0;

{* optional file with a listing of the results of the automatic bulk solvent optimization *}
{===>} sol_output="";

{* solvent mask parameters *}
{+ table: rows=1 "bulk solvent" cols=2 "probe radius (A) (usually set to 1)" "shrink radius (A) (usually set to 1)" +}
{===>} sol_rad=1.0;
{===>} sol_shrink=1.0;

{========================== atom selection ===========================}

{* select atoms to be included in statistics *}
{* this should include all conformations if multiple conformations are used *}
{===>} atom_select=(known and not hydrogen);

{* select fixed atoms *}
{* note: atoms at special positions are automatically fixed. So, 
   you don't have to explicitly fix them here. *}
{===>} atom_fixed=(none);

{* select atoms to be harmonically restrained during refinement *}
{===>} atom_harm=(none);

{* harmonic restraint constant - for harmonically restrained atoms *}
{===>} k_harmonic=10;

{* select atoms in alternate conformation 1 *}
{===>} conf_1=(none);

{* select atoms in alternate conformation 2 *}
{===>} conf_2=(none);

{* select atoms in alternate conformation 3 *}
{===>} conf_3=(none);

{* select atoms in alternate conformation 4 *}
{===>} conf_4=(none);

{* define main chain atoms - for B-factor analysis *}
{* note: atoms outside this selection will be considered to be
   side chain atoms *}
{===>} atom_main=(name ca or name n or name c or name o or name ot+);

{* additional restraints file *}
{* eg. auxiliary/dna-rna_restraints.def *}
{===>} restraints_infile="";

{======================= statistics parameters =======================}

{* number of bins for printing binwise statistics *}
{===>} print_bins=10;

{* threshhold for reported bond violations *}
{===>} bond_thresh=0.05;

{* threshhold for reported angle violations *}
{===>} angle_thresh=8.0;

{* threshhold for reported dihedral violations *}
{===>} dihe_thresh=60.0;

{* threshhold for reported improper violations *}
{===>} impr_thresh=3.0;

{* threshhold for reported close contacts *}
{===>} close=2.5;

{* lower resolution limit for coordinate error estimation *}
{===>} low_err_res=5.0;

{* refinement target *}
{+ list: mlf: maximum likelihood target using amplitudes
         mli: maximum likelihood target using intensities
        mlhl: maximum likelihood target using amplitudes
              and phase probability distribution
    residual: standard crystallographic residual
      vector: vector residual
       mixed: (1-fom)*residual + fom*vector
        e2e2: correlation coefficient using normalized E^2
        e1e1: correlation coefficient using normalized E
        f2f2: correlation coefficient using F^2
        f1f1: correlation coefficient using F +}
{+ choice: "mlf" "mli" "mlhl" "residual" "vector" "mixed"
           "e2e2" "e1e1" "f2f2" "f1f1" +}
{===>} reftarget="mlf";

{* number of bins for refinement target *}
{* this will be determined automatically if a negative value is given
   otherwise the specified number of bins will be used *}
{===>} target_bins=-1;

{* Wa weight for X-ray term *}
{* this is determined automatically if a negative value is given *}
{===>} wa=-1;

{* memory allocation for FFT calculation *}
{* this will be determined automatically if a negative value is given
   otherwise the specified number of words will be allocated *}
{===>} fft_memory=-1;

{=========================== output files ============================}

{* output statistics file *}
{===>} list_outfile="model_stats.list";

{* root for luzzati coordinate error plot files *}
{* files: "root".plot and "root"_cv.plot
   no plot file will be written if left blank *}
{===>} luzzati_error_plot="luzzati_error";

{* root for sigmaa coordinate error plot files *}
{* files: "root".plot and "root"_cv.plot
   no plot file will be written if left blank *}
{===>} sigmaa_error_plot="sigmaa_error";

{===========================================================================}
{        things below this line do not normally need to be changed          }
{===========================================================================}

 ) {- end block parameter definition -}

 checkversion 1.3

 evaluate ($log_level=quiet)

 if ( $log_level = verbose ) then
   set message=normal echo=on end
 else
   set message=off echo=off end
 end if

 if ( &BLANK%structure_infile = true ) then
 
    {- read topology files -}
    topology
     evaluate ($counter=1)
     evaluate ($done=false)
     while ( $done = false ) loop read
      if ( &exist_topology_infile_$counter = true ) then
         if ( &BLANK%topology_infile_$counter = false ) then
            @@&topology_infile_$counter
         end if
      else
        evaluate ($done=true)
      end if
      evaluate ($counter=$counter+1)
     end loop read
    end
    
    @CNS_XTALMODULE:mtfautogenerate (
                                  coordinate_infile=&coordinate_infile;
                                  convert=true;
                                  separate=true;
                                  atom_delete=(not known);
                                  hydrogen_flag=true;
                                  break_cutoff=2.5;
                                  disulphide_dist=3.0;
                                  carbo_dist=2.5;
                                  patch_infile=&patch_infile;
                                  O5_becomes="O";
                                 )

 else

   structure @&structure_infile end
   coordinates @&coordinate_infile

 end if

 {- read parameter files -}
 parameter
  evaluate ($counter=1)
  evaluate ($done=false)
  while ( $done = false ) loop read
   if ( &exist_parameter_infile_$counter = true ) then
      if ( &BLANK%parameter_infile_$counter = false ) then
         @@&parameter_infile_$counter
      end if
   else
     evaluate ($done=true)
   end if
   evaluate ($counter=$counter+1)
  end loop read
 end

 xray

   @CNS_XTALLIB:spacegroup.lib (sg=&sg;)

   a=&a b=&b c=&c  alpha=&alpha beta=&beta gamma=&gamma

   @CNS_XRAYLIB:scatter.lib

   evaluate ($counter=1)
   evaluate ($done=false)
   while ( $done = false ) loop read
    if ( &exist_reflection_infile_$counter = true ) then
      if ( &BLANK%reflection_infile_$counter = false ) then
       reflection
         @@&reflection_infile_$counter
       end
      end if
   else
     evaluate ($done=true)
   end if
   evaluate ($counter=$counter+1)
  end loop read

 end

 if ( &BLANK%anom_library = false ) then
   @@&anom_library
 else
   xray anomalous=? end
   if ( $result = true ) then
     display Warning: no anomalous library has been specified
     display          no anomalous contribution will used in refinement
   end if
 end if

 {- copy define parameters of optional arrays into symbols so 
    we can redefine them -}
    
 evaluate ($obs_i=&obs_i)
 evaluate ($obs_sigi=&obs_sigi)
 evaluate ($obs_w=&obs_w)
 xray
   @@CNS_XTALMODULE:checkrefinput (
                                  reftarget=&reftarget;
                                  obs_f=&obs_f;
                                  obs_sigf=&obs_sigf;
                                  test_set=&test_set;
                                  obs_pa=&obs_pa;
                                  obs_pb=&obs_pb;
                                  obs_pc=&obs_pc;
                                  obs_pd=&obs_pd;
                                  obs_phase=&obs_phase;
                                  obs_fom=&obs_fom;
                                  obs_w=$obs_w;
                                  obs_i=$obs_i;
                                  obs_sigi=$obs_sigi;
                                  )

   query name=fcalc domain=reciprocal end
   if ( $object_exist = false ) then
      declare name=fcalc domain=reciprocal type=complex end
   end if
   declare name=fbulk domain=reciprocal type=complex end

   do (fbulk=0) ( all ) 

   binresolution &low_res &high_res
   mapresolution &high_res

   if ( &obs_type = "intensity" ) then
     if ( &BLANK%obs_i = true ) then
       display  Error: observed intensity array is undefined
       display         aborting script
       abort
     end if
     evaluate ($reject_obs=&obs_i)
     evaluate ($reject_sig=&obs_sigi)
   else
     evaluate ($reject_obs=&obs_f)
     evaluate ($reject_sig=&obs_sigf)
   end if

   declare name=ref_active domain=reciprocal type=integer end
   declare name=tst_active domain=reciprocal type=integer end

   do (ref_active=0) ( all )
   do (ref_active=1) ( ( $STRIP%reject_sig # 0 ) and
                      ( &low_res >= d >= &high_res ) )

   statistics overall
     completeness
     selection=( ref_active=1 )
   end
   evaluate ($total_compl=$expression1)

   show sum(1) ( ref_active=1 )
   evaluate ($total_read=$select)
   evaluate ($total_theor=int(1./$total_compl * $total_read))

   show rms (amplitude($STRIP%reject_obs)) ( ref_active=1 )
   evaluate ($obs_high=$result*&obs_rms)
   show min (amplitude($STRIP%reject_obs)) ( ref_active=1 )
   evaluate ($obs_low=$result)

   do (ref_active=0) ( all )
   do (ref_active=1)
                  ( ( amplitude($STRIP%reject_obs) > &sigma_cut*$STRIP%reject_sig ) and
                    ( $STRIP%reject_sig # 0 ) and
                    ( $obs_low <= amplitude($STRIP%reject_obs) <= $obs_high ) and
                    ( &low_res >= d >= &high_res ) )

   do (tst_active=0) (all)
   if ( &BLANK%test_set = false ) then
     do (tst_active=1) (ref_active=1 and &STRIP%test_set=&test_flag)
   end if

   show sum(1) ( ref_active=1 and tst_active=0 )
   evaluate ($total_work=$select)
   show sum(1) ( ref_active=1 and tst_active=1 )
   evaluate ($total_test=$select)
   evaluate ($total_used=$total_work+$total_test)

   evaluate ($unobserved=$total_theor-$total_read)
   evaluate ($rejected=$total_read-$total_used)
   evaluate ($per_unobs=100*($unobserved/$total_theor))
   evaluate ($per_reject=100*($rejected/$total_theor))
   evaluate ($per_used=100*($total_used/$total_theor))
   evaluate ($per_work=100*($total_work/$total_theor))
   evaluate ($per_test=100*($total_test/$total_theor))

   associate fcalc ( &atom_select )

   tselection=( ref_active=1 )

   cvselection=( tst_active=1 )

   method=FFT          
   
 {- MODIFIED 2/15/06 -}
 end
 

 show min ( b ) ( &atom_select )
 evaluate ($b_min=$result)
 @@CNS_XTALMODULE:fft_parameter_check ( 
                             d_min=&high_res; 
                             b_min=$b_min;
                             grid=auto;
                             fft_memory=&fft_memory;
                             fft_grid=$fft_grid;   
                             fft_b_add=$fft_b_add; 
                             fft_elim=$fft_elim; 
                                      )
                            
 xray
 {- END MODIFICATION -}

   tolerance=0.0 lookup=false

   if ( &wa >= 0 ) then
      wa=&wa
   end if

 end                  

 if ( &BLANK%ncs_infile = false ) then
    inline @&ncs_infile
 end if

 if ( &BLANK%restraints_infile = false ) then
     @&restraints_infile
 end if

 do (store1=0) (all)

 evaluate ($nalt=1)
 evaluate ($alt=1)
 evaluate ($done=false)
 while ( $done = false ) loop nalt
   if ( &exist_conf_$alt = true ) then
     show sum(1) ( &conf_$alt )
     if ( $result > 0 ) then
       evaluate ($nalt=$nalt+1)
     end if
   else
     evaluate ($done=true)
     evaluate ($nalt=$nalt-1)
   end if
   evaluate ($alt=$alt+1)
 end loop nalt

 evaluate ($alt=1)
 while ( $alt <= $nalt ) loop alt
   do (store1=$alt) ( &conf_$alt )
   evaluate ($alt=$alt+1)
 end loop alt

 igroup
   interaction ( &atom_select and not(attr store1 > 0))
               ( &atom_select and not(attr store1 > 0))
   evaluate ($alt=1)
   while ( $alt <= $nalt ) loop alcs
     interaction ( &atom_select and ( attr store1 = $alt or attr store1 = 0 ))
                 ( &atom_select and ( attr store1 = $alt ))
     evaluate ($alt=$alt+1)
   end loop alcs
 end
 
 {- check isolated atoms and atoms at special positions and add to
    list of fixed atoms if needed - store1 will be used -}
 
 @CNS_XTALMODULE:setupfixed (
                           mode="minimization";
                           atom_select=&atom_select;
                           atom_fixed=&atom_fixed;
                           atom_total_fixed=store1;
                           atom_multiplicity=rmsd;
                           )
 
 fix selection=( store1 ) end


 fastnb grid end

 flags                                       
    include xref                   
   ?                                        
 end      
 if ( &BLANK%force_field_infile = true ) then
    flags                                       
       exclude elec pele include vdw pvdw                 
      ?                                        
    end      
 else
     @&force_field_infile
 end if                                   

 show sum(1) (&atom_harm)
 if ( $result > 0 ) then
   evaluate ($harmonic=true)
 else
   evaluate ($harmonic=false)
 end if

 xray
    predict
      mode=reciprocal
      to=fcalc
      selection=( ref_active=1 )
      atomselection=( &atom_select )
    end
    @@CNS_XTALMODULE:calculate_r (fobs=&STRIP%obs_f;
                                 fcalc=fcalc;
                                 fpart=fbulk;
                                 sel=(ref_active=1);
                                 sel_test=(tst_active=1);
                                 print=true;
                                 output=OUTPUT;
                                 r=$start_r;
                                 test_r=$start_test_r;)
 end


 {- BEGIN MODIFICATION -}
 @CNS_XTALMODULE:scale_and_solvent_grid_search (
                             bscale=&bscale;
                             sel=( ref_active=1 );
                             sel_test=( tst_active=1 );
                             atom_select=( &atom_select );
                             bulk_sol=&bulk_sol;
                             bulk_mask=&bulk_mask_infile;
                             bulk_atoms=( &atom_select );
                             
                             sol_auto=&sol_auto;
                             sol_k=&sol_k;
                             sol_b=&sol_b;
                             sol_rad=&sol_rad;
                             sol_shrink=&sol_shrink;
  
                             fcalc=fcalc;
                             obs_f=&STRIP%obs_f;
                             obs_sigf=&STRIP%obs_sigf;
                             obs_i=$STRIP%obs_i;
                             obs_sigi=$STRIP%obs_sigi;                             
                             fpart=fbulk;
                             
!
! Begin modification (6/28/06)                             
                             Baniso_11=$Baniso_11;
                             Baniso_22=$Baniso_22;
                             Baniso_33=$Baniso_33;
                             Baniso_12=$Baniso_12;
                             Baniso_13=$Baniso_13;
                             Baniso_23=$Baniso_23;
                             Biso=$Biso_model;
! End modification
! 
                             
                             sol_k_best=$sol_k_ref;
                             sol_b_best=$sol_b_ref;
			     solrad_best=$solrad_best;
			     shrink_best=$shrink_best;
                             
                             b=b;

                             low_b_flag=$low_b_flag;
                            
                             sol_output=&sol_output;
                             
                             )

 {- check the gridding again since the minimum B-factor may have changed -}
 show min ( b ) ( &atom_select )
 evaluate ($b_min=$result)
 @@CNS_XTALMODULE:fft_parameter_check ( 
                             d_min=&high_res; 
                             b_min=$b_min;
                             grid=auto;
                             fft_memory=&fft_memory;
                             fft_grid=$fft_grid;   
                             fft_b_add=$fft_b_add; 
                             fft_elim=$fft_elim; 
                                      )
 {- END MODIFICATION -}

 if ( $harmonic = true ) then
   do (refx=x) (all)
   do (refy=y) (all)
   do (refz=z) (all)
   do (harm=0) (all)
   do (harm=&k_harmonic) (&atom_harm)
   flags include harm end
 end if

 xray
   @@CNS_XTALMODULE:refinementtarget (target=&reftarget;
                                     sig_sigacv=0.07;
                                     mbins=&target_bins;
                                     fobs=&STRIP%obs_f;
                                     sigma=&STRIP%obs_sigf;
                                     weight=$STRIP%obs_w;
                                     iobs=$STRIP%obs_i;
                                     sigi=$STRIP%obs_sigi;
                                     test=tst_active;
                                     fcalc=fcalc;
                                     fpart=fbulk;
                                     pa=&STRIP%obs_pa;
                                     pb=&STRIP%obs_pb;
                                     pc=&STRIP%obs_pc;
                                     pd=&STRIP%obs_pd;
                                     phase=&STRIP%obs_phase;
                                     fom=&STRIP%obs_fom;
                                     sel=(ref_active=1);
                                     sel_test=(tst_active=1);
                                     statistics=true;)
 end

 if ( &wa < 0 ) then
   @@CNS_XTALMODULE:getweight (
                              selected=&atom_select;
                              fixed=(store1);
                             )
 end if

 xray
   @@CNS_XTALMODULE:definemonitor (monitor=&reftarget;
                                   fobs=&STRIP%obs_f;
                                   fcalc=fcalc;
                                   fpart=fbulk;
                                   phase=&STRIP%obs_phase;
                                   monitortype=$monitor_is;)
 end

 energy end
 evaluate ($curr_moni=$monitor) 
 evaluate ($test_curr_moni=$test_monitor)

 xray
    predict
      mode=reciprocal
      to=fcalc
      selection=( ref_active=1 )
      atomselection=( &atom_select )
    end
    @@CNS_XTALMODULE:calculate_r (fobs=&STRIP%obs_f;
                                 fcalc=fcalc;
                                 fpart=fbulk;
                                 sel=(ref_active=1);
                                 sel_test=(tst_active=1);
                                 print=true;
                                 output=OUTPUT;
                                 r=$full_r;
                                 test_r=$full_test_r;)
 end

 evaluate ($1=1)
 evaluate ($2=2)
 evaluate ($3=3)
 evaluate ($4=4)
 evaluate ($rmsd=RMSD)

 {- modification, allowing both NCS restraints and constraints, ATB, 12/20/08 -}
 
 if ( &BLANK%ncs_infile = false ) then

  eval ($ncs_strict=false)
  eval ($ncs_restrain=false)
  if ( &ncs_type = "strict" ) then
   eval ($ncs_strict=true) 
  elseif ( &ncs_type = "restrain" ) then
   eval ($ncs_restrain=true)
  elseif ( &ncs_type = "both" ) then
   eval ($ncs_strict=true) 
   eval ($ncs_restrain=true) 
  else
   display >>> Error: unknown NCS type, aborting
   abort
  end if

  if ($ncs_restrain = true) then
   ncs restrain ? end
   evaluate ($ngroup=1)
   evaluate ($done=false)
   while ( $done = false ) loop group
     if ( $exist_rot_$ngroup_$1_$1_$1 # true ) then
       evaluate ($done=true)
       evaluate ($ngroup=$ngroup-1)
     else
       evaluate ($ngroup=$ngroup+1)
     end if
   end loop group
  end if

  if ($ncs_strict = true) then
   ncs strict ? end
   evaluate ($num_op=1)
   evaluate ($done=false)
   while ( $done = false ) loop ncsop
     if ( $exist_ncsop_$num_op_$1_$1 # true ) then
       evaluate ($done=true)
       evaluate ($num_op=$num_op-1)
     else
       evaluate ($num_op=$num_op+1)
     end if
   end loop ncsop
   evaluate ($num_op_strict=$num_op)
  end if
 
  if ($ncs_restrain = true) then
   ncs restraint ? end
   evaluate ($group=1)
   while ($group <= $ngroup) loop group
     evaluate ($num_op=1)
     evaluate ($done=false)
     while ( $done = false ) loop ncsop
       if ( $exist_rot_$group_$num_op_$1_$1 # true ) then
         evaluate ($done=true)
         evaluate ($num_op=$num_op-1)
       else
         evaluate ($num_op=$num_op+1)
       end if
     end loop ncsop
     evaluate ($num_op_$group=$num_op)
     evaluate ($group=$group+1)
   end loop group
  end if
 
 end if
 {- end modification -}

 xray
   if ( $total_test > 0 ) then
     show sum (1) (amplitude(&STRIP%obs_f) > 0 and
                   &high_res <= d <= &low_res and tst_active=1)
     evaluate ($error_bins=int(max(10,$result/25)))
     evaluate ($error_bins=min($error_bins,50))
   else
     show sum (1) (amplitude(&STRIP%obs_f) > 0 and
                   &high_res <= d <= &low_res)
     evaluate ($error_bins=int(max(10,$result/250)))
     evaluate ($error_bins=min($error_bins,50))
   end if
 end

 xray
   if ( &BLANK%luzzati_error_plot = false ) then
     evaluate ($plotfile=&luzzati_error_plot + ".plot")
     evaluate ($mess_size=LONG)
   else
     evaluate ($plotfile=OUTPUT)
     evaluate ($mess_size=SHORT)
   end if
   @CNS_XTALMODULE:luzzaticoorderr ( bins=$error_bins;
                                     fobs=&STRIP%obs_f;
                                     fcalc=fcalc;
                                     fpart=fbulk;
                                     low=&low_err_res;
                                     high=&high_res;
                                     sel=(ref_active=1);
                                     disp=$plotfile;
                                     mess=$mess_size;
                                     esderr_luz=$esderr_luz; )

   if ( $total_test > 0 ) then
     if ( &BLANK%luzzati_error_plot = false ) then
       evaluate ($plotfile=&luzzati_error_plot + "_cv.plot")
       evaluate ($mess_size=LONG)
     else
       evaluate ($plotfile=OUTPUT)
       evaluate ($mess_size=SHORT)
     end if
     @CNS_XTALMODULE:luzzaticoorderr ( bins=$error_bins;
                                       fobs=&STRIP%obs_f;
                                       fcalc=fcalc;
                                       fpart=fbulk;
                                       low=&low_err_res;
                                       high=&high_res;
                                       sel=(ref_active=1 and tst_active=1);
                                       disp=$plotfile;
                                       mess=$mess_size;
                                       esderr_luz=$esderr_luz_cv; )
   end if
 end

 xray
   if ( &BLANK%sigmaa_error_plot = false ) then
     evaluate ($plotfile=&sigmaa_error_plot + ".plot")
     evaluate ($mess_size=LONG)
   else
     evaluate ($plotfile=OUTPUT)
     evaluate ($mess_size=SHORT)
   end if
   @CNS_XTALMODULE:sigmaacoorderr ( bins=$error_bins;
                                    fobs=&STRIP%obs_f;
                                    fcalc=fcalc;
                                    fpart=fbulk;
                                    sel=(ref_active=1 and d <= &low_err_res);
                                    disp=$plotfile;
                                    mess=$mess_size;
                                    esderr_sigmaa=$esderr_sigmaa; )

   if ( $total_test > 0 ) then
     if ( &BLANK%sigmaa_error_plot = false ) then
       evaluate ($plotfile=&sigmaa_error_plot + "_cv.plot")
       evaluate ($mess_size=LONG)
     else
       evaluate ($plotfile=OUTPUT)
       evaluate ($mess_size=SHORT)
     end if
     @CNS_XTALMODULE:sigmaacoorderr ( bins=$error_bins;
                                      fobs=&STRIP%obs_f;
                                      fcalc=fcalc;
                                      fpart=fbulk;
                                      sel=(ref_active=1 and tst_active=1 and
                                           d <= &low_err_res);
                                      disp=$plotfile;
                                      mess=$mess_size;
                                      esderr_sigmaa=$esderr_sigmaa_cv; )
   end if
 end

 print threshold=&bond_thresh bond 
 evaluate ($rmsd_bond=$result)
 evaluate ($viol_bond=$violations)

 print threshold=&angle_thresh angle 
 evaluate ($rmsd_angle=$result)
 evaluate ($viol_angle=$violations)

 print threshold=&dihe_thresh dihedral 
 evaluate ($rmsd_dihe=$result)
 evaluate ($viol_dihe=$violations)

 print threshold=&impr_thresh improper
 evaluate ($rmsd_impr=$result)
 evaluate ($viol_impr=$violations)

 xray
   optimize bfactors

     bmin=0
     bmax=500

     nstep=0
     drop=10.0

     bsigma=( &atom_select and &atom_main )=1.5
     bsigma=( &atom_select and not(&atom_main) )=2.0

     asigma=( &atom_select and &atom_main )=2.0
     asigma=( &atom_select and not(&atom_main) )=2.5

     rweight=-1
   end
 end

 show ave(b) (&atom_select)
 evaluate ($b_average=$result)
 show min(b) (&atom_select)
 evaluate ($b_min=$result)
 show max(b) (&atom_select)
 evaluate ($b_max=$result)

 set display=&list_outfile end

 display ==============================================================================
 display
 display >>> input coordinates: &STRIP%coordinate_infile

 evaluate ($counter=1)
 evaluate ($done=false)
 while ( $done = false ) loop read
  if ( &exist_parameter_infile_$counter = true ) then
    if ( &BLANK%parameter_infile_$counter = false ) then
      display >>> parameter file $counter  : &STRIP%parameter_infile_$counter
    end if
  else
   evaluate ($done=true)
  end if
  evaluate ($counter=$counter+1)
 end loop read
 
 if ( &BLANK%structure_infile = true ) then
   display >>> molecular structure file: automatic

   evaluate ($counter=1)
   evaluate ($done=false)
   while ( $done = false ) loop read
    if ( &exist_topology_infile_$counter = true ) then
      if ( &BLANK%topology_infile_$counter = false ) then
        display >>> topology file $counter  : &STRIP%topology_infile_$counter
      end if
    else
     evaluate ($done=true)
    end if
    evaluate ($counter=$counter+1)
   end loop read

   evaluate ($counter=1)
   evaluate ($done=false)
   while ( $done = false ) loop read
    if ( &exist_link_infile_$counter = true ) then
      if ( &BLANK%link_infile_$counter = false ) then
        display >>> linkage file $counter  : &STRIP%link_infile_$counter
      end if
    else
     evaluate ($done=true)
    end if
    evaluate ($counter=$counter+1)
   end loop read

   if ( &BLANK%patch_infile = false ) then
      display >>> custom patch file = &STRIP%patch_infile
   end if

 else
   display >>> molecular structure file: &STRIP%structure_infile
 end if

 if ( &BLANK%anom_library = false ) then
   display >>> anomalous f' f'' library: &STRIP%anom_library
 end if
 
 evaluate ($counter=1)
 evaluate ($done=false)
 while ( $done = false ) loop read
    if ( &exist_reflection_infile_$counter = true ) then
      if ( &BLANK%reflection_infile_$counter = false ) then
         display >>> reflection file $counter : &STRIP%reflection_infile_$counter       
      end if
   else
     evaluate ($done=true)
   end if
   evaluate ($counter=$counter+1)
 end loop read
 
 display >>> spacegroup: &STRIP%sg
 display >>> cell dimensions: a= &a b= &b c= &c alpha= &alpha beta= &beta gamma= &gamma
 xray wa=? end
 evaluate ($wa_print=$result)
 display >>> current wa= $wa_print for target= &STRIP%reftarget
 if ( &BLANK%restraints_infile = false ) then
   display >>> additional restraints file: &STRIP%restraints_infile
 end if
 if ( &BLANK%ncs_infile = false ) then
   display >>> ncs= &STRIP%ncs_type  ncs file= &STRIP%ncs_infile
 else
   display >>> ncs= none
 end if
 ! 
 ! Begin modification (6/28/06)
 if ( &bscale = "anisotropic" ) then
   display >>> Anisotropic B-factor tensor Ucart of atomic model without isotropic component :
   display >>>   B11=$Baniso_11[f8.3] B22=$Baniso_22[f8.3] B33=$Baniso_33[f8.3]
   display >>>   B12=$Baniso_12[f8.3] B13=$Baniso_13[f8.3] B23=$Baniso_23[f8.3]
   display >>> Isotropic component added to coordinate array B: $Biso_model[f8.3]
 elseif ( &bscale = "isotropic" ) then
   display >>> B-factor applied to coordinate array B: $Biso_model[f8.3]
 else
   display >>> initial B-factor correction: none
 end if
 ! End modification
 !
 if ( &bscale # "no" ) then
   if ( $low_b_flag = true ) then
     display >>> warning: B-correction gave atomic B-values less than zero
     display >>>          they have been reset to zero
   end if
 end if
 {- MODIFIED 5/18/05 -}
 if ( &bulk_sol = true ) then 
   display >>> bulk solvent: probe radius=$solrad_best, shrink value=$solrad_best
   display >>> bulk solvent: density level= $sol_k_ref e/A^3, B-factor= $sol_b_ref A^2
 else
   display >>> bulk solvent: false
 end if
 {- END MODIFICATION -}

 display
 display =================================== summary ==================================
 display

 display resolution range: &low_res - &high_res A 
 display   R-values:
 if ( $total_test > 0 ) then
   display   initial                                        r= $start_r[f6.4] free_r= $start_test_r[f6.4]
   display   after B-factor and/or bulk solvent correction  r= $full_r[f6.4] free_r= $full_test_r[f6.4]
 else
   display   initial                                        r= $start_r[f6.4]
   display   after B-factor and/or bulk solvent correction  r= $full_r[f6.4]
 end if
 display
 display   Monitor for target &reftarget is $monitor_is :
 if ( $total_test > 0 ) then
   display     working set= $curr_moni[f6.4]  test set= $test_curr_moni[f6.4]
 else
   display     working set= $curr_moni[f6.4]
 end if
 display
 display                 luzzati coordinate error (&low_err_res - &high_res A ): $esderr_luz[f6.2] A
 if ( $total_test > 0 ) then
   display cross-validated luzzati coordinate error (&low_err_res - &high_res A ): $esderr_luz_cv[f6.2] A
 end if
 if ( $esderr_sigmaa < 0 ) then
   display NB: Negative coordinate error estimates are sometimes obtained
   display     from sigmaa values - especially at the end of refinement.
   display     Manual interpretation of the sigmaa plot is suggested
 end if
 display                  sigmaa coordinate error (&low_err_res - &high_res A ): $esderr_sigmaa[f6.2] A
 if ( $total_test > 0 ) then
   if ( $esderr_sigmaa_cv < 0 ) then
     display NB: Negative coordinate error estimates are sometimes obtained
     display     from sigmaa values - especially at the end of refinement.
     display     Manual interpretation of the sigmaa plot is suggested
   end if
   display  cross-validated sigmaa coordinate error (&low_err_res - &high_res A ): $esderr_sigmaa_cv[f6.2] A
 end if
 display
 display rmsd bonds= $rmsd_bond[f8.6] with $viol_bond bond violations > &bond_thresh
 display rmsd angles= $rmsd_angle[f8.5] with $viol_angle angle violations >  &angle_thresh
 display rmsd dihedrals= $rmsd_dihe[f8.5] with $viol_dihe angle violations >  &dihe_thresh
 display rmsd improper= $rmsd_impr[f8.5] with $viol_impr angle violations >  &impr_thresh

 display
 display ================================== B-factors =================================
 display

 display average B-factor= $b_average
 display minimum B-factor= $b_min
 display maximum B-factor= $b_max
 display B rmsd for bonded mainchain atoms= $brms_bond_1[f6.3]
 if ( $exist_brms_bond_2 = true ) then
   display B rmsd for bonded sidechain atoms= $brms_bond_2[f6.3]
 end if
 display B rmsd for angle mainchain atoms= $brms_angl_1[f6.3]
 if ( $exist_brms_angl_2 = true ) then
   display B rmsd for angle sidechain atoms= $brms_angl_2[f6.3]
 end if

 display
 display ================================ diffraction data ============================
 display

 if ( &obs_type = "intensity" ) then
   display reflections with Iobs/sigma_I < &sigma_cut rejected
   display reflections with Iobs > &obs_rms * rms(Iobs) rejected
 else
   display reflections with |Fobs|/sigma_F < &sigma_cut rejected
   display reflections with |Fobs| > &obs_rms * rms(Fobs) rejected
 end if
 xray anomalous=? end
 if ( $result = true ) then
   display anomalous diffraction data was input
 end if
 
 {- MODIFIED 2/15/06 -}
 display fft gridding factor = $fft_grid, B factor offset = $fft_b_add A^2, Elimit = $fft_elim
 {- END MODIFICATION -}
 
 display theoretical total number of refl. in resol. range:    $total_theor[I6] ( 100.0 % )
 display number of unobserved reflections (no entry):          $unobserved[I6] ( $per_unobs[f5.1] % )
 display number of reflections rejected:                       $rejected[I6] ( $per_reject[f5.1] % )
 display total number of reflections used:                     $total_used[I6] ( $per_used[f5.1] % )
 display number of reflections in working set:                 $total_work[I6] ( $per_work[f5.1] % )
 display number of reflections in test set:                    $total_test[I6] ( $per_test[f5.1] % )

 display
 display =======> completeness
 if ( $total_test > 0 ) then
   display
   display  Test set (&STRIP%test_set = &test_flag):
   display

   xray
     mbins=&print_bins
     statistics
       completeness
       selection=(ref_active=1 and tst_active=1)
       output=&list_outfile
     end
   end
 end if

 display
 display  Working set:
 display

 xray
   mbins=&print_bins
   statistics
     completeness
     selection=(ref_active=1 and tst_active=0)
     output=&list_outfile
   end
 end

 display 
 display ================================= R-values ===================================
 display

 set print=&list_outfile end

 display
 display =======> R-values with |Fobs|/sigma cutoff= &sigma_cut
 if ( $total_test > 0 ) then
   display 
   display  Test set (&STRIP%test_set = &test_flag):
   display 
 
   xray
     mbins=&print_bins
     statistics
           (rvalue(&STRIP%obs_f,fcalc+fbulk))
           selection=(ref_active=1 and tst_active=1)
           output=&list_outfile
     end
   end
 end if

 display
 display  Working set:
 display

 xray
   mbins=&print_bins
   statistics
         (rvalue(&STRIP%obs_f,fcalc+fbulk))
         selection=(ref_active=1 and tst_active=0)
         output=&list_outfile
   end
 end

 set print=OUTPUT end

 display 
 display =============================== sigmaa-values ================================
 display

 xray
   declare name=x_eobs     domain=reciprocal type=real end
   declare name=x_ecalc    domain=reciprocal type=real end
   declare name=x_sigmaa   domain=reciprocal type=real end

   show sum (1) (tst_active=1)
   if ( $result > 0 ) then
     evaluate ($test_exist=true)
   else
     evaluate ($test_exist=false)
   end if

   if ( $test_exist=true ) then
     display sigmaa calculated using cross-validated data (test set)
   else
     display sigmaa calculated using all data
   end if

   mbins=?
   evaluate ($old_bins=$result)
   if ( &target_bins < 0 ) then
     if ( $test_exist = true ) then
       show sum (1) (tst_active=1)
     else
       show sum (1) (ref_active=1)
     end if
     evaluate ($x_mbins=int($result/50))
     if ( $test_exist = true ) then
       if ( $x_mbins <= 0 ) then
         display Error: there are less than 50 test set reflections
         display        please check the test set and the resolution limits
         abort
       end if
       if ( $x_mbins < 10 ) then
         display
         display Warning: there are less than 50 reflections per bin
         display          you may want to increase the size of the test set
         display
       end if
     end if
     evaluate ($x_mbins=max(10,$x_mbins))
     evaluate ($x_mbins=min($x_mbins,50))
     mbins=$x_mbins
   else
     evaluate ($x_mbin=&target_bins)
     mbins=&target_bins
   end if

   display number of bins for sigmaa calculation= $x_mbins

   if ( $test_exist=true ) then
     do (x_eobs=norm(amplitude(&STRIP%obs_f)))              ( tst_active=1 )
     do (x_ecalc=norm(amplitude(fcalc+fbulk)))              ( tst_active=1 )
     do (x_sigmaa=sigacv[sigma=0.07](x_eobs,x_ecalc))       ( tst_active=1 )
     do (x_sigmaa=distribute(x_sigmaa))                     ( ref_active=1 )
   else
     do (x_eobs=norm(amplitude(&STRIP%obs_f)))              ( ref_active=1 )
     do (x_ecalc=norm(amplitude(fcalc+fbulk)))              ( ref_active=1 )
     do (x_sigmaa=siga(x_eobs,x_ecalc))                     ( ref_active=1 )
   end if

   display

   statistics
     (x_sigmaa)
     selection=( ref_active=1 )
     output=&list_outfile
   end

   mbins=$old_bins

   undeclare name=x_eobs     domain=reciprocal end
   undeclare name=x_ecalc    domain=reciprocal end
   undeclare name=x_sigmaa   domain=reciprocal end

 end

 {- modification, allowing both NCS constraints and restraints, ATB, 12/20/08 -}
 if ( &BLANK%ncs_infile = false ) then
     
   if ($ncs_strict = true) then
   
      display 
      display ============= non-crystallographic symmetry constraints (NCS strict) ============
      display 
      evaluate ($ncsop=1)
      while ($ncsop <= $num_op_strict) loop ncsop
   
         display     NCS operator $ncsop ( 1 -> $ncsop ):
         display       matrix= $ncsop_$ncsop_$1_$1[f8.5] $ncsop_$ncsop_$1_$2[f8.5] $ncsop_$ncsop_$1_$3[f8.5]
         display               $ncsop_$ncsop_$2_$1[f8.5] $ncsop_$ncsop_$2_$2[f8.5] $ncsop_$ncsop_$2_$3[f8.5]
         display               $ncsop_$ncsop_$3_$1[f8.5] $ncsop_$ncsop_$3_$2[f8.5] $ncsop_$ncsop_$3_$3[f8.5]
         display       translation= $ncsop_$ncsop_$1_$4[f10.5] $ncsop_$ncsop_$2_$4[f10.5] $ncsop_$ncsop_$3_$4[f10.5]
         display       rms difference= $ncsop_$ncsop_$rmsd[f8.5]
      end loop ncsop
   end if
       
   if ($ncs_restrain = true) then
      display 
      display ============ non-crystallographic symmetry restraints (NCS restraints) ===========
      display 
      display >>> number of NCS groups= $ngroup
      evaluate ($group=1)
      while ($group <= $ngroup) loop group
        display >>> NCS group $group: number of NCS operators= $num_op_$group
        evaluate ($ncsop=1)
        while ($ncsop <= $num_op_$group) loop ncsop
       
         display     NCS operator $ncsop ( $ncsop -> 1 ):
         display       matrix= $rot_$group_$ncsop_$1_$1[f8.5] $rot_$group_$ncsop_$1_$2[f8.5] $rot_$group_$ncsop_$1_$3[f8.5]
         display               $rot_$group_$ncsop_$2_$1[f8.5] $rot_$group_$ncsop_$2_$2[f8.5] $rot_$group_$ncsop_$2_$3[f8.5]
         display               $rot_$group_$ncsop_$3_$1[f8.5] $rot_$group_$ncsop_$3_$2[f8.5] $rot_$group_$ncsop_$3_$3[f8.5]
         display       translation= $rot_$group_$ncsop_$1_$4[f10.5] $rot_$group_$ncsop_$2_$4[f10.5] $rot_$group_$ncsop_$3_$4[f10.5]
         display       rms difference= $rot_$group_$ncsop_$rmsd[f8.5]
         evaluate ($ncsop=$ncsop+1)
       end loop ncsop
       evaluate ($group=$group+1)
     end loop group
   end if
   
 end if
 
 {- end modification -}

 display 
 display ============================ non-trans peptides ==============================
 display 

 evaluate ($first=true)

 for $id in id ( &atom_select and name ca ) loop cis

   show (segid) (id $id)
   evaluate ($segid=$result)
   show (resid) (id $id)
   evaluate ($resid=$result)
   show (resname) (id $id)
   evaluate ($resname=$result)

   identity (store2) ( &atom_select and 
                     ( name c and bondedto ( name n and 
                       resid $resid and segid $segid ) ) )
   if ( $select = 1 ) then
     show element (store2) ( attribute store2 > 0 )
     evaluate ($id_prev=$result)
     show (segid) (id $id_prev)
     evaluate ($segid_prev=$result)
     show (resid) (id $id_prev)
     evaluate ($resid_prev=$result)
     show (resname) (id $id_prev)
     evaluate ($resname_prev=$result)
 
     evaluate ($result=180)  ! always define $result in case one or more of the selections below are zero
     pick dihedral
       (&atom_select and name ca and segid $segid_prev and resid $resid_prev)
       (&atom_select and name  c and segid $segid_prev and resid $resid_prev)
       (&atom_select and name  n and segid $segid and resid $resid)
       (&atom_select and name ca and segid $segid and resid $resid)
       geometry

     evaluate ($dihedral=mod($result+360,360))

     if ( $dihedral > 180 ) then
       evaluate ($dihedral=$dihedral-360)
     end if

     evaluate ($absdihedral=abs($dihedral))

     if ( $absdihedral < 15 ) then
       if ( $first = true ) then
         evaluate ($first=false)
       end if
       display cis-peptide: segid=$segid resid=$resid resname=$resname  
       display              current dihedral value= $dihedral[f8.3]
       display 
     elseif ( $absdihedral < 165 ) then
       if ( $first = true ) then
         evaluate ($first=false)
       end if
       display distorted peptide: segid= $segid resid= $resid resname=$resname
       display                    current dihedral value= $dihedral[f8.3]
       display >>>> this may require correction before starting refinement
       display 
     end if
   end if

 end loop cis

 if ( $first = true ) then
  display there are no distorted or cis- peptide planes
 end if

 display 
 display =============================== occupancies ==================================
 display 

 show sum(1) (&atom_select and (attr q = 0))
 if ( $result > 0 ) then
   display $result atoms have zero occupancy =>
   display
   for $id in id ( &atom_select and (attr q = 0) ) loop occ
     show (segid) (id $id)
     evaluate ($segid=$result)
     show (resid) (id $id)
     evaluate ($resid=$result)
     show (resname) (id $id)
     evaluate ($resname=$result)
     show (name) (id $id)
     evaluate ($name=$result)
     if ( $BLANK%segid = true ) then
       display  residue= $resname $resid  atom name= $name
     else
       display  segid= $segid  residue= $resname $resid  atom name= $name
     end if
   end loop occ
 else
   display no atoms have zero occupancy
 end if

 display 
 display ================================= geometry ===================================
 display 
 set print=&list_outfile end

 display =======> bond violations
 display
 print threshold=&bond_thresh bond 

 display
 display =======> angle violations
 display

 print threshold=&angle_thresh angle 

 display
 display =======> improper angle violations
 display

 print threshold=&impr_thresh improper

 display
 display =======> dihedral angle violations
 display

 print threshold=&dihe_thresh dihedral

 set print=OUTPUT end

 display 
 display ============================== close contacts ================================
 display 

 distance
   from=(&atom_select) to=(&atom_select)
   cutoff=&close
   output=&list_outfile
 end

 display 
 display ============================= crystal contacts ===============================
 display 

 flags exclude vdw end

 distance
   from=(&atom_select) to=(&atom_select)
   cutoff=3.5
   output=&list_outfile
 end


 display 
 display =============================== occupied volume ==============================
 display 

 xray
    declare name=mask domain=real end
    mask
      mode=vdw
      solrad=1.0
      shrink=1.0
      nshell=1
      to=mask
      selection=( &atom_select )
    end
 end
 evaluate ($percen_pack=$inside * 100 ) 

 display  unitcell volume occupied by selected atoms = $percen_pack[F8.2] %
 display 
 display ==============================================================================

 stop
'''
    out='model_stat.inp'
    fw=open(out, "w")
    fw.write ("%s" % inp)
    fw.close()

    return out

    
