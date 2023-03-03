#=================================================================
#  Calculate the Matthew coefficient and solvent content.
#  
#=================================================================
import os,sys,math
import util

##################################################################
def matthew_coeff_atom(file):
    '''calculate Matthew_coeff and solvent content: file is in pdb format;
       A simplyfied version : only use atom mass in cell.
    '''

    if not util.check_file(100, file):
        print('Error: file (%s) not exist' %file)
        return 0, 0

    fp=open(file, 'r')

    cell=[1,1,1,1,1,1]
    amass, nop, nmat, sym =0, 0, 1, 'X'

    for x in fp:
        if 'REMARK 290 ' in x[:12] and '555' in x and ',' in x[23:32]:
            nop=nop+1
        elif 'MTRIX3' in x[:6]  and '1' not in x[55:].strip() :
            nmat  = nmat +1
        elif 'CRYST1' in x[:6]:
            c=x[7:54].split()
            cell=[float(y) for y in c]
            sym=x[54:65].strip()

        elif ('ATOM' in x[:4]) : 
            occ=float(x[54:60])
            amass = amass + util.atom_mass(x[76:78].strip())*occ

    fp.close()

    cell_vol=cell_volume(cell)
    
    nsym=util.sg_nsym(sym)
    if nsym ==-1:
        print('Error: space group (%s) is not in the list (%s)' %(sym, file))
        nsym=nop

    matt, solv =   calc_matt(cell_vol, amass, nsym, nmat, 1) #by atom, occ

    return matt, solv  #only use mass of atom

##################################################################
def matthew_coeff(file):
    '''calculate Matthew_coeff and solven content: file is in pdb format;
    '''

    if not util.check_file(100, file):
        print('Error: file (%s) not exist' %file)
        return 0, 0

    fp=open(file, 'r')
    
    cell=[1,1,1,1,1,1]
    spt, nop, nmat, sym, res, atom, res= 1, 0, 1, 'X', [], [], []
    rmass, amass, armass = 0, 0, 0
    hetres, aname, rest = [],[], ''
    for x in fp:
        if 'REMARK 290 ' in x[:12] and '555' in x and ',' in x[23:32]:
            nop=nop+1
        elif 'SEQRES' in x[:6] :
            t=x[17:79].split()
            res.extend(t)
            
#        elif 'SPLIT' in x[:6] :
#            t=x[6:].split()
#            spt=len(t)
            
        elif 'MTRIX3' in x[:6]  and '1' not in x[55:].strip() :
            nmat  = nmat +1
        elif 'CRYST1' in x[:6]:
            c=x[7:54].split()
            cell=[float(y) for y in c]
            sym=x[54:65].strip()

        elif ('ATOM' in x[:4]): 
            
            atom.append(x)
            occ=float(x[54:60])
            amass = amass + util.atom_mass(x[76:78].strip())*occ
            #amass = amass + atom_mass(x[76:78].strip())
            t=x[17:27] #comp_ch_res_int
            aname.append(x[76:78].strip())
            
            if t!=rest :
              
                comp=t[:3].strip()
                restmp=util.residue_mass(comp)
                if restmp<1:
                    hetres.append(comp)
                else:
                    armass=armass + util.residue_mass(comp)
                    
                rest=t
                aname=[]
                

        elif 'ENDMDL' in x[:6]:
            break
    fp.close()

    cell_vol=cell_volume(cell)

    nsym=util.sg_nsym(sym)
    if nsym ==-1:
        print('Error: space group (%s) is not in the list (%s)' %(sym, file))
        nsym=nop
#----------

    for x in hetres: armass=armass + non_standard_res(x,atom)

    for x in res:
        resm=util.residue_mass(x)
        if resm<1:
            m1=non_standard_res(x,atom)
            rmass = rmass + m1
        else:
            rmass = rmass + resm
            
    amatt, asolv =   calc_matt(cell_vol, amass, nsym, nmat, spt) #by atom, occ 
    rmatt, rsolv =   calc_matt(cell_vol, rmass, nsym, nmat, spt) #by SEQRES
    armatt, arsolv = calc_matt(cell_vol, armass , nsym, nmat, spt) #residue 

    matt, solv=-1,-1
    
    if 1.8 < rmatt< 5 :
        matt, solv = rmatt, rsolv
    elif 1.8 < armatt< 5 :
        matt, solv = armatt, arsolv
    elif 1.7 < amatt< 5 :
        matt, solv = amatt, asolv
    else:
        matt, solv = armatt, arsolv
        #print('Warning: packing problem (%s), Matthew_coeff=%.2f; Solvent=%.2f' %(file, matt, solv))
    
          
    return matt, solv
##########################################################
def non_standard_res(res,atom):
    '''mass of none standard residues are calculated by each atoms
    in the atom list.
    '''

    amass, n =0, 0
    for i, x in enumerate(atom):
        if res == x[17:20].strip():
            if i>0 and x[22:26] != atom[i-1][22:26] and n>0 : break
            n=n+1
            amass = amass + util.atom_mass(x[76:78].strip())
     
    return amass

##########################################################
def calc_matt(cell_vol, mass, nsym, nmat, spt) :
    '''cell_vol: cell volume; mass: mass in ASU;
    nsym: symmetry operators;  nmat: number of matrix; spt: number of split
    '''
    
    if mass < 1  or cell_vol <2: return 0, 0
    if spt<=0 : spt=1
    if nsym<=0 : nsym=1
    
    mass = spt*nsym*mass
    if nmat>0 : mass=mass*nmat

    matt, solv = 0, 1
    if(cell_vol >10 and mass>10):
        matt = cell_vol/mass;
        solv =(1.0-1.239/matt)*100.0;

    return matt, solv

##########################################################
def cell_volume(cell):
    '''volume of cell.
    '''
       
    alpha = 3.14159 * cell[3]/180;
    beta  = 3.14159 * cell[4]/180;
    gamma = 3.14159 * cell[5]/180;

    cell_volume = cell[0] * cell[1] * cell[2] * \
math.sqrt(1.- math.cos(alpha)*math.cos(alpha) -  \
math.cos(beta)*math.cos(beta) - math.cos(gamma)*math.cos(gamma) \
     + 2.0 * math.cos(alpha) * math.cos(beta) * math.cos(gamma));
        
    return cell_volume;

