#===========================================================
# this module is to do the maps around the ligands/peptide.
#===========================================================

import os,sys, math,shutil
import config, util, tls
import cifparse as cif


##########################################################
def cut_map_around_ligand_peptide(dccfile, dic, mapfile_in, xyzfile_in):
    '''It generate a complete set for ligand (map, html, jmol).
    dccfile: the density file by dcc.
    dic: a directory to hold all the file for webpage (url below).
    mapfile_in: a input map file.
    xyzfile_in: a input coordinate file.
    '''
    
    print('Cutting the density maps for ligands/peptide')
    
    tmpxyz=xyzfile_in
    if util.is_cif(xyzfile_in): tmpxyz= cif.cif2pdb(xyzfile_in)
    pdbfile = os.path.basename(dic['pdbfile']) + '_new'
    if  pdbfile !=tmpxyz : shutil.copy(tmpxyz,pdbfile)     
    
    mapfile=os.path.basename(dic['pdbfile']) + '_2fofc.map'
    if dic['ligmapcif'] : mapfile=dic['xyzfile_orig'] + '_2fofc.map'
    shutil.move(mapfile_in,  mapfile)

    if dic['ligmapcif'] : #pre-parse the cif file.
        dic['cif']=1
        
        ciffile=dic['xyzfile_orig']

        flist=open(ciffile, 'r').readlines()
        cell_items,values = cif.cifparse(flist, '_cell.')
        cell=cif.get_rows(cell_items, values)
        dic['cell_items'], dic['lig_cell']=cell_items, cell 

        sym_items,values = cif.cifparse(flist, '_symmetry.')
        sym=cif.get_rows(sym_items, values)
        dic['sym_items'], dic['lig_sym']=sym_items, sym
        
        items,values = cif.cifparse(flist, '_atom_site.')
        comp=cif.parse_values(items,values,"_atom_site.auth_comp_id")
        asym=cif.parse_values(items,values,"_atom_site.auth_asym_id")
        seq=cif.parse_values(items,values,"_atom_site.auth_seq_id")
        alt=cif.parse_values(items,values,"_atom_site.label_alt_id")
        ins=cif.parse_values(items,values,"_atom_site.pdbx_PDB_ins_code");
        mod=cif.parse_values(items,values,"_atom_site.pdbx_PDB_model_num")
        row=cif.get_rows(items, values)
       
        dic['items'], dic['comp1'],dic['asym'],dic['seq']=items,comp,asym,seq
        dic['alt'],dic['ins'],dic['mod'],dic['row']=alt,ins,mod,row

    fw_itool=open('LIG_PEPTIDE.cif','w') #a cif file contains table, filenames
    fw_itool.write('data_lig_peptide\n')
    fw_itool.write('\n# A "!" will be given if the residue is bad with real_space_R.\n')
    fw_itool.write('\n# Criteria: (CC<0.7 and R>0.4) or CC<0.5 or R>0.5\n')

    url='http://sf-tool.wwpdb.org/users_data/dir_%s/' %dic['dir']
    #url=os.environ['THIS_SERVICE_URL__FIX_ME'] + '/users_data/dir_%s/' %dic['dir']

    ch_pep,chr_pep, ch_lig,chr_lig, ch_wat,chr_wat=tls.chain_res_range(pdbfile)

    ligpdb=non_poly_pdb(ch_pep, ch_lig, pdbfile) #get non-poly xyz file
    dcc=get_dcc(dccfile) #get a list for dcc of each residue

    if not dcc:
        util.perror('Warning: Failed to parse EDS values! No ligand/peptide maps will be generated. ')

    for k, v in ch_pep.items(): 
        if len(v)< 15: #length of peptide
            if not dic['sdsc_map'] : map_around_peptide(fw_itool,dic, mapfile, ligpdb, dcc, ch_pep, url)
            break

    if ch_lig : map_around_ligand(fw_itool,dic, mapfile, ligpdb, dcc, ch_lig, url)

    get_html_table_baddcc_general(mapfile, dcc) #for polymer/lig/peptide
    
    fw_itool.close()

    if dic['sdsc_map'] :
        arg='rm -f %s %s  %s  LIG_PEPTIDE.cif ' %(mapfile,  mapfile_in, xyzfile_in)
        arg=arg + ' %s_rcc_sum.cif.mtz  %s_2fofc.map_all.html ' %(dic['pdbfile'], dic['pdbfile'])
        os.system(arg)

#    util.delete_file(pdbfile) 
    return

##########################################################
def non_poly_pdb(ch_pep, ch_lig, pdbfile):
    '''get the xyz for peptides and ligands
    ch_lig={'ch':[...], };  ch_pep={'ch':[...], }
    '''

    ligpdb=[]
    pdb=open(pdbfile, 'r').readlines()
    
    for x in pdb: #
        if ('CRYST1' in x[:6] or 'SCALE' in x[:5]) :
            ligpdb.append(x)
        if ('ATOM' in x[:4] or 'HETATM' in x[:6]) : break
        
    
    for k, v in ch_pep.items(): #k=chain, v=list of resi number
        if len(v)>= 15: continue    #length of peptide
        for x in pdb:
            if (('ATOM' in x[:4] or 'HETATM' in x[:6]) and
                k==x[20:22].strip() and int(x[22:26]) in v):
                ligpdb.append(x)

    for k, v in ch_lig.items(): 
        for x in pdb:
            if (('ATOM' in x[:4] or 'HETATM' in x[:6]) and
                k==x[20:22].strip() and int(x[22:26]) in v):
                ligpdb.append(x)

    return ligpdb

##########################################################
def map_around_compound(mapfile, coord, compid):
    '''cut the ASU map to the residue level.
    compid:  model_compound_chainID_resnumber_alt_insertion.
    mapfile: the CCP4 map in ASU.
    coord:   the coordinate file for the map (in cif/pdb)
    '''

    if (not util.check_file(100, mapfile) or not util.check_file(100,coord)):
        err='Error: Either mapfile or coordinate file does not exist.'
        config.ERRLOG.append(err)
        print(err)
        return
    
    
    xyzlim,xyzcomp=find_xyzlim_compound(compid, coord) 
        
    if(len(xyzlim.strip())<2):
        err='Error: compound boundary in fraction not found. check with compound ID'
        config.ERRLOG.append(err)
        print(err)
        return

# below is cut map and get jmol
    t=compid.split('_')
    comp='_'.join([t[0], t[1], t[2]])
    mapout = comp + '_cut.map'
    maphtml = comp + '.html'
    
    mapscr=cut_map_bylimit(xyzlim)
    util.delete_file(mapout)
    arg = mapfile + ' ' + ' ' +  mapout
    command="chmod +x %s ; ./%s  " %(mapscr, mapscr) + arg 
    os.system(command)

    min,max,mean,sigma=map_info(mapfile)
    min1,max1,mean1,sigma1=map_info(mapout) 
    cont={'0.5':0.5, '0.7':0.7, '1.0':1.0, '1.5':1.5, '2.0':2.0} #contour 4  map in asu.
    cont1=cont #contour 4 sub map.

    scale=1.0
    if(float(sigma1)>0) :
        scale = float(sigma)/float(sigma1)
        for z in cont.keys(): cont1[z]=cont[z]*scale

    maphtml = get_html4jmol(comp, xyzcomp, mapout, cont1)

    return maphtml , mapout
     
##########################################################
def find_xyzlim_compound(compid, coord):
    '''find xyzlimit used by mapmask, and write the coord in cif or pdb format.
    compid: atom_group_id (model_compound_chainID_resnumber_alter_insertion)
    coord: the coordinate file
    idd = 0, cif format; =1, the pdb format
    '''

    comp='XXXX'
    t1=compid.split(':')
    for i,x in enumerate(t1):
        t=x.split('_')
        if i==0: comp='_'.join([t[0], t[1], t[2], t[3]])
        
        if len(t)!=6:
            print('Error: in group-id (%d). it should be (model_compound_chainID_resnumber_alter_insertion).' %(i+1))
            return '',''

    idd=util.is_cif(coord)
    xyzcomp=comp + '.pdb'
    if idd==1: xyzcomp=comp + '.cif'
    
    fw=open(xyzcomp, 'w')

    border=1 #extend a little to cover more density 
    xx,yy,zz=[],[],[]    
    if idd==1: #input cif format
        fw.write('data_xyzcomp\n#\n')

        flist=open(coord, 'r').readlines()
        items,values = cif.cifparse(flist, '_cell.')
        fw.write('\n#\n')
        for m, p in enumerate (items): fw.write("%s  %s\n" %(p, values[m]))

        cell=cif.get_cell(flist)
        
        items,values = cif.cifparse(flist, '_atom_site.')
        comp=cif.parse_values(items,values,"_atom_site.auth_comp_id");
        asym=cif.parse_values(items,values,"_atom_site.auth_asym_id");
        seq=cif.parse_values(items,values,"_atom_site.auth_seq_id");
        alt=cif.parse_values(items,values,"_atom_site.label_alt_id");
        ins=cif.parse_values(items,values,"_atom_site.pdbx_PDB_ins_code");
        x=cif.parse_values(items,values,"_atom_site.Cartn_x");
        y=cif.parse_values(items,values,"_atom_site.Cartn_y");
        z=cif.parse_values(items,values,"_atom_site.Cartn_z");
        model=cif.parse_values(items,values,"_atom_site.pdbx_PDB_model_num");

        if(not (alt and comp and ins and asym and seq and x and y and z)): 
            print('Error: not enough infor. extraced from (%s). Check ciftokens'%coord)
            sys.exit()
            
        fw.write('\n#\nloop_\n')
        for p in items: fw.write("%s\n" %p)
        row=cif.get_rows(items, values)

        for i in range(len(x)):
            alter, inst, mod ='.', '.', '1'
            if model and util.is_number(model[i]): mod=model[i]
            if alt and alt[i] != '?' : alter = alt[i]
            if ins and ins[i] != '?' : inst = ins[i]
            
            id1='_'.join([mod, comp[i], asym[i], seq[i], alter, inst])
            
            if id1 in compid:
                xx.append(float(x[i]))
                yy.append(float(y[i]))
                zz.append(float(z[i]))
                
                for m in row[i] : fw.write("%s " %m)
                fw.write('\n')

        
    else:  #pdb format
        fp=open(coord,'r')
        for x1 in fp:
                
            if ('CRYST1' in x1[:6] ):
                fw.write(x1)
                cell=[float(p) for p in x1[8:54].split()]
                
            elif ('ATOM' in x1[:4] or 'HETATM' in x1[:6] ):
                alt=x1[16:17]
                if alt.isspace() : alt='.'
                ins=x1[26:27]
                if ins.isspace() : ins='.'
                resname, chid, resnum = x1[17:20].strip(), x1[20:22].strip(), x1[22:26].strip()
                resid='_'.join([resname, chid, resnum, alt, ins])
                
                if resid in compid:
                    fw.write(x1) #only write the selected section 
                    xx.append(float(x1[30:38]))
                    yy.append(float(x1[38:46]))
                    zz.append(float(x1[46:54]))
        fp.close()

        
    if not xx or not yy or not zz :
        print('Error: %s can not be found in the coordinate. try a new id. ' %(compid))
        return '',''
        
    frac,orth=util.frac_orth_matrix(cell) #get matrix
    border=2.0
    xx_min,xx_max =min(xx)-border, max(xx)+border
    yy_min,yy_max =min(yy)-border, max(yy)+border
    zz_min,zz_max =min(zz)-border, max(zz)+border
        
    xf_min = util.matrix_prod(frac,[xx_min,yy_min,zz_min])
    xf_max = util.matrix_prod(frac,[xx_max,yy_max,zz_max])
                
    xyzlim='%.3f %.3f  %.3f %.3f  %.3f %.3f' %(xf_min[0],xf_max[0], xf_min[1],xf_max[1],xf_min[2],xf_max[2])

    fw.close()
    return xyzlim, xyzcomp

##########################################################
def parse_xyz_compound(compid, coord, idd):
    '''generate a file used by mapmask
    '''
    
    compdb=compid.split(':')[0] + '.pdb'
    fw=open(compdb,'w')
    
    if idd==0:
        tmpfile = cif.cif2pdb(coord)
        fp=open(tmpfile,'r')
        os.remove(tmpfile)
    else:
        fp=open(coord,'r')
        
    for x in fp:

        if ('CRYST1' in x[:6] or 'SCALE' in x[:5]):
            fw.write(x)
                
        elif ('ATOM' in x[:4] or 'HETATM' in x[:6] ):
            alt=x[16:17]
            if alt.isspace() : alt='.'
            ins=x[26:27]
            if ins.isspace() : ins='.'
            resname, chid, resnum = x[17:20].strip(), x[20:22].strip(), x[22:26].strip()
            resid='_'.join([resname, chid, resnum, alt, ins])
            if resid in compid: fw.write(x)
                    
                
    fp.close()        
    fw.close()
   
    return compdb

##########################################################
def remove_ligand(pdbfile):
    '''
    '''
    
    newpdb='%s_NOLIG' %pdbfile
    fr=open(pdbfile, 'r')
    fw=open(newpdb, 'w')
    
    ch_pep,chr_pep,ch_lig,chr_lig,ch_wat,chr_wat=tls.chain_res_range(pdbfile)

    for x in fr:
        if (('ATOM' in x[:4] or 'HETA' in x[:4] or 'ANISOU' in  x[:6]) ):
            ch=x[20:22].strip()
            nres=int(x[22:26])
            if ch in ch_lig and nres in ch_lig[ch]: continue
        fw.write(x)

    fw.close()
    fr.close()
#    shutil.move(newpdb,pdbfile)
    return newpdb
        
##########################################################
def get_html_table_baddcc_general(mapfile, dcc):
    '''
    '''
    
    if not dcc: return
    
    html_table = mapfile + '_all.html'
    fw=open(html_table, 'w')
    
    fw.write('<html xmlns="http://www.w3.org/1999/xhtml"> \n ')
    fw.write("<script language='Javascript' type='text/javascript'> </script>\n")
             
    html_table_head(fw, 'bad density correlation/real space R factor', 0)
    for x in dcc:
        
        if ((float(x[4])<0.7 and float(x[5])>0.4) or
            float(x[4])<0.5 or float(x[5])>0.5):
            tmp=x[3]
            if '.' in x[3]: tmp='_'
            idd=x[1] + x[0] + tmp
            ss=[x[2], idd, x[4], x[5], x[6], x[7]]
            s1=' '
            for y in ss: s1= s1 + '<td>%s</td>' %y
            all = '<tr>'  + s1 + '</tr>\n'
            fw.write(all)

            
    fw.write('</TABLE>\n </html>\n')

    ss= '<p><center><a href=\"javascript:window.history.back();\"><b>Go Back</b></a></a></center> \n </html>\n'
#    fw.write(ss)
    
    fw.close()
    return
    
##########################################################
def get_dcc(dccfile):
    '''put dccfile as list of list
    '''

    dcc=[]
    if(util.check_file(100, dccfile)==0): return dcc
    flist=open(dccfile, 'r').readlines()
    items,values = cif.cifparse(flist, '_pdbx_rscc_mapman.')
    nseq=cif.parse_values(items,values,"_pdbx_rscc_mapman.auth_seq_id");
    chid=cif.parse_values(items,values,"_pdbx_rscc_mapman.auth_asym_id");
    comp=cif.parse_values(items,values,"_pdbx_rscc_mapman.auth_comp_id");
    alt=cif.parse_values(items,values,"_pdbx_rscc_mapman.label_alt_id");
    ins=cif.parse_values(items,values,"_pdbx_rscc_mapman.label_ins_code");
    cc=cif.parse_values(items,values,"_pdbx_rscc_mapman.correlation");
    rsr=cif.parse_values(items,values,"_pdbx_rscc_mapman.real_space_R");
    zrsr=cif.parse_values(items,values,"_pdbx_rscc_mapman.real_space_Zscore");
    biso=cif.parse_values(items,values,"_pdbx_rscc_mapman.Biso_mean");
    occ=cif.parse_values(items,values,"_pdbx_rscc_mapman.occupancy_mean");
    #model=cif.parse_values(items,values,"_pdbx_rscc_mapman.model_id");
    pdbid=cif.parse_values(items,values,"_pdbx_rscc_mapman.pdb_id");
    if not items: return dcc 
    for i in range(len(chid)):
        a=[nseq[i], chid[i], comp[i], alt[i], cc[i], rsr[i], biso[i], occ[i],pdbid[i]]
        dcc.append(a)
    return dcc

##########################################################
def html_table_head(fw, title, idd):
    ''' title: name in caption;  idd==0, do not give map
    '''
    
    map1 = '<TH> electron desity map </TH>'
    if idd==0: map1 = ''
    fw.write('<br><TABLE BORDER="1">\n' )
    fw.write('<CAPTION><b><h3> Summary of %s</h3>  </b> </CAPTION> \n' %title )
    fw.write('<TR><TH> Name </TH> <TH> ID </TH> <TH>density <br> correlation </TH> \
    <TH> real_space_R </TH>  <TH> Biso </TH> <TH>Occ</TH> %s </TR> \n' %map1 )
    
##########################################################
def html_table_content(fw, idd, y, url, maphtml):
    '''Write a row in the table of HTML.
    idd: number, y: resid, id, rho, dcc, biso, occ
    '''
    warn=' . '

    if float(y[2])< -9.0:
        t='Warning: resid=%s_%s: No real_spaceR calculated (maybe 0 occupancy)\n' %(y[0], y[1])
        util.perror(t)
    elif ((float(y[2])<0.7 and float(y[3])>0.4) or float(y[2])<0.5 or float(y[3])>0.5):
        t='Warning: resid=%s_%s: bad real_spaceR (%s) or density correlation (%s)\n' %(y[0], y[1], y[3], y[2])
        config.ERRLOG.append(t)
        warn=' ! '
 
            
    ss=' '
    if len(y)<4: return
    for z in y:
        if ((float(y[2])<0.7 and float(y[3])>0.4) or
            float(y[2])<0.5 or float(y[3])>0.5):
            ss= ss + '<td><font color="red">%s</font></td> ' %z
        else:
            ss= ss + '<td>%s</td>' %z

    if idd ==0:
        map1='<A TARGET="_self" HREF="%s%s">view map</A>' %(url,maphtml)
        ss=ss+ '<td>%s</td>' %map1
    else:
        ss=ss+ '<td> </td> '

    all = '<tr>'  + ss + '</tr>\n'
    fw.write(all)

    return warn
#


##########################################################
def get_subpdb(pdb, k, v, idd):
    '''pick the atoms in pdb (a list) for chain k in residue number v (a list)
    '''

    peppdb='%s.pdb' %(idd)
    fwt = open(peppdb, 'w')

    natom=0
    for x in pdb:
        if ('CRYST1' in x[:6] or 'SCALE' in x[:5]):
            fwt.write(x)
        elif (('ATOM' in x[:6] or 'HETATM' in x[:6]) and
              k==x[20:22].strip() and int(x[22:26]) in v):
            fwt.write(x)
            natom=natom+1
       
    fwt.close()
    if natom <2 : util.delete_file(peppdb)  # ion 
    return natom, peppdb

##########################################################
def get_subcif(dic, k, v):
    '''k is chainID;  v is list of residue number
    return idd, natom, ciffile name
    '''
    '''
    ciffile=dic['pdbfile']
    
    flist=open(ciffile, 'r').readlines()

    cell_items,values = cif.cifparse(flist, '_cell.')
    cell=cif.get_rows(cell_items, values)

    sym_items,values = cif.cifparse(flist, '_symmetry.')
    sym=cif.get_rows(sym_items, values)

    
    items,values = cif.cifparse(flist, '_atom_site.')
    comp=cif.parse_values(items,values,"_atom_site.auth_comp_id")
    asym=cif.parse_values(items,values,"_atom_site.auth_asym_id")
    seq=cif.parse_values(items,values,"_atom_site.auth_seq_id")
    alt=cif.parse_values(items,values,"_atom_site.label_alt_id")
    ins=cif.parse_values(items,values,"_atom_site.pdbx_PDB_ins_code");
    mod=cif.parse_values(items,values,"_atom_site.pdbx_PDB_model_num") 
    row=cif.get_rows(items, values)
  
    '''
    cell_items, cell =dic['cell_items'], dic['lig_cell']
    sym_items, sym = dic['sym_items'], dic['lig_sym']
    items,comp,asym,seq,alt,ins,mod,row=dic['items'], dic['comp1'],dic['asym'],dic['seq'],dic['alt'],dic['ins'],dic['mod'],dic['row'] 

    idd, natom, atom='1_X_X_X__', 0, []

    for i in range(len(asym)):
        if asym and asym[i]==k and seq and int(seq[i]) in v:
            natom=natom+1
            #print(i, natom, k, v[0],asym[i], seq[i], comp[i])
            alter, inst = '',''
            if alt and alt[i] != '?' and alt[i] != '.' : alter=alt[i]
            if ins and ins[i] != '?' and ins[i] != '.' : inst=ins[i]
            atom.append(row[i])
            if natom==1 and mod and comp:
                
                idd='_'.join([mod[i], asym[i], comp[i], seq[i], inst, alter])

    
    subpep='%s.cif' %(idd) #write subcif name
    fw = open(subpep, 'w')

    fw.write('data_%s\n#\n' %idd)
    for i, p in enumerate(cell_items):
        if ' ' in cell[0][i]:
            fw.write("%s   '%s'\n" %(p, cell[0][i]))
        else:
            fw.write("%s   %s\n" %(p, cell[0][i]))
     
    fw.write('#\n')
    for i, p in enumerate(sym_items):
        if ' ' in sym[0][i]:
            fw.write("%s   '%s'\n" %(p, sym[0][i]))
        else:
            fw.write("%s   %s\n" %(p, sym[0][i]))

    
    fw.write('\n#\nloop_\n')
    for p in items: fw.write("%s\n" %p)
    for x in atom:
        for m in x:  fw.write("%s " %m)
        fw.write('\n')
  
    fw.close()
    if natom <2 : util.delete_file(subpep)  # ion 
    return idd, natom, subpep
    
##########################################################
def map_around_peptide(fw_itool,dic, mapfile, pdb, dcc, ch_pep, url):
    '''mapfile: the map in cell; pdb: a list;  dcc: a list of list; 
    ch_pep: a dict for peptide {'ch': [n1,n2 ..]}; url: the url for html
    '''
    
    if not ch_pep: return

    cif1='''
loop_
_dcc_peptide.id
_dcc_peptide.residue_name
_dcc_peptide.chain_id
_dcc_peptide.dcc_correlation 
_dcc_peptide.real_space_R
_dcc_peptide.Biso_mean 
_dcc_peptide.occupancy_mean
_dcc_peptide.warning
_dcc_peptide.file_name_map_html
_dcc_peptide.file_name_pdb
_dcc_peptide.file_name_map
_dcc_peptide.file_name_jmol
_dcc_peptide.full_map_sigma
_dcc_peptide.sub_map_sigma
'''
    fw_itool.write(cif1)
    
    html_table = mapfile + '_pep.html'
    fw=open(html_table, 'w')
    html_table_head(fw, 'Peptides (or nucleic acid)', 1)

    min, max,mean,sigma=map_info(mapfile)
    cont={'0.5':0.5, '0.7':0.7, '1.0':1.0, '1.5':1.5, '2.0':2.0} #contour 4  map in asu.
    cont1=cont #contour 4 sub map.
    
    pepid=0
    for k, v in ch_pep.items():  #k: chainID; v: a list of residue number
        if len(v)> 15: continue
        
        pep=get_table_value(pdb,k, v, dcc)
        idd=pep[0][1] #ligid_alt_chid_resn

        if dic['cif'] :
            idd, natom, peppdb=get_subcif(dic, k, v)
        else:
            natom, peppdb=get_subpdb(pdb, k, v, idd)
            
        if natom<2 : continue

        mapout = cut_map_around_xyz(mapfile, peppdb, idd)
        min1, max1,mean1,sigma1=map_info(mapout)
        print('%s: natom=%d: FullMap-sigma=%s: PepMap-sigma=%s'%(idd,natom,sigma,sigma1))

        scale=1.0
        if(float(sigma1)>0) :
            scale = float(sigma)/float(sigma1)
            for z in cont.keys(): cont1[z]=cont[z]*scale
        
        maphtml = get_html4jmol(idd, peppdb, mapout, cont1)
        

        pepid=pepid+1
        pp=pep[0][1]
        filename =' %s.html  %s.pdb  %s_cut.map  %s_com  %s %s' %(pp, pp, pp, pp, sigma, sigma1)
        if dic['cif'] :
            filename =' %s.html  %s.cif  %s_cut.map  %s_com  %s %s' %(idd, idd, idd, idd, sigma, sigma1)
            
        for i,y in enumerate(pep):
            warn=html_table_content(fw, i, y, url, maphtml) 
            ss= '%d '%pepid + ' '.join(y) + warn + filename + '\n'
            fw_itool.write('%s' %ss)
            
    fw.write('</TABLE>\n')
    fw.close()
    return
   

##########################################################
def map_around_ligand(fw_itool, dic, mapfile, pdb, dcc, ch_lig, url):
    '''mapfile: the map in cell; pdb: a list (non-poly); dcc: a list of list;
    ch_lig: a dict for ligand {'ch': [n1,n2 ..]}; url: the url for html
    '''
 
    if (not ch_lig) or (not dcc) : return

    cif1='''
loop_
_dcc_ligand.id
_dcc_ligand.residue_name
_dcc_ligand.chain_id
_dcc_ligand.dcc_correlation 
_dcc_ligand.real_space_R
_dcc_ligand.Biso_mean 
_dcc_ligand.occupancy_mean
_dcc_ligand.warning
_dcc_ligand.file_name_map_html
_dcc_ligand.file_name_pdb
_dcc_ligand.file_name_map
_dcc_ligand.file_name_jmol
_dcc_ligand.full_map_sigma
_dcc_ligand.sub_map_sigma
'''
    
    if not dic['sdsc_map'] : 
        fw_itool.write(cif1)
        html_table = mapfile + '_lig.html'
        fw=open(html_table, 'w')
        html_table_head(fw, 'Ligands', 1)
    
    min, max,mean,sigma=map_info(mapfile)
    cont={'0.5':0.5, '0.7':0.7, '1.0':1.0, '1.5':1.5, '2.0':2.0} #contour 4 map in asu.
    cont1={'0.5':0.5, '0.7':0.7, '1.0':1.0, '1.5':1.5, '2.0':2.0} #contour 4 sub map.
    contlist=['0.5','0.7', '1.0', '1.5', '2.0']
    lig_sdsc, level_sdsc, jmol_sdsc=[],[],[]

    nlig, ncov=0,0
    for k, v in ch_lig.items():  #k: chainID; v: a list of residue number
        
        nres_list = isolate_connect_ligand(pdb, k, v)     
        ligid=0
        for ii, x in enumerate(nres_list):
            
            pep=get_table_value(pdb,k, x, dcc) #look through dcc table natom>=2
            
            if not pep : continue
            
            if dic['cif'] :
                idd, natom, ligpdb=get_subcif(dic, k, x)
            else:
                idd=pep[0][1] #ligid_alt_chid_resn
                natom, ligpdb=get_subpdb(pdb, k, x, idd)

            if natom<2 : continue
            
            mapout = cut_map_around_xyz(mapfile, ligpdb, idd)
            min1,max1,mean1,sigma1=map_info(mapout)  #for 
            print('%s: natom=%d: FullMap-sigma=%s: LigMap-sigma=%s'%(idd,natom,sigma,sigma1)) 
            if len(x)>1: # exist of covelently bonded ligands
                ncov=ncov+1
                s1=[]
                for z in pep:  s1.append(z[1])
                ss='","'.join(s1)
                
                sss=' {"id":"composite_%d","ligands":["'%ncov +ss + '"]},' 
            else:
                sss=' {"id":"' + pep[0][1] + '","ligands":["' + pep[0][1] + '"]},'
        
            lig_sdsc.append([sss])
            scale=1.0
            if(float(sigma1)>0) :
                scale = float(sigma)/float(sigma1)
                for z in cont.keys(): cont1[z]=cont[z]*scale
            else:
                util.perror('Warning: Negative sigma scale, possible no map cut. Check needed.') 
        
            if dic['sdsc_map'] : 
                level, jmol=gen_ligmap_sdsc(idd, contlist, cont1, cont)
                level_sdsc.append(level)
                jmol_sdsc.append(jmol)
                continue

            maphtml = get_html4jmol(idd, ligpdb, mapout, cont1)
            
            ligid=ligid + 1
            filename =' %s.html  %s.pdb  %s_cut.map %s_com %s %s' %(idd,idd,idd,idd, sigma,sigma1)
            if dic['cif'] :
                filename =' %s.html  %s.cif  %s_cut.map %s_com %s %s' %(idd,idd,idd,idd,sigma,sigma1)
            
            for i,y in enumerate(pep):
                warn=html_table_content(fw, i, y, url, maphtml)
                ss= '%d '%ligid + ' '.join(y) + warn + filename + '\n'
                fw_itool.write('%s' %ss)
                nlig=nlig+1

    if not dic['sdsc_map'] and nlig==0 : fw_itool.write('? ? ? ? ? ? ? ? ? ? ? ? ? ? \n')
    
    if dic['sdsc_map'] >0:
        fw=open('ERF_table.json', 'w')
        if len(level_sdsc)<=0:
            fw.close()
            return
        fw.write('{\n  "components":[\n')
        write_sdsc_map(fw,lig_sdsc)
        fw.write('  ],\n')

        fw.write('\n  "ligmap":[\n')
        write_sdsc_map(fw,level_sdsc)
        fw.write('  ],\n')

        fw.write('\n  "contour_level":[\n')
        write_sdsc_map(fw,jmol_sdsc)
        fw.write('  ]\n}\n')
        
        fw.close()
    
    else:
        fw.write('</TABLE>\n')
        fw.close()

##########################################################
def write_sdsc_map(fw,alist):
    '''remove the last comma
    '''
    all_1=[]
    for  x in alist :
        for y in x:
            all_1.append(y.rstrip())
            
    all_1[-1]=all_1[-1][0:-1]   
    for x in all_1 : fw.write('%s\n'%x)
    
    
##########################################################
def isolate_connect_ligand(pdb, k, v):
    '''separate isolated and connected ligands
    pdb: the list of coordinate (all atoms)
    k: the chainID;  v: the list of residue numbers. 
    '''
    
    if len(v)==1: return [v]

    idd = '%s_%s_all' %(k, v[0])
    natom, pdb_lig=get_subpdb(pdb, k, v, idd)

    tmp=[v[0]]
    for i,x in enumerate(v):
        if i==0 : continue
        n1, n2 = i-1, i
        nc=connect(pdb_lig, n1, n2, k, v)
        if not nc: tmp.append(99999)
        tmp.append(x)
        
    ss=''    
    for x in tmp: ss = ss + ' %d' %x
    t1=(ss.split('99999'))
    nres_list=[]
    for x in t1:
        tt=x.split()
        nres_list.append([int(i) for i in tt])
        
    util.delete_file(pdb_lig)
#    print(nres_list)
    return nres_list

##########################################################
def get_table_value(pdb, ch, nres, dcc):
    '''dcc is a list of list; nres is a list of residue number
    '''

    
    pep=[]
    for x in nres:  #residue number
        for y in dcc:
            if ch == y[1] and x == int(y[0]):
                
                natom=0  #
                for z in pdb:
                    if (('ATOM' in z[:6] or 'HETATM' in z[:6]) and
                        ch==z[20:22].strip() and x==int(z[22:26])):
                        natom=natom+1
                if natom<2: continue
                
                s='%s_%s_%s_%d' %(y[2], y[3],y[1], x)
                pep.append([y[2], s, y[4],y[5],y[6],y[7]])
                break

    return pep

##########################################################
def connect(pdb_ligfile, n1, n2, k, dic):
    '''check if residue n1 and n2 is connected. k is chainID
    pdb_lig is a pdb file
    
    '''
    
    cutoff=1.6
    
    if not os.path.exists(pdb_ligfile) or  os.path.getsize(pdb_ligfile)<10 : return 
#    if util.check_file(30, pdb_ligfile)==0: return
    
    pdb_lig=open(pdb_ligfile, 'r').readlines()
    data1, data2 = [], []
    for x in pdb_lig:
        if k==x[20:22].strip() :
            if dic[n1]==int(x[22:26]):
                data1.append(x)
            elif dic[n2]==int(x[22:26]):
                data2.append(x)
    conn=0            
    if data1 and data2:
        for x in data1:
            if ('ATOM' not in x[:6] and 'HETATM' not in x[:6]): continue
            x1,y1,z1=float(x[28:38]), float(x[38:46]), float(x[46:54])
            
            for y in data2:
                if ('ATOM' not in y[:6] and 'HETATM' not in y[:6]): continue
                x2,y2,z2=float(y[28:38]), float(y[38:46]), float(y[46:54])
                
                d=math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
                if d<cutoff:
                    conn=1
                    break
            if conn==1 : break
#    print( conn, n1, n2, k, dic )          
    return conn
        
    
##########################################################
def cut_map_around_xyz(mapfile, peppdb_in, pepid):
    '''using several CCP4 until to cut map around the selected molecule
    color: [x800080] [xc00000] [xb0b0b0]
    '''
    
    mapout = pepid + '_cut.map'
    util.delete_file(mapout)

    
    '''
    xyzlim,xyzcomp=find_xyzlim_compound(pepid, peppdb_in)  
    mapscr=cut_map_bylimit(xyzlim)
    arg = mapfile + ' ' + ' ' +  mapout
    command="chmod +x %s ; ./%s  " %(mapscr, mapscr) + arg 
    os.system(command)
    
    return  mapout
    '''
    
    peppdb = peppdb_in
    if util.is_cif(peppdb_in): peppdb=cif.cif2pdb(peppdb_in)
    
    mapscr=cut_map_scr() 
    arg = mapfile + ' ' + peppdb + ' ' +  mapout
    command="chmod +x %s ; ./%s  " %(mapscr, mapscr) + arg 
    os.system(command)

    if util.is_cif(peppdb_in):
        util.delete_file('%s.PDB' %peppdb_in)
        peppdb = peppdb_in

    return mapout
            
            
##############################################################
def gen_ligmap_sdsc(idd, contlist, cont, cont1):
    
    color={'0.5':'[xb0b0b0]','0.7':'[xFFFF30]',
           '1.0':'[x00AB24]','1.5':'[x3050F8]','2.0':'[xFFFF30]'} 
    level,jmol=[],[]
    comma=','
    n=len(contlist)
    for i,x in enumerate(contlist):
        
        ss= '    {"id":"%s", "name":"%s_cut.map", ' %(idd, idd) + \
            '"actual_contour_level_cut":%.1f, "actual_contour_level_full":%.1f}%s' %(cont[x], cont1[x], comma)
          
        
        level.append(ss)
        s1='''    {"id":"%s", "script":"isosurface s_one color %s sigma %.1f within 2.0 {*} \
'%s_cut.map' mesh nofill;"}%s''' %(idd, color[x], cont[x], idd, comma)
        
        jmol.append(s1)

    return level, jmol

##############################################################
def get_html4jmol(ligid, ligpdb,  mapout, cont):
    '''Each ligand has a .com file
    The generated HTML file is all is needed for display.
    The Jmol package has to be specified. (borrow Jmol.js from somewhere,
    if not in the package!)
    
    color red[x800080] [xd000d0] [xc00000] pink[xb0b0b0] green[x008000]
    http://jmol.sourceforge.net/jscolors/#color_Sr
    
    '''

    
    jmol_pth='../..'  #put source : /net/wwwdev/auto-check/html/Jmol-13/
#    jmol_pth='/net/techusers/hyang/program'
    maphtml = ligid + '.html'
    jmol_com_name = ligid + '_com'

    html="""
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
        "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">

<head>
  <script type="text/javascript" src="%s/Jmol-13/Jmol.js"></script>
  
  <script language='Javascript' type='text/javascript'>
  // if second parameter is true, use signed applet
  jmolInitialize("%s/Jmol-13/", window.location.protocol == "file:");
  </script>

	<title>Displaying Electron Density Map by Jmol</title>

</head>
<body bgcolor="#d8d8d8" >


<br>
<!-- JMOL 0 ==================== -->
<table align="left" hspace="10"><tr><td>
<script type="text/javascript">
	jmolApplet(600, "script %s", "0");

jmolBr();
jmolLink("isosurface s_one color [xb0b0b0] sigma %.2f within 2.0 {*} '%s' mesh nofill;", "Sigma cutoff 0.5");
jmolBr();
jmolLink("isosurface s_one color [xFFFF30] sigma %.2f within 2.0 {*} '%s' mesh nofill;", "Sigma cutoff 0.7");
jmolBr();
jmolLink("isosurface s_one color [x00AB24] sigma %.2f within 2.0 {*} '%s' mesh nofill;", "Sigma cutoff 1.0");
jmolBr();
jmolLink("isosurface s_one color [x3050F8] sigma %.2f within 2.0 {*} '%s' mesh nofill;", "Sigma cutoff 1.5");
jmolBr();
jmolLink("isosurface s_one color [xFFFF30] sigma %.2f within 2.0 {*} '%s' mesh nofill;", "Sigma cutoff 2.0");
jmolBr();        
</script>
</td><td> &nbsp; &nbsp; </td>

</tr></table>
<!-- =========================== -->

<p><center><a href=\"javascript:window.history.back();\"><b>Go Back</b></a></a></center>

</body>
</html>
""" %(jmol_pth,jmol_pth, jmol_com_name, cont['0.5'],mapout, cont['0.7'],mapout,
      cont['1.0'], mapout, cont['1.5'],mapout, cont['2.0'],mapout)
      

    fw=open(maphtml, "w")
    fw.write ("%s" % html )
    fw.close()

    mapscr="""
load %s; spacefill off; rotate x 90;

isosurface s_one color [x00AB24]  within 2.0 {*} "%s" mesh nofill;
#isosurface cutoff 0.2
set echo e1 [0 35];
#echo cutoff 0.2
echo "2mFo-DFc"
color echo magenta
""" %(ligpdb,  mapout)  

    fw=open(jmol_com_name, "w")
    fw.write ("%s" %mapscr )
    fw.close()

    return maphtml
    
    
#####################################################################
def cut_map_scr():
    
                      
    csh_script="""#!/bin/csh -f

######################################################################
# This script is used to cut ASU CCP4 map to a small size around a 
# selected part of the model 
# (created 10/31/2011)
#=================== Usage ========================== 
# mapcut mapfile pdbfile  mapout
######################################################################

set mapfile=$1
set pdbfile=$2
set mapout=$3
set maptmp=${1}_tmp

if( $#argv == 2 ) then
    set mapout="${mapfile}_cut"
endif


#cut map around the selected volume (pdbfile) 
# Border can not be smaller that the value set for the big map (>4.)
mapmask MAPIN $mapfile  xyzin $pdbfile  MAPOUT  $mapout <<EOF >& /dev/null
Border 4.0
SCALE FACTOR 1.0 0.0
MODE mapin
EOF

"""

    mapscr="mapcut_TMP.csh"
    fw=open(mapscr, "w")
    fw.write ("%s" % csh_script )
    fw.close()

    return mapscr

#####################################################################
def cut_map_bylimit(xyzlim):
    
    csh_script="""#!/bin/csh -f

######################################################################
# This script is used to cut ASU CCP4 map to a small size around a 
# selected xyz limit in fraction (2013-02-03)
#
#---------------------------- Usage ----------------------------------
# mapcut mapfile  mapout
#######################################################################

set mapfile=$1
set mapout=$2
set maptmp=${1}_tmp

if( $#argv == 1 ) then
    set mapout="${mapfile}_cut"
endif

#cut map around the selected box. No Border should be given!!
mapmask MAPIN $mapfile  MAPOUT  $mapout <<EOF >& /dev/null
XYZLIM  %s
SCALE FACTOR 1.0 0.0
MODE mapin
EOF


""" %(xyzlim)

    mapscr="mapmask_TMP.csh"
    fw=open(mapscr, "w")
    fw.write ("%s" % csh_script )
    fw.close()

    return mapscr

##########################################################
def map_info(mapfile):
    '''get the min, max, mean, sigma from the map
    '''

    min, max, mean, sigma='-1', '-1','-1','-1'
    log=mapfile + '_header'
    scr=mapfile + '.sh'
    
    arg='mapdump mapin %s <<eof >%s \neof\n' %(mapfile, log)
    
    fw = open(scr, 'w')
    fw.write(arg)
    fw.close()
    os.system('chmod +x %s ; ./%s' %(scr,scr))


    if not util.check_file(10,log): return  min, max, mean, sigma
    fp=open(log,'r')
    for x in fp:
        if 'Minimum density ...' in x:
            min=x.rstrip().split( ' ')[-1]
        elif 'Maximum density ...' in x:
            max=x.rstrip().split( ' ')[-1]
        elif '    Mean density ..' in x:
            mean=x.rstrip().split( ' ')[-1]
            
        elif 'deviation from mean density ..' in x:
            sigma=x.rstrip().split( ' ')[-1]

    util.delete_file(scr, log)
            
    return  min, max, mean, sigma        
