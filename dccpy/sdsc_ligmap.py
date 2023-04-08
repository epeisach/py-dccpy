#!/usr/bin/env python

import os
import math
import sys
import ligand
import parse
import util
import cifparse as cif
#


def process(*files):

    arg = files[0]
    narg = len(arg)
    xyz = ''
    map2fofc = ''
    for k in range(narg):
        if (arg[k].upper() == "-PDB" or arg[k].upper() == "-CIF"):
            xyz = arg[k + 1]
        elif arg[k].upper() == "-MAP":
            map2fofc = arg[k + 1]

    dic = {}
    dic['pdbfile'] = xyz
    _pdb = cif.cif2pdb(xyz)  # noqa: F841
    pdbfile = xyz + '.PDB'

    sdsc_ligand(dic, map2fofc, pdbfile)


######################################################
def sdsc_ligand(dic, mapfile, pdbfile):
    '''This is the module to cut the ligand maps for SDSC. The PRD and ligand
    are parsed from CIF only. The convalently bonded ligands are determined
    by the distance based on the consequent residue number sequence. This
    could mis-identify the covalent bond. For example.4aqd (or 5d4q)
    'NAG_A_611', 'NAG_A_612', 'FUL_A_613' should be connected, but only 611&612
    is identified. It is ok at the moment

    mapfile: the input mapfile (must have more boudary than the cutted residues)
    pdbfile: The pdbfile
    dic: a dictionary, it should contains the cif xyz file.
    '''

    pdblist = open(pdbfile, 'r').readlines()

    ciffile = dic['xyzfile_orig']
    if not util.is_cif(ciffile):
        util.perror('Warning: Input file %s is not in CIF format. Not map cut!' % ciffile)
        return

    xyzlist = open(ciffile, 'r').readlines()
    nonploy = parse.parse_cif(xyzlist, 'pdbx_nonpoly_scheme')  # for ligands
    conn = parse.parse_cif(xyzlist, 'struct_conn')  # for link of ligands
    prd = parse.parse_cif(xyzlist, 'pdbx_molecule')  # for PRD
    poly = parse.parse_cif(xyzlist, 'pdbx_poly_seq_scheme')  # for PRD

    connect = []  # test covalent bonded ligand
    if conn:
        nn = len(conn['id'])
        for i in range(nn):
            if 'covale' not in conn['id'][i]:
                continue
            res1 = '%s_%s_%s' % (conn['comp1'][i], conn['asym1'][i], conn['nseq1'][i])
            res2 = '%s_%s_%s' % (conn['comp2'][i], conn['asym2'][i], conn['nseq2'][i])
            connect.append([res1, res2])

    lig_ch = {}  # All the ligands
    res_pdb = {}
    if nonploy:
        nn = len(nonploy['monid'])
        for i in range(nn):
            if 'MG' == nonploy['monid'][i] or 'CL' == nonploy['monid'] \
               or 'HOH' in nonploy['monid'][i] or 'DOD' in nonploy['monid'][i]:
                continue
            comp, ch, nres = nonploy['monid'][i], nonploy['chid'][i], nonploy['nseq'][i]
            res = '%s_._%s_%s' % (comp, ch, nres)
            res_pdb[res] = residue_xyz(comp, ch, nres, pdblist)
            if ch not in lig_ch.keys():
                lig_ch[ch] = []
            lig_ch[ch].append(res)

    lig_group = []  # group the connected ligands
    for ch, lig in lig_ch.iteritems():
        tmp = get_lig_group(lig, res_pdb)
        if tmp:
            lig_group.extend(tmp)

    prd_all = []  # for PRD molecules
    if prd:
        nn = len(prd['ins_id'])
        for i in range(nn):
            asym = prd['asm_id'][i]
            if asym not in poly['asym']:
                continue
            prd_tmp = []
            for j, x in enumerate(poly['asym']):
                if asym == x:
                    comp, ch, nres = poly['monid'][j], poly['chid'][j], poly['nseq'][j]
                    if '?' in comp:
                        continue  # not in pdb
                    res = '%s_._%s_%s' % (comp, ch, nres)
                    res_pdb[res] = residue_xyz(comp, ch, nres, pdblist)
                    prd_tmp.append(res)

            if prd_tmp:
                prd_all.append(prd_tmp)
    if prd_all:
        lig_group.extend(prd_all)

    # print '\nfinal=', lig_group

    scale_card = []
    for x in pdblist:
        if ('CRYST1' in x[:6] or 'SCALE' in x[:5]):
            scale_card.append(x)
        elif ('ATOM' in x[:4] or 'HETATM' in x[:6]):
            break

    mini_map(mapfile, pdbfile, lig_group, res_pdb, scale_card)

    arg = 'rm -f %s %s.PDB  %s_rcc_sum.cif.mtz  ' % (mapfile, dic['pdbfile'], dic['pdbfile'])
    if (dic['verb'] == 0):
        os.system(arg)


##########################################################
def get_minipdb(x, scale, res_pdb):
    '''generate the protein for the non-poly
    x: the ligand or connected or the PRD
    scale: the scale card in the PDB
    res_pdb
    '''

    natom = 0
    idd = x[0]
    pdb = '%s.pdb' % idd
    fw = open(pdb, 'w')
    for z in scale:
        fw.write(z)

    for z in x:
        for k in res_pdb[z]:
            natom = natom + 1
            fw.write(k)

    fw.close()
    return natom, pdb


##########################################################
def residue_xyz(comp, ch, nres, pdblist):
    xyz = []
    for i, x in enumerate(pdblist):
        if ('ATOM' not in x[:6] and 'HETATM' not in x[:6]):
            continue

        if xyz and ch == x[20:22] != pdblist[i - 1][20:22]:
            break
        elif comp in x[17:20] and ch == x[20:22].strip() and nres == x[22:26].strip():
            xyz.append(x)

    return xyz


##########################################################
def test_connect(res1, res2, data1, data2):

    # print res1,res2,len(data1), len(data2)
    cutoff = 1.6
    conn = 0
    if data1 and data2:
        for x in data1:
            if ('ATOM' not in x[:6] and 'HETATM' not in x[:6]):
                continue
            x1, y1, z1 = float(x[28:38]), float(x[38:46]), float(x[46:54])

            for y in data2:
                if ('ATOM' not in y[:6] and 'HETATM' not in y[:6]):
                    continue
                x2, y2, z2 = float(y[28:38]), float(y[38:46]), float(y[46:54])

                d = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
                # print res1,res2,d
                if d < cutoff:
                    # print res1,res2,d
                    conn = 1
                    break
            if conn == 1:
                break
    # print( conn, n1, n2, k, dic )
    return conn


##########################################################
def get_lig_group(liglist, res_pdb):
    '''group the ligand if bonded [[L1, L2], [L3],,].

    liglist:  the list of ligand
    res_pdb:  the dictionary to contain xyz. res_pdb[res:[,,,]]
    '''

    if len(liglist) == 1:
        return [liglist]

    tmp = [liglist[0]]
    for i, x in enumerate(liglist):
        if i == 0:
            continue
        xyz1 = res_pdb[liglist[i - 1]]
        xyz2 = res_pdb[liglist[i]]

        nc = test_connect(liglist[i - 1], liglist[i], xyz1, xyz2)
        if not nc:
            tmp.append(' 99999')
        tmp.append(x)

    ss = ' '.join(tmp)

    t1 = (ss.split(' 99999'))
    nres_list = []
    for x in t1:
        tt = x.split()
        nres_list.append([y for y in tt])

    return nres_list


##########################################################
def mini_map(mapfile, pdbfile, nonpoly, res_pdb, scale_card):
    '''mapfile: the map in asu;
    pdbfile: a list (non-poly);
    nonpoly: [[comp_alt_ch_nres, ..], [..]]
    '''

    cont = {'0.5': 0.5, '0.7': 0.7, '1.0': 1.0, '1.5': 1.5, '2.0': 2.0}  # contour 4 map asu.
    cont1 = {'0.5': 0.5, '0.7': 0.7, '1.0': 1.0, '1.5': 1.5, '2.0': 2.0}  # contour 4 sub map.
    contlist = ['0.5', '0.7', '1.0', '1.5', '2.0']

    min, max, mean, sigma = ligand.map_info(mapfile)  # big map info
    lig_sdsc, level_sdsc, jmol_sdsc = [], [], []

    ncov = 0
    for x in nonpoly:  # [[..], [..] ..] a list of ligand,prd
        idd = x[0]
        natom, ligpdb = get_minipdb(x, scale_card, res_pdb)
        if natom < 2:
            os.system('rm -f %s' % ligpdb)
            continue
        mapout = ligand.cut_map_around_xyz(mapfile, ligpdb, idd)
        min1, max1, mean1, sigma1 = ligand.map_info(mapout)  # sub map info
        print('%s: natom=%d: FullMap-sigma=%s: LigMap-sigma=%s' % (idd, natom, sigma, sigma1))

        if len(x) > 1:  # exist of covalently bonded ligands, pr PRD
            ncov = ncov + 1
            ss = '","'.join(x)
            sss = ' {"id":"composite_%d","ligands":["' % ncov + ss + '"]},'
        else:
            sss = ' {"id":"' + x[0] + '","ligands":["' + x[0] + '"]},'

        lig_sdsc.append([sss])
        scale = 1.0
        if (float(sigma1) > 0):  # re-scale the contour level to match the big map
            scale = float(sigma) / float(sigma1)
            for z in cont.keys():
                cont1[z] = cont[z] * scale
        else:
            util.perror('Error: Mapcut is wrong. No scale is applied. (check needed).')

        level, jmol = ligand.gen_ligmap_sdsc(idd, contlist, cont1, cont)
        level_sdsc.append(level)
        jmol_sdsc.append(jmol)

    fw = open('ERF_table.json', 'w')
    if len(level_sdsc) <= 0:
        fw.close()
        return
    fw.write('{\n  "components":[\n')
    ligand.write_sdsc_map(fw, lig_sdsc)
    fw.write('  ],\n')

    fw.write('\n  "ligmap":[\n')
    ligand.write_sdsc_map(fw, level_sdsc)
    fw.write('  ],\n')

    fw.write('\n  "contour_level":[\n')
    ligand.write_sdsc_map(fw, jmol_sdsc)
    fw.write('  ]\n}\n')

    fw.close()


##########################################################
if __name__ == '__main__':
    process(sys.argv)
