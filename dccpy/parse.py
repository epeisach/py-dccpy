import os
import util
import cifparse as cif


# ===========================================================
# this module is to parse all the data values and put them into dict
# ===========================================================

def parse_cif(flist, table):

    dic = {}

    ciftable = '_%s.' % table
    items, values = cif.cifparse(flist, ciftable)
    if ciftable == '_pdbx_nonpoly_scheme.':
        asym = cif.parse_values(items, values, '_pdbx_nonpoly_scheme.asym_id')
        monid = cif.parse_values(items, values, '_pdbx_nonpoly_scheme.pdb_mon_id')
        nseq = cif.parse_values(items, values, '_pdbx_nonpoly_scheme.pdb_seq_num')
        chid = cif.parse_values(items, values, '_pdbx_nonpoly_scheme.pdb_strand_id')
        ins = cif.parse_values(items, values, '_pdbx_nonpoly_scheme.pdb_ins_code')
        if monid and nseq and asym:
            dic = {'monid': monid, 'nseq': nseq, 'asym': asym, 'chid': chid, 'ins': ins}

    elif ciftable == '_struct_conn.':
        id = cif.parse_values(items, values, '_struct_conn.id')  # pylint: disable=redefined-builtin
        asym1 = cif.parse_values(items, values, '_struct_conn.ptnr1_auth_asym_id')
        comp1 = cif.parse_values(items, values, '_struct_conn.ptnr1_auth_comp_id')
        nseq1 = cif.parse_values(items, values, '_struct_conn.ptnr1_auth_seq_id')
        alt1 = cif.parse_values(items, values, '_struct_conn.pdbx_ptnr1_label_alt_id')
        ins1 = cif.parse_values(items, values, '_struct_conn.pdbx_ptnr1_PDB_ins_code')

        asym2 = cif.parse_values(items, values, '_struct_conn.ptnr2_auth_asym_id')
        comp2 = cif.parse_values(items, values, '_struct_conn.ptnr2_auth_comp_id')
        nseq2 = cif.parse_values(items, values, '_struct_conn.ptnr2_auth_seq_id')
        alt2 = cif.parse_values(items, values, '_struct_conn.pdbx_ptnr2_label_alt_id')
        ins2 = cif.parse_values(items, values, '_struct_conn.pdbx_ptnr2_PDB_ins_code')

        if (id and asym1 and comp1 and nseq1 and asym2 and comp2 and nseq2):
            dic = {'id': id, 'asym1': asym1, 'comp1': comp1, 'nseq1': nseq1, 'alt1': alt1, 'ins1': ins1,
                   'asym2': asym2, 'comp2': comp2, 'nseq2': nseq2, 'alt2': alt2, 'ins2': ins2}

    elif ciftable == '_pdbx_molecule.':
        ins_id = cif.parse_values(items, values, '_pdbx_molecule.instance_id')
        prd_id = cif.parse_values(items, values, '_pdbx_molecule.prd_id')
        asm_id = cif.parse_values(items, values, '_pdbx_molecule.asym_id')

        if (ins_id and prd_id and asm_id):
            dic = {'ins_id': ins_id, 'prd_id': prd_id, 'asm_id': asm_id}

    elif ciftable == '_pdbx_poly_seq_scheme.':
        asym = cif.parse_values(items, values, '_pdbx_poly_seq_scheme.asym_id')
        nseq = cif.parse_values(items, values, '_pdbx_poly_seq_scheme.pdb_seq_num')
        monid = cif.parse_values(items, values, '_pdbx_poly_seq_scheme.pdb_mon_id')
        chid = cif.parse_values(items, values, '_pdbx_poly_seq_scheme.pdb_strand_id')
        ins = cif.parse_values(items, values, '_pdbx_poly_seq_scheme.pdb_ins_code')
        if asym and monid:
            dic = {'monid': monid, 'nseq': nseq, 'asym': asym, 'chid': chid, 'ins': ins}

    return dic


##########################################################
def refmac_log(logfile, dic):
    '''update dictionary (dic) from filename. opt is an option
    '''
    if not util.check_file(100, logfile):
        return
    fp = open(logfile, "r")
    dic['prog'] = 'REFMAC'
    dic['twin_op'], dic['twin_fr'] = '', ''
    ndom = 0
    for x in fp:
        if 'Resolution limits' in x[:24]:
            t = util.get_value_after_id(x[28:], '=').split()

            dic['resl'] = t[0].strip()
            dic['resh'] = t[1].strip()

        elif 'Partial structure    1:' in x[:24]:
            dic['k_sol'] = util.str_between_id(x, '=', ',')
            dic['b_sol'] = util.get_value_after_id(x[33:], '=')

        elif 'Number of used reflections' in x[:30]:
            dic['nref'] = util.get_value_after_id(x, '=')

        elif 'Percentage observed' in x[:21]:
            dic['comp'] = util.get_value_after_id(x, '=')

        elif 'Percentage of free reflections' in x[:31]:
            dic['pfree'] = util.get_value_after_id(x, '=')

        elif 'Number of "free" reflections' in x:
            dic['nfree'] = util.get_value_after_id(x, 'tions')

        elif 'Overall R factor' in x[:17]:
            dic['rfact'] = util.get_value_after_id(x, '=')

        elif 'Free R factor' in x[:15]:
            dic['rfree'] = util.get_value_after_id(x, '=')

        elif 'Overall correlation coefficient' in x[:33]:
            dic['fcc'] = util.get_value_after_id(x, '=')

        elif 'Free correlation coefficient' in x[:30]:
            dic['rfcc'] = util.get_value_after_id(x, '=')

        elif 'Cruickshanks DPI for coordinate error' in x:
            dic['dpi_xyz'] = util.get_value_after_id(x, '=')

        elif 'DPI based on free R factor' in x:
            dic['dpi_fr'] = util.get_value_after_id(x, '=')

        elif 'Overall figure of merit' in x:
            dic['fom'] = util.get_value_after_id(x, '=')

        elif 'Twin operator:' in x[:16] and 'Fraction =' in x:
            ndom = ndom + 1
            op = util.str_between_id(x, 'or:', ': Fr').strip().replace(' ', '')
            dic['twin_op'] = dic['twin_op'] + '%d: %s ' % (ndom, op)

        elif 'Twin fractions  ' in x[:16]:
            t = x.strip().split('=')[1].split()
            f = ['%d: %s' % (i + 1, y) for i, y in enumerate(t)]

            if len(f) > 1:
                dic['twin_fr'] = ' '.join(f)

    fp.close()


##########################################################
def data_from_refmac_cif(inpfile, dd):
    if not os.path.exists(inpfile) or os.path.getsize(inpfile) < 100:
        return

    fr = open(inpfile, "r")
    for ln in fr:

        if '_refine.ls_d_res_high ' in ln:  # from refmac
            dd['resh'] = ln.split()[1]

        elif '_refine.ls_d_res_low' in ln:
            dd['resl'] = ln.split()[1]

        elif '_refine.ls_number_reflns_R_work' in ln:
            dd['nref'] = ln.split()[1]

        elif '_refine.ls_number_reflns_R_free' in ln:
            dd['nfree'] = ln.split()[1]

        elif '_refine.ls_percent_reflns_R_free' in ln:
            dd['pfree'] = ln.split()[1]

        elif '_refine.ls_R_factor_R_work' in ln:
            dd['rfact'] = ln.split()[1]

        elif '_refine.ls_R_factor_R_free ' in ln:
            dd['rfree'] = ln.split()[1]

        elif '_refine.ls_R_factor_R_all ' in ln:
            dd['rall'] = ln.split()[1]

        elif '_refine.ls_percent_reflns_obs ' in ln:
            dd['comp'] = ln.split()[1]

        elif '_refine.correlation_coeff_Fo_to_Fc ' in ln:
            dd['fcc'] = ln.split()[1]

        elif "_pdbx_rscc_refmac_overall.ls_d_res_high" in ln:  # refmac
            dd['resh'] = util.get_value_after_id(ln, "high")

        elif "_pdbx_rscc_mapman_overall.correlation " in ln:
            dd['dcorr'] = util.get_value_after_id(ln, "elation")

        elif "_pdbx_rscc_mapman_overall.real_space_R " in ln:
            dd['realr'] = util.get_value_after_id(ln, "space_R")

    if dd['prog'].upper() == 'REFMAC' and dd['nfree'].isdigit():
        tmp = int(dd['nref']) + int(dd['nfree'])
        dd['nref'] = '%s' % tmp

    fr.close()


##########################################################
def sfcheck_log(inpfile, dd, detail):
    '''extract calculated data from file and put them into dd
    '''

    if not util.check_file(100, inpfile):
        return

    dd['prog'] = 'SFCHECK'
    if detail == 'NOTLS':
        dd['detail'] = 'without TLS correction'
    else:
        dd['detail'] = 'with    TLS correction'

    fr = open(inpfile, "r")
    for ln in fr:
        if "Resolution /input data/ :" in ln:  # sfcheck
            tmp = util.get_value_after_id(ln, ':')
            dd['resl'], dd['resh'] = [x.strip() for x in tmp.split('-')]

        elif 'Completeness :' in ln:
            tmp = util.get_value_after_id(ln, ':')
            dd['comp'] = tmp.split()[0]

        elif 'R-factor           :' in ln:
            dd['rall'] = util.get_value_after_id(ln, ':')
            dd['rfact'] = dd['rall']

        elif 'Correlation factor :' in ln:
            dd['fcc'] = util.get_value_after_id(ln, ':')

        # elif 'N_refls        :' in ln:
        #    dd['nref']= util.get_value_after_id(ln, ':')

        elif 'Number of reflections :' in ln and 'in_file ' in ln:
            dd['nref'] = util.get_value_after_id(ln, ':').split('(')[0]

        elif 'Nfree,Nrest        :' in ln:
            dd['nfree'] = util.get_value_after_id(ln, ':').split()[0]

        elif 'Rfree,Rrest        :' in ln:
            tmp = util.get_value_after_id(ln, ':')
            dd['rfree'], t = [x.strip() for x in tmp.split()]
            if (util.is_number(t)):
                dd['rfact'] = t
        elif 'Boverall /estimated by Wilson plot/     : ' in ln:
            dd['bwilson'] = util.get_value_after_id(ln, ':')

        elif 'PDB_code :' in ln:
            dd['pdbid'] = util.get_value_after_id(ln, ':')

    fr.close()


##########################################################
def sfcheck_cif(inpfile, dd, detail):  # pylint: disable=unused-argument

    if not util.check_file(100, inpfile):
        return

    fr = open(inpfile, "r")
    for ln in fr:
        if "_pdbx_map_overall.RSR" in ln:  # sfcheck,
            dd['realr'] = util.get_value_after_id(ln, "RSR")

        elif "_pdbx_map_overall.RSCC" in ln:
            dd['dcorr'] = util.get_value_after_id(ln, "RSCC")
    fr.close()


##########################################################
def model_vs_data_log(infile, dic1, dic):
    '''extract values and put them into dic (should be often updated)
    '''

    if not util.check_file(100, infile):
        return

    dic['prog_c'] = 'PHENIX'

    fr = open(infile, 'r')
    for ln in fr:
        if 'Stereochemistry' in ln and 'overall:' in ln:
            for y in fr:
                t = y.split(':')[1].split()

                if 'bonds            :' in y:
                    dic['bond'], dic['bond_max'] = t[0], t[1]
                elif 'angles           :' in y:
                    dic['angle'], dic['angle_max'] = t[0], t[1]
                elif 'dihedrals        :' in y:
                    dic['dihedral'], dic['dihedral_max'] = t[0], t[1]
                elif 'chirality        :' in y:
                    dic['chirality'], dic['chirality_max'] = t[0], t[1]
                elif 'planarity        :' in y:
                    dic['planarity'], dic['planarity_max'] = t[0], t[1]
                elif 'non-bonded (min) :' in y:
                    dic['non_bond'] = t[0].strip()
                    break

        elif 'Stereochemistry' in ln and 'ligands:' in ln:
            for y in fr:
                t = y.split(':')[1].split()
                if 'bonds            :' in y:
                    dic['bond_lig'], dic['bond_lig_max'] = t[0], t[1]
                elif 'angles           :' in y:
                    dic['angle_lig'], dic['angle_lig_max'] = t[0], t[1]
                elif 'non-bonded (min) :' in y:
                    break

        elif 'Molprobity statistics:' in ln:
            for y in fr:
                if 'outliers :' in y:
                    dic['outlier_num'] = util.str_between_id(y, ':', '(')
                    dic['outlier'] = util.str_between_id(y, '(', '%')
                elif '  allowed  :' in y:
                    dic['allowed_num'] = util.str_between_id(y, ':', '(')
                    dic['allowed'] = util.str_between_id(y, '(', '%')
                elif '  favored  :' in y:
                    dic['favored_num'] = util.str_between_id(y, ':', '(')
                    dic['favored'] = util.str_between_id(y, '(', '%')
                elif ' Rotamer outliers' in y:
                    dic['rot_outlier_num'] = util.str_between_id(y, ':', '(')
                    dic['rot_outlier'] = util.str_between_id(y, '(', '%')
                elif 'Cbeta deviations' in y:
                    dic['cbeta'] = y.strip().split(':')[1].split()[0]
                elif 'All-atom clashscore ' in y:
                    dic['clashscore'] = y.strip().split(':')[1].split()[0]
                elif 'Overall score ' in y:
                    dic['oscore'] = y.strip().split(':')[1].split()[0]
                    break

        elif 'mFo-DFc map:' in ln:
            for y in fr:
                if ' >  3 sigma:' in y:
                    dic['fofc_3sig_p'] = y.strip().split(':')[1]
                elif ' >  6 sigma:' in y:
                    dic['fofc_6sig_p'] = y.strip().split(':')[1]
                elif '< -3 sigma:' in y:
                    dic['fofc_3sig_n'] = y.strip().split(':')[1]
                elif '< -6 sigma:' in y:
                    dic['fofc_6sig_n'] = y.strip().split(':')[1]
                    break

        elif 'bulk_solvent_(k_sol,b_sol' in ln:
            tmp = ln.split(':')[1].strip().split()
            dic1['k_sol'] = tmp[0]
            dic1['b_sol'] = tmp[1]

        elif 'high_resolution         :' in ln or 'high_resolution                      :' in ln:
            dic1['resh'] = ln.split(':')[1].strip()

        elif 'low_resolution          :' in ln or 'low_resolution                       :' in ln:
            dic1['resl'] = ln.split(':')[1].strip()

        elif 'completeness_in_range' in ln and ':' in ln:
            tmp = ln.split(':')[1].strip()
            if float(tmp) < 1.1:
                dic1['comp'] = "%.3f" % (float(tmp) * 100)
            else:
                dic1['comp'] = tmp

        elif 'wilson_b                :' in ln or 'wilson_b                             :' in ln:
            dic1['bwilson'] = ln.split(':')[1].strip()

        elif ' number_of_reflections   :' in ln or 'number_of_reflections                :' in ln:
            dic1['nref'] = ln.split(':')[1].strip()

        elif 'test_set_size           :' in ln or 'test_set_size                        :' in ln:
            tmp = ln.split(':')[1].strip()
            if '?' not in dic1['nref']:
                dic1['nfree'] = "%d" % (float(tmp) * float(dic1['nref']))

        elif 'twinned                 :' in ln:
            twin = ln.split(':')[1].strip()
            if 'False' not in twin:
                dic1['twin'] = 'Y'

        elif 'r_work(re-computed)                :' in ln:
            dic1['rfact'] = ln.split(':')[1].strip()

        elif 'r_free(re-computed)                :' in ln:
            t = ln.split(':')[1].strip()
            if not util.is_number(t):
                t = '?'
            dic1['rfree'] = t

        elif 'TLS             : ' in ln:
            tls_rep = ln.split(':')[1].strip()
            if 'False' not in tls_rep:
                dic1['tls_c'] = 'Y'
        elif 'program_name    :' in ln:
            dic1['prog_r'] = ln.split(':')[1].strip()

        elif 'r_work          :' in ln:
            dic1['rfact_r'] = ln.split(':')[1].strip()

        elif 'r_free          :' in ln:
            dic1['rfree_r'] = ln.split(':')[1].strip()

        elif 'high_resolution :' in ln:
            dic1['resh_r'] = ln.split(':')[1].strip()

        elif 'low_resolution  :' in ln:
            dic1['resl_r'] = ln.split(':')[1].strip()

        elif 'r_work_cutoff :' in ln:
            dic1['r_cutoff'] = util.float_after_id(ln, ':')
        elif 'r_free_cutoff :' in ln:
            dic1['rf_cutoff'] = util.float_after_id(ln, ':')

        elif 'PDB_string' in ln:  # for density correlation(atom list)
            dic1['dcc'] = get_dcc_avg(fr, 1)

        elif 'resseq resname' in ln:  # for density correlation(residue list)
            dic1['dcc'] = get_dcc_avg(fr, 2)

    # print(dic['dcc'])

    fr.close()


##########################################################
def model_vs_data_v13_log(infile, dic1, dic):
    '''extract values and put them into dic (should be often updated)
    '''

    if not util.check_file(100, infile):
        return

    dic['prog_c'] = 'PHENIX'

    fr = open(infile, 'r')
    for ln in fr:

        if ln[0:6] == 'Data: ':
            for y in fr:
                t = y.split(':')[1].split()

                if 'Completeness in resolution range' in y:
                    tmp = t[0].strip()
                    if float(tmp) < 1.1:
                        dic1['comp'] = "%.3f" % (float(tmp) * 100)
                    else:
                        dic1['comp'] = tmp

                elif 'Resolution range:' in y:
                    dic1['resh'] = t[1]
                    dic1['resl'] = t[0]
                elif 'Number of Miller indices:' in y:
                    dic1['nref'] = t[0]
                elif 'Wavelength' in y:
                    # Last one
                    break

        elif ln[0:15] == 'Model vs Data: ':
            intable = False
            nnfree = 0
            for y in fr:
                t = y.split(':')

                if "Resolution    Compl Nwork" in y:
                    intable = True
                    continue

                if intable:
                    t = y.split()

                    if len(t) == 0:
                        if nnfree:
                            dic1['nfree'] = nnfree
                        intable = False
                        continue
                    else:
                        nnfree += int(t[3])
                        continue

                if 'r_work:' in y:
                    dic1['rfact'] = t[1].strip()
                elif 'r_free:' in y:
                    tmp = t[1].strip()
                    if not util.is_number(tmp):
                        tmp = '?'
                    dic1['rfree'] = tmp
                elif 'Number of F-obs outliers' in y:
                    # End of category
                    break

        elif ln[0:43] == 'Information extracted from PDB file header:':
            for y in fr:
                t = y.split(':')

                if 'program_name' in y:
                    dic1['prog_r'] = t[1].strip()
                elif 'r_work' in y:
                    dic1['rfact_r'] = t[1].strip()
                elif 'r_free' in y:
                    dic1['rfree_r'] = t[1].strip()
                elif 'high_resolution' in y:
                    dic1['resh_r'] = t[1].strip()
                elif 'low_resolution' in y:
                    dic1['resl_r'] = t[1].strip()
                elif 'exptl_method' in y:
                    # End of category
                    break

        elif ln[0:44] == 'After applying resolution and sigma cutoffs:':
            for y in fr:
                t = y.split(':')

                if 'r_work' in y:
                    dic1['r_cutoff'] = t[1].strip()
                elif 'r_free' in y:
                    dic1['rf_cutoff'] = t[1].strip()
                elif 'Number of F-obs outliers' in y:
                    # End of category
                    break

    # print(dic['dcc'])

    fr.close()


##########################################################
def get_dcc_avg(fr, id):  # pylint: disable=redefined-builtin
    ''' get the average Biso, Occ, CC for atom list
    Reorganize the cif list
    '''

    list_in = []
    for x in fr:
        if 'standard' in x:
            break
        list_in.append(x)

    n = len(list_in)
    tmp = []
    if id == 1:
        b, cc, occ, m = 0, 0, 0, 0
        for i in range(n):
            m = m + 1
            res = list_in[i][14:17]
            ch = list_in[i][18:19]
            if not ch.strip():
                ch = '?'
            if util.is_number(list_in[i][19:23]):
                nr = int(list_in[i][19:23])
            ns = 34
            if 'segid=' in list_in[i]:
                ns = 47

            if util.is_number(list_in[i][ns:ns + 5]):
                occ = occ + float(list_in[i][ns:ns + 5])
            if util.is_number(list_in[i][ns + 5:ns + 12]):
                b = b + float(list_in[i][ns + 5:ns + 12])
            if util.is_number(list_in[i][ns + 12:ns + 20]):
                cc = cc + float(list_in[i][ns + 12:ns + 20])

            if i == n - 1 or list_in[i][19:23] != list_in[i + 1][19:23]:
                d = [nr, ch, res, cc / m, b / m, occ / m]
                tmp.append(d)
                b, cc, occ, m = 0, 0, 0, 0

    elif id == 2:
        for i in range(n):

            res = list_in[i][31:34]
            ch = list_in[i][18:19]
            if not ch.strip():
                ch = '?'
            nr = int(list_in[i][20:25])
            occ = float(list_in[i][35:40])
            b = float(list_in[i][40:47])
            cc = float(list_in[i][47:55])

            d = [nr, ch, res, cc, b, occ]
            tmp.append(d)

    return tmp


##########################################################
def xtriage_log(file, dic):
    ''' parse data from the output of xtriage'''

    if util.check_file(500, file) == 0:
        return

    fp = open(file, 'r')

    dic['twin_ph'] = '?'
    out1 = [0] * 20
    spg, twin = '?', 'Y'
    isigi_bin = []
    for ln in fp:

        if 'Space group:' in ln[:12]:
            spg = ln[12:].split('(')[0].strip()

        elif 'Resolution range:' in ln[:]:
            out1[0] = float(ln[17:].split()[1])

        elif '  Number of centrics  :' in ln[:23]:
            out1[1] = int(ln[23:].split()[0])

        elif '  Number of acentrics :' in ln[:23]:
            out1[2] = int(ln[23:].split()[0])

        elif 'p_value(height) ' in ln and ' : ' in ln:
            dic['pvalue'] = util.get_value_after_id(ln, ":")

        elif '[MaxAnisoB-MinAnisoB]' in ln[:47]:
            dic['aniso'] = util.get_value_after_id(ln, ":")

        elif 'max. difference between axes =' in ln:
            dic['aniso'] = util.get_value_after_id(ln, "=")

        elif 'The combined Z-score of ' in ln[:25]:
            t = util.get_value_after_id(ln, "Z-score of")
            dic['zscore'] = t.split()[0]

        elif 'No ice ring related problems detected.' in ln:
            dic['ice_ring'] = '0'

        elif 'There were ' in ln[:20] and 'ice ring related warnings.' in ln[10:]:
            dic['ice_ring'] = ln.split()[3]

        elif ' Patterson function reveals a significant off-origin' in ln:
            dic['stran'] = 'Y'

        elif ' No significant pseudotranslation is detected' in ln:
            dic['stran'] = 'N'

        elif 'Acentric reflections' in ln:
            n = 0
            for tmp in fp:
                if len(tmp.strip()) < 5:
                    continue
                n = n + 1
                if '   <I^2>/<I>^2    :' in tmp:
                    out1[3] = float(tmp[19:].split()[0])
                elif '   <F>^2/<F^2>    :' in tmp:
                    out1[4] = float(tmp[19:].split()[0])
                elif '   <|E^2 - 1|>    :' in tmp:
                    out1[5] = float(tmp[19:].split()[0])
                elif n >= 3:
                    break

        elif 'Centric reflections' in ln:
            n = 0
            for tmp in fp:
                if len(tmp.strip()) < 5:
                    continue
                n = n + 1
                if '   <I^2>/<I>^2    :' in tmp:
                    out1[6] = float(tmp[19:].split()[0])
                elif '   <F>^2/<F^2>    :' in tmp:
                    out1[7] = float(tmp[19:].split()[0])
                elif '   <|E^2 - 1|>    :' in tmp:
                    out1[8] = float(tmp[19:].split()[0])
                elif n >= 3:
                    break

        elif '  Mean |L|   :' in ln[:14]:
            out1[9] = float(ln[14:].split()[0])

        elif '  Mean  L^2  :' in ln[:14]:
            out1[10] = float(ln[14:].split()[0])

        elif '  Multivariate Z score L-test:' in ln[:30]:
            out1[11] = float(ln[30:].split()[0])

        elif 'ML estimate of overall B value of ' in ln:
            n = 0
            for tmp in fp:
                n = n + 1
                if 'A**' in tmp:
                    dic['bwilson'] = tmp.split()[0]
                if n >= 3:
                    dic['bwilson_s'] = tmp.split()[0]
                    break

        elif 'Statistics independent of twin laws' in ln:

            n = 0
            for tmp in fp:
                if len(tmp.strip()) < 5:
                    continue
                n = n + 1
                if ' <I^2>/<I>^2 :' in tmp:
                    dic['I2'] = tmp.split(':')[1].strip().split()[0]
                elif ' <F>^2/<F^2> : ' in tmp:
                    dic['F2'] = tmp.split(':')[1].strip().split()[0]
                elif ' <|E^2-1|>   :' in tmp:
                    dic['E2'] = tmp.split(':')[1].strip().split()[0]
                elif '<|L|>, <L^2>:' in tmp:
                    ttmp = tmp.split(':')[1]
                    dic['L'] = ttmp.split(',')[0].strip()
                    dic['L2'] = ttmp.split(',')[1].strip()

                elif '<|L|>       :' in tmp:
                    t = tmp.split(':')[1]
                    dic['L'] = t.split()[0]
                elif '<L^2>       :' in tmp:
                    t = tmp.split(':')[1]
                    dic['L2'] = t.split()[0]

                elif 'Multivariate Z score L-test:' in tmp:
                    dic['Z'] = tmp.split(':')[1].strip()
                elif n >= 6:
                    break

        elif ('No (pseudo)merohedral twin laws were found' in ln
              or (' No twinning is suspected.' in ln
                  and util.is_number(dic['twin_r']) and float(dic['twin_r']) >= 0.3)):
            twin = 'N'
            dic['twin_ph'] = 'N'

        elif '<I/sigma_I>;_Signal_to_noise' in ln[:41]:
            n = 0
            for tmp in fp:
                n = n + 1
                if n > 2 and '$$' in tmp[:4]:
                    break
                t = tmp.strip().split()
                if len(t) == 2:
                    isigi_bin.append(t)

        elif ' type | R obs. | Britton alpha | H alpha | ML alpha |' in ln:
            n = 0
            for tmp in fp:
                n = n + 1
                if n > 4:
                    break
                t = tmp.strip().split('|')

                if len(t) == 8:
                    dic['twin_op_ph'] = t[1].strip()
                    dic['twin_type'] = t[2].strip()
                    dic['twin_r'] = t[3].strip()
                    dic['twin_fr_ph'] = t[6].strip()

    if (dic['twin_ph'] != 'N' and util.is_number(dic['twin_r']) and float(dic['twin_r']) < 0.3):
        dic['twin_ph'] = 'Y'

    if len(isigi_bin) > 3:
        dic['i_sigi_resh'] = isigi_bin[-1][1]
        if (util.is_number(isigi_bin[0][1]) and util.is_number(isigi_bin[1][1])):
            a = float(isigi_bin[0][1]) - float(isigi_bin[1][1])
            dic['i_sigi_slop'] = '%.1f' % a

    fp.close()
    return spg, twin, out1


##########################################################
def model_stat(file, logfile, dic):
    ''' parse data from the output of cns (model_stats)'''

    dic['prog'] = 'CNS'
    if util.check_file(100, file):
        fp = open(file, 'r')
    elif util.check_file(100, logfile):
        fp = open(logfile, 'r')
    else:
        util.perror('Error: CNS failed!')
        return

    for x in fp:
        if 'resolution range:' in x[:18] :
            t = x.split(':')[1].split()
            dic['resl'], dic['resh'] = t[0], t[2]

        elif ' after B-factor and/or bulk solvent correction  r=' in x:
            t = x.split('r=', 1)[1].split()
            dic['rfact'], dic['rfree'] = t[0], t[2]

        elif '>>> bulk solvent: density level=' in x:
            t = x.split('level=', 1)[1].split()
            dic['k_sol'], dic['b_sol'] = t[0], t[3]

        elif 'total number of reflections used: ' in x:
            t = x.split('used:', 1)[1].split()
            dic['nref'], dic['comp'] = t[0], t[2]

        elif 'number of reflections in test set:' in x:
            t = x.split('set:', 1)[1].split()
            dic['nfree'], dic['pfree'] = t[0], t[2]

        # below from stat.log

        elif ' low_res= ' in x :
            dic['resl'] = x.split('low_res=')[1].strip().split(';')[0]

        elif ' high_res= ' in x :
            dic['resh'] = x.split('high_res=')[1].strip().split(';')[0]

        elif ' REFLection> NREFlection= ' in x:
            dic['nref'] = x.split('NREFlection=')[1].strip()

        elif 'best sol_k = ' in x and '.' in x[12:20]:
            t = x.split('=', 1)[1].strip()
            dic['k_sol'] = t

        elif 'best sol_b = ' in x and '.' in x[12:20]:
            t = x.split('=', 1)[1].strip()
            dic['b_sol'] = t

        elif 'Overall R-value for test set:' in x:
            t = x.split(':', 1)[1]
            dic['rfree'] = t.strip()

        elif 'Overall R-value for working set:' in x:
            t = x.split(':', 1)[1]
            dic['rfact'] = t.strip()

    if '?' in dic['rfact']:
        util.perror('Error: CNS failed!')
    fp.close()


##########################################################
def log_from_shelx(shelx_log, dic):
    '''parse stat from the log of shelx
    '''

    if util.check_file(100, shelx_log) == 0:
        util.perror('Error: SHELX failed!')
        return

    dic['prog'] = 'SHELX'

    fp = open(shelx_log, 'r')
    for x in fp:
        if ' R1 = ' in x and 'unique reflections' in x:
            t = x.split('=', 1)[1].split()
            dic['rfact'], dic['nref'] = t[0], t[2]

    fp.close()


##########################################################
def pointless_log(logfile, dic):
    '''parse stat from the log of pointless

    '''

    if util.check_file(100, logfile) == 0:
        return

    for x in open(logfile, 'r').readlines():
        if ('<L^2>:' in x[:11] or '< L^2 >' in x[:11]) and 'untwinned' in x:
            dic['L2_pt'] = x.split(':')[1].lstrip().split()[0]
        elif 'Best Solution: ' in x[:17] and 'space group' in x:
            dic['spg_pt'] = x.split('group')[1].strip()


##########################################################
def parse_table(alist, table, id):  # pylint: disable=redefined-builtin
    '''parse the table in the DCC:  flist is a list
    id=0, input file; id=1. input list
    '''

    if id == 0:  # a file
        flist = open(alist, 'r').readlines()
    else:
        flist = alist

    items, values = cif.cifparse(flist, table)
    rows = cif.get_rows(items, values)

    return items, rows


##########################################################
def parse_dcc(flist, table):
    '''parse the table in the DCC:  flist is a list
    '''

    dcc = []
    if table == '_pdbx_rscc_mapman.' :
        items, values = cif.cifparse(flist, '_pdbx_rscc_mapman.')

        nseq = cif.parse_values(items, values, "_pdbx_rscc_mapman.auth_seq_id")
        chid = cif.parse_values(items, values, "_pdbx_rscc_mapman.auth_asym_id")
        comp = cif.parse_values(items, values, "_pdbx_rscc_mapman.auth_comp_id")
        alt = cif.parse_values(items, values, "_pdbx_rscc_mapman.label_alt_id")
        ins = cif.parse_values(items, values, "_pdbx_rscc_mapman.label_ins_code")
        cc = cif.parse_values(items, values, "_pdbx_rscc_mapman.correlation")
        rsr = cif.parse_values(items, values, "_pdbx_rscc_mapman.real_space_R")
        zrsr = cif.parse_values(items, values, "_pdbx_rscc_mapman.real_space_Zscore")
        biso = cif.parse_values(items, values, "_pdbx_rscc_mapman.Biso_mean")
        occ = cif.parse_values(items, values, "_pdbx_rscc_mapman.occupancy_mean")
        modid = cif.parse_values(items, values, "_pdbx_rscc_mapman.model_id")
        pdbid = cif.parse_values(items, values, "_pdbx_rscc_mapman.pdb_id")

        if not items:
            return dcc
        for i in range(len(chid)):
            if (modid and int(modid[i]) > 1):
                break
            a = [nseq[i], chid[i], comp[i], alt[i], cc[i],
                 rsr[i], biso[i], occ[i], pdbid[i], zrsr[i], ins[i]]
            dcc.append(a)

    elif table == '_pdbx_map.':

        items, values = cif.cifparse(flist, '_pdbx_map.')

        nseq = cif.parse_values(items, values, "_pdbx_map.auth_seq_id")
        chid = cif.parse_values(items, values, "_pdbx_map.auth_asym_id")
        comp = cif.parse_values(items, values, "_pdbx_map.auth_comp_id")

        biso = cif.parse_values(items, values, "_pdbx_map.Biso_mean_overall")
        cc = cif.parse_values(items, values, "_pdbx_map.density_correlation_overall")
        rsr = cif.parse_values(items, values, "_pdbx_map.real_space_R_overall")
        zobs = cif.parse_values(items, values, "_pdbx_map.ZOBS_overall")
        zdiff = cif.parse_values(items, values, "_pdbx_map.ZDIFF_overall")
        zdplus = cif.parse_values(items, values, "_pdbx_map.ZDplus_overall")
        zdminus = cif.parse_values(items, values, "_pdbx_map.ZDminus_overall")

        cc_m = cif.parse_values(items, values, "_pdbx_map.density_correlation_main_chain")
        rsr_m = cif.parse_values(items, values, "_pdbx_map.real_space_R_main_chain")
        zobs_m = cif.parse_values(items, values, "_pdbx_map.ZOBS_main_chain")
        zdiff_m = cif.parse_values(items, values, "_pdbx_map.ZDIFF_main_chain")

        cc_s = cif.parse_values(items, values, "_pdbx_map.density_correlation_side_chain")
        rsr_s = cif.parse_values(items, values, "_pdbx_map.real_space_R_side_chain")
        zobs_s = cif.parse_values(items, values, "_pdbx_map.ZOBS_side_chain")
        zdiff_s = cif.parse_values(items, values, "_pdbx_map.ZDIFF_side_chain")

        if not items:
            return dcc
        for i in range(len(chid)):
            if zdiff and not zdiff_s:
                a = [comp[i], chid[i], nseq[i], biso[i], cc[i], rsr[i], zobs[i], zdiff[i], zdplus[i], zdminus[i]]
            elif zdiff and zdiff_s:
                a = [comp[i], chid[i], nseq[i], biso[i], cc[i], rsr[i], zobs[i], zdiff[i], zdplus[i], zdminus[i],
                     cc_m[i], rsr_m[i], zobs_m[i], zdiff_m[i], cc_s[i], rsr_s[i], zobs_s[i], zdiff_s[i]]
            dcc.append(a)

    return dcc  # nseq,chid,comp,alt,cc,rsr,biso,occ,pdid


##########################################################
def get_rsr_database(res, rfact, file):  # pylint: disable=unused-argument
    '''parse the mean and dev of the RsR for the res range for each
    residue.
    '''

    if not util.check_file(30, file):
        print('Warning: no RsR data base is found. ')

    flist = open(file, 'r').readlines()
    items, values = cif.cifparse(flist, '_rsr_shell.')  # a loop

    resname = cif.parse_values(items, values, "_rsr_shell.residue")
    resh = cif.parse_values(items, values, "_rsr_shell.d_res_high")
    resl = cif.parse_values(items, values, "_rsr_shell.d_res_low")
    # rfh = cif.parse_values(items, values, "_rsr_shell.rfact_high")
    # rfl = cif.parse_values(items, values, "_rsr_shell.rfact_low")
    sst = cif.parse_values(items, values, "_rsr_shell.secondary_structure")
    mean_all = cif.parse_values(items, values, "_rsr_shell.mean_all")
    dev_all = cif.parse_values(items, values, "_rsr_shell.deviation_all")
    num_all = cif.parse_values(items, values, "_rsr_shell.number_all")
    mean_fil = cif.parse_values(items, values, "_rsr_shell.mean_filter")
    dev_fil = cif.parse_values(items, values, "_rsr_shell.deviation_filter")
    num_fil = cif.parse_values(items, values, "_rsr_shell.number_filter")

    resn = {}
    for i, _x in enumerate(resname):
        # if float(resh[i]) <=res <= float(resl[i]):
        id = '%s_%s' % (resname[i], sst[i])  # pylint: disable=redefined-builtin
        # if float(resh[i]) <= res <= float(resl[i]) and float(rfl[i]) <= rfact <= float(rfh[i]):
        if float(resh[i]) <= res <= float(resl[i]):
            t = [float(mean_all[i]), float(dev_all[i]), int(num_all[i]),
                 float(mean_fil[i]), float(dev_fil[i]), int(num_fil[i])]
            if id not in resn:
                resn[id] = t

    return resn


##########################################################
def parse_edstat(out):
    '''parse the output file generated by EDSTAT. use 1st line items as key.
    '''

    edst = {}

    if util.check_file(100, out) == 0:
        return edst

    flist = open(out, 'r').readlines()
    first = []
    for x in flist:

        if 'RT ' in x and 'RN ' in x and 'CI ' in x:
            first = x.split()
            for y in first:
                edst[y] = []
        else:
            t = x.split()
            for i, y in enumerate(t):
                edst[first[i]].append(y)
    return edst


##########################################################
