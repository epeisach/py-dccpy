import util
import config
import prog
import parse


##########################################################
def map_detail(id):  # pylint: disable=redefined-builtin
    """a simple explain about the statistics."""

    detail_edstat = """
#_pdbx_edstat.details
#;Using Tickle's analysis (Acta Cryst. D. 2012) . Instruction for the cif columns

#RSCC (density correlation), RSR (real space R).
#RSZD : real space difference density Z score (defined as Delta_rho/sigma(Delta_rho))
#       The max value of (RSZD- and RSZD+). It is related to model accuracy.
#RSZO (real space observed density Z score (defined as <rho_obs>/sigma(Delta_rho))
#Biso_mean: scaled/weighted B factor by Tickle.

#RSRZ : Zscore of RSR by HW (extre outlier removed using IQR,k=2.5).
#weighted_RSRZ: Zscore of the weighted RSR=RSR/RSCC by HW.
#RSZO_Zscore: Zscore calculated from the value RSZO by HW.
#;
"""

    if id == "edstat":
        return detail_edstat


##########################################################
def edstat2cif(out):
    """The initial assignment. Only for tempary use (used by parse_dcc)"""

    dd = parse.parse_edstat(out)

    head1 = """
#
#The stats are from Tickle's anal.
#RSCCZ-> the Zscore calculated from correlation;
#RSDO=mean_density/(sigma of difference map); It is precision related.
#RSDZ=max of (RSZD- and RSZD+) ; It is accuracy related.
#RSRD-+: sets of negative and positive values of DIFF
#
loop_
 _pdbx_map.auth_comp_id
 _pdbx_map.auth_asym_id
 _pdbx_map.auth_seq_id
 _pdbx_map.Biso_mean_overall
 _pdbx_map.real_space_R_overall
 _pdbx_map.density_correlation_overall
 _pdbx_map.ZCC_overall
 _pdbx_map.ZOBS_overall
 _pdbx_map.ZDIFF_overall
 _pdbx_map.ZDplus_overall
 _pdbx_map.ZDminus_overall
"""
    head2 = """ _pdbx_map.Biso_mean_main_chain
 _pdbx_map.real_space_R_main_chain
 _pdbx_map.density_correlation_main_chain
 _pdbx_map.ZCC_main_chain
 _pdbx_map.ZOBS_main_chain
 _pdbx_map.ZDIFF_main_chain
 _pdbx_map.Biso_mean_side_chain
 _pdbx_map.real_space_R_side_chain
 _pdbx_map.density_correlation_side_chain
 _pdbx_map.ZCC_side_chain
 _pdbx_map.ZOBS_side_chain
 _pdbx_map.ZDIFF_side_chain
"""
    nres = len(dd.keys())

    outcif = out + ".cif"
    fw = open(outcif, "w")
    fw.write(head1)
    if nres > 35:
        fw.write(head2)

    for i, x in enumerate(dd["RT"]):
        idd = "%3s %2s %4s " % (x, dd["CI"][i], dd["RN"][i])
        if "CCSa" in dd.keys() and nres > 35:
            zda = zdiff(dd["ZD-a"][i], dd["ZD+a"][i])
            zdm = zdiff(dd["ZD-m"][i], dd["ZD+m"][i])
            zds = zdiff(dd["ZD-s"][i], dd["ZD+s"][i])

            # Log file label fixed in newer versions
            if "ZCCa" in dd:
                zcca = dd["ZCCa"][i]
                zccm = dd["ZCCm"][i]
                zccs = dd["ZCCs"][i]
            else:
                zcca = dd["ZCCPa"][i]
                zccm = dd["ZCCPm"][i]
                zccs = dd["ZCCPs"][i]
            ss = " %5s %5s %5s %4s %4s %4s %4s %4s  %5s %5s %5s %4s %4s %4s   %5s %5s %5s %4s %4s %4s" % (
                dd["BAa"][i],
                dd["Ra"][i],
                dd["CCSa"][i],
                zcca,
                dd["ZOa"][i],
                zda,
                dd["ZD+a"][i],
                dd["ZD-a"][i],
                dd["BAm"][i],
                dd["Rm"][i],
                dd["CCSm"][i],
                zccm,
                dd["ZOm"][i],
                zdm,
                dd["BAs"][i],
                dd["Rs"][i],
                dd["CCSs"][i],
                zccs,
                dd["ZOs"][i],
                zds,
            )

        elif "CCSa" not in dd.keys() or nres < 35:  # this is actually the -a
            zdm = zdiff(dd["ZD-m"][i], dd["ZD+m"][i])
            ss = " %5s %5s %5s %5s %5s %4s %4s %4s " % (dd["BAm"][i], dd["Rm"][i], dd["CCSm"][i], dd["ZCCm"][i], dd["ZOm"][i], zdm, dd["ZD+m"][i], dd["ZD-m"][i])

        tt = "%s %s\n" % (idd, ss)
        fw.write(tt)
    fw.close()
    return outcif


##########################################################
def zdiff(a, b):
    if util.is_number(a) and util.is_number(b) and abs(float(a)) > abs(float(b)):
        za = a
    else:
        za = b
    return za


##########################################################
def rsr_edstat(fw, dic):
    """write a new cif token for all the rsr"""

    head = """
#
loop_
 _pdbx_map.id
 _pdbx_map.auth_comp_id
 _pdbx_map.auth_asym_id
 _pdbx_map.auth_seq_id
 _pdbx_map.label_ins_code
 _pdbx_map.quality_indicator
 _pdbx_map.Biso_Zscore
 _pdbx_map.Biso_mean
 _pdbx_map.RSCC
 _pdbx_map.RSR
 _pdbx_map.RSZO
 _pdbx_map.RSZD
 _pdbx_map.RSZD_plus
 _pdbx_map.RSZD_minus
 _pdbx_map.RSRZ
 _pdbx_map.weighted_RSRZ
 _pdbx_map.RSZO_Zscore
 _pdbx_map.RSCC_main_chain
 _pdbx_map.RSR_main_chain
 _pdbx_map.RSZO_main_chain
 _pdbx_map.RSZD_main_chain
 _pdbx_map.RSCC_side_chain
 _pdbx_map.RSR_side_chain
 _pdbx_map.RSZO_side_chain
 _pdbx_map.RSZD_side_chain
 _pdbx_map.RSCC_phosphate_group
 _pdbx_map.RSR_phosphate_group
 _pdbx_map.RSZO_phosphate_group
 _pdbx_map.RSZD_phosphate_group
"""

    detail = map_detail("edstat")
    fw.write(detail)
    fw.write(head)

    cifrefm, cifrefm1, cifrefm2, cifrefm3 = dic["cifrefm"], dic["cifrefm1"], dic["cifrefm2"], dic["cifrefm3"]
    flist = open(cifrefm, "r").readlines()
    dcc = parse.parse_dcc(flist, "_pdbx_map.")

    ids = []
    for x in dcc:
        ch = x[1]
        if ch == "_":
            ch = "."
        id = "%s_%s_%s" % (x[0], ch, x[2])  # pylint: disable=redefined-builtin
        ids.append(id)

    dic, dic1, dic2, dic3 = {}, {}, {}, {}
    dic = get_dic_edstat(cifrefm, ids)
    dic1 = get_dic_edstat(cifrefm1, ids)
    dic2 = get_dic_edstat(cifrefm2, ids)
    dic3 = get_dic_edstat(cifrefm3, ids)

    wrsr_avg_w, wrsr_dev_w = get_avg_dev(ids, dic, -1, "hoh_wrsr")  # for HOH weighted (rsr/cc)
    wrsr_avg, wrsr_dev = get_avg_dev(ids, dic, -1, "wrsr")  # for weighted (rsr/cc)
    rsr_avg, rsr_dev = get_avg_dev(ids, dic, 2, "rsr ")
    zo_avg, zo_dev = get_avg_dev(ids, dic, 3, "zo  ")
    b_avg, b_dev = get_avg_dev(ids, dic, 0, "Bfactor ")
    nid = 0
    prob = []

    for x in ids:  # dic[x] ; 0=B, 1,CC; 2, RSR; 3, ZO; 4,ZD; 5-8 main chain; 9-12 side chain
        y = x.split("_")
        ins = "."
        nres = y[2]

        if not util.is_number(nres[-1]):
            nres = y[2][:-1]
            ins = nres[-1]
        resn = y[0]
        idd = "%3s %2s %4s %s " % (y[0], y[1], nres, ins)
        b, cc, rsr, zo, zd, zdp, zdm = dic[x][0], dic[x][1], dic[x][2], dic[x][3], dic[x][4], dic[x][5], dic[x][6]
        cc_m, rsr_m, zo_m, zd_m = dic1[x][1], dic1[x][2], dic1[x][3], dic1[x][4]
        cc_s, rsr_s, zo_s, zd_s = dic2[x][1], dic2[x][2], dic2[x][3], dic2[x][4]
        cc_p, rsr_p, zo_p, zd_p = dic3[x][1], dic3[x][2], dic3[x][3], dic3[x][4]

        if resn in util.protein() and len(dic[x]) > 11:
            cc_m, rsr_m, zo_m, zd_m = dic[x][5], dic[x][6], dic[x][7], dic[x][8]
            cc_s, rsr_s, zo_s, zd_s = dic[x][9], dic[x][10], dic[x][11], dic[x][12]

        zwrsr, zrsr, z_zo, z_b = ".", ".", ".", "."
        if util.is_number(cc) and abs(float(cc)) > 0 and util.is_number(rsr):
            wrsr = float(rsr) / float(cc)
            if "HOH" in resn or "DOD" in resn:
                zwrsr = get_z(wrsr, wrsr_avg_w, wrsr_dev_w)
            else:
                zwrsr = get_z(wrsr, wrsr_avg, wrsr_dev)
                zrsr = get_z(rsr, rsr_avg, rsr_dev)
                z_zo = get_z(zo, zo_avg, zo_dev)
                z_b = get_z(b, b_avg, b_dev)

        select = select_bad_residue(zrsr, zwrsr, z_zo, zd, idd)
        nid = nid + 1
        t = "%3d %s %2s  %5s %5s %5s %4s  %4s %4s %4s %4s %4s %4s %4s   %5s %5s %4s %4s   %5s %5s %4s %4s  %s %s %s %s \n" % (
            nid,
            idd,
            select,
            z_b,
            b,
            cc,
            rsr,
            zo,
            zd,
            zdp,
            zdm,
            zrsr,
            zwrsr,
            z_zo,
            cc_m,
            rsr_m,
            zo_m,
            zd_m,
            cc_s,
            rsr_s,
            zo_s,
            zd_s,
            cc_p,
            rsr_p,
            zo_p,
            zd_p,
        )
        ss = t.replace("n/a", "  .")

        if "." not in select:
            prob.append(ss)
            continue
        fw.write(ss)
    for y in prob:
        fw.write(y)


######################################################################
def select_bad_residue(rsrz, wrsrz, z_zo, zdiff, idd):  # pylint: disable=redefined-outer-name
    rsr, zd = "", ""

    if "HOH " in idd or "DOD " in idd:
        if util.is_number(wrsrz) and float(wrsrz) > 3.0:
            rsr = "w"
    else:
        if util.is_number(rsrz) and util.is_number(wrsrz) and util.is_number(z_zo):
            if (float(rsrz) > 2.0 and float(wrsrz) > 2.0) or float(rsrz) > 2.5 or float(wrsrz) > 2.5 or float(z_zo) < -2.5:
                rsr = "w"

    if util.is_number(zdiff) and abs(float(zdiff)) >= 3.0:
        zd = "d"

    ss = "%s%s" % (rsr, zd)
    if not ss:
        ss = "."

    return ss


######################################################################
def get_avg_dev(ids, dic, col, name):
    datat = []
    for x in ids:
        if "hoh_wrsr" in name:
            if "HOH" not in x:
                continue
            cc, rsr = dic[x][1], dic[x][2]
            if util.is_number(cc) and abs(float(cc)) > 0.01 and util.is_number(rsr):
                wrsr = float(rsr) / float(cc)
                datat.append([wrsr])
        else:
            if "HOH" in x:
                continue
            if col < 0:  # weighted RsR
                cc, rsr = dic[x][1], dic[x][2]
                if util.is_number(cc) and abs(float(cc)) > 0.01 and util.is_number(rsr):
                    wrsr = float(rsr) / float(cc)
                    datat.append([wrsr])
            else:
                v = dic[x][col]
                if util.is_number(v):
                    datat.append([float(v)])

    # b1, b2, q2 = util.limits_by_iqr(datat, 0, scale=3.0) #remove extreme values
    # data, outlier=filter_data(datat, [b1,b2], 0)
    # avg, dev, mini, maxi = util.mean_dev(data, 0)
    # print '%s: avg=%.3f, dev=%.3f, min=%.3f, max=%.3f, b1=%.3f, b2=%.3f, Nall=%d, Nfilt=%d' %(name,avg, dev, mini, maxi,b1, b2, len(datat), len(data))

    avg, dev, mini, maxi = util.mean_dev(datat, 0)
    print("%-9s: avg=%.3f, dev=%.3f, min=%.3f, max=%.3f, Nall=%d" % (name, avg, dev, mini, maxi, len(datat)))

    return avg, dev


#######################################################################
def get_z(val, avg, sig):
    z_zo = "."
    if sig == 0:
        return z_zo
    if util.is_number(val):
        z = (float(val) - avg) / sig
        z_zo = "%.1f" % z

    return z_zo


#####################################################################
def filter_data(data_in, shell, col):
    """filter the data according to the shell & col values."""

    data_new, outlier = [], []
    for x in data_in:
        if x[col] < shell[0] or x[col] > shell[-1]:
            outlier.append([x[0], x[col]])
        else:
            data_new.append(x)

    return data_new, outlier


##########################################################
def get_dic_edstat(cifrefm, ids):
    dic = {}
    for x in ids:
        dic[x] = [".", ".", ".", ".", ".", ".", ".", "."]

    if not util.check_file(100, cifrefm):
        return dic

    flist = open(cifrefm, "r").readlines()
    dcc = parse.parse_dcc(flist, "_pdbx_map.")

    for x in dcc:
        id = "%s_%s_%s" % (x[0], x[1], x[2])  # pylint: disable=redefined-builtin
        if len(x) < 11:
            dic[id] = [x[3], x[4], x[5], x[6], x[7], x[8], x[9]]  # for whole residue
        else:
            dic[id] = [x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17]]
    return dic


##########################################################
def rsr_edsall(fw, dic):
    """write a new cif token for all the rsr"""

    head = """
# By REFMAC/Mapman
loop_
 _pdbx_map.id
 _pdbx_map.auth_seq_id
 _pdbx_map.auth_asym_id
 _pdbx_map.auth_comp_id
 _pdbx_map.label_alt_id
 _pdbx_map.label_ins_code
 _pdbx_map.quality_indicator
 _pdbx_map.protein_motif
 _pdbx_map.Biso_Zscore
 _pdbx_map.ASA
 _pdbx_map.RSCC
 _pdbx_map.RSR
 _pdbx_map.Biso_mean
 _pdbx_map.occupancy_mean
 _pdbx_map.RSRZ
 _pdbx_map.RSCC_main_chain
 _pdbx_map.RSR_main_chain
 _pdbx_map.Biso_mean_main_chain
 _pdbx_map.RSRZ_main_chain
 _pdbx_map.RSCC_side_chain
 _pdbx_map.RSR_side_chain
 _pdbx_map.Biso_mean_side_chain
 _pdbx_map.RSRZ_side_chain
 _pdbx_map.RSCC_phosphate_group
 _pdbx_map.RSR_phosphate_group
 _pdbx_map.Biso_mean_phosphate_group
 _pdbx_map.RSRZ_phosphate_group
"""
    fw.write(head)
    cifrefm, cifrefm1, cifrefm2, cifrefm3 = dic["cifrefm"], dic["cifrefm1"], dic["cifrefm2"], dic["cifrefm3"]
    pdb = dic["pdbfile_sig"]
    asa_mean, helix, sstr = {}, [], []
    if dic["add_asa"]:
        asa_mean = prog.calc_asa_areaimol(pdb, 0)
    if dic["add_motif"]:
        helix, sstr = prog.run_motif(pdb)

    flist = open(cifrefm, "r").readlines()
    dcc = parse.parse_dcc(flist, "_pdbx_rscc_mapman.")
    ids, bfact = [], []
    for x in dcc:
        id = "%s_%s_%s_%s_%s" % (x[0], x[1], x[2], x[3], x[10])  # pylint: disable=redefined-builtin
        ids.append(id)

        if float(x[7]) == 0 or float(x[6]) > 990 or "HOH" in id:
            continue
        b = float(x[6]) * float(x[7])  # get weighted B factors
        bfact.append(b)

    dic0, dic1, dic2, dic3 = {}, {}, {}, {}
    dic0 = get_dic(cifrefm, ids)
    dic1 = get_dic(cifrefm1, ids)
    dic2 = get_dic(cifrefm2, ids)
    dic3 = get_dic(cifrefm3, ids)

    #    b1, b2, q2 = util.limits_by_iqr(bfact, -1, scale=4) #remove extre values
    #    data, outlier=filter_data(bfact, [b1,b2], 0)
    b_avg, b_dev, _minb, _maxb = util.mean_dev(bfact, -1)
    #    print 'b_avg, b_dev ', b_avg, b_dev, minb, maxb, b1, b2, q2

    prob = []
    nres, nprob = 0, 0
    for i, x in enumerate(ids):  #
        y = x.split("_")
        resseq, ch, resname, alt, ins = y[0], y[1], y[2], y[3], y[4]
        idd = "%4s %2s %3s %s %s " % (resseq, ch, resname, alt, ins)

        asa = "."
        asa_id = "%s_%s_%s" % (resname, ch, resseq)

        if asa_mean and asa_id in asa_mean:
            asa = asa_mean[asa_id]
        motif_res = "."
        if dic["add_motif"]:
            motif_res = test_motif(helix, sstr, resname, ch, int(resseq))
        bf = float(dic0[x][2]) * float(dic0[x][3])
        z_b = get_z(bf, b_avg, b_dev)
        select = select_bad_residue_eds(asa_id, dic0[x][4], dic1[x][4], dic2[x][4])
        t = "%3d %s %s %s %5s %5s %5s %5s %6s %5s %5s   %5s %5s %6s %5s   %5s %5s %6s %5s   %s %s %s %s\n" % (
            i + 1,
            idd,
            select,
            motif_res,
            z_b,
            asa,
            dic0[x][0],
            dic0[x][1],
            dic0[x][2],
            dic0[x][3],
            dic0[x][4],
            dic1[x][0],
            dic1[x][1],
            dic1[x][2],
            dic1[x][4],
            dic2[x][0],
            dic2[x][1],
            dic2[x][2],
            dic2[x][4],
            dic3[x][0],
            dic3[x][1],
            dic3[x][2],
            dic3[x][4],
        )

        if "HOH" not in asa_id and "DOD" not in asa_id:
            nres = nres + 1
        if "." not in select:
            prob.append(t)
            continue
        fw.write(t)
    for y in prob:
        if "HOH" not in y and "DOD" not in y:
            nprob = nprob + 1
        fw.write(y)

    dic["rsrz_p"] = "%5.2f" % (100 * nprob / float(nres))

    if dic["verb"] == 0:
        util.delete_file("TMP_PROMOTIF_FILE.*", "phipsi.mat", "%s_asa_out" % pdb)


######################################################################
def select_bad_residue_eds(idd, rsrz, rsrz_m, rsrz_s):
    rsr = "."

    if "HOH" in idd or "DOD" in idd:
        if util.is_number(rsrz) and float(rsrz) > 2.5:
            rsr = "w"
    else:
        if util.is_number(rsrz) and float(rsrz) > 2.0:
            rsr = "w"
        elif util.is_number(rsrz) and float(rsrz) < 2.0 and util.is_number(rsrz_s) and float(rsrz_s) > 2.0:
            rsr = "s"
        if util.is_number(rsrz) and float(rsrz) > 2.0 and util.is_number(rsrz_m) and float(rsrz_m) < 2.0:
            rsr = "s"

    return rsr


##########################################################
def get_dic(cifrefm, ids):
    dic = {}
    init_dic(dic, ids)
    if not util.check_file(100, cifrefm):
        return dic

    flist = open(cifrefm, "r").readlines()
    dcc = parse.parse_dcc(flist, "_pdbx_rscc_mapman.")
    for x in dcc:
        id = "%s_%s_%s_%s_%s" % (x[0], x[1], x[2], x[3], x[10])  # pylint: disable=redefined-builtin
        dic[id] = [x[4], x[5], x[6], x[7], x[9]]

    return dic


##########################################################
def init_dic(dic, ids):
    for x in ids:
        dic[x] = [".", ".", ".", ".", "."]


##########################################################
def dcc_residual(fout, inpfile):
    """write the DCC file to fout"""

    if not util.check_file(300, inpfile):
        return

    fr = open(inpfile, "r")
    for line in fr:  # if id==1, 3 : #sfcheck  cif files
        # if (id==3 or id==4): line=change_cif_tag(line)
        if line.lstrip()[0:5] == "data_":
            continue
        fout.write(line)

    fr.close()


##########################################################
def phenix_cc2cif(ccfile, dic):
    """parse ccfile to mmcif (three tables)"""

    print("Calculating the density correlations..")
    # if util.check_file(100, ccfile) == 0 : return ''

    fr = open(ccfile, "r")

    radius, cc_all, cc_model, cc_data = "?", "?", "?", []
    nres, nres_w, res_no, res_vweak, res_weak = "?", "?", "?", "?", "?"
    for x in fr:
        if "Radius used for map calculation:" in x:
            radius = x.split(":")[1].strip().split()[0]
        elif "Overall map correlation:" in x:
            cc_all = x.split(":")[1].strip()
        elif "Map correlation in region of model:" in x:
            cc_model = x.split(":")[1].strip()
        elif "RESIDUE       CC       ALL     MAIN     SIDE" in x:
            for y in fr:
                if len(y.strip()) < 20 or "========" in y:
                    break

                note = "."
                ln = y[:50].split()
                if len(y) > 50 and "DENSITY" in y:
                    note = y[50:].strip()

                ln.extend([note])
                cc_data.append(ln)
            for y in fr:
                if "Total residues:" in y:
                    nres = y.split(":")[1].strip()
                elif "Residues with some weak density:" in y:
                    nres_w = y.split(":")[1].strip()
                elif "Residues out of density:" in y:
                    res_no = y.split(":")[1].strip()
                elif "Residues in very weak density:" in y:
                    res_vweak = y.split(":")[1].strip()
                elif "Residues in weak density:" in y:
                    res_weak = y.split(":")[1].strip()

    outcc = "phenix_cc.cif"  # write the file
    fw = open(outcc, "w")

    cifhead1 = dcc_ciftoken("detail_phenix", dic)
    fw.write(cifhead1)

    fw.write("#\n#Overall properties from phenix:\n")
    fw.write("_pdbx_rscc_phenix_overall.radius  %s\n" % radius)
    fw.write("_pdbx_rscc_phenix_overall.correlation_all  %s\n" % cc_all)
    fw.write("_pdbx_rscc_phenix_overall.correlation_model  %s\n" % cc_model)
    fw.write("_pdbx_rscc_phenix_overall.residue_number  %s\n" % nres)
    fw.write("_pdbx_rscc_phenix_overall.residue_some_weak  %s\n" % nres_w)
    fw.write("_pdbx_rscc_phenix_overall.residue_weak_density  %s\n" % res_weak)
    fw.write("_pdbx_rscc_phenix_overall.residue_very_weak_density  %s\n" % res_vweak)
    fw.write("_pdbx_rscc_phenix_overall.residue_no_density  %s\n" % res_no)

    if len(cc_data) > 0:
        cifhead2 = dcc_ciftoken("pdbx_rscc_phenix", dic)
        fw.write(cifhead2)
        for x in cc_data:
            ss = '%3s %5s %2s . . %4s  %3s  %3s  %3s  "%s" \n' % (x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7])
            fw.write(ss)

        cifhead3 = dcc_ciftoken("pdbx_rscc_phenix_prob", dic)
        fw.write(cifhead3)
        for x in cc_data:
            if len(x[7]) < 3:
                continue
            ss = '%3s %5s %2s . . %4s  %3s  %3s  %3s  "%s" \n' % (x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7])
            fw.write(ss)

    fr.close()
    fw.close()


##########################################################
def pdbch(pdbfile):
    as2ch = {}
    if util.check_file(200, pdbfile) == 0:
        return as2ch
    fp = open(pdbfile, "r")
    mm = 0
    for x in fp:
        if "ATOM" in x[:4] or "HETA" in x[:4]:
            chid = x[20:22].strip()
            segid = x[72:76].strip()
            if chid not in as2ch.keys():
                as2ch[chid] = segid
            if len(segid) > 2:
                mm = 1
    fp.close()
    if mm == 1:
        return as2ch
    else:
        return {}


##########################################################
def eds2cif(edsfile, dic):
    """convert the EDS output file (edsfile) to mmcif format."""

    # asym2ch=pdbch(dic['pdb_new'])
    if util.check_file(200, edsfile) == 0:
        return ""
    fr = open(edsfile, "r")

    if dic["rsrd"]:
        # file = config.PATH3 + 'rsrw_stat_all.cif'
        file = config.PATH3 + "rsr_stat_all.cif"

        resh, rfact = float(dic["resh"]), float(dic["rfact"])
        rsrd = parse.get_rsr_database(resh, rfact, file)
        helix, sstr = prog.run_motif(dic["pdbfile"])

    dcc_all, sdcc_all, rsr_all, srsr_all = "", "", "", ""
    for x in fr:
        if "Average RSCC for these   :" in x:
            dcc_all = util.get_value_after_id(x, ":")
            dic["dcorr"] = dcc_all
        elif "St. dev. RSCC" in x:
            sdcc_all = util.get_value_after_id(x, ":")
            dic["dcorr_d"] = sdcc_all

        elif "Average RSR for these    :" in x:
            rsr_all = util.get_value_after_id(x, ":")
            dic["realr"] = rsr_all
        elif "St. dev. RSR " in x:
            srsr_all = util.get_value_after_id(x, ":")
            dic["realr_d"] = srsr_all

    nres = 0
    fr.seek(0)
    rscc, rscc_p = [], []

    inv_map = {}
    if "pdbfile_map" in dic:
        # Create reverse mapping for replaced PDB id with longer
        inv_map = {v: k for k, v in dic["pdbfile_map"].items()}

    for x in fr:
        if "[" in x and "]" in x and len(x) > 127:
            if "residue_id" in x:
                continue
            (resname, ch, resseq, ins, dcc, rsr, biso, ac, occ) = (x[8:11], x[11:13].strip(), x[13:17], x[17:18], x[19:26], x[26:32], x[32:40], x[59:60], x[94:100])

            # strip spaces
            if "pdbfile_map" in dic and resname.strip() in inv_map:
                resname = inv_map[resname.strip()]

            if not ch:
                ch = "."
                util.perror("Warning:  Missing Chain-ID, a dot (.) is given for %s" % x[7:19])

            if ac == "F":
                ac = "."
            if ins == " ":
                ins = "."
            zscore, zscore_fil, zscore_all = -9.99, -9.99, -9.99
            num_fil, num_all = -9.99, -9.99
            rsrw = -9.99
            motif = "?"

            nres = nres + 1

            if float(dcc) != 0 and float(dcc) > -9.0:
                rsrw = float(rsr) / float(dcc)

            if rsr_all and srsr_all and util.is_number(rsr) and util.is_number(rsr_all) and float(rsr) > -9.9 and abs(float(srsr_all)) > 0:
                if dic["rsrd"]:
                    motif = test_motif(helix, sstr, resname, ch, int(resseq))
                    resna = "%s_%s" % (resname, motif)
                    if resna in rsrd:
                        mean_all, dev_all, num_all = rsrd[resna][0], rsrd[resna][1], rsrd[resna][2]
                        mean_fil, dev_fil, num_fil = rsrd[resna][3], rsrd[resna][4], rsrd[resna][5]
                        if dev_all > 0:
                            zscore_all = (float(rsr) - mean_all) / dev_all
                        if dev_fil > 0:
                            zscore_fil = (float(rsr) - mean_fil) / dev_fil

                zscore = (float(rsr) - float(rsr_all)) / float(srsr_all)

            # dcc_v, rsr_v = float(dcc), float(rsr)
            #  criteria=((dcc_v <0.75 and rsr_v>0.3) or (dcc_v <0.75 and abs(zscore) >1.5 and rsr_v>0.2 ) or rsrw>0.35 or  (abs(zscore) >3 and rsr_v>0.3))
            criteria = zscore > 2.0
            # if asym2ch and ch in asym2ch.keys() : ch = asym2ch[ch]
            if dic["rsrd"]:
                arg = "%4s %s %s %4s %s %s %s %s %7.2f %s %s  %7.3f %7.2f %6d %7.2f %6d %s\n" % (
                    dic["pdbid"],
                    ch,
                    resname,
                    resseq,
                    ac,
                    ins,
                    dcc,
                    rsr,
                    zscore,
                    biso,
                    occ,
                    rsrw,
                    zscore_all,
                    num_all,
                    zscore_fil,
                    num_fil,
                    motif,
                )

                rscc.append(arg)
                if criteria:
                    rscc_p.append(arg)
            else:
                arg = "%4s %s %s %4s %s %s %s %s %7.2f %s %s \n" % (dic["pdbid"], ch, resname, resseq, ac, ins, dcc, rsr, zscore, biso, occ)
                rscc.append(arg)
                if criteria:
                    rscc_p.append(arg)

    fr.close()

    # --
    if len(rscc) == 0:
        print("Error: No RsR and correlation are generated!")
        return ""

    outeds = edsfile + ".cif"  # write the file
    if util.check_file(10, outeds):
        util.delete_file(outeds)

    fw = open(outeds, "w")

    cifhead1 = dcc_ciftoken("detail", dic)
    cifhead2 = dcc_ciftoken("pdbx_rscc_mapman", dic)
    fw.write(cifhead1)
    fw.write(cifhead2)

    for i, x in enumerate(rscc):
        fw.write("%s %3d\n" % (x.strip(), i + 1))

    if len(rscc_p) > 0:
        cifhead3 = dcc_ciftoken("pdbx_rscc_prob", dic)
        fw.write(cifhead3)
        for x in rscc_p:
            fw.write(x)

    if dcc_all and sdcc_all and rsr_all and srsr_all:
        fw.write("#\n#Overall properties from mapman:\n")
        fw.write("_pdbx_rscc_mapman_overall.pdbid  %s\n" % dic["pdbid"])
        fw.write("_pdbx_rscc_mapman_overall.correlation  %s\n" % dcc_all)
        fw.write("_pdbx_rscc_mapman_overall.correlation_sigma  %s\n" % sdcc_all)
        fw.write("_pdbx_rscc_mapman_overall.real_space_R  %s\n" % rsr_all)
        fw.write("_pdbx_rscc_mapman_overall.real_space_R_sigma  %s\n" % srsr_all)

    fw.close()

    return outeds


##########################################################
def test_motif(helix, sstr, resname, ch, resseq):
    """ """

    if not helix and not sstr:
        return "."
    resn = resname.strip()
    if resn not in util.protein():
        return "."

    id = "T"  # pylint: disable=redefined-builtin
    for x in helix:
        if x[0] == ch and x[1] <= resseq <= x[2]:
            return "H"

    for x in sstr:
        if x[0] == ch and x[1] <= resseq <= x[2]:
            return "S"

    return id


##########################################################
def dcc_ciftoken(id, dic):  # pylint: disable=redefined-builtin
    """contains the cif tokens for mapman and phenix.
    dic contains the cif tokens.
    """

    if id == "detail":
        cifhead = (
            """data_eds_rscc
_pdbx_dcc_mapman.pdbid  %s
_pdbx_dcc_mapman.details
;Items below are the local density correlation using mapman and refmac(Dcc).
correlation: Dcc=(<xy>-<x><y>)/[sqrt(<x**2>-<x>**2)*sqrt(<y**2>-<y>**2)]

Real spaceR: RSR = sum(|x-y|)/sum(|x+y|) sum over all grid around residue
             x=Do (observed density 2mFo-dFc); y=Dc (caculated denisty Fc)

real_space_Zscore:  RSRZ = (RSR-<RSR>)/sigma
where <RSR> is calculated from all in ASU, and sigma is by all the grid points.

Biso_mean:  occupancy-weighted average B = (SUM B*Q)/(SUM Q)
occupancy_mean:  the average occupancy of each residue = S_occ / Nuniq
;
"""
            % dic["pdbid"]
        )
        return cifhead

    if id == "pdbx_non_poly":
        cifhead1 = """
#
loop_
 _pdbx_non_poly.pdb_id
 _pdbx_non_poly.auth_asym_id
 _pdbx_non_poly.auth_comp_id
 _pdbx_non_poly.auth_seq_id
 _pdbx_non_poly.label_alt_id
 _pdbx_non_poly.correlation
 _pdbx_non_poly.real_space_R
 _pdbx_non_poly.real_space_Zscore
 _pdbx_non_poly.LLDF
 _pdbx_non_poly.Biso_mean
 _pdbx_non_poly.occupancy_mean
"""
        return cifhead1

    if id == "pdbx_rscc_mapman":
        cifhead2 = ""

        cifhead1 = """
#
loop_
 _pdbx_rscc_mapman.pdb_id
 _pdbx_rscc_mapman.auth_asym_id
 _pdbx_rscc_mapman.auth_comp_id
 _pdbx_rscc_mapman.auth_seq_id
 _pdbx_rscc_mapman.label_alt_id
 _pdbx_rscc_mapman.label_ins_code
 _pdbx_rscc_mapman.correlation
 _pdbx_rscc_mapman.real_space_R
 _pdbx_rscc_mapman.real_space_Zscore
 _pdbx_rscc_mapman.Biso_mean
 _pdbx_rscc_mapman.occupancy_mean
 _pdbx_rscc_mapman.id
"""
        if dic["rsrd"]:
            cifhead2 = """ _pdbx_rscc_mapman.RsR_over_correlation
 _pdbx_rscc_mapman.secondary_structure
"""
        return cifhead1 + cifhead2

    if id == "pdbx_rscc_prob":
        cifhead2 = ""

        cifhead1 = """
#
loop_
 _pdbx_rscc_mapman_prob.pdb_id
 _pdbx_rscc_mapman_prob.auth_asym_id
 _pdbx_rscc_mapman_prob.auth_comp_id
 _pdbx_rscc_mapman_prob.auth_seq_id
 _pdbx_rscc_mapman_prob.label_alt_id
 _pdbx_rscc_mapman_prob.label_ins_code
 _pdbx_rscc_mapman_prob.correlation
 _pdbx_rscc_mapman_prob.real_space_R
 _pdbx_rscc_mapman_prob.real_space_Zscore
 _pdbx_rscc_mapman_prob.Biso_mean
 _pdbx_rscc_mapman_prob.occupancy_mean
"""
        if dic["rsrd"]:
            cifhead2 = """ _pdbx_rscc_mapman_prob.RsR_over_correlation
 _pdbx_rscc_mapman_prob.secondary_structure
"""

        return cifhead1 + cifhead2

    # -----------------below for phenix ------
    if id == "detail_phenix":
        cifhead = """data_phenix_rscc

_pdbx_dcc_phenix.details
; Some Definitions:
"out of density":  RHO < 2 SD below 1/2 mean density for this atom/group
"very weak density":  RHO < 1 SD below 1/2 mean density for this atom/group
"weak density":  RHO < 1/2 mean density for this atom/group
"acceptable density":  No atom is < 1 SD below 1/2 mean for that atom type and
no group is < 1/2 mean for that group

"DENSITY" is the electron density in the map, normalized so that
 the mean and SD in the region of the model are zero and 1.0, respectively.
;
"""
        return cifhead

    if id == "pdbx_rscc_phenix":
        cifhead1 = """
#
loop_
 _pdbx_rscc_phenix.auth_comp_id
 _pdbx_rscc_phenix.auth_seq_id
 _pdbx_rscc_phenix.auth_asym_id
 _pdbx_rscc_phenix.label_alt_id
 _pdbx_rscc_phenix.label_ins_code
 _pdbx_rscc_phenix.correlation
 _pdbx_rscc_phenix.mean_sigma_all
 _pdbx_rscc_phenix.mean_sigma_mainchain
 _pdbx_rscc_phenix.mean_sigma_sidechain
 _pdbx_rscc_phenix.details
"""
        return cifhead1

    if id == "pdbx_rscc_phenix_prob":
        cifhead1 = """
#
loop_
 _pdbx_rscc_phenix_prob.auth_comp_id
 _pdbx_rscc_phenix_prob.auth_seq_id
 _pdbx_rscc_phenix_prob.auth_asym_id
 _pdbx_rscc_phenix_prob.label_alt_id
 _pdbx_rscc_phenix_prob.label_ins_code
 _pdbx_rscc_phenix_prob.correlation
 _pdbx_rscc_phenix_prob.mean_sigma_all
 _pdbx_rscc_phenix_prob.mean_sigma_mainchain
 _pdbx_rscc_phenix_prob.mean_sigma_sidechain
 _pdbx_rscc_phenix_prob.details
"""

        return cifhead1


# ######################################################### config.VERSION
def summary1(fout, dic):
    fout.write("\n\n####################Overall values#######################\n#\n")
    fout.write("_pdbx_density.DCC_version    '%s' \n" % config.VERSION)
    fout.write("_pdbx_density.pdbid     %s \n" % dic["pdbid"])
    fout.write("_pdbx_density.pdbtype  '%s' \n" % dic["pdbhead"])
    fout.write("_pdbx_density.unit_cell  '%s' \n" % (" ".join(dic["cell"].split())))
    fout.write("_pdbx_density.space_group_name_H-M                  '%s' \n" % (dic["spg"]))
    fout.write("_pdbx_density.space_group_pointless                 '%s' \n" % dic["spg_pt"])
    fout.write("_pdbx_density.ls_d_res_high                         %s \n" % dic["resh"])
    fout.write("_pdbx_density.ls_d_res_high_sf                      %s \n" % dic["resh_sf"])
    fout.write("_pdbx_density.ls_d_res_low_sf                       %s \n" % dic["resl_sf"])
    fout.write("_pdbx_density.R_value_R_work                        %s \n" % dic["rfact"])
    fout.write("_pdbx_density.R_value_R_free                        %s \n" % dic["rfree"])
    fout.write("_pdbx_density.working_set_count                     %s \n" % dic["nref"])
    fout.write("_pdbx_density.free_set_count                        %s \n" % dic["nfree"])
    fout.write("_pdbx_density.occupancy_min                         %s \n" % dic["occ_min"])
    fout.write("_pdbx_density.occupancy_max                         %s \n" % dic["occ_max"])
    fout.write("_pdbx_density.occupancy_mean                        %s \n" % dic["occ_avg"])
    fout.write("_pdbx_density.Biso_min                              %s \n" % dic["bmin"])
    fout.write("_pdbx_density.Biso_max                              %s \n" % dic["bmax"])
    fout.write("_pdbx_density.Biso_mean                             %s \n" % dic["biso"])

    fout.write("_pdbx_density.B_wilson                              %s \n" % dic["bwilson"])
    fout.write("_pdbx_density.B_wilson_scale                        %s \n" % dic["bwilson_s"])
    fout.write("_pdbx_density.mean_I2_over_mean_I_square            %s \n" % dic["I2"])
    fout.write("_pdbx_density.mean_F_square_over_mean_F2            %s \n" % dic["F2"])
    fout.write("_pdbx_density.mean_E2_1_abs                         %s \n" % dic["E2"])
    fout.write("_pdbx_density.Padilla-Yeates_L_mean                 %s \n" % dic["L"])
    fout.write("_pdbx_density.Padilla-Yeates_L2_mean                %s \n" % dic["L2"])
    fout.write("_pdbx_density.Padilla-Yeates_L2_mean_pointless      %s \n" % dic["L2_pt"])
    fout.write("_pdbx_density.Z_score_L_test                        %s \n" % dic["Z"])
    fout.write("_pdbx_density.twin_type                             %s \n" % dic["twin_type"])
    fout.write("_pdbx_density.twin_operator_xtriage                 %s \n" % dic["twin_op_ph"])
    fout.write("_pdbx_density.twin_fraction_xtriage                 %s \n" % dic["twin_fr_ph"])
    fout.write("_pdbx_density.twin_Rfactor                          %s \n" % dic["twin_r"])

    fout.write("_pdbx_density.I_over_sigI_resh                      %s \n" % dic["i_sigi_resh"])
    fout.write("_pdbx_density.I_over_sigI_diff                      %s \n" % dic["i_sigi_slop"])
    fout.write("_pdbx_density.I_over_sigI_mean                      %s \n" % dic["i_sigi_mean"])
    fout.write("_pdbx_density.ice_ring                              %s \n" % dic["ice_ring"])
    fout.write("_pdbx_density.anisotropy                            %s \n" % dic["aniso"])
    fout.write("_pdbx_density.Z-score                               %s \n" % dic["zscore"])
    fout.write("_pdbx_density.prob_peak_value                       %s \n" % dic["pvalue"])
    fout.write("_pdbx_density.translational_pseudo_symmetry         %s \n" % dic["stran"])

    fout.write("_pdbx_density.wavelength                            %s \n" % dic["wavelength"])

    if "b_sol" in dic.keys():
        fout.write("_pdbx_density.B_solvent                             %s \n" % dic["b_sol"])
    if "k_sol" in dic.keys():
        fout.write("_pdbx_density.K_solvent                             %s \n" % dic["k_sol"])

    fout.write("_pdbx_density.TLS_refinement_reported               %s \n" % dic["tls_report"])
    fout.write("_pdbx_density.partial_B_value_correction_attempted  %s \n" % dic["tls_attempt"])
    fout.write("_pdbx_density.partial_B_value_correction_success    %s \n" % dic["tls_ok"])
    fout.write("_pdbx_density.reflection_status_archived            %s \n" % dic["status"])
    fout.write("_pdbx_density.reflection_status_used                %s \n" % dic["status_used"])
    fout.write("_pdbx_density.iso_B_value_type                      %s \n" % dic["isob"])
    fout.write("_pdbx_density.reflns_twin                           %s \n" % dic["twin"])
    fout.write("_pdbx_density.twin_by_xtriage                       %s \n" % dic["twin_ph"])
    fout.write("_pdbx_density.twin_operator    '%s ' \n" % dic["twin_op"])
    fout.write("_pdbx_density.twin_fraction    '%s ' \n" % dic["twin_fr"])

    fout.write("_pdbx_density.tls_group_number                      %d \n" % dic["tls"])
    fout.write("_pdbx_density.ncs_group_number                      %d \n" % dic["ncs"])
    fout.write("_pdbx_density.mtrix_number                          %d \n" % dic["mtr"])
    fout.write("_pdbx_density.Matthew_coeff                         %s \n" % dic["matt"])
    fout.write("_pdbx_density.solvent_content                       %s \n" % dic["solvent"])
    fout.write("_pdbx_density.Cruickshank_dpi_xyz                   %s \n" % dic["dpi_xyz"])
    fout.write("_pdbx_density.dpi_free_R                            %s \n" % dic["dpi_fr"])
    fout.write("_pdbx_density.fom                                   %s \n" % dic["fom"])
    fout.write("_pdbx_density.correlation_overall                   %s \n" % dic["dcorr"])
    fout.write("_pdbx_density.real_space_R_overall                  %s \n" % dic["realr"])
    fout.write("_pdbx_density.rsrz_percent                          %s \n" % dic["rsrz_p"])

    fout.write("_pdbx_density.mFo-DFc-3sigma_positive               %s \n" % dic["fofc_3sig_p"])
    fout.write("_pdbx_density.mFo-DFc-6sigma_positive               %s \n" % dic["fofc_6sig_p"])
    fout.write("_pdbx_density.mFo-DFc-3sigma_negative               %s \n" % dic["fofc_3sig_n"])
    fout.write("_pdbx_density.mFo-DFc-6sigma_negative               %s \n" % dic["fofc_6sig_n"])

    if util.is_number(dic["biso"]) and util.is_number(dic["bwilson"]):
        n = float(dic["biso"]) - float(dic["bwilson"])
        dic["bm_bw"] = "%.3f" % n
        fout.write("_pdbx_density.Bmean-Bwilson                        %s \n" % dic["bm_bw"])

    if util.is_number(dic["rfact"]) and util.is_number(dic["rfree"]):
        n = float(dic["rfree"]) - float(dic["rfact"])
        dic["rf_rw"] = "%.4f" % n
        fout.write("_pdbx_density.Rfree-Rwork                          %s \n" % dic["rf_rw"])

    fout.write("_pdbx_density.error                      \n;\n")

    for x in config.ERRLOG:
        fout.write(x.strip())
        fout.write("\n")
    fout.write("\n;\n")

    if "diags" in dic.keys():  # to generate a log file for D&A
        f = dic["diags"]
        fw = open(f, "w")
        if len(config.ERRLOG):
            for x in config.ERRLOG:
                fw.write(x.strip())
                fw.write("\n")
        else:
            fw.write("No Warning/Error were found.\n")
        fw.close()


# ######################################################### config.VERSION
def summary2(fout, dic):
    """This function exports the geometry quality by Molprobity through phenix"""

    fout.write("\n#\n")
    fout.write("_pdbx_geometry.pdbid                            %s \n" % dic["pdbid"])
    fout.write("_pdbx_geometry.Ramachandran_outlier_percent     %s \n" % dic["outlier"])
    fout.write("_pdbx_geometry.Ramachandran_outlier_number      %s \n" % dic["outlier_num"])
    fout.write("_pdbx_geometry.Ramachandran_allowed_percent     %s \n" % dic["allowed"])
    fout.write("_pdbx_geometry.Ramachandran_allowed_number      %s \n" % dic["allowed_num"])
    fout.write("_pdbx_geometry.Ramachandran_favored_percent     %s \n" % dic["favored"])
    fout.write("_pdbx_geometry.Ramachandran_favored_number      %s \n" % dic["favored_num"])
    fout.write("_pdbx_geometry.rotamer_outliers_percent         %s \n" % dic["rot_outlier"])
    fout.write("_pdbx_geometry.rotamer_outliers_number          %s \n" % dic["rot_outlier_num"])
    fout.write("_pdbx_geometry.cbeta_deviations                 %s \n" % dic["cbeta"])
    fout.write("_pdbx_geometry.all_atom_clashscore              %s \n" % dic["clashscore"])
    fout.write("_pdbx_geometry.overall_score                    %s \n" % dic["oscore"])
    fout.write("_pdbx_geometry.bond_overall_rms                 %s \n" % dic["bond"])
    fout.write("_pdbx_geometry.bond_overall_max                 %s \n" % dic["bond_max"])
    fout.write("_pdbx_geometry.bond_ligand_rms                  %s \n" % dic["bond_lig"])
    fout.write("_pdbx_geometry.bond_ligand_max                  %s \n" % dic["bond_lig_max"])
    fout.write("_pdbx_geometry.angle_overall_rms                %s \n" % dic["angle"])
    fout.write("_pdbx_geometry.angle_overall_max                %s \n" % dic["angle_max"])
    fout.write("_pdbx_geometry.angle_ligand_rms                 %s \n" % dic["angle_lig"])
    fout.write("_pdbx_geometry.angle_ligand_max                 %s \n" % dic["angle_lig_max"])
    fout.write("_pdbx_geometry.dihedral_overall_rms             %s \n" % dic["dihedral"])
    fout.write("_pdbx_geometry.dihedral_overall_max             %s \n" % dic["dihedral_max"])
    fout.write("_pdbx_geometry.chirality_overall_rms            %s \n" % dic["chirality"])
    fout.write("_pdbx_geometry.chirality_overall_max            %s \n" % dic["chirality_max"])
    fout.write("_pdbx_geometry.planarity_overall_rms            %s \n" % dic["planarity"])
    fout.write("_pdbx_geometry.planarity_overall_max            %s \n" % dic["planarity_max"])
    fout.write("_pdbx_geometry.non-bonded_rms                   %s \n" % dic["non_bond"])


##########################################################
def write_software(fw, dic):
    """write the wrapped software dic['software'] to the output (fw) of DCC."""

    items = """\n#\nloop_
_pdbx_dcc_software.ordinal
_pdbx_dcc_software.name
_pdbx_dcc_software.version
_pdbx_dcc_software.description
"""

    if not dic["software"]:
        print("Warning: No software program is reported!")
        return
    else:
        fw.write("%s" % items)

    n = 0
    for x in dic["software"]:
        vers = x[1].strip().replace(" ", "")
        if "refinement program reported" in x[2]:
            continue
        n = n + 1
        ss = '%d  %-10s  %-8s "%s"' % (n, x[0].strip(), vers, x[2].strip())
        fw.write("%s \n" % ss)


##########################################################
def write_dcc_sf(fw, dic):
    """write the sf statistics (through sf_convert)"""

    pdbid = "_pdbx_dcc_sf.id"
    section_id = "_pdbx_dcc_sf.section_id"
    sf_nobs = "_pdbx_dcc_sf.number_reflns_obs"
    sf_nwork = "_pdbx_dcc_sf.number_reflns_R_work"
    sf_nfree = "_pdbx_dcc_sf.number_reflns_R_free"
    sf_pfree = "_pdbx_dcc_sf.percent_reflns_R_free"
    sf_nfpair = "_pdbx_dcc_sf.number_reflns_Friedel_pair_F"
    sf_nfpair_p = "_pdbx_dcc_sf.number_reflns_F_plus"
    sf_nfpair_m = "_pdbx_dcc_sf.number_reflns_F_minus"
    sf_nfpair_pm = "_pdbx_dcc_sf.number_reflns_F_all"
    sf_nipair = "_pdbx_dcc_sf.number_reflns_Friedel_pair_I"
    sf_nipair_p = "_pdbx_dcc_sf.number_reflns_I_plus"
    sf_nipair_m = "_pdbx_dcc_sf.number_reflns_I_minus"
    sf_nipair_pm = "_pdbx_dcc_sf.number_reflns_I_all"

    fw.write("\n#\n")

    if dic["pdbid"]:
        fw.write("%-43s %s \n" % (pdbid, dic["pdbid"]))
    if dic["pdbid"]:
        fw.write("%-43s r%ssf \n" % (section_id, dic["pdbid"]))
    if dic["sf_nobs"]:
        fw.write("%-43s %s \n" % (sf_nobs, dic["sf_nobs"]))
    if dic["sf_nwork"]:
        fw.write("%-43s %s \n" % (sf_nwork, dic["sf_nwork"]))
    if dic["sf_nfree"]:
        fw.write("%-43s %s \n" % (sf_nfree, dic["sf_nfree"]))
    if dic["sf_pfree"]:
        fw.write("%-43s %s \n" % (sf_pfree, dic["sf_pfree"]))
    if dic["sf_nfpair"]:
        fw.write("%-43s %s \n" % (sf_nfpair, dic["sf_nfpair"]))
    if dic["sf_nfpair_p"]:
        fw.write("%-43s %s \n" % (sf_nfpair_p, dic["sf_nfpair_p"]))
    if dic["sf_nfpair_m"]:
        fw.write("%-43s %s \n" % (sf_nfpair_m, dic["sf_nfpair_m"]))
    if dic["sf_nfpair_pm"]:
        fw.write("%-43s %s \n" % (sf_nfpair_pm, dic["sf_nfpair_pm"]))
    if dic["sf_nipair"]:
        fw.write("%-43s %s \n" % (sf_nipair, dic["sf_nipair"]))
    if dic["sf_nipair_p"]:
        fw.write("%-43s %s \n" % (sf_nipair_p, dic["sf_nipair_p"]))
    if dic["sf_nipair_m"]:
        fw.write("%-43s %s \n" % (sf_nipair_m, dic["sf_nipair_m"]))
    if dic["sf_nipair_pm"]:
        fw.write("%-43s %s \n" % (sf_nipair_pm, dic["sf_nipair_pm"]))
    fw.write("\n#\n")


##########################################################
def final_cif_item(fout, dic, ordinal, dic_vers):
    """write the final items in a loop"""

    vers = "?"
    for x in dic_vers:
        if ordinal == 1:
            if "coordinate" in x[2] and dic["prog"] in x[0]:
                vers = x[1].replace(" ", "")
        elif ordinal > 1:
            if dic["prog"] in x[0]:
                vers = x[1].replace(" ", "")

    final = """#
# the final items
#
loop_
_pdbx_density_corr.ordinal
_pdbx_density_corr.program
_pdbx_density_corr.ls_d_res_high
_pdbx_density_corr.ls_d_res_low
_pdbx_density_corr.ls_R_factor_R_all
_pdbx_density_corr.ls_R_factor_R_work
_pdbx_density_corr.ls_R_factor_R_free
_pdbx_density_corr.ls_number_reflns_obs
_pdbx_density_corr.ls_percent_reflns_obs
_pdbx_density_corr.ls_number_reflns_R_free
_pdbx_density_corr.correlation_coeff_Fo_to_Fc
_pdbx_density_corr.real_space_R
_pdbx_density_corr.correlation
_pdbx_density_corr.program_version
_pdbx_density_corr.details
"""
    if len(dic) == 0:
        fout.write(final)
    else:
        rsr = dic["realr"]
        rscc = dic["dcorr"]
        if ordinal == 1:
            rsr, rscc = "?", "?"
        fout.write(
            "%d %8s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %8s '%s'\n"
            % (
                ordinal,
                dic["prog"],
                dic["resh"],
                dic["resl"],
                dic["rall"],
                dic["rfact"],
                dic["rfree"],
                dic["nref"],
                dic["comp"],
                dic["nfree"],
                dic["fcc"],
                rsr,
                rscc,
                vers,
                dic["detail"],
            )
        )


# #########################################################config.VERSION
