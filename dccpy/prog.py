import os
import sys
import shutil
import util
import config
import parse
import cns_inp
import dcc_calc as main
import cifparse as cif
import write_cif


# ===========================================================
# this module is a interface for all the external packages
# ===========================================================


##########################################################
def sf_convertor(sffile, pdbfile, type):  # pylint: disable=redefined-builtin
    """convert to any type of formats"""

    print("\nConverting SF file to %s format ... " % type)
    out = sffile + "_" + type
    if not util.check_file(500, pdbfile):
        arg = "%ssf_convert -o %s -sf %s -out %s>/dev/null " % (config.PATH1, type, sffile, out)
    else:
        arg = "%ssf_convert -o %s -sf %s -pdb %s -out %s>/dev/null " % (config.PATH1, type, sffile, pdbfile, out)
        # print 'arg=',sffile, pdbfile, arg
    os.system(arg)
    return out


##########################################################
def calc_fem(pdbfile, sffile):
    """calc rsr/cc from feature enhanced map (phenix)"""

    util.delete_file("fem.mtz")
    print("Calculating the feature enhanced map by phenix.fem")
    # arg='phenix.fem  %s  %s  ' %(pdbfile, sffile)
    os.system("phenix.fem  %s  %s >/dev/null " % (pdbfile, sffile))
    fem_mtz = "fem.mtz"
    return fem_mtz


##########################################################
def calc_fc(pdbfile, sffile):
    """calc FC/PHIC"""

    sfout = "%s_fc.mtz" % pdbfile

    if util.is_cif(sffile):
        mtz = sf_convertor(sffile, pdbfile, "mtz")
    else:
        mtz = sffile

    arg = """#!/bin/sh

sfall xyzin %s hklin %s hklout %s <<eof-sfalli >/dev/null
        TITLE Fc and Phic calculation
        MODE SFCALC HKLIN XYZIN
        LABIN FP=FP SIGFP=SIGFP
        LABOUT FC=FC  PHIC=PHIC
        END
eof-sfall
    """ % (
        pdbfile,
        mtz,
        sfout,
    )

    scr_name = "add_fc2mtz.sh"
    fw = open(scr_name, "w")
    fw.write(arg)
    fw.close()
    command = "chmod +x %s ; ./%s  " % (scr_name, scr_name)
    os.system(command)
    return sfout


##########################################################
def cat_mtz(mtz1, mtz2, mtz1_lab, mtz2_lab, sfout):
    """cat two mtzs with given labels for mtz1, mtz2"""

    lab1 = ["E%d=%s" % (i + 1, x) for i, x in enumerate(mtz1_lab)]
    lab2 = ["E%d=%s" % (i + 1, x) for i, x in enumerate(mtz2_lab)]

    itm1 = " ".join(lab1)
    itm2 = " ".join(lab2)

    arg = """#!/bin/csh -f

cad  hklin1 %s hklin2 %s hklout %s <<eof >/dev/null
LABIN FILE 1 %s
LABIN FILE 2 %s
eof
    """ % (
        mtz1,
        mtz2,
        sfout,
        itm1,
        itm2,
    )

    scr_name = "cad_mtz.csh"
    fw = open(scr_name, "w")
    fw.write(arg)
    fw.close()
    command = "chmod +x %s ; ./%s  " % (scr_name, scr_name)
    os.system(command)


##########################################################
def merge_mtz(mtz1, mtz2, sfout):
    arg = """#!/bin/csh -f

mtzutils hklin1 %s hklin2 %s hklout %s <<eof >/dev/null
eof
    """ % (
        mtz1,
        mtz2,
        sfout,
    )

    scr_name = "merge_mtz.csh"
    fw = open(scr_name, "w")
    fw.write(arg)
    fw.close()
    command = "chmod +x %s ; ./%s  " % (scr_name, scr_name)
    os.system(command)


##########################################################
def run_buster(pdbfile, sffile, dic):
    """type refine to get all the options
    nthreads = 6 optimize   (UsePdbchk="no")
    """
    if "BDG_home" not in os.environ:
        print("\nError: No env is found. To run Buster program, please source 'setup.csh or setup.sh'.\n")
        return " "
    res = " "
    if util.is_number(dic["resh"]) and util.is_number(dic["resl"]):
        res = " -R %s %s " % (dic["resl"], dic["resh"])
    dir = "buster_%d" % os.getgid()  # pylint: disable=redefined-builtin
    arg = "refine -p %s  -m %s -d %s -nbig 1 -nsmall 0 -nthreads  6  %s  StopOnGellySanityCheckError=no " % (pdbfile, sffile, dir, res)
    os.system(arg)


##########################################################
def run_ccp4(pdbfile, sffile, type1):  # pylint: disable=unused-argument
    """run sub-programs of CCP4"""
    outfile = "ccp4__%s.log" % type1
    util.delete_file(outfile)
    print("\nDoing %s ..." % type1)

    if type1 == "ctruncate":
        arg = 'ctruncate -hklin %s -amplitudes -colin "/*/*/[FP,SIGFP]" -hklout ctruncate-SF.mtz>& %s' % (sffile, outfile)
        os.system(arg)

    elif type1 == "pointless":
        arg = "pointless hklin %s > %s" % (sffile, outfile)
        os.system(arg)
        os.system('grep "Best Solution" %s ' % outfile)

    return outfile


##########################################################
def calc_asa_areaimol(file, id):  # pylint: disable=redefined-builtin
    """calculate solvent accessible area"""

    print("Calulating ASA by areaimol (%s)..." % file)

    outf = file + "_asa_out"
    arg = """
areaimol XYZIN %s  XYZOUT %s  <<eof-area >/dev/null
VERB      ! Verbose output
#SMODE IMOL
#PROBE 1.5  #defaut 1.4A
#OUTPUT  RESIDUE  ! Output  pseudo-pdb file,  total ASA in residue
OUTPUT ATOM ! Output pseudo-pdb file, area for each atom
END
eof-area

    """ % (
        file,
        outf,
    )

    os.system("%s " % arg)
    mean_asa = get_resid_percent_asa(outf, id)

    return mean_asa


##########################################################
def get_resid_percent_asa(outf, id):  # pylint: disable=redefined-builtin
    """get the percentage of the ASA for each residue"""

    asa_perc = {}
    resid = {}
    keys = []

    if not util.check_file(10, outf):
        return asa_perc

    fp = open(outf, "r")
    for x in fp:
        if "ATOM" not in x[:4] and "HETATM" not in x[:6]:
            continue
        comp, ch, nres = x[17:20].strip(), x[20:22].strip(), x[22:26].strip()
        key = "%s_%s_%s" % (comp, ch, nres)
        if key not in resid:
            resid[key] = []
        if key in resid:
            if key not in keys:
                keys.append(key)
            resid[key].append(float(x[60:67]))

    for key in keys:
        val = sum(resid[key]) / len(resid[key])
        asa_perc[key] = "%5.2f" % val

    fp.close()
    if id == 1:
        out = "%s_mean_asa" % id
        fw = open(out, "w")
        for x in keys:
            ss = "%15s  %s\n" % (x, asa_perc[x])
            fw.write(ss)
        fw.close()
    return asa_perc


##########################################################
def run_motif(pdbfile_in):
    """calculate the helix and strand for protein."""

    print("\nRuning motif (%s)..." % (pdbfile_in))
    base = "TMP_PROMOTIF_FILE"
    newpdb = base + ".pdb"
    shutil.copy(pdbfile_in, newpdb)
    arg = "%smotif %s >/dev/null" % (config.PATH1, newpdb)
    os.system(arg)
    hf = base + ".hlx"
    sf = base + ".str"
    # print 'arg=', arg, hf, sf

    helix = []
    if util.check_file(10, hf):
        for x in open(hf, "r").readlines():
            t = x.split()
            if len(t) > 14 and util.is_number(t[0]) and util.is_number(t[3]) and util.is_number(t[5]):
                helix.append([t[2], int(t[3]), int(t[5])])
            if "Interactions:" in t:
                break

    sstr = []
    if util.check_file(10, sf):
        for x in open(sf, "r").readlines():
            t = x.split()

            if len(t) > 14 and util.is_number(t[0]) and util.is_number(t[2]) and util.is_number(t[4]):
                sstr.append([t[1], int(t[2]), int(t[4])])

    # util.delete_file('TMP_PROMOTIF_FILE.???')
    return helix, sstr


##########################################################
def run_sfcheck(pdbfile_in, sffile, dic, detail):
    print("\nRuning SFCHECK (%s, %s )..." % (pdbfile_in, sffile))
    ciffile = "extract_map_sfcheck.mmcif"
    log = "sfcheck.log"

    idd = util.is_cif(pdbfile_in)

    pdbfile = pdbfile_in
    if idd:
        pdbfile = cif.cif2pdb(pdbfile_in)

    if util.is_number(dic["resh"]):
        arg = "%ssf_convert  -o mmcif  -sf %s -cut %s >/dev/null" % (config.PATH1, sffile, dic["resh"])
        os.system(arg)
        util.delete_file("SF_4_validate.cif.mmcif")

    util.delete_file(ciffile, log)
    arg = "%ssfcheck -pdb %s -sf %s > /dev/null" % (config.PATH1, pdbfile, sffile)
    os.system(arg)

    dic1 = main.initialize_dic1()
    if not util.check_file(300, ciffile, log):
        print("Error: SFCHECK failed.(%s)\n" % pdbfile)
        return dic1, ""

    parse.sfcheck_log(log, dic1, detail)
    parse.sfcheck_cif(ciffile, dic1, detail)
    cif_new = ciffile + "_" + detail

    util.move(ciffile, cif_new)
    if idd:
        util.delete_file(pdbfile)
    return dic1, cif_new


##########################################################
def run_phenix(dic, dic1, type1):
    """use the static model_vs_data in phenix (Default params:)
        f_obs_label = None
        r_free_flags_label = None
        scattering_table = wk1995 it1992 *n_gaussian neutron
        map = None
        high_resolution = None
        comprehensive = False
        dump_result_object_as_pickle = False
        ignore_giant_models_and_datasets = True
        skip_twin_detection = False
        unmerged_data = None
        unmerged_labels = None
        n_bins = 20

      refinement.main.scattering_table=electron
    ---- implement CC ---
    1. phenix.model_vs_data  pdb sf map="2mFo-DFc"  >log
    2. phenix.get_cc_mtz_pdb pdb mtz labin='FP=2mFobs-DFmodel PHIB=P2mFobs-DFmodel'
    3. parse cc.log
    """

    main.check_env(dic)
    vers = util.get_phenix_version_from_dic(dic)
    # print "Version is ", vers

    pdbfile = dic["pdbfile"]

    if dic.get("xyzfile_cif", False):
        pdbfile = dic["xyzfile_orig"]
        # print("Using original PDBx/mmCIF file %s" % pdbfile)

    if (dic["sf_f"] and not dic["sf_sigf"]) or (dic["sf_i"] and not dic["sf_sigi"]):
        sffile = dic["sffile"]  # (sigma added)
    else:
        sffile = dic["sf_xray"]
    if "neut" in type1 and (dic["exp"] == "xn" or dic["exp"] == "nx"):
        sffile = dic["sf_neut"]

    mtz = sf_convertor(sffile, pdbfile, "mtz")

    label = ""
    if dic["sf_f"]:
        label = "FP,SIGFP"
    elif dic["sf_i"]:
        label = "I,SIGI"
    elif dic["anom_f"]:
        label = "F(+),SIGF(+),F(-),SIGF(-)"
    elif dic["anom_i"]:
        label = "I(+),SIGI(+),I(-),SIGI(-)"

    if not label:
        util.perror("Error: Can not find proper labels (F/I/F+F-/I+I-).Stopped!")
        sys.exit()

    outfile = "phenix__%s.log" % type1
    util.delete_file(outfile)
    print("\nDoing PHENIX calculation for (%s) ..." % type1)

    if type1 == "xtriage":
        arg = "phenix.xtriage %s scaling.input.xray_data.obs_labels='%s' log=%s >/dev/null" % (mtz, label, outfile)
        os.system(arg)
        parse.xtriage_log(outfile, dic1)

    elif type1 == "xray":
        if vers >= 1013:
            run_phenix_v13_xray(dic, dic1, pdbfile, mtz, label, outfile)
        else:
            dic1["prog"] = "PHENIX"

            resol = " "
            if "?" not in dic["resh"] and "n" not in dic["exp"]:
                resol = "high_resolution=%s" % dic["resh"]

            lib = " "
            if dic["lib"]:
                lib = dic["lib"]
            maps = ""
            if dic["rsrd"]:
                maps = ' map="2mFo-DFc" '

            arg = "phenix.model_vs_data %s %s %s %s %s f_obs_label='%s' > %s " % (pdbfile, mtz, lib, resol, maps, label, outfile)
            print("arg_x=%s" % arg)
            os.system(arg)
            # '''
            # if not util.check_file(30,outfile) :
            #     print 'Using mmtbx.model_vs_data for validation ..'
            #     arg="mmtbx.model_vs_data %s %s %s %s %s f_obs_label='%s' > %s " %(pdbfile,mtz,lib,resol,maps,label,outfile)
            #     print 'arg_x=', arg
            #     os.system(arg)
            # '''
            if not util.check_file(30, outfile):
                print("Initial trial of phenix is not successful. Generating new ligand and trying again")

                arg = "phenix.elbow --do-all %s --output=elbow__LIG >&/dev/null" % pdbfile
                os.system(arg)
                print("Trying phenix again ...")
                arg = "phenix.model_vs_data %s %s %s %s elbow__LIG.cif f_obs_label='%s' > %s " % (pdbfile, mtz, resol, maps, label, outfile)
                os.system(arg)
                if dic["verb"] == 0:
                    os.system("rm -f elbow__LIG.*")

            parse.model_vs_data_log(outfile, dic1, dic)

            mtzfile = sffile.split(".")[0] + "_2mFobs-DFmodel_map_coeffs.mtz"
            if dic["rsrd"] and util.check_file(30, mtzfile):  # get density correlation
                dir_tmp = "PH_CC_TMP.DIR"
                cc_log = "phenix_cc.log"
                if os.path.dirname(dir_tmp):
                    os.system("rm -rf %s " % dir_tmp)
                arg = (
                    'phenix.get_cc_mtz_pdb  pdb_in="%s"  mtz_in="%s" \
                labin="FP=2mFobs-DFmodel PHIB=P2mFobs-DFmodel" resolution=%s \
                temp_dir = "%s" > %s '
                    % (pdbfile, mtzfile, dic["resh"], dir_tmp, cc_log)
                )

                print("Doing density correlations ..")
                os.system(arg)

                ccfile = "%s/cc.log" % dir_tmp  # file contains CC
                if util.check_file(100, ccfile):
                    write_cif.phenix_cc2cif(ccfile, dic)
                    os.system("rm -rf %s " % dir_tmp)

    elif type1 == "neut":
        if vers >= 1013:
            run_phenix_v13_neutron(dic, dic1, pdbfile, mtz, label, outfile)
        else:
            dic1["prog"] = "PHENIX"
            arg = "phenix.model_vs_data %s %s scattering_table=neutron f_obs_label='%s' >%s " % (pdbfile, mtz, label, outfile)
            # print 'arg_n=',arg
            os.system(arg)

            if not util.check_file(30, outfile):
                print("Initial trial of phenix is not successful. Generating new ligand and trying again")

                arg = "phenix.elbow --do-all %s --output=elbow__LIG >&/dev/null" % pdbfile
                os.system(arg)
                print("Tring phenix again ...")
                arg = "phenix.model_vs_data %s %s scattering_table=neutron elbow__LIG.cif  f_obs_label='%s' >%s " % (pdbfile, mtz, label, outfile)
                os.system(arg)
                if dic["verb"] == 0:
                    os.system("rm -f elbow__LIG.*")

            if not util.check_file(300, outfile):  # try the older version
                arg = "phenix.model_vs_data %s %s  --scattering-table=neutron >%s " % (pdbfile, sffile, outfile)
                os.system(arg)

            parse.model_vs_data_log(outfile, dic1, dic)

    if not dic["verb"]:
        util.delete_file(mtz)
    return outfile


##########################################################
def run_phenix_v13_xray(dic, dic1, pdbfile, mtz, label, outfile):
    """Runs xray phenix command. Needs to split out RSRD calculations as phenix.map_vs_model does not produce in newer phenix versions"""

    dic1["prog"] = "PHENIX"

    resol = " "
    if "?" not in dic["resh"] and "n" not in dic["exp"]:
        resol = "high_resolution=%s" % dic["resh"]

    arg = "phenix.model_vs_data %s %s %s f_obs_label='%s' > %s " % (pdbfile, mtz, resol, label, outfile)
    # print 'arg_x=',arg
    os.system(arg)

    parse.model_vs_data_v13_log(outfile, dic1, dic)
    # rsrd not supported


##########################################################
def run_phenix_v13_neutron(dic, dic1, pdbfile, mtz, label, outfile):
    """Runs neutron phenix command. Needs to split out RSRD calculations as phenix.map_vs_model does not produce in newer phenix versions"""

    dic1["prog"] = "PHENIX"

    arg = "phenix.model_vs_data %s %s scattering_table=neutron f_obs_label='%s' >%s " % (pdbfile, mtz, label, outfile)
    # print 'arg_x=',arg
    os.system(arg)

    parse.model_vs_data_v13_log(outfile, dic1, dic)

    # rsrd not supported


##########################################################
def get_matthew_coeff(pdbfile):
    """do solvent"""
    matth_coeff, solvent = -1.0, -1.0
    outf = pdbfile + "_matth_coeff"

    util.delete_file(outf)
    os.system("rwcontents XYZIN %s </dev/null>%s" % (pdbfile, outf))

    if not util.check_file(30, outf):
        return matth_coeff, solvent

    fr = open(outf, "r")
    for line in fr:
        if "The Matthews Coefficient is :" in line:
            matth_coeff = util.float_after_id(line, ":")
        elif " Assuming protein density is 1.34, the" in line:
            solvent = util.float_after_id(line, ":")
    fr.close()
    util.delete_file(outf)

    return matth_coeff, solvent


##########################################################
def grep_items(file, pattern):
    """grep items from the pattern and put them  in a list"""
    item = os.popen('grep "%s" %s ' % (pattern, file)).read().split("\n")

    return item


##########################################################
def phenix_refine(pdb_new, sf_new):
    if "PHENIX" not in os.environ:
        print("\nWarning: To run phenix program, please source 'phenix_env'.\n")
        return " "

    prefix = "PHENIX_REF_TMP"

    arg = "phenix.refine  %s %s stop_for_unknowns=False refinement.output.prefix=%s >/dev/null" % (pdb_new, sf_new, prefix)
    os.system(arg)

    for x in ("_data.mtz", "_001.eff", "_001.geo", "_001.mtz", "_001.log", "_002.def"):
        util.delete_file("%s%s" % (prefix, x))

    return prefix + "_001.pdb"


##########################################################
def refmac_refine(pdb_new, sf_new, dic):
    """do 4 cycle of refinement by refmac"""

    mtz = sf_convertor(sf_new, pdb_new, "mtz")

    pdb_ref = pdb_new + "_refmac0"
    refmac_cif = "NATIVE.refmac"
    log = pdb_new + "refmac0.log"

    util.delete_file(log, pdb_ref, refmac_cif)

    if dic["tls"] > 0:
        mtz3, pdb_ref, log, _scr = run_refmac(pdb_new, mtz, dic, 1, " -rest -ncyc -tls ")
    else:
        mtz3, pdb_ref, log, _scr = run_refmac(pdb_new, mtz, dic, 1, " -rest -ncyc ")

    # arg='auto -refine refmac -pdb %s -sf %s -ncyc 2 > /dev/null ' %(pdb_new, sffile_new)
    # os.system(arg)
    if dic["verb"] == 0:
        util.delete_file(mtz3)

    return refmac_cif, log, pdb_ref


##########################################################
def do_extra_check(cns_pdb):
    """only for cns"""

    if not util.check_file(10, cns_pdb):
        return [], ""
    fp = open(cns_pdb, "r")
    lig = []
    for x in fp:
        res = x[17:20]
        if res not in util.residue() and "9999.0009999.0" in x:
            if res not in lig:
                lig.append(res)

    newpdb = cns_pdb + ".NEWPDB"
    fw = open(newpdb, "w")
    fp.seek(0)
    for x in fp:
        res = x[17:20]
        if res in lig:
            continue
        fw.write(x)
    fw.close()
    fp.close()
    return lig, newpdb


##########################################################
def run_cns(dic, pdb_new, cns_sf):
    """do model_stat to validate the structure. the file names are given"""

    util.delete_file("model_stats.list", "model_stats.log")

    print("Doing cns (%s, %s)..." % (pdb_new, cns_sf))

    pdb_tmp = pdb_new
    if util.is_cif(pdb_new):
        pdb_tmp = cif.cif2pdb(pdb_new)

    # pdb_1 = fix_dna_rna_cns(pdb_tmp)
    (cns_pdb, cns_mtf) = generate_cns_mtf(pdb_tmp)

    lig, pdbnew = do_extra_check(cns_pdb)
    if len(lig) > 0:
        (cns_pdb, cns_mtf) = generate_cns_mtf(pdbnew)

    outfile, logfile = run_model_stat(dic, cns_sf, cns_pdb, cns_mtf)

    if dic["verb"] == 0:
        util.delete_file("luzzati_error.plot luzzati_error_cv.plot sigmaa_error.plot")
        util.delete_file("sigmaa_error_cv.plot generate.pdb generate.mtf ")
        util.delete_file("LIGAND.top LIGAND.param  ", pdbnew)
        if dic["cif"]:
            util.delete_file(pdb_tmp)
        util.delete_file("generate.inp  model_stat.inp generate.log SF_4_validate.cif_cns ")

    return outfile, logfile


##########################################################
def fix_dna_rna_cns(pdb_inp):
    fp = open(pdb_inp, "r")
    pdb_out = pdb_inp + "_cns"
    fw = open(pdb_out, "w")

    for x in fp:
        if "ATOM" in x[:4] or "HETATM" in x[:6] or "ANISOU" in x[:6]:
            res = x[17:20]

            if res == "  C" or res == " DC":
                y = x[:17] + "CYT" + x[20:]
                print("%s %s" % (res, y))
            elif res == "  A" or res == " DA":
                y = x[:17] + "ADT" + x[20:]
            elif res == "  T" or res == " DT":
                y = x[:17] + "THY" + x[20:]

            elif res == "  G" or res == " DG":
                y = x[:17] + "GUA" + x[20:]

            elif res == "  U":
                y = x[:17] + "URA" + x[20:]
            else:
                y = x
            fw.write(y)
        else:
            fw.write(x)
    fw.close()
    fp.close()
    return pdb_out


##########################################################
def symmetry_to_cns(spg):
    """convert space group to CNS type."""

    sym = spg.split()
    if len(sym) > 4 or len(sym) < 2:
        print("Error: symmetry %s is not good for CNS" % spg)
        return spg

    t1 = []
    for i, x in enumerate(sym):
        if i == 0:
            t = x
        elif len(x) == 2:
            t = x[0] + "(" + x[1] + ")"
        else:
            t = x
        t1.append(t)

    symmetry = "".join(t1)

    if symmetry[1] == "1" and symmetry[-1] == "1" and len(symmetry) > 2:
        symmetry = symmetry[0] + symmetry[2:-1]
    print('Note: symmetry "%s" changed to CNS style "%s".' % (spg, symmetry))

    # print('symmetry=', spg, sym, symmetry)

    return symmetry


##########################################################
def run_model_stat(dic, cns_sf, cns_pdb, cns_mtf):
    cell = [float(x) for x in dic["cell"].split()]

    spg = dic["spg"].strip()
    symmetry = symmetry_to_cns(spg)

    inp = cns_inp.model_stat()
    inp_new = inp + "new"

    fp = open(inp, "r")
    fw = open(inp_new, "w")
    for x in fp:
        if "{===>} structure_infile=" in x:
            fw.write('{===>} structure_infile= "%s";\n' % cns_mtf)

        elif "{===>} coordinate_infile=" in x:
            fw.write('{===>} coordinate_infile= "%s";\n' % cns_pdb)

        elif "{===>} parameter_infile_5=" in x:
            if util.check_file(100, "LIGAND.param"):
                fw.write('{===>} parameter_infile_5= "LIGAND.param";\n')
            else:
                fw.write('{===>} parameter_infile_5= "";\n')

        elif "reflection_infile_1=" in x:
            fw.write('{===>} reflection_infile_1= "%s";\n' % cns_sf)

        elif "{===>} sg=" in x:
            fw.write('{===>} sg="%s";\n' % symmetry)

        elif "{===>} a=" in x:
            fw.write("{===>} a=%.2f;\n" % cell[0])

        elif "{===>} b=" in x:
            fw.write("{===>} b=%.2f;\n" % cell[1])

        elif "{===>} c=" in x:
            fw.write("{===>} c=%.2f;\n" % cell[2])

        elif "{===>} alpha=" in x:
            fw.write("{===>} alpha=%.2f;\n" % cell[3])

        elif "{===>} beta=" in x:
            fw.write("{===>} beta=%.2f;\n" % cell[4])

        elif "{===>} gamma=" in x:
            fw.write("{===>} gamma=%.2f;\n" % cell[5])

        elif "{===>} low_res=" in x:
            fw.write("{===>} low_res= %s;\n" % dic["resl"])

        elif "{===>} high_res=" in x:
            fw.write("{===>} high_res= %s;\n" % dic["resh"])

        else:
            fw.write(x)

    fp.close()
    fw.close()
    util.move(inp_new, inp)

    arg = "cns_solve < model_stat.inp  > model_stats.log"
    os.system(arg)

    return "model_stats.list", "model_stats.log"


##########################################################
def generate_cns_mtf(pdb_new):
    geninp = cns_inp.generate()
    inp_new = geninp + "new"

    param, top = cns_param_top(pdb_new)
    cns_mtf, cns_pdb = "generate.mtf", "generate.pdb"

    fp = open(geninp, "r")
    fw = open(inp_new, "w")
    for x in fp:
        if "{===>} coordinate_infile=" in x:
            s = '{===>} coordinate_infile= "%s";\n' % pdb_new
            fw.write(s)
        elif "{===>} structure_outfile=" in x:
            fw.write('{===>} structure_outfile="%s"; \n' % cns_mtf)

        elif "{===>} coordinate_outfile=" in x:
            fw.write('{===>} coordinate_outfile="%s"; \n' % cns_pdb)

        elif "{===>} ligand_topology_infile=" in x:
            if util.check_file(5, top):
                fw.write('{===>} ligand_topology_infile="%s"; \n' % top)
            else:
                fw.write('{===>} ligand_topology_infile=""; \n')
        elif "{===>} ligand_parameter_infile=" in x:
            if util.check_file(5, param):
                fw.write('{===>} ligand_parameter_infile="%s"; \n' % param)
            else:
                fw.write('{===>} ligand_parameter_infile=""; \n')
        else:
            fw.write(x)

    fp.close()
    fw.close()

    util.move(inp_new, geninp)

    arg = "cns_solve < generate.inp > generate.log"
    os.system(arg)

    return cns_pdb, cns_mtf


##########################################################
def cns_param_top(pdb_new):
    param, top = "LIGAND.param", "LIGAND.top"

    lig = unique_ligand(pdb_new)

    if not lig:
        return param, top

    fw1 = open(param, "w")
    fw2 = open(top, "w")
    print("list of unique ligands:", lig.keys())
    for x in lig.keys():
        ligpdb = "%s.pdb" % x.strip()
        fw = open(ligpdb, "w")
        for y in lig[x]:
            fw.write(y)
        fw.close()
        run_xplo2d(ligpdb)
        util.delete_file("%s_xplo2d %s_clean %s_min.inp " % (ligpdb, ligpdb, x.lower().strip()))

        ligname = x.lower().strip()
        parfile, topfile = "%s.par" % ligname, "%s.top" % ligname

        if util.check_file(10, "%s.par" % ligname):
            for y in open(parfile, "r").readlines():
                fw1.write(y)

        if util.check_file(10, "%s.top" % ligname):
            for y in open(topfile, "r").readlines():
                fw2.write(y)

        util.delete_file(parfile, topfile, ligpdb)

    fw1.close()
    fw2.close()
    return param, top


##########################################################
def unique_ligand(pdb_new):
    """generate a dic: unique ligand with coordinates"""

    lig = {}
    if util.check_file(200, pdb_new) == 0:
        print("Error, Can not open pdb file in unique_ligand!")
        return lig

    dd = util.chain_res_atom(pdb_new)  # lig{'chain': {'res':[atom list]}
    for x in dd.keys():
        # for y in sorted(dd[x].iterkeys()) :
        for y in dd[x]:
            res = dd[x][y][0][17:20]
            if res not in util.residue() and len(dd[x][y]) > 1:
                lig[res] = dd[x][y]
                #  print(x, y, res,  dd[x][y][0],len(dd[x][y]) )
    return lig


##########################################################
def run_xplo2d(ligpdb):
    inp = cns_inp.xplo2d()
    inp_new = ligpdb + ".inp"

    fw = open(inp_new, "w")
    fp = open(inp, "r").readlines()
    for x in fp:
        if "INSERT" in x:
            fw.write("%s\n" % ligpdb)
        else:
            fw.write(x)

    fw.close()
    lx_xplo2d = util.program_other_bin("lx_xplo2d")
    arg = "%s < %s >/dev/null  " % (lx_xplo2d, inp_new)
    os.system(arg)
    util.delete_file(inp, inp_new)


##########################################################
def run_shelx(pdb_new, shelx_sf):
    """run shelxH (the large version of shelxH)"""

    print("Running shelx (%s, %s)..." % (pdb_new, shelx_sf))

    _ins_name, base1 = run_shelxpro(pdb_new)

    shelx_log = "shelx_%s.log" % base1
    util.move(shelx_sf, base1 + ".hkl")

    shelxh = util.program_other_bin("shelxh")
    arg = "%s %s >%s" % (shelxh, base1, shelx_log)
    os.system(arg)

    return shelx_log


##########################################################
def run_shelxpro(pdb_new):
    """use pdb file to get the .ins in order to do shelx refinement"""

    base1 = os.path.splitext(os.path.basename(pdb_new))[0] + "_base"
    ins_name = base1 + ".ins"
    scr = base1 + ".sh"

    fw = open(scr, "w")  # do this as the interactive mode!
    fw.write("#!/bin/sh\n")
    shelxpro = util.program_other_bin("shelxpro")
    fw.write("%s  %s <<eof\nI\n\n%s.ins\n%s\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nQ\n\n\nQ\n\n\neof\n" % (shelxpro, base1, base1, pdb_new))
    fw.close()
    arg = "chmod +x %s ; ./%s >/dev/null " % (scr, scr)
    os.system(arg)

    return ins_name, base1


##########################################################
def find_xyzlimt(extend, coord):
    """extend: extend the coordinate. return the xyzlimit in fract"""

    cell, xx, yy, zz = [], [], [], []
    xyzlim = ""

    if not util.check_file(100, coord):
        return xyzlim
    fp = open(coord, "r")
    for x1 in fp:
        if "CRYST1" in x1[:6]:
            cell = [float(p) for p in x1[8:54].split()]
        elif "ATOM" in x1[:4] or "HETATM" in x1[:6]:
            xx.append(float(x1[30:38]))
            yy.append(float(x1[38:46]))
            zz.append(float(x1[46:54]))
        elif "ENDMDL" in x1[:6]:
            break

    fp.close()

    if not xx or not yy or not zz:
        print("Error: residue can not be found in the coordinate.")
        return xyzlim

    frac, _orth = util.frac_orth_matrix(cell)  # get matrix
    xx_min, xx_max = min(xx) - extend, max(xx) + extend
    yy_min, yy_max = min(yy) - extend, max(yy) + extend
    zz_min, zz_max = min(zz) - extend, max(zz) + extend

    xf_min = util.matrix_prod(frac, [xx_min, yy_min, zz_min])
    xf_max = util.matrix_prod(frac, [xx_max, yy_max, zz_max])

    xyzlim = "%.3f %.3f  %.3f %.3f  %.3f %.3f" % (xf_min[0], xf_max[0], xf_min[1], xf_max[1], xf_min[2], xf_max[2])
    return xyzlim


##########################################################
def mtz2map(mtz, pdbfile, mapout, maptype, dic):
    """run fft to get the map: the mtz must be from refmac at the moment!"""

    if maptype == "2FO_FC":
        f1, phi = "FWT", "PHWT"
    elif maptype == "FO_FC":
        f1, phi = "DELFWT", "PHDELWT"
    elif maptype == "FC":
        f1, phi = "FC", "PHIC"
    elif maptype == "FEM":
        f1, phi = "FEM", "PHIFEM"

    size = ""
    if dic["mapsize"] and util.is_number(dic["resh"]):
        r = 2.3 + 0.15 * float(dic["resh"]) ** (3 / 2.0)
        if r > 3:
            r = 3.0
        size = "grid samp %.3f" % r

    scr = """#!/bin/tcsh  -f

##################### mtz2map #####################
echo "Calculating %s map using %s %s ..."

fft  HKLIN  %s  MAPOUT  TMPMAP_BY_MTZ.map  <<end >/dev/null
title mtz2map
#grid samp 4.5
%s
 labin F1=%s  PHI=%s
end

# Extend the map to cover the volume of a model (about 0.5A)
mapmask MAPIN  TMPMAP_BY_MTZ.map  XYZIN  %s   MAPOUT %s <<end >/dev/null
BORDER 0.5
end

rm -f TMPMAP_BY_MTZ.map

""" % (
        maptype,
        f1,
        phi,
        mtz,
        size,
        f1,
        phi,
        pdbfile,
        mapout,
    )

    scr_name = "get_mtz2map.csh"
    fw = open(scr_name, "w")
    fw.write(scr)
    fw.close()
    command = "chmod +x %s ; ./%s  " % (scr_name, scr_name)
    os.system(command)


########################################################################
def eds_radius(resh):
    """The atomic radius used in the mapman. resh: the resolution"""

    radius = 0.7 + (resh - 0.6) / 3.0
    if resh < 0.6:
        radius = 0.7
    elif resh > 3.0:
        radius = resh / 2.0

    return radius


########################################################################
def run_mapman(pdbfile, map1, map2, resh, edsout):
    """Calculating density correlation and real space R factors by mapman.
    map1: the FC map;  map2: the 2Fo_Fc map;
    resh: the resolution to determine the distance around the atoms.
    """

    print("Calculating density correlation and real space R factors by mapman ...")
    radius = eds_radius(resh)

    mapman = util.program_other_bin("lx_mapman")

    script = "run_mapman.csh"
    fw = open(script, "w")

    csh_script = """#!/bin/tcsh  -f

#setenv MAPSIZE 250000000
%s -b  <<EOF >& /dev/null
Read m1 %s ccp4
Read m2 %s ccp4
rsfit m1 m2 %s  %s  %.3f
quit
EOF
    """ % (
        mapman,
        map1,
        map2,
        pdbfile,
        edsout,
        radius,
    )

    fw.write("%s" % csh_script)
    fw.close()

    command = "chmod +x %s ; ./%s " % (script, script)
    os.system(command)

    if util.check_file(500, edsout) == 0:
        modify_eds_script(script)
        command = "chmod +x %s ; ./%s  " % (script, script)
        os.system(command)


##########################################################
def get_contact(chid, nseq, xyzfile, radius):
    """get contact of chid/nseq with other residues in 5A sphere"""

    ncont_log = "NCONT_out.log"
    script = "NCONT_out.csh"
    fw = open(script, "w")

    sh_script = """#!/bin/sh
ncont xyzin %s << eof >%s
source %s/%d
target *
maxdist %.2f
cell 1
eof
    """ % (
        xyzfile,
        ncont_log,
        chid,
        nseq,
        radius,
    )

    fw.write("%s" % sh_script)
    fw.close()

    command = "chmod +x %s ; ./%s " % (script, script)
    os.system(command)

    # util.delete_file(script)

    return ncont_log


##########################################################
def get_mtzfc(pdbfile, mtz_in):
    """get the FC/PHIC from the pdffile and sf file.
    output mtz
    """

    mtzout = mtz_in + "_FC.mtz"
    script = "calc_fc_mtz.csh"
    fw = open(script, "w")

    sh_script = """#!/bin/sh
sfall  HKLIN %s  HKLOUT %s  XYZIN %s << END-sfall >/dev/null
TITL Phasing on initial PA structure (refined twice)
#GRID 52 64 76
MODE SFCALC XYZIN HKLIN
#RESO 30 2.1
#BINS  60
#RSCB 5.0 2.1
SFSG 1
LABI FP=FP SIGFP=SIGFP FREE=FREE
LABO FC=FC PHIC=PHIC
END
END-sfall
    """ % (
        mtz_in,
        mtzout,
        pdbfile,
    )

    fw.write("%s" % sh_script)
    fw.close()

    command = "chmod +x %s ; ./%s " % (script, script)
    os.system(command)

    util.delete_file(script)

    return mtzout


##########################################################
def run_eds_omit(mtzo, pdbfile, dic, dic3):
    """mtzo: optimized phase for 2fo-fc; FC must be calculated with full xyz!"""

    mtz_val = "SF_4_validate.cif_mtz"

    mtz3 = get_mtzfc(pdbfile, mtz_val)
    # mtz3,pdb3,log3,scr=run_refmac(pdbfile, mtz_val, dic, 2, '  ')
    # mtz3 has FC PHIC with the ligand

    map_fc = "CCP4_FC.map"
    mtz2map(mtz3, pdbfile, map_fc, "FC", dic)

    map_2fofc = "CCP4_2FOFC.map"
    mtz2map(mtzo, pdbfile, map_2fofc, "2FO_FC", dic)

    map_fofc = "CCP4_FOFC.map"
    mtz2map(mtzo, pdbfile, map_fofc, "FO_FC", dic)

    edsout = "EDS_OMIT.OUT"
    resh = float(dic["resh"])
    run_mapman(pdbfile, map_fc, map_2fofc, resh, edsout)

    dic3["pdbfile"] = dic["pdbfile"]
    dic3["pdb_new"] = pdbfile
    dic3["pdbid"] = dic["pdbid"]
    dic3["rsrd"] = dic["rsrd"]
    cifrefm1 = write_cif.eds2cif(edsout, dic3)  # get edsout in cif
    if not dic["verb"]:
        # util.delete_file(pdb3,log3,scr )
        util.delete_file("run_mapman.csh", "get_mtz2map.csh", "calc_fc_mtz.csh", mtz3)

    return cifrefm1, edsout


##########################################################
def run_eds_fem(mtzfile, pdbfile, dic, dic3):
    """execute the script to run EDS/MAPS. option: extra options
    dic3 is to be updated after running eds2cif
    """

    if util.check_file(500, mtzfile, pdbfile) == 0:
        util.perror("Error: Refmac failed to generate mtzfile (for EDS).")
        return "", ""

    edsout = pdbfile + "_eds_out"
    util.delete_file(edsout)
    arg = pdbfile + " " + mtzfile + " "

    if dic["xyzlim"]:  # do mapmask by xylim in fract coord
        print("Using xyzlimit for mapmask ...")
        resh = 2.0
        if util.is_number(dic["resh"]):
            resh = float(dic["resh"])
        script = "mtz2eds_xyzlim.csh"
        xyzlim = find_xyzlimt(4.1, pdbfile)  # extend 4.1
        gen_mtz2eds_script_xyzlim(xyzlim, resh, script)
    else:
        if dic["map_dsn6"]:
            arg = arg + " -map_dsn6 "
        if dic["map"]:
            arg = arg + " -map "
        if dic["noeds"]:
            arg = arg + " -noeds "
        script = "mtz2eds.csh"
        # gen_mtz2eds_script(script,'FEM', 'PHIFEM')  # fem
        gen_mtz2eds_script(script, "2mFoDFc", "PHI2mFoDFc")  # ph

    command = "chmod +x %s ; ./%s  %s " % (script, script, arg)
    os.system(command)

    # print(command)
    dic3["pdbfile"] = dic["pdbfile"]
    dic3["pdb_new"] = pdbfile
    dic3["pdbid"] = dic["pdbid"]
    dic3["rsrd"] = dic["rsrd"]
    cifrefm1 = write_cif.eds2cif(edsout, dic3)  # get edsout in cif
    if not dic["verb"]:
        util.delete_file("lx_mapman.log", "mtz2eds_xyzlim.csh")

    return cifrefm1, edsout


##########################################################
def modify_eds_script(script):
    script1 = open(script, "r").readlines()
    fw = open(script, "w")
    for x in script1:
        t = x.lstrip()
        if len(t) > 0 and t[0] == "#" and "MAPSIZE" in t:
            fw.write(t[1:])
        else:
            fw.write(x)
    fw.close()


##########################################################
def run_eds(mtzfile, pdbfile, dic, dic3):
    """execute the script to run EDS/MAPS. option: extra options
    dic3 is to be updated after running eds2cif
    """

    if util.check_file(500, mtzfile, pdbfile) == 0:
        util.perror("Error: Refmac failed to generate mtzfile (for EDS).")
        return "", ""

    edsout = pdbfile + "_eds_out"
    util.delete_file(edsout)
    arg = pdbfile + " " + mtzfile + " "

    if dic["xyzlim"]:  # do mapmask by xylim in fract coord
        print("Using xyzlimit for mapmask ...")
        resh = 2.0
        if util.is_number(dic["resh"]):
            resh = float(dic["resh"])
        script = "mtz2eds_xyzlim.csh"
        xyzlim = find_xyzlimt(4.1, pdbfile)  # extend 4.1
        gen_mtz2eds_script_xyzlim(xyzlim, resh, script)
    else:
        if dic["map_dsn6"]:
            arg = arg + " -map_dsn6 "
        if dic["map"]:
            arg = arg + " -map "
        if dic["noeds"]:
            arg = arg + " -noeds "
        script = "mtz2eds.csh"
        gen_mtz2eds_script(script, "FWT", "PHWT")
    command = "chmod +x %s ; ./%s  %s " % (script, script, arg)
    os.system(command)

    if util.check_file(500, edsout) == 0:
        modify_eds_script(script)
        command = "chmod +x %s ; ./%s  %s " % (script, script, arg)
        os.system(command)

        # print(command)

    dic3["pdbfile"] = dic["pdbfile"]
    dic3["pdb_new"] = pdbfile
    dic3["pdbid"] = dic["pdbid"]
    dic3["rsrd"] = dic["rsrd"]
    cifrefm1 = write_cif.eds2cif(edsout, dic3)  # get edsout in cif
    if not dic["verb"]:
        util.delete_file("lx_mapman.log", "mtz2eds_xyzlim.csh")

    util.get_software_version(edsout, "MAPMAN", dic)
    return cifrefm1, edsout


##########################################################
# def run_fixmap(mtzfile):
#     ''' not used yet.
#     '''

#     map_2fofc = 'CCP4_2FOFC_fixed.map'
#     map_fofc = 'CCP4_FOFC_fixed.map'
#     fix_mtz = '%s_fixed' % mtzfile

#     arg = '''#!/bin/tcsh
# # Fix up the map coefficients: FLABEL specifies the label for Fobs &
# # sigma(Fobs) (defaults are F/SIGF or FOSC/SIGFOSC).  Here, 'in.mtz'
# # is the output reflection file from the refinement program in MTZ
# # format.

# mtzfix  FLABEL FP  HKLIN %s  HKLOUT %s > mtzfix.log
# if($?) exit $?


# # Compute the 2mFo-DFc map; specify the correct labels for
# # the F and phi columns: 'FWT' & 'PHWT' should work for Refmac.
# # Note that EDSTATS needs only 1 asymmetric unit (but will also work
# # with more).  Grid sampling must be at least 4.

# echo 'labi F1=FWT PHI=PHWT\nxyzl asu\ngrid samp 4.5'  | fft  \
# HKLIN %s  MAPOUT %s
# if($?) exit $?

# # Compute the 2(mFo-DFc) map; again need to specify the right labels.

# echo 'labi F1=DELFWT PHI=PHDELWT\nxyzl asu\ngrid samp 4.5'  | fft  \
# HKLIN %s MAPOUT %s
# if($?) exit $?


# # Main- & side-chain atom statistics, using chains A & I only & writing
# # PDB file with per-atom Zdiff metrics.
# echo mole=P,resl=50,resh=1.75,main=atom,side=atom  | edstats  \
# XYZIN pdb3str.ent  MAPIN1  fo.map  MAPIN2 df.map  XYZOUT out.pdb  \
# OUT stats.out
# if($?) exit $?

# ''' % (mtzfile, fix_mtz, fix_mtz, map_2fofc, fix_mtz, map_fofc)


##########################################################
def run_edstat(mtzfile, pdbfile, dic, dic3):
    """Method by Tickle (diff map, Acta Cryst D 2012))"""

    if util.check_file(500, mtzfile, pdbfile) == 0:
        util.perror("Error: Refmac failed to generate mtzfile (for EDSTAT).")
        return

    pdb_main, pdb_side, pdb_phos = split_pdb_edstat(pdbfile)

    dic["cifrefm"] = run_ed(mtzfile, pdbfile, dic, dic3)
    dic["cifrefm1"] = run_ed(mtzfile, pdb_main, dic, dic3)
    dic["cifrefm2"] = run_ed(mtzfile, pdb_side, dic, dic3)
    dic["cifrefm3"] = run_ed(mtzfile, pdb_phos, dic, dic3)


##########################################################
def run_ed(mtzfile, pdbfile, dic, dic3):
    """Method by Tickle (diff map)"""

    if not util.check_file(100, pdbfile):
        return ""

    out1 = pdbfile + "_edstat.out"
    log1 = pdbfile + "_edstat.log"
    dic3["pdb_new"] = pdbfile
    eddir = ""
    if dic["eddir"] != "?":
        eddir = " -tmpdir %s " % dic["eddir"]
    # Take edstat from CCP4 - which should be in path, unless override
    if "EDSTATSBINDIR" in os.environ:
        arg = "perl %s/edstats.pl -exec %s" % (os.environ["EDSTATSBINDIR"], os.environ["EDSTATSBINDIR"])
    else:
        arg = "perl %s/edstats.pl" % os.environ["CBIN"]
    arg = arg + " -hklin %s -xyzin %s  -output %s -log %s %s >/dev/null" % (mtzfile, pdbfile, out1, log1, eddir)
    os.system(arg)

    edstatcif = write_cif.edstat2cif(out1)

    util.get_software_version(log1, "EDSTAT", dic)

    if not dic["verb"]:
        util.delete_file(out1, log1)

    return edstatcif


##########################################################
def run_edsall(mtzfile, pdbfile, dic, dic3):
    """make the xyz be three files (main, sidechain for protein, base,sugar,phos for rna)
    dic3 is to be updated after running eds2cif
    """

    if util.check_file(500, mtzfile, pdbfile) == 0:
        util.perror("Error: Refmac failed to generate mtzfile (for EDS).")
        return

    dic3["pdb_new"] = pdbfile
    pdb1, pdb2, pdb3 = split_pdb(pdbfile)

    map_fc = "CCP4_FC.map"
    mtz2map(mtzfile, pdbfile, map_fc, "FC", dic)

    map_2fofc = "CCP4_2FOFC.map"
    mtz2map(mtzfile, pdbfile, map_2fofc, "2FO_FC", dic)

    edsout, edsout1, edsout2, edsout3 = "EDStmp.OUT", "EDStmp1.OUT", "EDStmp2.OUT", "EDStmp3.OUT"
    resh = float(dic["resh"])
    run_mapman(pdbfile, map_fc, map_2fofc, resh, edsout)
    dic["cifrefm"] = write_cif.eds2cif(edsout, dic3)

    if util.check_file(200, pdb1):
        run_mapman(pdb1, map_fc, map_2fofc, resh, edsout1)
        dic["cifrefm1"] = write_cif.eds2cif(edsout1, dic3)

    if util.check_file(200, pdb2):
        run_mapman(pdb2, map_fc, map_2fofc, resh, edsout2)
        dic["cifrefm2"] = write_cif.eds2cif(edsout2, dic3)

    if util.check_file(200, pdb3):
        run_mapman(pdb3, map_fc, map_2fofc, resh, edsout3)
        dic["cifrefm3"] = write_cif.eds2cif(edsout3, dic3)

    dic3["pdbfile"] = dic["pdbfile"]
    dic3["pdbid"] = dic["pdbid"]
    dic3["rsrd"] = dic["rsrd"]
    if not dic["verb"]:
        util.delete_file("run_mapman.csh", "get_mtz2map.csh")

    util.delete_file(map_fc, map_2fofc, pdb1, pdb2, pdb3, edsout, edsout1, edsout2, edsout3)


##########################################################
def split_pdb_edstat(pdb):
    """make the xyz be three files (base,sugar,phos for rna)"""

    record = util.chain_res_atom(pdb)
    na = 0
    for ch in record.keys():
        nres = record[ch].keys()
        etype = pdb_type_of_chain(nres, record[ch])
        if etype == "na":
            na = 1
            break

    if na == 0:
        return "", "", ""  # no NA involved!

    fp = open(pdb, "r")
    pdb_main = getpdb_na(fp, "sugar")
    fp.seek(0)
    pdb_side = getpdb_na(fp, "side")
    fp.seek(0)
    pdb_phos = getpdb_na(fp, "phos")
    fp.close()

    return pdb_main, pdb_side, pdb_phos


##########################################################
def getpdb_na(fp, id):  # pylint: disable=redefined-builtin
    """get the NA pdb file (base, sugar, phos)"""

    pdbo = "NEW_%s.pdb" % id
    fw = open(pdbo, "w")
    for x in fp:
        # print x
        if not ("ATOM" in x[:4] or "HETA" in x[:4] or "TER " in x[:4] or x[:6] in ["CRYST1", "EXPDTA"] or "SCALE" in x[:5]):
            continue
        if x[:6] in ["CRYST1", "EXPDTA"] or "SCALE" in x[:5]:
            fw.write(x)
            continue

        atom = x[12:16]
        comp = x[17:20]
        if id == "sugar" and comp in util.dna_rna():
            if not (" C1'" in atom or " C2'" in atom or " C3'" in atom or " O4'" in atom or " C4'" in atom or " C5'" in atom or " O2'" in atom):  # sugar
                t = "".join([x[:54], "  0.00", x[60:]])
                fw.write(t)
            else:
                fw.write(x)
        elif id == "phos" and comp in util.dna_rna():
            if not (" P  " in atom or " OP1" in atom or " OP2" in atom or " O3'" in atom or " OP3" in atom or " O5'" in atom):  # phosphate
                t = "".join([x[:54], "  0.00", x[60:]])
                fw.write(t)
            else:
                fw.write(x)

        elif id == " side" and comp in util.dna_rna():
            if (
                " P  " in atom
                or " OP1" in atom
                or " OP2" in atom
                or " O3'" in atom
                or " OP3" in atom
                or " O5'" in atom
                or " C1'" in atom
                or " C2'" in atom
                or " C3'" in atom
                or " O4'" in atom
                or " C4'" in atom
                or " C5'" in atom
                or " O2'" in atom
            ):  # phosphate
                t = "".join([x[:54], "  0.00", x[60:]])
                fw.write(t)
            else:
                fw.write(x)

    fw.close()
    return pdbo


##########################################################
def split_pdb(pdb):
    """make the xyz be three files (main, sidechain for protein, base,sugar,phos for rna)"""

    out1, out2, out3 = "%s_1" % pdb, "%s_2" % pdb, "%s_3" % pdb
    fw1, fw2, fw3 = open(out1, "w"), open(out2, "w"), open(out3, "w")

    fw3.write("REMARK      residue number (n) the atom O3' is put to (n+1)\n")
    record = util.chain_res_atom(pdb)

    for ch in record.keys():
        nres = sorted(list(record[ch].keys()))  # newer version of python
        etype = pdb_type_of_chain(nres, record[ch])
        phos = []
        nr_len = len(nres)
        for i, n in enumerate(nres):
            atoms = record[ch][n]
            comp = atoms[0][17:20]
            if "HOH" in comp or "DOD" in comp:
                continue
            if comp not in util.protein() and comp not in util.dna_rna():
                continue
            for x in atoms:
                if "ATOM" not in x[:4] and "HETA" not in x[:4]:
                    continue
                atom = x[12:16]
                if etype == "prot":
                    if " N  " in atom or " CA " in atom or " C  " in atom or " O  " in atom or " CB " in atom:
                        fw1.write(x)  # main chain
                    else:
                        fw2.write(x)  # sidechain

                elif etype == "na":
                    if " P  " in atom or " OP1" in atom or " OP2" in atom or " OP3" in atom or " O5'" in atom:  # phosphate
                        phos.append(x)
                    elif i < nr_len - 2 and " O3'" in atom:
                        nr = nres[i + 1]
                        bln = record[ch][nr][0]
                        t = x[:23] + bln[23:28] + x[28:]
                        phos.append(t)

                    elif " C1'" in atom or " C2'" in atom or " C3'" in atom or " O4'" in atom or " C4'" in atom or " C5'" in atom or " O2'" in atom:  # sugar
                        fw1.write(x)
                    else:
                        fw2.write(x)  # base

        phos_order = util.sort_column(phos, 22, 27)
        for z in phos_order:
            fw3.write("".join(z))
    fw1.close()
    fw2.close()
    fw3.close()

    return out1, out2, out3


##########################################################
def pdb_type_of_chain(nres, record):
    """identify the type of chain"""
    etype = "other"
    for y in nres:
        for x in record[y]:
            comp = x[17:20]
            if comp in util.protein():
                etype = "prot"
                break
            elif comp in util.dna_rna():
                etype = "na"
                break
            else:
                etype = "other"
        if etype == "prot" or etype == "na":
            break
    return etype


##########################################################
def run_refmac(pdbfile, sffile, dic, idd, option):
    """execute the script to run REFMAC. option: extra options
    update dic,  return mtz, pdb, log, scr
    """

    if dic["cif"] > 0:  # the cif version depend on current path
        if pdbfile[0] != "." and pdbfile[0] != "/":
            pdbfile = "./" + pdbfile
        if sffile[0] != "." and sffile[0] != "/":
            sffile = "./" + sffile

    arg = pdbfile + " " + sffile + " -mtz " + option

    if dic["refmac_tls"]:
        arg = arg + " -tls "
    if "-scale" not in option and dic["scale"] > 0:
        arg = arg + " -scale "
    if "-rest" not in option and dic["rest"] > 0:
        arg = arg + " -rest "
    # if '-reso'  not in option and dic['reso']>0 :  arg = arg + ' -reso '
    arg = arg + " -reso "
    if "-wave" not in option and dic["wave"] > 0:
        arg = arg + " -wave "
    if "-nofree" not in option and dic["status"] == "N":
        arg = arg + " -nofree "
    if "-twin" not in option and (dic["refmactwin"] > 0 or dic["twin"] == "Y"):
        arg = arg + " -twin "
        print("Note: validation is set in twin mode!\n")

    pdbout = pdbfile + "_refmac0"
    logout = pdbfile + "_refmac0_log"
    mtzout = pdbfile + "_refmac0.mtz"
    util.delete_file(pdbout, logout, mtzout)

    script = "run_refmac.csh"
    if dic["cif"] > 0:  # only test the cif version
        if "REFMAC_CIF_BINARY" not in os.environ:
            print("Error: please point REFMAC_CIF_BINARY to refmac in full path")
            sys.exit()
        else:
            refmac = os.environ["REFMAC_CIF_BINARY"]
            gen_refmac_script(script, refmac)
            print("Testing the cif version of %s." % refmac)
    else:
        gen_refmac_script(script, "refmac5")

    command = "chmod +x %s ; ./%s  %s " % (script, script, arg)
    # print 'comand=', (command)

    os.system(command)

    mtzo, pdbo, logo = "%s_%d.mtz" % (pdbfile, idd), "%s_%d.pdb" % (pdbfile, idd), "%s_%d.log" % (pdbfile, idd)

    util.move(mtzout, mtzo)
    util.move(pdbout, pdbo)
    util.move(logout, logo)

    util.get_software_version(pdbo, "PDBCAL", dic)
    return mtzo, pdbo, logo, script


##########################################################
def run_refmac_more(pdb, mtz, dic, dici):
    """do extra validations, change resolution, scale, twin, wave"""

    rcut = 0.015  # 1.0% difference tolerance
    if dic["omitmap"]:
        rcut = 0.03  # 3.0% difference tolerance
    if util.is_number(dic["resh"]) and float(dic["resh"]) <= 1.5:
        rcut = 0.010

    dicout = main.initialize_dic1()
    if not (util.is_number(dic["rfact"]) and util.is_number(dici["rfact"])) or float(dici["rfact"]) - float(dic["rfact"]) < rcut:
        return "", "", "", "", dicout

    dics, mtzs, pdbs, logs = [], [], [], []
    option = [" -scale -reso ", " -rest -wave -reso ", " -twin -reso ", " -scale -twin -reso "]
    # option=[' -scale  ', ' -rest -wave  ', ' -twin  ', ' -scale -twin  ']
    detail = ["Use LSSC EXPE", "Use WAVE+REST", "Use TWIN ", "Use EXPE+TWIN "]
    # option=[' -twin ',  ' -scale ', ' -scale -twin ', ' -rest -wave ', ' -reso ']

    ntry = len(option)
    for i, x in enumerate(option):
        dicx = main.initialize_dic1()
        if "-rest " in detail[i]:
            if (util.is_number(dic["resh"]) and float(dic["resh"]) > 2.5) or dic["model"] == "Y":
                continue

        dicx["detail"] = detail[i]

        id = ntry + 1 + i  # pylint: disable=redefined-builtin
        mtzo, pdbo, logo, scr = run_refmac(pdb, mtz, dic, id, x)
        parse.refmac_log(logo, dicx)

        dics.append(dicx)
        mtzs.append(mtzo)
        pdbs.append(pdbo)
        logs.append(logo)

        print("Running refmac with option (%s).." % detail[i])
        print("Rfact: reported=%s, calculated=%s, fom=%s " % (dic["rfact"], dicx["rfact"], dicx["fom"]))
        if (
            util.is_number(dicx["fom"])
            and float(dicx["fom"]) > 0.6
            and util.is_number(dicx["rfact"])
            and util.is_number(dic["rfact"])
            and float(dicx["rfact"]) - float(dic["rfact"]) < rcut
        ):
            if "-twin " in option[i] and float(dici["rfact"]) - float(dic["rfact"]) > 0.05 and "2:" in dic["twin_fr"]:
                util.perror("Warning: Crystal may be twinned. Please double check!")
            break

    nn = main.best_of_solutions(dics)
    dicout = dics[nn]

    for i in range(len(mtzs)):
        if i == nn:
            continue
        util.delete_file(mtzs[i], pdbs[i], logs[i])

    return mtzs[nn], pdbs[nn], logs[nn], scr, dicout


########################################################################
def gen_refmac_script(script_name, refmac):
    """generate a script to run refmac only"""

    csh_script = """#!/bin/csh  -f


###########################################################################
# NOTE: this script is generated by DCC. It can be run independently
#
#  This script is to use multi-programs to calculate statistics (such as
#  R_work, Rfree, real space R factors, by refmac)
#
#  Below are the programs used
#  1. sf_convert: to convert SF file from mmCIF to mtz format
#  2. refmac: to calculated statistics
#  3. cad:  sort all reflection in asymmetric unit.
#  4. cif2mtz: optional. if sf_convert is not installed.
#
############################### Usage ######################################
#
#  usage: "refmac.csh  pdbfile  sffile (add options)"
#
#############################################################################

if( $#argv < 1 ) then
    echo "Please input a PDB file and a SF file in mmcif format. "
    echo "Usage: ./refmac.csh  pdbfile  sffile "
    echo "If add (-twin), do validation in twin mode."
    echo "If add (-scale), use 'SCALe type BULK  LSSC  ANISO  EXPE' "
    echo "If add (-reso), use the high resolution from the PDB file. "
    echo "If add (-rest), do zero cycle with restraint(default, unrestraint)."
    echo "If add (-wave), use wavelength to get correct scattering factor."
    echo "If add (-mtz), input sf is mtz, else, input other formats. "
    echo "If add (-nofree), do not use the free set. (for legacy files). "
    echo "If add (-ncyc), do 4 cycle refinement."
    echo "If add (-tls), include tls (obtained by tlsextract)."
    echo "If add (-addu), export full B factor "
    echo "If add (-isot), Indiv. isot. B-factor refine, else,default by xyz"
    echo "If add (-all0), all 0 cycle(NCYC, TLSC). update full B by TLS group"
    exit()
endif

############ identify pdb and sffile ###########
if ( $#argv >= 2 ) then
    set tmp = `egrep "^CRYST1|_atom_site.occupancy" -m 1 $1 |wc -c`
    if ($tmp >= 6) then
        set pdbfile = $1
        set sffile  = $2
    else
        set pdbfile = $2
        set sffile  = $1
    endif
else
    echo "Too few arguments. input files:  pdbfile  sffile."
    exit()
endif

############# set options ############
set twin=""
set tls=0
set tls_addu=""
set scale=""
set reso=0
set rest="UNrestraint"
set wave=0
set ncyc=0
set mtz=0
set nofree=0
set isot=0
set all0=0  #only three options (-tls -all0 -addu)
foreach arg  ($argv)
    if ( "$arg" == "-twin" ) then
      #  set twin="twin  FilterLevel 0.10 "
        set twin="twin   "
    else if ( "$arg" == "-scale" ) then
        set scale="SCALe  type BULK    LSSC   ANISO   EXPE "
    else if ( "$arg" == "-reso" ) then
        set reso=1
    else if ( "$arg" == "-rest" ) then
        set rest="REstraint "
    else if ( "$arg" == "-wave" ) then
        set wave=1
    else if ( "$arg" == "-mtz" ) then
        set mtz=1
    else if ( "$arg" == "-nofree" ) then
        set nofree=1
    else if ( "$arg" == "-ncyc" ) then
        set ncyc=4
    else if ( "$arg" == "-tls" ) then
        set tls=1
    else if ( "$arg" == "-addu" ) then
        set tls_addu="tlsout addu"  #export Full Bfactors
    else if ( "$arg" == "-isot" ) then
        set isot=1
    else if ( "$arg" == "-all0" ) then
        set all0=1
    endif
end


############ convert cif to mtz format if not mtz  ###########
if ($mtz>0) then
    set sfmtz = $sffile
else
    echo "Converting SF file to MTZ format ..."
    set sfmtz = "${sffile}.mtz"
    %ssf_convert -o mtz -pdb $pdbfile -sf $sffile -out $sfmtz >/dev/null
    set mtzsize=`wc -c $sfmtz | awk '{print $1}'`
    if( ! -e $sfmtz || $mtzsize < 1000 ) then
        set cel = `grep "^CRYST1 " -m 1 $pdbfile | cut -c 8-54`
        set sym = `grep "^CRYST1 " -m 1 $pdbfile | cut -c 55-66`
        set tmp_mtz = "${sffile}__tmp"

cif2mtz HKLIN $sffile HKLOUT $tmp_mtz << end > /dev/null
TITLE testing
symmetry '$sym'
cell $cel
PNAME mtz2cif
DNAME data_1
XNAME cryst_1
end

cad  HKLIN1 $tmp_mtz HKLOUT  $sfmtz <<end > /dev/null
TITLE Sorting HKL with cad
monitor BRIEF
labin file 1 ALL
sort H K L
end

        if( -e  $tmp_mtz ) rm -f $tmp_mtz

        set mtzsize1 = `wc -c $sfmtz | awk '{print $1}'`
        if( ! -e $sfmtz || $mtzsize1 < 2000 ) then  # mtz2cif failed.
            echo "Error: SF ($sffile) can not be converted to mtz format.";
            exit()
        endif

    endif
endif

############################ Run REFMAC5 Below ##########################
# This script is to do unrestrained refinement for ZERO cycles.
# The purpose is to generate Fc and weighted 2Fo-Fc (FWT).
# A new MTZ file contains FWT and phase information.
# use 'scale type bulk' for Babinet's bulk solvent correction.
# use 'solvent yes/no' or 'default' for mask based bulk solvent.
#########################################################################

set refmac_mtz_out="${pdbfile}_refmac0.mtz"

set anis="bref MIXEd"
if($isot>0) set anis="bref isot"

set free_inp = "FREE=FREE"
if($nofree>0) set free_inp=""

set reso_inp = "#"
if($reso>0) then
    set tmp1 = `grep "^REMARK.*RESOLUTION.*LOW.*ANGSTROMS" -m 1 $pdbfile |cut -c49-55 `
    set tmp2 = `grep "^REMARK.*RESOLUTION.*HIGH.*ANGSTROMS" -m 1 $pdbfile |cut -c49-55 `
#    if ($tmp1 != "NULL" && $tmp2 != "NULL") set reso_inp = "reso $tmp1  $tmp2 "
    if ($tmp2 != "NULL") set reso_inp = "reso   $tmp2 "
endif

set wave_inp = "#"
if ($wave > 0) then
    set wavelength = "0.9795"
    set tmp2=`grep "^REMARK.*WAVELENGTH OR RANGE.*(A).*:" -m 1 $pdbfile |cut -c45-53 |awk '{print $1}' `
    set test=`echo $tmp2 | egrep "^[0-9]" | wc -c `
    if ( $test >0 ) set wavelength = "$tmp2"
    set wave_inp = "anomalous wavelength $wavelength "
endif


set tls_info = ""
set ntls = ""
if($tls>0) then
    set tls_inp = "${pdbfile}_tls.inp"
    tlsextract XYZIN $pdbfile TLSOUT $tls_inp >/dev/null
    set size=`wc -l $tls_inp  | awk '{print $1}'`
    if( -e $tls_inp ||  $size > 6 ) set tls_info = "TLSIN $tls_inp "
    if($all0>0) then
        set ncyc = 0
        set ntls = "tlsc 0"
    endif

endif

set hydr = "MAKE HYDRogens yes"
if ("$rest" == "REstraint ") then
   set hydr = "MAKE HYDRogens all"
endif

echo " "
#echo "Doing REFMAC refinement($ncyc cycle, $rest) ($pdbfile, $sfmtz)..."
echo "Doing REFMAC refinement($ncyc cycle, $rest)..."
%s  XYZIN  $pdbfile XYZOUT ${pdbfile}_refmac0 HKLIN  $sfmtz $tls_info HKLOUT  $refmac_mtz_out  << end > ${pdbfile}_refmac0_log

MAKE CHECK NONE
MAKE NEWLigand Noexit
$hydr
labin  FP=FP SIGFP=SIGFP  $free_inp
#labin  IP=I SIGIP=SIGI  $free_inp
#labin  F+=F(+)  SIGF+=SIGF(+)  F-=F(-)  SIGF-=SIGF(-)   $free_inp

labout  FC=FC FWT=FWT PHIC=PHIC PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM
#REFI tlsc 10 TYPE $rest  RESI MLKF  $anis
REFI $ntls TYPE $rest  RESI MLKF  $anis
BLIM 0.001 1000
$scale
$twin
$reso_inp
$wave_inp
$tls_addu
NCYC $ncyc

#Anomalous atom SE f' f"
#anomalous wavelength 0.9683
#WEIGHT MATRIX 0.35
#SCALe LSSC FIXBulk BBULk 200
#SCALe TYPE SIMP
#SPHEricity 2.0
#BFACtor SET 20.0
#mapcalculate shar # regularised map sharpening, v5.6
#ncsr local       # automatic and local ncs  v5.6
#ridg dist sigm 0.05 # jelly body restraints  v5.6
#SOLVENT VDWProb 1.4 IONProb 0.8 RSHRink 0.8

PNAME TEMP
DNAME NATIVE
RSIZE 80
USECWD
end

# check for success
set mtzsize=`wc -c $refmac_mtz_out | awk '{print $1}'`
if( ! -e $refmac_mtz_out || $mtzsize < 1000 ) then
    echo "Error: REFMAC failed! (pdb=$pdbfile, sf=$sffile)!";
    exit()
endif


#---------------- End of script ----------------
""" % (
        config.PATH1,
        refmac,
    )

    fw = open(script_name, "w")
    fw.write("%s" % csh_script)
    fw.close()


########################################################################
def gen_mtz2eds_script_xyzlim(xyzlim, resh, script):
    mapman = util.program_other_bin("lx_mapman")

    csh_script = """#!/bin/csh  -f

#--------------created 2014-01-16 ---------------------------------
#  A new script to handle huge structures by mapmask using xyzlimit
#  mtz2eds-xyzlim.csh
#------------------------------------------------------------------

#------------------- set output files -------------------
set pdbfile = $1
set refmac_mtz_out = $2

set map_fc = "${pdbfile}_fc.map"
set map_fofc = "${pdbfile}_fofc.map"
set map_2fofc = "${pdbfile}_2fofc.map"
set eds_outtmp="${pdbfile}_eds_out"

#------------------- generate FC map -------------------
set tmp_map1 = "${pdbfile}_tmp_fc.map"

echo "Calculating Fc map using FC/PHIC in the mtz file ..."
fft  HKLIN  $refmac_mtz_out  MAPOUT  $tmp_map1 <<end >& /dev/null
title FC
#grid samp 4.5
 labin F1=FC PHI=PHIC
end

# Extend the FC map to cover the volume of a model
mapmask MAPIN  $tmp_map1 MAPOUT $map_fc  <<end >& /dev/null
xyzlimit  %s
end

#------------------- generate 2Fo-FC map -------------------
set tmp_map2 = "${pdbfile}_tmp_fo.map"

echo "Calculating 2mFo-DFc map using FWT/PHWT in the mtz file ..."
fft  HKLIN  $refmac_mtz_out  MAPOUT $tmp_map2 <<end >& /dev/null
title 2FO-1FC
#grid samp 4.5
 labin F1=FWT  PHI=PHWT
end

# Extend the 2Fo-Fc map to cover the volume of a model
mapmask MAPIN  $tmp_map2 MAPOUT $map_2fofc  <<end >& /dev/null
xyzlimit  %s
end

#------------------- generate Fo-FC map -------------------
set tmp_map3 = "${pdbfile}_tmp_fofc.map"

echo "Calculating Fo-Fc map using FWT/PHWT in the mtz file ..."
fft  HKLIN  $refmac_mtz_out  MAPOUT $tmp_map3 <<end >& /dev/null
title FO-FC
#grid samp 4.5
 labin F1=DELFWT  PHI=PHDELWT
end

# Extend the 2Fo-Fc map to cover the volume of a model
mapmask MAPIN  $tmp_map3 MAPOUT $map_fofc  <<end >& /dev/null
xyzlimit  %s
end


#------------------- check generated maps -------------------


set size1=`wc -c $map_fc | awk '{print $1}'`
set size2=`wc -c $map_2fofc | awk '{print $1}'`
if( ! -e $map_fc || ! -e $map_2fofc || $size1 < 1000 || $size2 < 1000 ) then
    echo "Error: maps were not generated. (No EDS calculation performed!)";
    exit()
endif

#------------------- Doing EDS calculation below-------------------
echo "Calculating density correlation and real space R factors by mapman ..."
set resh = "%.3f"
set radius = ` echo $resh |  awk '{ if( $resh < 0.6 ) { print 0.7} else if \
($resh>3.0) {print $resh/2.0} else { print 0.7 + ($resh - 0.6)/3.0}}'`

#echo " resolution=$resh ;  radius for EDS=$radius "

# start to run
setenv MAPSIZE 80000000
%s -b  <<EOF >& lx_mapman.log
Read m1 $map_fc ccp4
Read m2 $map_2fofc ccp4
rsfit m1 m2 $pdbfile  $eds_outtmp  $radius
quit
EOF

rm -f $tmp_map1 $tmp_map2 $tmp_map3 $map_fc


set size1=`wc -c  $eds_outtmp | awk '{print $1}'`
if( ! -e $eds_outtmp || $size1 < 1000 ) then
    echo "Error: EDS calculation failed!";
#    exit()
endif
#-------end of script -----

""" % (
        xyzlim,
        xyzlim,
        xyzlim,
        resh,
        mapman,
    )

    fw = open(script, "w")
    fw.write("%s" % csh_script)
    fw.close()


########################################################################
def gen_mtz2eds_script(script_name, ampl, phase):
    """generate a script to run maps and do EDS only"""

    mapman = util.program_other_bin("lx_mapman")
    # print 'mapman=', mapman

    csh_script = """#!/bin/csh  -f


#########################created 2009-01-25################################
#  This script is to use multi-programs to calculate maps,
#  density correlation... using the mtz file produced by refmac)
#
#  Below are the programs used
#  1. fft: to generate Fc and 2Fo-Fc maps using the new MTZ file
#  2. mapmask: to cover the volume of a model (radius = 4.1 angstrom)
#  3. lx_mapman: calculate real space R factors & density correlation
#     It also converts map format from CCP4 to BRIX/DSN6 (O type maps)
#
############################### Usage ######################################
#
#  density_eds.csh   pdbfile  mtzfile
#
#  The following options (as below) can be added:
# -map :  output maps  (2Fo-Fc & Fo-Fc) in CCP4 formats.
# -map_dsn6 : output maps (2Fo-Fc & Fo-Fc) in BRIX/DSN6 formats.
# -noeds : do not calculate EDS
#
#############################################################################

##########################################################
# things below are for map and DCC calculations
##########################################################

if( $#argv < 1 ) then
    echo "Input a PDB file and a mtz file with phase generated by REFMAC"
    echo "Usage: ./mtz2eds.csh  pdbfile  mtzfile "
    echo "If add (-map), output maps (2Fo-Fc & Fo-Fc) in CCP4 formats"
    echo "If add (-map_dsn6), output maps(2Fo-Fc & Fo-Fc)in BRIX/DSN6 formats"
    echo "If add (-noeds), do not calculate EDS. Only maps"
    exit()
endif

set map = 0
set map_dsn6 = 0
set noeds = 0
foreach arg  ($argv)
    if ( "$arg" == "-map" ) then
        set map = 1
    else if ( "$arg" == "-map_dsn6" ) then
        set map_dsn6 = 1
    else if ( "$arg" == "-noeds" ) then
        set noeds = 1
    endif
end

############ identify pdb and sffile ###########
if ( $#argv >= 2 ) then
    set tmp = `egrep "^CRYST1|_atom_site.occupancy" -m 1 $1 |wc -c`
    if ($tmp >= 6) then
        set pdbfile = $1
        set refmac_mtz_out = $2
    else
        set pdbfile = $2
        set refmac_mtz_out  = $1
    endif
else
    echo "Too few arguments. input files:  pdbfile  mtzfile."
    exit()
endif


set iscif = `grep "_atom_site.occupancy"  -m 1 $pdbfile |wc -c`
set pdbnew = "${pdbfile}.PDB"
set resh = "1.5"
if($iscif >0 ) then
    echo "Note: Input coordinate file ($pdbfile) is in mmcif format."
    set resh = `grep "_refine.ls_d_res_high" -m 1 $pdbfile |awk '{print $2}' `
    %sdcc -cif2pdb $pdbfile >/dev/null
else
    set resh = `grep "^REMARK.*RESOLUTION.*HIGH.*ANGSTROMS" -m 1 $pdbfile |cut -c49-55 `
    cp -f $pdbfile $pdbnew
endif


##################### set output files #################
set map_fc = "${pdbfile}_fc.map"
set map_2fofc = "${pdbfile}_2fofc.map"
set eds_outtmp="${pdbfile}_eds_out"

set map_fofc = "${pdbfile}_fofc.map"
set map_fofc_dsn6 = "${pdbfile}_fofc_dsn6.map"
set map_2fofc_dsn6 = "${pdbfile}_2fofc_dsn6.map"


#####################generate 2Fo-Fc map #####################
set tmp_map2 = "${pdbfile}_maptmp2"

echo "Calculating 2mFo-DFc map using FWT/PHWT by REFMAC ..."
fft  HKLIN  $refmac_mtz_out   MAPOUT $tmp_map2 <<end >/dev/null
title 2FO-1FC
#grid samp 4.5
 labin F1=%s  PHI=%s
end

# Extend the 2Fo-Fc map to cover the volume of a model
mapmask MAPIN  $tmp_map2 MAPOUT $map_2fofc XYZIN  $pdbnew <<end >/dev/null
BORDER 4.1
end

rm -f $tmp_map2

#####################generate Fo-Fc map #####################
#echo "Calculating Fo-Fc map using DELFWT/PHDELWT by REFMAC ..."

set tmp_map = "${pdbfile}_maptmp"

#echo "Calculating  mFo-DFc map using DELFWT/PHDELWT by REFMAC ..."
fft  HKLIN  $refmac_mtz_out   MAPOUT  $tmp_map  <<end >/dev/null
title FO-FC
#grid samp 4.5
 labin F1=DELFWT  PHI=PHDELWT
end

# Extend the Fo-Fc map to cover the volume of a model
mapmask MAPIN  $tmp_map MAPOUT $map_fofc XYZIN  $pdbnew <<end >/dev/null
BORDER 4.1
end

rm -f $tmp_map


#--------------------------------------------------------

if( $noeds ) then
    echo "Note: No EDS calculation. Only maps are generated.";
    rm -f  $pdbnew
    exit()
endif

##############################################################
# Doing EDS calculation below
##############################################################

#####################generate FC map #####################
set tmp_map1 = "${pdbfile}_maptmp1"

echo "Calculating Fc map using FC/PHIC by REFMAC ..."
fft  HKLIN  $refmac_mtz_out  MAPOUT  $tmp_map1 <<end >/dev/null
title FC
#grid samp 4.5
 labin F1=FC PHI=PHIC
end

# Extend the FC map to cover the volume of a model (about 4.1A)
mapmask MAPIN  $tmp_map1 MAPOUT $map_fc XYZIN  $pdbnew <<end >/dev/null
BORDER 4.1
end

rm -f $tmp_map1

###################### check generated maps #####################

set size1=`wc -c $map_fc | awk '{print $1}'`
set size2=`wc -c $map_2fofc | awk '{print $1}'`
if( ! -e $map_fc || ! -e $map_2fofc || $size1 < 1000 || $size2 < 1000 ) then
    echo "Error: No maps generated for ($pdbfile). (stop).";
    exit()
endif

####################### EDS calculation ##########################

echo "Calculating density correlation and real space R factors by mapman ..."


set radius = ` echo $resh | awk '{ if( $resh < 0.6 ) { print 0.7} else if ($resh>3.0) {print $resh/2.0} else { print 0.7 + ($resh - 0.6)/3.0}}'`

#echo " resolution=$resh ;  radius=$radius "

# start to run
#setenv MAPSIZE 250000000
%s -b  <<EOF >& /dev/null
Read m1 $map_fc ccp4
Read m2 $map_2fofc ccp4
rsfit m1 m2 $pdbnew  $eds_outtmp  $radius
quit
EOF

rm -f $map_fc    #no need anymore
set size1=`wc -c  $eds_outtmp | awk '{print $1}'`
if( ! -e $eds_outtmp || $size1 < 1000 ) then
    echo "Error: EDS calculation was not successful!";
    exit()
endif

######################################################################
# Below is for map conversion
######################################################################
if ( $map == 0 && $map_dsn6 == 0 ) then
    rm -f  $map_2fofc $map_fofc $pdbnew
    exit()
endif

################################################
if ( $map == 1 ) then
    rm -f  $pdbnew
    exit()
endif


#############Convert CCP4 map to O style (BRIX/DSN6) ############
#setenv MAPSIZE 250000000
setenv NUMMAP 1
%s -b << EOF >& /dev/null
Read m1 $map_fofc CCP4
Write m1 $map_fofc_dsn6  DSN6
quit
EOF

%s -b << EOF >& /dev/null
Read m1 $map_2fofc CCP4
Write m1 $map_2fofc_dsn6  DSN6
quit
EOF

set size2=`wc -c $map_2fofc_dsn6 | awk '{print $1}'`
if( ! -e $map_2fofc_dsn6 ||  $size2 < 2000 ) then
    echo "Error: map conversion to DSN6 was not successful! (check size)";
    exit()
endif

rm -f   $map_fofc  $map_2fofc  $pdbnew
gzip -f  $map_fofc_dsn6  $map_2fofc_dsn6

############## End of script ############
""" % (
        config.PATH1,
        ampl,
        phase,
        mapman,
        mapman,
        mapman,
    )

    fw = open(script_name, "w")
    fw.write("%s" % csh_script)
    fw.close()


########################################################################
def run_pointless(mtz, dic):
    out_pointless = run_ccp4(dic["pdbfile"], mtz, "pointless")
    parse.pointless_log(out_pointless, dic)
    util.get_software_version(out_pointless, "POINTLESS", dic)
    if not dic["verb"]:
        util.delete_file(out_pointless)


########################################################################
