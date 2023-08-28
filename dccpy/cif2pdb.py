##################################################
# Convert cif2pdb - separated out from cifparse.py
##################################################

import sys
from cifparse import parse_values, is_a_number, get_cell, cifparse, get_symm, get_uij, get_prog, get_xyz


##########################################################
def cif2pdb(ciffile):
    """convert the cif to simple pdb file. (good enough for validation)"""

    fname_out, _mapfile = cif2pdb_ext(ciffile)
    return fname_out


def cif2pdb_ext(ciffile):
    """convert the cif to simple pdb file. (good enough for validation).  For CCD > 5 characters return mapping file"""

    pdbf = ciffile + ".PDB"
    fw = open(pdbf, "w")

    flist = open(ciffile, "r").readlines()
    method = _add_header_pdb(fw, flist)
    _add_remark3(fw, flist, method)
    _add_twin_pdb(fw, flist)
    _add_tls_pdb(fw, flist)
    _add_remark200(fw, flist)
    _add_remark285(fw, flist)
    _add_ncs(fw, flist)

    c = get_cell(flist)
    s = get_symm(flist)

    ss = "CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f %-10s \n" % (c[0], c[1], c[2], c[3], c[4], c[5], s)
    fw.write(ss)
    _add_scale_pdb(fw, flist)

    # for anisou records
    _unatm, uatom, uasym, useq, ucomp, _uins, _ualt, u11, u22, u33, u12, u13, u23 = get_uij(flist)

    # for xyz
    group, natm, atom, asym, seq, comp, ins, alt, x, y, z, occ, biso, model, symbol = get_xyz(flist)
    if not (atom and comp and asym and seq and x and y and z and occ and biso and symbol):
        print("Error: Cif (%s) to pdb conversion failed. Missing the key cif tokens on atom_site." % ciffile)
        sys.exit()

    # Create map for CCD if needed

    # Get list of CCDs used in this file
    hasextended = False

    usedccd = set()
    # In case incoming file has a special CCD - so we cannot use
    for idx in range(len(comp)):
        ccd = comp[idx]
        if ccd not in usedccd:
            usedccd.add(ccd)
            if len(ccd) > 3:
                hasextended = True

    ccdmap = {}
    if hasextended:
        nextid = 0
        # Map
        for idx in range(len(comp)):
            ccd = comp[idx]
            if len(ccd) > 3:
                if ccd not in ccdmap:
                    # Create lookup
                    haveext = False
                    while (not haveext):
                        extccdname = "%02d" % nextid
                        if extccdname in usedccd:
                            # Used already from incoming model
                            nextid += 1
                            continue
                        haveext = True
                        ccdmap[ccd] = extccdname
                    if not haveext:
                        print("Error: Cif (%s) to pdb conversion failed. Could not map extended CCD." % ciffile)
                        sys.exit()
                comp[idx] = ccdmap[ccd]

    ######

    nmodel = 1
    if model and is_a_number(model[-1]) and int(model[-1]) > 1:
        nmodel = int(model[-1])

    asym2ch = _asym2chain(asym)

    nline = len(x)
    if not group:
        group = ["ATOM" for i in x]
    if not alt:
        alt = ["." for i in x]
    if not ins:
        ins = ["." for i in x]
    if not natm:
        natm = ["%d" % i for i in range(nline)]

    n = 0
    for i in range(nline):
        nat = int(natm[i])
        nres = int(seq[i])
        if nat > 99999:
            nat = 99999
        if nres > 9999:
            nres = 9999

        if len(symbol[i].strip()) > 1:
            tmp1 = "%s    " % atom[i].strip()
        else:
            tmp1 = " %s    " % atom[i].strip()
        if len(atom[i].strip()) == 4:
            tmp1 = "%s    " % atom[i]

        atomname = tmp1[:4]
        if nmodel > 1:
            if int(model[i]) >= 180:
                break  # if too many atom, refmac breaks.
            if i == 0:
                fw.write("MODEL %8s \n" % model[i])
            elif 0 < i < nline - 1 and model[i - 1] != model[i]:
                fw.write("ENDMDL  \n")
                fw.write("MODEL %8s \n" % model[i])

        inst = ins[i]
        if ins[i] == "." or ins[i] == "?":
            inst = " "

        alter = alt[i]
        if alt[i] == "." or alt[i] == "?":
            alter = " "
        v = [float(xx) for xx in (x[i], y[i], z[i], occ[i], biso[i])]
        asym2 = asym[i]
        if asym2ch:
            asym2 = asym2ch[asym[i]]
        space = " "
        if v[4] < 1000:
            ss = "%-6s%5d %4s%1s%3s%2s%4d%1s   %8.3f%8.3f%8.3f%6.3f%6.2f%6s%-4s%2s  \n" % (
                group[i],
                nat,
                atomname,
                alter,
                comp[i],
                asym2,
                nres,
                inst,
                v[0],
                v[1],
                v[2],
                v[3],
                v[4],
                space,
                asym[i],
                symbol[i],
            )
        else:
            ss = "%-6s%5d %4s%1s%3s%2s%4d%1s   %8.3f%8.3f%8.3f%6.3f%6.1f%6s%-4s%2s  \n" % (
                group[i],
                nat,
                atomname,
                alter,
                comp[i],
                asym2,
                nres,
                inst,
                v[0],
                v[1],
                v[2],
                v[3],
                v[4],
                space,
                asym[i],
                symbol[i],
            )

        fw.write(ss)

        if n == len(uatom):
            continue
        if uatom and ucomp and uasym and useq:
            # print(n, ucomp[n], comp[i], uatom[n], atom[i],uasym[n], asym[i])
            if n < len(ucomp) and ucomp[n] == comp[i] and uatom[n] == atom[i] and uasym[n] == asym[i] and useq[n] == seq[i]:
                u = [10000 * float(xx) for xx in (u11[n], u22[n], u33[n], u12[n], u13[n], u23[n])]
                # uu='%7.0f%7.0f%7.0f%7.0f%7.0f%7.0f' %(u[0], u[1],u[2],u[3],u[4],u[5])
                ss = "ANISOU%5d %4s%1s%3s%2s%4d%1s %7.0f%7.0f%7.0f%7.0f%7.0f%7.0f%2s%-4s%2s  \n" % (
                    nat,
                    atomname,
                    alter,
                    comp[i],
                    asym2,
                    nres,
                    inst,
                    u[0],
                    u[1],
                    u[2],
                    u[3],
                    u[4],
                    u[5],
                    space,
                    asym[i],
                    symbol[i],
                )
                fw.write(ss)
                n = n + 1

    if nmodel > 1:
        fw.write("ENDMDL  \n")
    fw.write("END\n")
    fw.close()

    return pdbf, ccdmap


##########################################################
def _add_twin_pdb(fw, flist):
    """Add twin records to PDB, if exist"""

    items, values = cifparse(flist, "_pdbx_reflns_twin.")  # a loop

    if not items:
        return

    domain_id = parse_values(items, values, "_pdbx_reflns_twin.domain_id")
    operator = parse_values(items, values, "_pdbx_reflns_twin.operator")
    fraction = parse_values(items, values, "_pdbx_reflns_twin.fraction")

    ndomain = len(domain_id)
    fw.write("REMARK   3  TWIN DETAILS \n")
    fw.write("REMARK   3   NUMBER OF TWIN DOMAINS  :   %s \n" % ndomain)

    for i in range(ndomain):
        if domain_id:
            fw.write("REMARK   3      TWIN DOMAIN   :  %s\n" % (domain_id[i]))
        if operator:
            fw.write("REMARK   3      TWIN OPERATOR :  %s\n" % (operator[i]))
        if fraction:
            fw.write("REMARK   3      TWIN FRACTION :  %s\n" % (fraction[i]))

    fw.write("REMARK   3    \n")


##########################################################
def _add_tls_pdb(fw, flist):
    """Add TLS records to PDB, if exist"""

    prog, _version = get_prog(flist)

    items, values = cifparse(flist, "_ccp4_refine_tls.")  # a loop
    if items:
        tab = "_ccp4_refine_tls"
    else:
        items, values = cifparse(flist, "_pdbx_refine_tls.")
        tab = "_pdbx_refine_tls"

    # refid = parse_values(items, values, '%s.pdbx_refine_id' % tab)
    tlsid = parse_values(items, values, "%s.id" % tab)
    # method = parse_values(items, values, '%s.method' % tab)
    xo = parse_values(items, values, "%s.origin_x" % tab)
    yo = parse_values(items, values, "%s.origin_y" % tab)
    zo = parse_values(items, values, "%s.origin_z" % tab)
    t11 = parse_values(items, values, "%s.T[1][1]" % tab)
    t22 = parse_values(items, values, "%s.T[2][2]" % tab)
    t33 = parse_values(items, values, "%s.T[3][3]" % tab)
    t12 = parse_values(items, values, "%s.T[1][2]" % tab)
    t13 = parse_values(items, values, "%s.T[1][3]" % tab)
    t23 = parse_values(items, values, "%s.T[2][3]" % tab)
    l11 = parse_values(items, values, "%s.L[1][1]" % tab)
    l22 = parse_values(items, values, "%s.L[2][2]" % tab)
    l33 = parse_values(items, values, "%s.L[3][3]" % tab)
    l12 = parse_values(items, values, "%s.L[1][2]" % tab)
    l13 = parse_values(items, values, "%s.L[1][3]" % tab)
    l23 = parse_values(items, values, "%s.L[2][3]" % tab)
    s11 = parse_values(items, values, "%s.S[1][1]" % tab)
    s12 = parse_values(items, values, "%s.S[1][2]" % tab)
    s13 = parse_values(items, values, "%s.S[1][3]" % tab)
    s21 = parse_values(items, values, "%s.S[2][1]" % tab)
    s22 = parse_values(items, values, "%s.S[2][2]" % tab)
    s23 = parse_values(items, values, "%s.S[2][3]" % tab)
    s31 = parse_values(items, values, "%s.S[3][1]" % tab)
    s32 = parse_values(items, values, "%s.S[3][2]" % tab)
    s33 = parse_values(items, values, "%s.S[3][3]" % tab)

    items, values = cifparse(flist, "%s_group." % tab)  # a loop
    # refid1 = parse_values(items, values, '%s_group.pdbx_refine_id' % tab)
    # groupid = parse_values(items, values, '%s_group.id' % tab)
    tlsid1 = parse_values(items, values, "%s_group.refine_tls_id" % tab)
    asymid1 = parse_values(items, values, "%s_group.beg_auth_asym_id" % tab)
    seqid1 = parse_values(items, values, "%s_group.beg_auth_seq_id" % tab)
    asymid2 = parse_values(items, values, "%s_group.end_auth_asym_id" % tab)
    seqid2 = parse_values(items, values, "%s_group.end_auth_seq_id" % tab)
    detail = parse_values(items, values, "%s_group.selection_details" % tab)

    if not (tlsid and tlsid1):
        return

    bpart1 = 0
    for x in flist:
        if "ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY" in x:
            bpart1 = 1
            break

    if not (is_a_number(tlsid[-1]) and is_a_number(tlsid1[-1]) and tlsid[-1] == tlsid1[-1]):
        print("Error: Number of TLS in %s & %s_group differ." % (tab, tab))
        return

    fw.write("REMARK   3  TLS DETAILS. \n")
    fw.write("REMARK   3   NUMBER OF TLS GROUPS: %s \n" % tlsid[-1])
    if bpart1:
        fw.write("REMARK   3   ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY\n")

    if "REFMAC" in prog:
        fw.write("REMARK   3  \n")
    else:
        fw.write("REMARK   3   ORIGIN: CENTER OF MASS \n")

    nlen = len(tlsid1)
    n = 0
    tmp = []
    for i in range(len(tlsid1)):
        if n >= len(tlsid):
            print("Error: Number of TLS group in _pdbx_refine_tls not consistent with _pdbx_refine_tls_group.")
            break

        if "REFMAC" not in prog and i < nlen - 1:
            if tlsid1[i] == tlsid1[i + 1]:
                continue
        elif "REFMAC" in prog and i < nlen - 1:
            if tlsid1[i] == tlsid1[i + 1]:
                tmp.append([asymid1[i], seqid1[i], asymid2[i], seqid2[i]])
                continue

        fw.write("REMARK   3   TLS GROUP : %2d \n" % (n + 1))
        ntmp = len(tmp) + 1

        if "REFMAC" in prog:
            fw.write("REMARK   3    NUMBER OF COMPONENTS GROUP : %d\n" % ntmp)
            fw.write("REMARK   3    COMPONENTS        C SSSEQI   TO  C SSSEQI \n")
            tmp.append([asymid1[i], seqid1[i], asymid2[i], seqid2[i]])
            for y in tmp:
                fw.write("REMARK   3    RESIDUE RANGE :   %s %6s    %s %6s \n" % (y[0], y[1], y[2], y[3]))
            tmp = []
        elif detail:
            details = detail[i].replace("\n", " ")
            _write_select_phenix_buster_pdb(fw, details, prog)
        elif "PHENIX" in prog or "BUSTER" in prog:
            print("Error: problem in atom selection for the TLS group (%d)" % n)
        else:
            print("Warning: structure refined by program (%s), but TLS not supported by DCC." % prog)

        if is_a_number(xo[n]) and is_a_number(yo[n]) and is_a_number(zo[n]):
            if "BUSTER" in prog:
                orig = "%9.4f%10.4f%10.4f " % (float(xo[n]), float(yo[n]), float(zo[n]))

            else:
                orig = "%8.3f %8.3f %8.3f " % (float(xo[n]), float(yo[n]), float(zo[n]))
        else:
            orig = "   NULL     NULL     NULL "

        fw.write("REMARK   3    ORIGIN FOR THE GROUP (A): %s \n" % orig)
        fw.write("REMARK   3    T TENSOR \n")
        fw.write("REMARK   3      T11: %8s T22: %8s \n" % (t11[n], t22[n]))
        fw.write("REMARK   3      T33: %8s T12: %8s \n" % (t33[n], t12[n]))
        fw.write("REMARK   3      T13: %8s T23: %8s \n" % (t13[n], t23[n]))

        fw.write("REMARK   3    L TENSOR \n")
        fw.write("REMARK   3      L11: %8s L22: %8s \n" % (l11[n], l22[n]))
        fw.write("REMARK   3      L33: %8s L12: %8s \n" % (l33[n], l12[n]))
        fw.write("REMARK   3      L13: %8s L23: %8s \n" % (l13[n], l23[n]))

        fw.write("REMARK   3    S TENSOR \n")
        fw.write("REMARK   3      S11: %8s S12: %8s S13: %8s \n" % (s11[n], s12[n], s13[n]))
        fw.write("REMARK   3      S21: %8s S22: %8s S23: %8s \n" % (s21[n], s22[n], s23[n]))
        fw.write("REMARK   3      S31: %8s S32: %8s S33: %8s \n" % (s31[n], s32[n], s33[n]))
        n = n + 1

    fw.write("REMARK   3    \n")
    fw.write("REMARK   3  BULK SOLVENT MODELLING. \n")
    # print(items)


##########################################################
def _write_select_phenix_buster_pdb(fw, detail, prog):
    select = "REMARK   3    SELECTION:"
    select1 = "REMARK   3             :"
    if "BUSTER" in prog:
        select = "REMARK   3    SET :"
        select1 = "REMARK   3        :"

    if len(detail) > 55:
        t = detail.split()
        s, ss = "", ""
        for j, x in enumerate(t):
            if j < len(t) - 1 and len(s) >= 46:
                ss = ss + s + "\n%s " % select1
                s = ""
            elif j == len(t) - 1:
                ss = ss + s + x
                break
            s = s + x + " "
            # print(s, ss)

        fw.write("%s %s\n" % (select, ss))
    else:
        fw.write("%s %s\n" % (select, detail))


##########################################################
def _add_header_pdb(fw, flist):
    """Add simple stuff for the head"""

    items, values = cifparse(flist, "_database_2.")
    id = parse_values(items, values, "_database_2.database_id")  # pylint: disable=redefined-builtin
    code = parse_values(items, values, "_database_2.database_code")
    pdbid = ["XXXX"]
    if id and code:
        for i, x in enumerate(code):
            if id[i] == "PDB":
                pdbid = [x.strip()]
                break

    name = parse_values(items, values, "_struct_keywords.pdbx_keywords")
    if not name:
        name = [""]

    items, values = cifparse(flist, "_pdbx_database_related.")
    ids = parse_values(items, values, "_pdbx_database_related.db_id")
    type = parse_values(items, values, "_pdbx_database_related.content_type")  # pylint: disable=redefined-builtin
    pdbids = ""
    if ids and type:
        for i, x in enumerate(type):
            if x == "split":
                pdbids = pdbids + ids[i] + " "

    fw.write("HEADER    %-39s %16s    \n" % (name[0].upper(), pdbid[0].upper()))
    if len(pdbids) > 2:
        fw.write("SPLIT      %s  \n" % pdbids)  # legacy

    items, values = cifparse(flist, "_exptl.")
    method = parse_values(items, values, "_exptl.method")
    ss_method = ""
    if method:
        for x in method:
            ss_method = ss_method + x.strip() + "; "
    if ss_method.strip():
        ss = "EXPDTA    %-s\n" % ss_method
        fw.write(ss)
    return ss_method


##########################################################
def _add_remark3(fw, flist, method):
    """Add simple stuff for REMARK3"""

    items, values = cifparse(flist, "_refine.")
    refid = parse_values(items, values, "_refine.pdbx_refine_id")

    nobs = parse_values(items, values, "_refine.ls_number_reflns_obs")
    if not nobs or "?" in nobs:
        nobs = ["NULL", "NULL"]
    resh = parse_values(items, values, "_refine.ls_d_res_high")
    if not resh or "?" in resh:
        resh = ["NULL", "NULL"]
    resl = parse_values(items, values, "_refine.ls_d_res_low")
    if not resl or "?" in resl:
        resl = ["NULL", "NULL"]
    comp = parse_values(items, values, "_refine.ls_percent_reflns_obs")
    if not comp or "?" in comp:
        comp = ["NULL", "NULL"]
    robs = parse_values(items, values, "_refine.ls_R_factor_obs")
    if not robs or "?" in robs:
        robs = ["NULL", "NULL"]
    rwork = parse_values(items, values, "_refine.ls_R_factor_R_work")
    if not rwork or "?" in rwork:
        rwork = ["NULL", "NULL"]
    if rwork[0] == "NULL" and robs[0] != "NULL":
        rwork = robs
    rfree = parse_values(items, values, "_refine.ls_R_factor_R_free")
    if not rfree or "?" in rfree:
        rfree = ["NULL", "NULL"]
    compfr = parse_values(items, values, "_refine.ls_percent_reflns_R_free")
    if not compfr or "?" in compfr:
        compfr = ["NULL", "NULL"]
    nfree = parse_values(items, values, "_refine.ls_number_reflns_R_free")
    if not nfree or "?" in nfree:
        nfree = ["NULL", "NULL"]

    items, values = cifparse(flist, "_database_PDB_remark.")
    text_id = parse_values(items, values, "_database_PDB_remark.id")
    text = parse_values(items, values, "_database_PDB_remark.text")
    bres = ""
    for i, x in enumerate(text_id):
        if x == "3":
            if "U VALUES      : RESIDUAL ONLY" in text[i]:
                bres = "U VALUES      : RESIDUAL ONLY"

    prog, version = get_prog(flist)
    idd = ["REMARK   3 ", "REMARK   3 "]
    if refid and len(refid) == 2 and "X-RAY " in refid[0]:
        idd[0] = "REMARK   3  X-RAY DATA. "
        idd[1] = "REMARK   3  NEUTRON DATA. "
    elif refid and len(refid) == 2 and "NEUTRON " in refid[0]:
        idd[0] = "REMARK   3  NEUTRON DATA. "
        idd[1] = "REMARK   3  X-RAY DATA. "

    wds = """REMARK   2
REMARK   3
REMARK   3 REFINEMENT.
REMARK   3   PROGRAM     : %s  %s
REMARK   3
""" % (
        prog,
        version,
    )
    fw.write(wds)

    wds = """REMARK   3
%s
REMARK   3
REMARK   3  REFINEMENT TARGET : NULL
REMARK   3
REMARK   3  DATA USED IN REFINEMENT.
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : %s
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : %s
REMARK   3   COMPLETENESS (WORKING+TEST)   (%%) : %s
REMARK   3   NUMBER OF REFLECTIONS             : %s
REMARK   3
REMARK   3  FIT TO DATA USED IN REFINEMENT.
REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT
REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM
REMARK   3   R VALUE            (WORKING SET) : %s
REMARK   3   FREE R VALUE                     : %s
REMARK   3   FREE R VALUE TEST SET SIZE   (%%) : %s
REMARK   3   FREE R VALUE TEST SET COUNT      : %s
REMARK   3  %s
REMARK   3
""" % (
        idd[0],
        resh[0],
        resl[0],
        comp[0],
        nobs[0],
        rwork[0],
        rfree[0],
        compfr[0],
        nfree[0],
        bres,
    )
    fw.write(wds)

    if "X-RAY" in method and "NEUTRON" in method and len(resh) > 1:
        wds = """REMARK   3
%s
REMARK   3
REMARK   3  REFINEMENT TARGET : NULL
REMARK   3
REMARK   3  DATA USED IN REFINEMENT.
REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : %s
REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : %s
REMARK   3   COMPLETENESS (WORKING+TEST)   (%%) : %s
REMARK   3   NUMBER OF REFLECTIONS             : %s
REMARK   3
REMARK   3  FIT TO DATA USED IN REFINEMENT.
REMARK   3   CROSS-VALIDATION METHOD          : THROUGHOUT
REMARK   3   FREE R VALUE TEST SET SELECTION  : RANDOM
REMARK   3   R VALUE            (WORKING SET) : %s
REMARK   3   FREE R VALUE                     : %s
REMARK   3   FREE R VALUE TEST SET SIZE   (%%) : %s
REMARK   3   FREE R VALUE TEST SET COUNT      : %s
REMARK   3  %s
REMARK   3
""" % (
            idd[1],
            resh[1],
            resl[1],
            comp[1],
            nobs[1],
            rwork[1],
            rfree[1],
            compfr[1],
            nfree[1],
            bres,
        )

        fw.write(wds)


##########################################################
def _add_remark285(fw, flist):
    """ """
    # for 285
    items, values = cifparse(flist, "_pdbx_database_remark.")
    text = parse_values(items, values, "_pdbx_database_remark.text")
    if text and "FOLLOWING TRANSFORMATION MATRIX" in text[0]:
        t = text[0].split("\n")
        for x in t:
            if "REMARK 285" not in x and len(x) < 72:
                fw.write("REMARK 285 %-s\n" % x)

        fw.write("REMARK 285 \n")

    # for 350
    items, values = cifparse(flist, "_pdbx_struct_oper_list.")
    if items:
        btype = parse_values(items, values, "_pdbx_struct_oper_list.type")
        b11 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[1][1]")
        b12 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[1][2]")
        b13 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[1][3]")
        t1 = parse_values(items, values, "_pdbx_struct_oper_list.vector[1]")

        b21 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[2][1]")
        b22 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[2][2]")
        b23 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[2][3]")
        t2 = parse_values(items, values, "_pdbx_struct_oper_list.vector[2]")

        b31 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[3][1]")
        b32 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[3][2]")
        b33 = parse_values(items, values, "_pdbx_struct_oper_list.matrix[3][3]")
        t3 = parse_values(items, values, "_pdbx_struct_oper_list.vector[3]")

        if not (b11 and b12 and b13 and t1 and b21 and b22 and b23 and t2 and b31 and b32 and b33 and t3):
            return
        n = 0
        for i in range(len(b11)):
            if btype and "transform " in btype[i]:
                continue
            n = n + 1
            v = [float(x) for x in (b11[i], b12[i], b13[i], t1[i]) if is_a_number(x)]
            if len(v) == 4:
                fw.write("REMARK 350   BIOMT1%4d %9.6f %9.6f %9.6f %14.5f \n" % (n, v[0], v[1], v[2], v[3]))
            v = [float(x) for x in (b21[i], b22[i], b23[i], t2[i]) if is_a_number(x)]
            if len(v) == 4:
                fw.write("REMARK 350   BIOMT2%4d %9.6f %9.6f %9.6f %14.5f \n" % (n, v[0], v[1], v[2], v[3]))
            v = [float(x) for x in (b31[i], b32[i], b33[i], t3[i]) if is_a_number(x)]
            if len(v) == 4:
                fw.write("REMARK 350   BIOMT3%4d %9.6f %9.6f %9.6f %14.5f \n" % (n, v[0], v[1], v[2], v[3]))


##########################################################
def _add_scale_pdb(fw, flist):
    """some programs need it."""

    items, values = cifparse(flist, "_atom_sites.")
    if items:
        b11 = parse_values(items, values, "_atom_sites.fract_transf_matrix[1][1]")
        b12 = parse_values(items, values, "_atom_sites.fract_transf_matrix[1][2]")
        b13 = parse_values(items, values, "_atom_sites.fract_transf_matrix[1][3]")
        t1 = parse_values(items, values, "_atom_sites.fract_transf_vector[1]")

        b21 = parse_values(items, values, "_atom_sites.fract_transf_matrix[2][1]")
        b22 = parse_values(items, values, "_atom_sites.fract_transf_matrix[2][2]")
        b23 = parse_values(items, values, "_atom_sites.fract_transf_matrix[2][3]")
        t2 = parse_values(items, values, "_atom_sites.fract_transf_vector[2]")

        b31 = parse_values(items, values, "_atom_sites.fract_transf_matrix[3][1]")
        b32 = parse_values(items, values, "_atom_sites.fract_transf_matrix[3][2]")
        b33 = parse_values(items, values, "_atom_sites.fract_transf_matrix[3][3]")
        t3 = parse_values(items, values, "_atom_sites.fract_transf_vector[3]")

        if not (b11 and b12 and b13 and t1 and b21 and b22 and b23 and t2 and b31 and b32 and b33 and t3):
            return

        for i in range(len(b11)):
            v = [float(x) for x in (b11[i], b12[i], b13[i], t1[i]) if is_a_number(x)]
            if not v:
                continue
            fw.write("SCALE1     %9.6f %9.6f %9.6f %14.5f \n" % (v[0], v[1], v[2], v[3]))
            v = [float(x) for x in (b21[i], b22[i], b23[i], t2[i])]
            fw.write("SCALE2     %9.6f %9.6f %9.6f %14.5f \n" % (v[0], v[1], v[2], v[3]))
            v = [float(x) for x in (b31[i], b32[i], b33[i], t3[i])]
            fw.write("SCALE3     %9.6f %9.6f %9.6f %14.5f \n" % (v[0], v[1], v[2], v[3]))


##########################################################
def _add_ncs(fw, flist):
    """mainly for virus"""

    items, values = cifparse(flist, "_struct_ncs_oper.")
    if items:
        id = parse_values(items, values, "_struct_ncs_oper.id")  # pylint: disable=redefined-builtin
        code = parse_values(items, values, "_struct_ncs_oper.code")
        b11 = parse_values(items, values, "_struct_ncs_oper.matrix[1][1]")
        b12 = parse_values(items, values, "_struct_ncs_oper.matrix[1][2]")
        b13 = parse_values(items, values, "_struct_ncs_oper.matrix[1][3]")

        b21 = parse_values(items, values, "_struct_ncs_oper.matrix[2][1]")
        b22 = parse_values(items, values, "_struct_ncs_oper.matrix[2][2]")
        b23 = parse_values(items, values, "_struct_ncs_oper.matrix[2][3]")

        b31 = parse_values(items, values, "_struct_ncs_oper.matrix[3][1]")
        b32 = parse_values(items, values, "_struct_ncs_oper.matrix[3][2]")
        b33 = parse_values(items, values, "_struct_ncs_oper.matrix[3][3]")

        t1 = parse_values(items, values, "_struct_ncs_oper.vector[1]")
        t2 = parse_values(items, values, "_struct_ncs_oper.vector[2]")
        t3 = parse_values(items, values, "_struct_ncs_oper.vector[3]")

        if not (id and code and b11 and b12 and b13 and t1 and b21 and b22 and b23 and t2 and b31 and b32 and b33 and t3):
            return
        for i in range(len(b11)):
            idd = " "
            if code[i] == "given":
                idd = "1"
            v = [float(x) for x in (b11[i], b12[i], b13[i], t1[i]) if is_a_number(x)]
            if len(v) != 4:
                print("Error: matrix can not be parsed or not a numerical value.")
                continue
            fw.write("MTRIX1%4d %9.6f %9.6f %9.6f %14.5f %4s \n" % (i + 1, v[0], v[1], v[2], v[3], idd))
            v = [float(x) for x in (b21[i], b22[i], b23[i], t2[i])]
            fw.write("MTRIX2%4d %9.6f %9.6f %9.6f %14.5f %4s \n" % (i + 1, v[0], v[1], v[2], v[3], idd))
            v = [float(x) for x in (b31[i], b32[i], b33[i], t3[i])]
            fw.write("MTRIX3%4d %9.6f %9.6f %9.6f %14.5f %4s \n" % (i + 1, v[0], v[1], v[2], v[3], idd))


##########################################################
def _add_remark200(fw, flist):
    """only for source"""
    items, values = cifparse(flist, "_diffrn_source.")
    syc = parse_values(items, values, "_diffrn_source.pdbx_synchrotron_y_n")
    if not syc:
        syc = parse_values(items, values, "_diffrn_source.ndb_synchrotron_y_n")
    if not syc:
        syc = ["NULL"]
    site = parse_values(items, values, "_diffrn_source.pdbx_synchrotron_site")
    if not site:
        site = parse_values(items, values, "_diffrn_source.ndb_synchrotron_site")
    if not site:
        site = ["NULL"]
    beam = parse_values(items, values, "_diffrn_source.pdbx_synchrotron_beamline")
    if not beam:
        beam = parse_values(items, values, "_diffrn_source.pdbx_synchrotron_beamline")
    if not beam:
        beam = ["NULL"]
    wave = parse_values(items, values, "_diffrn_source.pdbx_wavelength_list")
    if not wave:
        wave = parse_values(items, values, "_diffrn_source.rcsb_wavelength_list")
    if not wave:
        wave = ["NULL"]

    wds = """REMARK 200
REMARK 200  SYNCHROTRON              (Y/N) : %s
REMARK 200  RADIATION SOURCE               : %s
REMARK 200  BEAMLINE                       : %s
REMARK 200  X-RAY GENERATOR MODEL          : NULL
REMARK 200  MONOCHROMATIC OR LAUE    (M/L) : M
REMARK 200  WAVELENGTH OR RANGE        (A) : %s
REMARK 200
""" % (
        syc[0],
        site[0].strip(),
        beam[0].strip(),
        wave[0].strip(),
    )

    fw.write(wds)


##########################################################
def _asym2chain(asym):
    """If asym has length>2, assign it to 2 letters"""

    ch = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789"

    as2ch = {}
    # print set(asym)
    uch = set(asym)
    id = 0  # pylint: disable=redefined-builtin
    for x in uch:
        if len(x) > 2:
            id = 1
            break
    nch = len(uch)
    if id == 0 or nch > 35 * 35:
        return as2ch

    n = 0
    ch2 = []
    for x in ch:
        for y in ch:
            n = n + 1
            xy = "%s%s" % (x, y)
            ch2.append(xy)
            if n > nch + 1:
                break
        if n > nch + 1:
            break

    for i, x in enumerate(uch):
        as2ch[x] = ch2[i]

    return as2ch
