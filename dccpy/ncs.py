# ===========================================================
# this module is to check all the NCS related problems
# (created 2013-09-24)
# ===========================================================
import util
import cifparse as cif


def check_ncs(file):
    '''full check of the NCS groups (file in pdb or cif)
    '''

    if util.is_cif(file) :
        check_ncs_cif(file)
    else:
        # print 'Please input a cif file.'
        return
        # check_ncs_pdb(file)


####################################################################
def check_ncs_cif(file):
    '''check all errors in the cif file:
    check the four cif tables agaist the scheme
    '''

    flist = open(file, 'r').readlines()

    ncs = ncs_from_head(flist)

    items, values = cif.cifparse(flist, '_struct_ncs_ens.')
    ens_id = cif.parse_values(items, values, '_struct_ncs_ens.id')

    items, values = cif.cifparse(flist, '_refine_ls_restr_ncs.')
    res_ord = cif.parse_values(items, values, '_refine_ls_restr_ncs.pdbx_ordinal')
    res_ref = cif.parse_values(items, values, '_refine_ls_restr_ncs.pdbx_refine_id')
    res_ens = cif.parse_values(items, values, '_refine_ls_restr_ncs.pdbx_ens_id')
    res_dom = cif.parse_values(items, values, '_refine_ls_restr_ncs.dom_id')
    res_typ = cif.parse_values(items, values, '_refine_ls_restr_ncs.pdbx_type')
    res_asy = cif.parse_values(items, values, '_refine_ls_restr_ncs.pdbx_auth_asym_id')
    res_num = cif.parse_values(items, values, '_refine_ls_restr_ncs.pdbx_number')
    res_rms = cif.parse_values(items, values, '_refine_ls_restr_ncs.rms_dev_position')

    if '?' in res_rms:
        res_rms = cif.parse_values(items, values, '_refine_ls_restr_ncs.pdbx_rms')
    # res_wgh = cif.parse_values(items, values, '_refine_ls_restr_ncs.weight_position')

    items, values = cif.cifparse(flist, '_struct_ncs_dom.')
    dom_ens = cif.parse_values(items, values, '_struct_ncs_dom.pdbx_ens_id')
    dom_id = cif.parse_values(items, values, '_struct_ncs_dom.id')
    # dom_all = cif.parse_values(items, values, '_struct_ncs_dom.details')

    items, values = cif.cifparse(flist, '_struct_ncs_dom_lim.')
    lim_ens = cif.parse_values(items, values, '_struct_ncs_dom_lim.pdbx_ens_id')
    lim_dom = cif.parse_values(items, values, '_struct_ncs_dom_lim.dom_id')
    lim_com = cif.parse_values(items, values, '_struct_ncs_dom_lim.pdbx_component_id')
    lim_basy = cif.parse_values(items, values, '_struct_ncs_dom_lim.beg_auth_asym_id')
    lim_bseq = cif.parse_values(items, values, '_struct_ncs_dom_lim.beg_auth_seq_id')
    lim_easy = cif.parse_values(items, values, '_struct_ncs_dom_lim.end_auth_asym_id')
    lim_eseq = cif.parse_values(items, values, '_struct_ncs_dom_lim.end_auth_seq_id')
    # lim_all = cif.parse_values(items, values, '_struct_ncs_dom_lim.selection_details')

    if len(ncs):
        if not (res_ens or dom_ens or lim_ens or ens_id):
            util.perror('Warning: NCS records exist, but no cif tables for it.')
            return
    else:
        if not len(ens_id):
            return

    if not (res_ens or res_dom or res_ord or res_ref or res_asy or res_typ):
        util.perror('Warning: No cif table (refine_ls_restr_ncs) or missing key items for NCS.')
    if not (dom_ens or dom_id):
        util.perror('Warning: No cif table (struct_ncs_dom) or missing key items for NCS.')
    if not (lim_ens or lim_dom or lim_com):
        util.perror('Warning: No cif table (struct_ncs_dom_lim) or missing key items for NCS.')
    if not (ens_id):
        util.perror('Warning: No cif table (struct_ncs_ens) for NCS.')

    chain = get_chain_seq(flist)  # chain is dic for int

    def tmp(list1, list2):  # put ensemble=key and doms=[] to a dic
        tmpd = {}
        for j, y in enumerate(list1):
            if y not in tmpd.keys():
                tmpd[y] = []
            tmpd[y].append(list2[j])
        return tmpd

    dom = tmp(dom_ens, dom_id)
    lim = tmp(lim_ens, lim_dom)
    res = tmp(res_ens, res_dom)

    def tmp1(list1, list2, s2):
        for n in list1:
            if n not in list2:
                util.perror('Error: NCS ID (%s) not in table (%s).' % (n, s2))
    if dom:
        tmp1(ens_id, dom.keys(), 'struct_ncs_dom')
    if lim:
        tmp1(ens_id, lim.keys(), 'struct_ncs_dom_lim')
    if res:
        tmp1(ens_id, res.keys(), 'refine_ls_restr_ncs')

    # print chain.keys(), chain, dom, lim, res

    if (lim_basy and lim_easy and lim_bseq and lim_eseq):
        check_lim(lim_ens, chain, lim_basy, lim_easy, lim_bseq, lim_eseq)

    check_res(res_ens, chain, res_ref, res_asy, res_typ, res_num, res_rms)


####################################################################
def check_lim(lim_ens, chain, lim_basy, lim_easy, lim_bseq, lim_eseq):
    '''check struct_ncs_dom_lim
    '''

    tab = 'struct_ncs_dom_lim'
    for i, _x in enumerate(lim_ens):  #
        ch1, ch2, ns1, ns2 = lim_basy[i], lim_easy[i], lim_bseq[i], lim_eseq[i]
        if ch1 != ch2:
            util.perror('Warning: two chainID (%s & %s) are not the same. look at row %2d in %s.' % (ch1, ch2, i + 1, tab))

        if ch1 not in chain.keys():
            util.perror('Warning: chainID (%s) not in coordinate. look at row %2d in %s.' % (ch1, i + 1, tab))

        if ((ch1 in chain.keys() and int(ns1) not in chain[ch1])):
            util.perror('Warning: residue number (%s_%s) not in coordinate. look at row %2d in %s.' % (ch1, ns1, i + 1, tab))

        if ((ch2 in chain.keys() and int(ns2) not in chain[ch2])):
            util.perror('Warning: residue number (%s_%s) not in coordinate. look at row %2d in %s.' % (ch2, ns2, i + 1, tab))


####################################################################
def check_res(res_ens, chain, res_ref, res_asy, res_type, res_num, res_rms):
    '''refine_ls_restr_ncs
    '''

    tab = 'refine_ls_restr_ncs'
    for i, _x in enumerate(res_ens):
        if res_ref and 'X-RAY DIFFRACTION' not in res_ref[i]:
            util.perror('Warning: refine_id must be "X-RAY DIFFRACTION". look at row %2d in %s.' % (i + 1, tab))

        if res_asy and res_asy[i] not in chain.keys():
            util.perror('Warning: chainID (%s) not in coordinate. look at row %2d in %s.' % (res_asy[i], i + 1, tab))

        if res_num and util.is_number(res_num[i]) and int(res_num[i]) == 0:
            util.perror('Warning: Restrain number (or atom pairs)=0. look at row %2d in %s.' % (i + 1, tab))

        if (res_type and ('POSITIONAL' not in res_type[i].upper()
                          and 'TORSIONAL' not in res_type[i].upper()
                          and 'LOCAL' not in res_type[i].upper()
                          and 'INTERATOMIC DISTANCE' not in res_type[i].upper()
                          and 'THERMAL' not in res_type[i].upper())):
            util.perror('Warning: wrong restrain type (%s). look at row %2d in %s.' % (res_type[i], i + 1, tab))

        if res_rms:
            rms = res_rms[i]
            if util.is_number(rms):
                if ((res_type and 'POSITIONAL' in res_type[i] and float(rms) > 2)
                        or (res_type and 'TORTIONAL' in res_type[i] and float(rms) > 15)):
                    util.perror('Warning: RMSD (%s) is large. look at row %2d in %s.' % (rms, i + 1, tab))
            else:
                util.perror('Warning: RMSD (%s) is not a number. look at row %2d in %s. ' % (rms, i + 1, tab))


####################################################################
def get_chain_seq(flist):
    '''return a dic (each chain has a list nres).
    '''

    chain = {}

    items, values = cif.cifparse(flist, '_atom_site.')
    asym = cif.parse_values(items, values, "_atom_site.auth_asym_id")
    seq = cif.parse_values(items, values, "_atom_site.auth_seq_id")

    if not (asym or seq):
        return chain

    ch_old, ns_old = '', ''

    for x in asym:
        ch = x
        if ch != ch_old:
            chain[ch] = []
            ch_old = ch

    for i, x in enumerate(seq):
        ns = x
        ch = asym[i]

        if ns != ns_old and util.is_number(ns):
            chain[ch].append(int(ns))
            ns_old = ns

    # print chain
    return chain


####################################################################
def ncs_from_head(flist):
    '''parse all the NCS info from the PDB-like header
    '''

    n, ncs_pos , ncs = 0, [], []
    for i, x in enumerate(flist):
        if n == 0 and ' NCS GROUP :' in x:
            n = i
        elif n > 0 and ' NCS GROUP :' in flist[i + 1]:
            ncs.append(flist[n:i])
            n = 0
        elif n > 0 and (' Histogram of ' in x or ';' in x or 'OTHER REFINEMENT' in x):
            ncs.append(flist[n:i - 1])
            n = 0
            break

    if len(ncs_pos) == 1:
        print('Warning: no ending mark of NCS is found')
        n1 = ncs_pos[0]
        ncs = [n1, len(flist)]

    # for x in ncs: print x
    return ncs
