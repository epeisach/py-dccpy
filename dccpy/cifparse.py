##########################################################
# A cif parser in python (2013-01-06 ), HY
##########################################################

import os
import sys


##########################################################
def cifparse(flist, cate):
    """parse cif tokens. flist is a list of content; cate is the category"""

    items, values = [], []
    # remove the empty string and front end space
    line = [x.strip() for x in flist if not x.isspace()]

    m, loop, start = 0, 0, 0
    nline, clen = len(line), len(cate)
    for i in range(0, nline):  # look for the table
        if line[i][:clen] == cate:
            if i > 0 and line[i - 1][:5] == "loop_":
                loop = 1
            start = i
            m = m + 1
            break  # find start position and loop or nonloop

    if m == 0:
        # print("Error: cif table (%s) is not found." %cate)
        return items, values

    if loop == 0:
        for i in range(start, nline):
            val = line[i].split()
            m = len(val)
            if m > 1 and val[0][:clen] == cate:  # ciftoken & value in same line
                ss = ""
                if val[1][0] == "'":
                    ss = _string_between_char(line[i], "'", "'")
                elif val[1][0] == '"':
                    ss = _string_between_char(line[i], '"', '"')
                else:
                    ss = val[1]

                items.append(val[0])
                values.append(ss)

            elif m == 1 and i < nline - 1 and val[0][:clen] == cate:
                ss = ""
                if line[i + 1][0] == ";":  # go through lines until hit ; or
                    ssnew = []
                    n = 0
                    for j in range(i + 1, nline):
                        if n > 0 and line[j][0] == ";":
                            start = j
                            break
                        ssnew.append(line[j])
                        n = n + 1
                        if line[j] == "loop_" or line[j] == "#" or line[j][:clen] == cate:
                            print('Error: can not find the second ";" in (%s). Stop here!' % cate)
                            start = j
                            break

                    ss = "\n".join(ssnew)
                    if ss[0] == ";":
                        ss = ss[1:]

                elif line[i + 1][0] == "'" or line[i + 1][0] == '"':
                    cc = line[i + 1][0]
                    ss = line[i + 1][1:]
                    if ss[-1] != cc:
                        print("Error: in table (%s). End of string is not (%s)" % (cate, cc))
                    else:
                        ss = ss[:-1]
                elif cate not in line[i + 1][:clen]:
                    ss = line[i + 1]

                items.append(val[0])
                values.append(ss)

            elif line[i][0] == "#" or line[i][:5] == "loop_" or (line[i][0] == "_" and cate not in line[i][:clen]):
                break

    else:  # parse looped items
        for i in range(start, nline):
            val = line[i].split()

            m = len(val)
            if m == 1 and val[0][:clen] == cate:  # get items
                items.append(val[0])
            elif (len(items) > 0 and i < nline - 1 and cate not in line[i + 1]) or i == nline - 1:
                start = i
                break

        j = start
        while j < nline:  # through list from start position
            if line[j][0] == ";":
                m, ss = _lines_between_char(j, line, ";")
                values.append(ss)
                j = m

                continue
            elif line[j][0] == "#" or line[j][:clen] == "loop_" or (line[j][0] == "_" and "." in line[j]):
                break
            else:
                tmp = _clean_str(line[j])

                m1 = 0
                while m1 < 800:
                    m1 = m1 + 1

                    if tmp[0] == "'" or tmp[0] == '"':
                        cc = tmp[0]
                        m, ss = _string_after_char(0, tmp, cc)
                        tmp = tmp[m:].lstrip()
                    else:
                        m, ss = _string_after_char(0, tmp, " ")
                        tmp = tmp[m:].lstrip()

                    values.append(ss)
                    if m == 0 or not tmp:
                        break
            j = j + 1

    # ## finshed parsing
    nitem, nval = len(items), len(values)
    n = nval % nitem
    nrow = nval / nitem

    if nval % nitem != 0:
        print("Parsing (%s) is NOT succesful!  Row=%d : nitem=%d : loop=%d " % (cate, nrow, nitem, loop))

    return items, values


##########################################################
def get_rows(items, values):
    """items: all the cif items for the table;  values: all the values in a big row.
    return the row of list [[...] ...[...]]
    """

    row, nval, nitem = [], len(values), len(items)
    for i in range(nval):
        if (i + 1) % nitem == 0:
            row.append(values[i + 1 - nitem : i + 1])

    return row


##########################################################
def parse_values(items, values, item):
    m = -1
    for i, x in enumerate(items):
        if x == item:
            m = i
            break
    if m == -1:
        # print("Error: the item (%s) is not in the cif file." %item)
        return []

    nitem = len(items)
    nrow = len(values) // nitem
    val = []
    for i in range(nrow):
        val.append(values[m + i * nitem])

    return val


##########################################################
#  Below are the applications
##########################################################
def _cif_nonloop_format(rows, items):  # pylint: disable=unused-argument
    """get the maximum length in the items"""

    fmt = []  # the writing format
    for x in items:
        n = len(x)
        fmt.append(n)
    return max(fmt)


##########################################################
def _cif_loop_format(rows):
    """the columns are formated (the max length is determined for each column)"""

    fmt = []  # the writing format
    nrow, ncol = len(rows), len(rows[0])
    for i in range(ncol):  # put column format in a list
        tmp = []
        for j in range(nrow):
            if " " in rows[j][i] or "'" in rows[j][i] or '"' in rows[j][i]:
                n = len(rows[j][i]) + 2  # add two " or ' in the item
            else:
                n = len(rows[j][i])
            tmp.append(n)
        fmt.append(max(tmp))

    return fmt


##########################################################
def _string_between_char(ss, c1, c2):
    """return the string between c1 & c2; ss is the input."""

    return ss[ss.find(c1) + 1 : ss.rfind(c2)]


##########################################################
def _string_after_char(k, ss, cc):
    """get one string starting from k and return index and the new string"""

    end = k
    ssnew = ss[k:]
    nss = len(ss)
    if cc.isspace():
        n = 0
        for i in range(k, nss):
            if n > 0 and (ss[i].isspace() or i == nss - 1):
                if ss[i].isspace():
                    ssnew = ss[k:i]
                elif i == nss - 1:  # the last one
                    ssnew = ss[k : i + 1]
                end = i + 1
                break
            elif not ss[i].isspace():
                n = n + 1

    else:
        n = 0
        for i in range(k + 1, nss):
            if n > 0 and ((ss[i] == cc and i < nss - 1 and ss[i + 1] == " ") or i == nss - 1):
                ssnew = ss[k + 1 : i]
                end = i + 1
                break
            elif not ss[i].isspace():
                n = n + 1
    return end - k, ssnew


##########################################################
def _lines_between_char(k, lines, cc):
    """k: the starting position; return cc lines and end position"""

    ss, ssnew = "", []
    n = 0
    for i in range(k, len(lines)):
        if n > 0 and lines[i] == cc:
            end = i + 1
            break
        ssnew.append(lines[i])
        n = n + 1
    ss = "\n".join(ssnew)
    if ss[0] == cc:
        ss = ss[1:]

    return end, ss


##########################################################
def _clean_str(ss):
    ssnew = ss.split()
    return " ".join(ssnew)


##########################################################
def cif2cif(file):
    """parse cif and rewrite in new cif"""

    if not (os.path.exists(file) and os.path.getsize(file) > 0):
        print("Error: file (%s) is empty or zero size")
        return ""

    blockid = "data_xxxx\n"
    nfile = file + "_new.cif"
    fw = open(nfile, "w")

    flist = open(file, "r").readlines()
    nblock = []
    for i, x in enumerate(flist):
        if "data_" in x.lstrip()[:5] and len(x.strip().split("data_")[1]) > 0:
            nblock.append(i)

    if nblock:  # has data block name
        nblock.append(len(flist))
        for i in range(1, len(nblock)):
            alist = flist[nblock[i - 1] : nblock[i]]
            blockid = flist[nblock[i - 1]]
            _cif2cif_single_block(fw, blockid, alist)
    else:  # no data block name
        alist = flist[0 : len(flist)]
        _cif2cif_single_block(fw, blockid, alist)

    return nfile


##########################################################
def _cif2cif_single_block(fw, block, flist):
    """parse cif and rewrite in new cif. block: the line with block id.; flist: a list"""

    fw.write(block)
    cif = _cif_table_items(flist)
    for x in cif:
        items, values = cifparse(flist, x[1])
        rows = get_rows(items, values)

        n = len(values) / len(items)
        if n == 1:  # non-looped
            _write_cif_non_loop(fw, items, values, rows)
        elif n > 1:  # looped
            _write_cif_loop(fw, items, rows)


##########################################################
def _write_cif_non_loop(fw, items, values, rows):
    fmt = _cif_nonloop_format(rows, items)
    fw.write("#\n")
    for i, p in enumerate(items):
        v = values[i].strip()
        if not v:
            v = "?"
        if "\n" in v:
            fw.write("%s  \n" % p)
            fw.write(";%s\n;\n" % v)
        elif " " in v or '"' in v or "'" in v:
            fw.write(p.ljust(fmt + 1))
            if "'" in v:
                t = '"%s"' % v
            else:
                t = "'%s'" % v
            fw.write("   %s\n" % t)
        else:
            fw.write(p.ljust(fmt + 1))
            fw.write("   %s\n" % v)


##########################################################
def _write_cif_loop(fw, items, rows):
    fmt = _cif_loop_format(rows)
    fw.write("\n#\nloop_\n")
    for p in items:
        fw.write("%s\n" % p)
    for row in rows:  # (left formated)
        for ii, zz in enumerate(row):
            yy = zz.strip()
            if not yy:
                yy = "?"
            if len(yy) > 60:
                fw.write("\n;%s\n;\n" % yy)
            else:
                if " " in yy or '"' in yy or "'" in yy:  # related to fmt
                    if "'" in yy:
                        t = '"%s"' % yy
                    else:
                        t = "'%s'" % yy
                    fw.write(t.rjust(fmt[ii] + 1))
                else:
                    fw.write(yy.rjust(fmt[ii] + 1))
        fw.write("\n")


##########################################################
def _cif_table_items(flist_all):
    """get all the block_category_items from the file. Return a list.
    list=[[block, table, [items..], [block, table, [items..]...]
    """

    flist = []
    for x in flist_all:
        y = x.strip()
        if len(y) > 5 and "_" in y[0] or "loop_" in y[:5] or "data_" in y[:5]:
            flist.append(y)

    # print('file= ', file, len(flist), flist)

    cif = []
    nlen = len(flist)
    st = 0
    block = "X"
    while st < nlen:
        x = flist[st].strip()
        if len(x) > 5 and "data_" in x[:5]:
            block = x[5:].split()[0]

        n, cate, item = check_item(x)
        if n > 0:  # potentional cif item
            cates, items = "", []
            items.append(item)
            if len(items) == 1:
                cates = cate
            if st == nlen - 1:
                cif.append([block, cates, items])

            for j in range(st + 1, nlen):
                y = flist[j].strip()
                if len(y) > 5 and "data_" in y[:5]:
                    block = y[5:].split()[0]
                m, cate, item = check_item(y)

                if len(items) > 0 and ("loop_" in y[:5] or (m > 0 and cate != cates)):
                    st = j - 1
                    cif.append([block, cates, items])
                    cates, items = "", []
                    break
                elif (m > 0 and cate == cates) and j == nlen - 1:
                    st = j
                    items.append(item)
                    cif.append([block, cates, items])
                    cates, items = "", []
                    break
                if m > 0:
                    items.append(item)
        st = st + 1
    return cif


##########################################################
def check_item(y):
    """y: a one line string: check if it is a possible cif item (n>0)"""

    n = 0
    ss = y.split()[0]
    nlen = len(ss)
    cate, item = "", ""
    if nlen > 4 and ss[0] == "_" and "." in ss[1:]:
        m = ss.find(".") + 1
        cate = ss[:m]
        item = ss[m:]
        if m > 1 and "." not in ss[m + 1 :] and len(cate) > 2 and len(item) > 0:
            n = 1

    return n, cate, item


##########################################################
def get_cell(flist):
    """flist is a list of the cif file
    return cell as a list of float!
    """

    cell = [0, 0, 0, 0, 0, 0]
    items, values = cifparse(flist, "_cell.")
    if not len(items):
        return cell
    a = parse_values(items, values, "_cell.length_a")
    b = parse_values(items, values, "_cell.length_b")
    c = parse_values(items, values, "_cell.length_c")
    alpha = parse_values(items, values, "_cell.angle_alpha")
    beta = parse_values(items, values, "_cell.angle_beta")
    gamma = parse_values(items, values, "_cell.angle_gamma")

    if not (a and b and c and alpha and beta and gamma):
        print("Warning: cells not extracted. Check ciftokens")

    # print(a , b , c , alpha , beta , gamma)
    for i, x in enumerate([a, b, c, alpha, beta, gamma]):
        if len(x) == 0 or not is_a_number(x[0]):
            print("Error: cell has wrong (%s) values" % x)
            continue
        cell[i] = float(x[0].strip())

    return cell


##########################################################
def get_symm(flist):
    """get symmetry"""

    spg = ""
    items, values = cifparse(flist, "_symmetry.")
    symm = parse_values(items, values, "_symmetry.space_group_name_H-M")
    if symm:
        spg = symm[0].replace("'", "").replace('"', "").strip()
    # else:
    #    print('Warning: space group not extracted. Check ciftokens')

    return spg


##########################################################
def get_xyz(flist):
    """From list of lines of cif data, extact coordinates"""
    items, values = cifparse(flist, "_atom_site.")  # a loop
    group = parse_values(items, values, "_atom_site.group_PDB")
    natm = parse_values(items, values, "_atom_site.id")
    symbol = parse_values(items, values, "_atom_site.type_symbol")
    if not symbol:
        symbol = parse_values(items, values, "_atom_site.atom_type_symbol")
    atom = parse_values(items, values, "_atom_site.label_atom_id")
    asym = parse_values(items, values, "_atom_site.auth_asym_id")
    comp = parse_values(items, values, "_atom_site.label_comp_id")
    seq = parse_values(items, values, "_atom_site.auth_seq_id")
    biso = parse_values(items, values, "_atom_site.B_iso_or_equiv")
    occ = parse_values(items, values, "_atom_site.occupancy")
    x = parse_values(items, values, "_atom_site.Cartn_x")
    y = parse_values(items, values, "_atom_site.Cartn_y")
    z = parse_values(items, values, "_atom_site.Cartn_z")
    ins = parse_values(items, values, "_atom_site.pdbx_PDB_ins_code")
    if not ins:
        ins = parse_values(items, values, "_atom_site.ndb_ins_code")
    alt = parse_values(items, values, "_atom_site.label_alt_id")
    model = parse_values(items, values, "_atom_site.pdbx_PDB_model_num")
    if not model:
        model = parse_values(items, values, "_atom_site.pdbx_model")
    if not model:
        model = parse_values(items, values, "_atom_site.ndb_model")

    return group, natm, atom, asym, seq, comp, ins, alt, x, y, z, occ, biso, model, symbol


##########################################################
def get_uij(flist):
    """pars the anisou records"""

    items, values = cifparse(flist, "_atom_site_anisotrop.")  # a loop

    unatm = parse_values(items, values, "_atom_site_anisotrop.id")

    uatom = parse_values(items, values, "_atom_site_anisotrop.pdbx_auth_atom_id")
    if not uatom:
        uatom = parse_values(items, values, "_atom_site_anisotrop.ndb_label_atom_id")
        if not uatom:
            uatom = parse_values(items, values, "_atom_site_anisotrop.ndb_PDB_atom_name")

    uasym = parse_values(items, values, "_atom_site_anisotrop.pdbx_auth_asym_id")
    if not uasym:
        uasym = parse_values(items, values, "_atom_site_anisotrop.ndb_PDB_strand_id")
        if not uasym:
            uasym = parse_values(items, values, "_atom_site_anisotrop.ndb_auth_asym_id")

    useq = parse_values(items, values, "_atom_site_anisotrop.pdbx_auth_seq_id")
    if not useq:
        useq = parse_values(items, values, "_atom_site_anisotrop.ndb_auth_seq_id")
        if not useq:
            useq = parse_values(items, values, "_atom_site_anisotrop.ndb_PDB_residue_no")

    ucomp = parse_values(items, values, "_atom_site_anisotrop.pdbx_auth_comp_id")
    if not ucomp:
        ucomp = parse_values(items, values, "_atom_site_anisotrop.ndb_label_comp_id")
        if not ucomp:
            ucomp = parse_values(items, values, "_atom_site_anisotrop.ndb_PDB_residue_name")

    uins = parse_values(items, values, "_atom_site_anisotrop.pdbx_PDB_ins_code")
    if not uins:
        uins = parse_values(items, values, "_atom_site_anisotrop.ndb_label_ins_code")

    ualt = parse_values(items, values, "_atom_site_anisotrop.pdbx_label_alt_id")
    if not ualt:
        ualt = parse_values(items, values, "_atom_site_anisotrop.ndb_label_alt_id")

    u11 = parse_values(items, values, "_atom_site_anisotrop.U[1][1]")
    u22 = parse_values(items, values, "_atom_site_anisotrop.U[2][2]")
    u33 = parse_values(items, values, "_atom_site_anisotrop.U[3][3]")
    u12 = parse_values(items, values, "_atom_site_anisotrop.U[1][2]")
    u13 = parse_values(items, values, "_atom_site_anisotrop.U[1][3]")
    u23 = parse_values(items, values, "_atom_site_anisotrop.U[2][3]")

    return unatm, uatom, uasym, useq, ucomp, uins, ualt, u11, u22, u33, u12, u13, u23


##########################################################
def get_prog(flist):
    """get program and version"""

    prog, vers = "?", "?"

    items, values = cifparse(flist, "_software.")
    p = parse_values(items, values, "_software.name")
    d = parse_values(items, values, "_software.classification")
    v = parse_values(items, values, "_software.version")
    if p and d:
        for i, x in enumerate(d):
            if "refinement" in x:
                prog = p[i].strip()
                vers = v[i].strip().replace(" ", "")
                if "PHENIX" in prog.upper() and "REFINE" not in prog.upper():
                    vers = "(phenix.refine: %s)" % vers
                break

    else:
        items, values = cifparse(flist, "_computing.")
        p = parse_values(items, values, "_computing.structure_refinement")
        if p:
            prog = p[0]

    return prog, vers


##########################################################
def is_a_number(s):
    try:
        float(s)
        return True
    except ValueError:
        # print('Error: %s is not a number.' %s)
        return False


##########################################################
def usage():
    content = """

-------------------------------------------------------------------
 A program for fast parsing mmcif files. (HY,2013-01-10)
-------------------------------------------------------------------
      How to use the cif parser (see the main below)


-------------------------------------------------------------------
              Defination of valid mmcif format:
1. for non-looped items, below are the valid syntaxes
_citation.title "The Structure and Evolution of the Major Capsid Protein"
or
_citation.title 'The Structure and Evolution of the Major Capsid Protein'
or
_citation.title
;The Structure and Evolution of
the Major Capsid Protein
;

2. for looped items, below are the valid syntaxes
loop_
_citation_author.citation_id
_citation_author.name
_citation_author.ordinal
primary '???, N.' 1
primary '???, A.'     2

or
primary
"???, N." 1
primary
'???, A.'     2

or
primary
;???, N.
;
1
primary
'???, A.'     2

-------------------------------------------------------------------
"""
    print(content)
    sys.exit()


#     '''
# ##########################################################
# if __name__ == '__main__':

#     files=sys.argv[1]
#     print('file=%s' %files)
#     flist=open(files, 'r').readlines()
#     items,values = cifparse(flist, '_cell.')
#     v1=parse_values(items,values, '_cell.entry_id')
#     v2=parse_values(items,values, '_cell.length_c')
#     v3=parse_values(items,values, '_cell.angle_beta')
#     print(v1,v2,v3)

#     items,values = cifparse(flist, '_database_PDB_rev.')
#     v1=parse_values(items,values, '_database_PDB_rev.num')
#     v2=parse_values(items,values, '_database_PDB_rev.date')
#     v3=parse_values(items,values, '_database_PDB_rev.status')

#     print(v1,v2,v3)

# ##########################################################

#     '''
