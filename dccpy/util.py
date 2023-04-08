# ===========================================================
# this module contains the utility functions
# ===========================================================


import os
import sys
import math
import shutil
import re
import config
import cifparse as cif


##########################################################
def gsort(data, col, idd):
    '''A general function to sort a list (col=-1) or a list of list (col>=0)
    idd : 1 for sorting from small to large; -1 from large to small
    '''

    if col < 0:  # a single list [..]
        if idd >= 0:
            data.sort()
        else:
            data.sort(reverse=True)

    else:  # a list of list [[..], [..]]
        if idd >= 0:
            data.sort(key=lambda y: y[col])
        else:
            data.sort(key=lambda y: y[col], reverse=True)


######################################################################
def limits_by_iqr(alist, col, scale):
    '''Using the quartile method (IQR) to select the outliers
    col: the column of [[col,..], [...]]
    scale : a scale factor (1.5, 2.2, 3.0). scale=2.22 is equivalent to 3 standard deviations
    '''

    if not alist:
        return 0, 0, 0

    nt = len(alist)
    n1 = int(nt / 4.)
    n2 = int(nt / 2.)  # the median
    n3 = int(nt * 3 / 4.)

    # print 'numbers=',  nt, n1,n2,n3

    if col < 0:
        gsort(alist, -1, 0)
        q1, q2, q3 = float(alist[n1]), float(alist[n2]), float(alist[n3])
    else:
        alist.sort(key=lambda y: y[col])  # sort fp by the col
        q1, q2, q3 = float(alist[n1][col]), float(alist[n2][col]), float(alist[n3][col])

    iqr = q3 - q1
    b1 = q1 - iqr * scale
    b2 = q3 + iqr * scale

    return b1, b2, q2


##########################################################
def mean_dev(data_in, col):
    '''get the min, max, average and deviation
    if col=-1, data_in is a list of data, else a list of list!
    '''
    if len(data_in) == 0:
        return 0, 0, 0, 0
    data = []
    if col < 0:
        data = data_in
    else:
        for x in data_in:
            data.append(float(x[col]))
    # print(data)
    mini, maxi = min(data), max(data)
    n = len(data)
    if n == 0:
        return 0, 0, 0, 0
    avg = sum(data) / float(n)
    s1 = 0
    for x in data:
        a2 = (x - avg) * (x - avg)
        s1 = s1 + a2
    dev = math.sqrt(s1 / float(n))

    return avg, dev, mini, maxi


##########################################################
def chain_res_atom(pdb):
    ''' parse the pdb  into chain-residnum-atom list. a dictionary
    dd[x][y]: atom list, x:chainID, y:residue_number_in_string,
    dd={ch:{resn:[atom_lines]}}
    '''

    if check_file(300, pdb) == 0:
        return

    fp = open(pdb, 'r')

    id = 0
    dd = {}
    for x in fp:
        if not ('ATOM' in x[:4] or 'HETATM' in x[:6] or len(x.strip()) < 50):
            continue
        ch = x[20:22].strip()
        if ch == ' ':
            id = 1

        if ch not in dd.keys() and ch != ' ':
            dd[ch] = {}
        res = x[22:27]  # include inserted
        if ch != ' ' and ch in dd.keys() and res not in dd[ch].keys():
            dd[ch][res] = []
        if ch != ' ':
            dd[ch][res].append(x)

    if id == 1:
        print('Warning: file (%s) has no ChainID.' % pdb)
    return dd


##########################################################
def delete_file(*files):
    for x in files:
        os.system('rm -f ' + x)


##########################################################
def myreplace(before, after, ss):
    '''before: string to be replaced by after. ss is the string
    '''
    return re.sub(before, after, ss)


##########################################################
def myreplace_case(before, after, ss):
    '''before: string to be replaced by after. ss is the string
    '''
    beforein = '(?i)%s' % before
    return re.sub(beforein, after, ss)


##########################################################
def check_file(size, *files):

    if len(files) == 0:
        print('Error: File size is not specified. check_file(size, file)')
        return 0

    n = 1
    for f in files:
        if not os.path.exists(f) or os.path.getsize(f) < size:
            # print('Error: file (%s) does not exist (or file size=0).' %f)
            n = 0
            break
    return n


##########################################################
def str_after_id(line, id):
    """ get string after a given id """

    if id not in line:
        print('Warning: %s not in string (%s)' % (id, line))
        return line

    n = line.index(id)
    value = line[n + len(id):].strip()
    if len(value) <= 0 or value.upper() == 'NULL':
        value = "?"
    return value


##########################################################
def float_after_id(line, id):
    """ get float after a given id """

    if id not in line:
        print('Warning: %s not in string.' % (id))
        return 9999.0

    n = line.index(id)
    li = line[n + 1:].strip().split()
    if len(li) <= 0:
        print('Error: No value after id (%s).' % (id))
        return 9999.0
    else:
        if is_number(li[0]):
            return float(li[0])
        else:
            print('Error: %s is not a numeric.' % li[0])
            return 9999.0


##########################################################
def int_after_id(line, id):
    """ get int after a given id """

    if id not in line:
        print('Warning: %s not in string.' % (id))
        return 9999

    n = line.index(id)
    li = line[n + 1:].strip().split()
    if len(li) <= 0:
        print('Error: No value after id (%s).' % (id))
        return 9999
    else:
        if is_number(li[0]):
            return int(li[0])
        else:
            print('Error: %s is not a numeric.' % li[0])
            return 9999


##########################################################
def str_between_id(line, id1, id2):
    """ get string between the two ids """

    if id1 not in line or id2 not in line :
        print('Warning: either %s or %s not in string.' % (id1, id2))
        return "?"
    n1, n2 = line.index(id1) + len(id1), line.index(id2)
    return line[n1:n2].strip()


##########################################################
def float_between_id(line, id1, id2):
    """ get float value between the two ids """

    if id1 not in line or id2 not in line:
        print('Warning: either %s or %s not in string' % (id1, id2))
        return 9999.
    n1, n2 = line.index(id1), line.index(id2)
    return float(line[n1 + 1:n2])


##########################################################
def int_between_id(line, id1, id2):
    """ get int value between the two ids """

    if id1 not in line or id2 not in line :
        print('Warning: either %s or %s not in string' % (id1, id2))
        return 9999
    n1, n2 = line.index(id1), line.index(id2)
    return int(line[n1 + 1:n2])


##########################################################
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        # print('Error: %s is not a number.' %s)
        return False


##########################################################
def is_digit(s):

    n = 0
    if len(s) == 0:
        return 0
    if s[0:1] == '-':
        if s[1:].isdigit():
            n = 1
    else:
        if s.isdigit():
            n = 1

    return n


##########################################################
def get_value_after_id(line, id):
    """ get value after a given id """

    if id not in line:
        return "?"
    li = line.split(id)
    value = li[1].strip()
    if len(value) <= 0 or value.upper() == 'NULL':
        value = "?"
    return value


##########################################################
def is_bigcif(file):
    '''big cif has chainID more than one charactors
    return 0: not cif
    '''
    n = 0
    flist = open(file, 'r').readlines()
    items, values = cif.cifparse(flist, '_atom_site.')  # a loop
    asym = cif.parse_values(items, values, "_atom_site.auth_asym_id")
    if asym :
        for x in asym:
            if len(x) > 1:
                n = 2
                print('Note: Cif file (%s) has chainID >1 charactor.' % file)
                break

    return n


##########################################################
def is_cif(file):
    '''fast detect if it is a mmcif file (n>0)
    '''

    n, m = 0, 0
    if not os.path.exists(file):
        return n

    fp = open(file, 'r')
    for x in fp:
        t = x.strip()
        if len(t) < 4 or t[0] == '#':
            continue
        m = m + 1
        if m > 100:
            break

        if (('data_' in t[:5] and len(t) > 5) or ('loop_' in t[:5] and len(t) == 5)):
            n = 1
            break
        elif (t[0] == '_' and '.' in t):
            n1, cate, item = cif.check_item(t)
            if n1 > 0:
                n = 1
                break

        elif ('HEADER  ' in x[:8] or 'COMPND  ' in x[:8] or 'CRYST1 ' in x[:7]
              or 'ATOM  ' in x[:6]):
            break

    fp.close()
    return n


##########################################################
def value_between_string(ln, s1, s2):
    ''' get the string between two strings (s1,s2) in a line(case sensitive)
    '''
    if s1 not in ln or s2 not in ln:
        print('Warning: string (%s or %s) not in %s' % (s1, s2, ln))
        return ''
    n = ln.index(s1) + len(s1)
    t1 = ln[n:]
    n = t1.index(s2)
    return t1[:n].strip()


##########################################################
def residue():
    '''non-ligand residues
    '''

    res = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'CYS', 'MET', 'PRO', 'PHE', 'TRP', 'TYR',
           'HIS', 'LYS', 'ARG', 'ASP', 'GLU', 'ASN', 'GLN', 'THR', 'SER', 'MSE',
           'HOH', 'DOD',
           '  A', '  G', '  C', '  T', ' DA', ' DG', ' DC', ' DT', ' DI', '  U', '  I',
           'ADE', 'THY', 'GUA', 'CYT', 'URA']

    return res


##########################################################
def protein():
    res = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'CYS', 'MET', 'PRO', 'PHE', 'TRP', 'TYR',
           'HIS', 'LYS', 'ARG', 'ASP', 'GLU', 'ASN', 'GLN', 'THR', 'SER', 'MSE']

    return res


##########################################################
def dnarna():
    res = ['A', 'G', 'C', 'T', 'DA', 'DG', 'DC', 'DT', 'DI', 'U', 'I',
           'ADE', 'THY', 'GUA', 'CYT', 'URA']

    return res


##########################################################
def dna_rna():
    res = ['  A', '  G', '  C', '  T', ' DA', ' DG', ' DC', ' DT', ' DI', '  U', '  I',
           'ADE', 'THY', 'GUA', 'CYT', 'URA']

    return res


##########################################################
def perror(info):
    '''print error messages
    '''
    if info not in config.ERRLOG:
        config.ERRLOG.append(info)
        print(info.strip())


##########################################################
def space_group_crystal_match(spg):
    '''
    The dic contains all the chiral space groups and crystal system.
    http://www.ruppweb.org/Xray/comp/space_instr.htm
    '''
    spg_cryst = {

        # TRICLINIC
        "A1": 1,
        "P1": 1,
        "P1-": 1,
        "P-1": 1,

        # MONOCLINIC

        "A121": 2,
        "A2": 2,
        "B112": 20,  # C is unique axis. alpha, beta, =90
        "B2": 20,  # C is unique axis. alpha, beta, =90
        "C2": 2,
        "C121": 2,
        "C21": 2,
        "C1211": 2,
        "C2(A112)": 2,
        "I121": 2,
        "I1211": 2,
        "I2": 2,
        "I21": 2,
        "P2": 2,
        "P121": 2,
        "P112": 20,
        "P21": 2,
        "P1211": 2,
        "P1121": 20,
        "P21(C)": 2,

        # ORTHORHOMBIC

        "P222": 3,
        "P2221": 3,
        "P21212": 3,
        "P212121": 3,
        "P22121": 3,
        "P21221": 3,
        "P21212A": 3,
        "B2212": 3,
        "C222": 3,
        "C2221": 3,
        "I222": 3,
        "I212121": 3,
        "F222": 3,

        # TETRAGONAL

        "P4": 4,
        "P41": 4,
        "P42": 4,
        "P43": 4,
        "P422": 4,
        "P4212": 4,
        "P4122": 4,
        "P41212": 4,
        "P4222": 4,
        "P42212": 4,
        "P4322": 4,
        "P43212": 4,
        "I4": 4,
        "I41": 4,
        "I422": 4,
        "I4122": 4,
        "I-42D": 4,

        # TRIGONAL   #hexagonal axis a=b gamma=120

        "P3": 5,
        "P31": 5,
        "P32": 5,
        "P312": 5,
        "P321": 5,
        "P3112": 5,
        "P3121": 5,
        "P3212": 5,
        "P3221": 5,
        "R3": 50,  # (rhombohedral axis, a=b=c & alpha=beta=gamma)
        "R32": 50,  # (rhombohedral axis, a=b=c & alpha=beta=gamma)
        "H3": 5,
        "H32": 5,

        # HEXAGONAL

        "P6": 6,
        "P61": 6,
        "P62": 6,
        "P63": 6,
        "P64": 6,
        "P65": 6,
        "P622": 6,
        "P6122": 6,
        "P6222": 6,
        "P6322": 6,
        "P6422": 6,
        "P6522": 6,

        # CUBIC

        "C4212": 7,  # ?
        "F422": 7,  # ?

        "P23": 7,
        "F23": 7,
        "I23": 7,
        "P213": 7,
        "I213": 7,

        "P432": 7,
        "P4132": 7,
        "P4232": 7,
        "P4332": 7,
        "F432": 7,
        "F4132": 7,
        "I432": 7,
        "I4132": 7,

        # others
        "P121/c1": 2,
        "P121/n1": 2,
        "I41/a": 4,
        "I-4C2": 4,
        "I-42d": 4
    }

    tmp = spg.replace(' ', '')
    if tmp in spg_cryst.keys():
        return spg_cryst[tmp]
    else:
        print('Warning: The space group (%s) is not in the list' % spg)
        return 0


##########################################################
def sg_nsym(sg):
    '''contains space group names and the number of operators.
    '''

    sym = {
        "A 1": 2,
        "A 1 2 1": 4,
        "A 2": 4,
        "B 1 1 2": 4,
        "B 2": 4,
        "B 2 21 2": 8,
        "C 2": 4,
        "C 1 2 1": 4,
        "C 21": 4,
        "C 1 21 1": 4,
        "C 2(A 112)": 4,
        "C 2 2 2": 8,
        "C 2 2 21": 8,
        "C 4 21 2": 16,
        "F 2 2 2": 16,
        "F 2 3": 48,
        "F 4 2 2": 32,
        "F 4 3 2": 96,
        "F 41 3 2": 96,
        "I 1 2 1": 4,
        "I 1 21 1": 4,
        "I 2": 4,
        "I 2 2 2": 8,
        "I 2 3": 24,
        "I 21": 4,
        "I 21 3": 24,
        "I 4": 8,
        "I 21 21 21": 8,
        "I 4 2 2": 16,
        "I 4 3 2": 48,
        "I 41": 8,
        "I 41 2 2": 16,
        "I 41 3 2": 48,
        "P 1": 1,
        "P -1": 2,
        "P 2": 2,
        "P 1 2 1": 2,
        "P 1 1 2": 2,
        "P 2 2 2": 4,
        "P 2 3": 12,
        "P 2 2 21": 4,
        "P 2 21 21": 4,
        "P 2 21 2": 4,
        "P 21 2 2": 4,
        "P 21": 2,
        "P 1 21 1": 2,
        "P 1 1 21": 2,
        "P 21(C)": 2,
        "P 21 2 21": 4,
        "P 21 3": 12,
        "P 21 21 2": 4,
        "P 21 21 2 A": 4,
        "P 21 21 21": 4,
        "P 3": 3,
        "P 3 1 2": 6,
        "P 3 2 1": 6,
        "P 31": 3,
        "P 31 1 2": 6,
        "P 31 2 1": 6,
        "P 32": 3,
        "P 32 1 2": 6,
        "P 32 2 1": 6,
        "P 4": 4,
        "P 4 2 2": 8,
        "P 4 3 2": 24,
        "P 4 21 2": 8,
        "P 41": 4,
        "P 41 2 2": 8,
        "P 41 3 2": 24,
        "P 41 21 2": 8,
        "P 42": 4,
        "P 42 2 2": 8,
        "P 42 3 2": 24,
        "P 42 21 2": 8,
        "P 43": 4,
        "P 43 2 2": 8,
        "P 43 3 2": 24,
        "P 43 21 2": 8,
        "P 6": 6,
        "P 6 2 2": 12,
        "P 61": 6,
        "P 61 2 2": 12,
        "P 62": 6,
        "P 62 2 2": 12,
        "P 63": 6,
        "P 63 2 2": 12,
        "P 64": 6,
        "P 64 2 2": 12,
        "P 65": 6,
        "P 65 2 2": 12,
        "H 3": 9,
        "R 3": 3,
        "H 3 2": 18,
        "R 3 2": 6 ,
        "I 41/A": 16 ,
        "P 1 21/C 1": 4 ,
        "P 1 21/N 1": 4 ,
        "I -4 C 2": 16,
        "I -4 2 d": 16,
        "P b c n": 8
    }

    if sg in sym:
        nsym = sym[sg]
    else:
        nsym = -1
    return nsym


##########################################################
def atom_mass(atom):
    '''selected mass of atom (only for estimating matthew coeff.)
    '''
    at = {
        "H": 1.0079,
        "D": 1.0079,
        "C": 12.011,
        "N": 14.0067,
        "O": 15.9994,
        "F": 18.9984,
        "MG": 24.305,
        "P": 30.97376,
        "S": 32.06,
        "ZN": 65.38,
        "SE": 78.96,
        "BR": 79.904,
        "X": 40
    }
    mass = 40
    if atom in at:
        mass = at[atom]

    return mass


##########################################################
def residue_mass(resname):
    '''get mass of residue;
       AAD: the residue mass of modified aa; DD is modified NA.
    '''
    resd = {
        "GLY": 57.05,
        "ALA": 71.08,
        "VAL": 99.13,
        "LEU": 113.16,
        "ILE": 113.16,
        "SER": 87.08,
        "THR": 101.10,
        "CYS": 103.14,
        "PRO": 97.12,
        "PHE": 147.18,
        "TYR": 163.18,
        "TRP": 186.21,
        "HIS": 138.15,
        "ASP": 110.05,
        "ASN": 114.10,
        "GLU": 122.06,
        "GLN": 128.13,
        "MET": 131.19,
        "LYS": 129.18,
        "ARG": 157.19,
        "DA": 328.20,
        "DT": 303.19,
        "DG": 344.20,
        "DC": 304.18,
        "A": 328.20,
        "T": 303.19,
        "G": 344.20,
        "C": 304.18,
        "U": 305.16,
        "MG": 24.305,
        "ZN": 65.38
    }
    mass = 0
    if resname in resd:
        mass = resd[resname]
        # elif len(resname)<3:
        # mass=310

    return mass


##########################################################

def get_file_by_pdbid(pdbid_in, idd):
    '''pdbid_in : 4 char PDBID;   idd='cifid|pdbid|cif|pdb|sf'
    '''

    pth = '/net/ftp_tree_v5/ftp-v5.0/pdb/data/structures/divided/'
    # www_path = "http://www.rcsb.org/pdb/files"

    pdbid = pdbid_in.lower()
    if len(pdbid) != 4:
        print('Error: PDBID (%s) is not 4 characters. ' % pdbid)
        sys.exit()

    hash = pdbid[1:3]

    cif = '%s/mmCIF/%s/%s.cif.gz' % (pth, hash, pdbid)
    pdb = '%s/pdb/%s/pdb%s.ent.gz' % (pth, hash, pdbid)
    sf = '%s/structure_factors/%s/r%ssf.ent.gz' % (pth, hash, pdbid)

    sffile = "r" + pdbid + "sf.ent"
    pdbfile = "pdb" + pdbid + ".ent"
    ciffile = pdbid + ".cif"

    if idd == 'pdb':
        os.system("zcat  %s > %s " % (pdb, pdbfile))
        if not check_file(10, pdbfile):
            os.system("zcat  %s > %s " % (cif, pdbfile))  # for large entry
        return pdbfile, ''
    elif idd == 'sf':
        os.system("zcat  %s > %s " % (sf, sffile))
        return '', sffile
    elif idd == 'cif':
        os.system("zcat  %s > %s " % (cif, ciffile))
        return ciffile, ''
    elif idd == 'cifid':
        os.system("zcat  %s > %s " % (cif, ciffile))
        os.system("zcat  %s > %s " % (sf, sffile))
        return ciffile, sffile
    elif idd == 'pdbid':
        os.system("zcat  %s > %s " % (pdb, pdbfile))
        os.system("zcat  %s > %s " % (sf, sffile))
        return pdbfile, sffile
    else:
        print("Error: wrong file format. ")
        sys.exit()


##########################################################
def frac_orth_matrix(cell):
    '''get fractional and orthogonal matrix to apply xyz
    cell is a list containing of six floats
    '''

    frac = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    orth = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    d2r = 3.141592654 / 180
    a = cell[0]
    b = cell[1]
    c = cell[2]

    sa = math.sin(cell[3] * d2r)
    sb = math.sin(cell[4] * d2r)
    sg = math.sin(cell[5] * d2r)

    ca = math.cos(cell[3] * d2r)
    cb = math.cos(cell[4] * d2r)
    cg = math.cos(cell[5] * d2r)

    vol = sg * sb * sa
    frac[0][0] = 1.0 / a
    frac[0][1] = (-cg) / (a * sg)
    frac[0][2] = (cg * ca - cb) / (a * vol * sg)
    frac[1][0] = 0.0
    frac[1][1] = 1.0 / (b * sg)
    frac[1][2] = (cg * cb - ca) / (b * vol * sg)
    frac[2][0] = 0.0
    frac[2][1] = 0.0
    frac[2][2] = (sg) / (c * vol)

    orth[0][0] = a
    orth[0][1] = b * cg
    orth[0][2] = c * cb
    orth[1][0] = 0.0
    orth[1][1] = b * sg
    orth[1][2] = c * (ca - cb * cg) / sg
    orth[2][0] = 0.0
    orth[2][1] = 0.0
    orth[2][2] = c * sb * sa

    return frac, orth


##########################################################
def matrix_prod(f, x):
    '''f: a list of list [[],[],[]];  x: a list []

    '''

    xo = [0, 0, 0]

    xo[0] = x[0] * f[0][0] + x[1] * f[0][1] + x[2] * f[0][2]
    xo[1] = x[0] * f[1][0] + x[1] * f[1][1] + x[2] * f[1][2]
    xo[2] = x[0] * f[2][0] + x[1] * f[2][1] + x[2] * f[2][2]

    return [xo[0], xo[1], xo[2]]


##########################################################
def matrix_prod_(m1, m2):
    '''input matrix m1=m1Xn ;  m2=nXm2 ::  mo = m1Xm2 is the output
    '''
    mo = [[0 for row in range(len(m1))] for col in range(len(m2[0]))]

    for i in range(len(m1)):
        for j in range(len(m2[0])):
            for k in range(len(m2)):
                mo[i][j] += m1[i][k] * m2[k][j]

    return mo


##########################################################
def move(finp, fout):
    if os.path.exists(finp):
        shutil.move(finp, fout)


##########################################################
def position_cif_table(cif, table):
    '''a simple way to get the start-end position of the ciflist ( cif )
    '''

    n1, n2 , m, ntab = 0, 0, len(cif), len(table)
    for i in range(0, m):
        if table in cif[i].strip()[:ntab]:
            n1 = i
            if i > 0 and 'loop_' in cif[i - 1]:
                n1 = i - 1
            for j in range(i, m):
                x = cif[j].lstrip()
                if ('#' in x[:1] or 'loop_' in x[:5] or ('_' in x[:1] and table not in x) or j == m - 1):
                    n2 = j
                    break
            if n2 > n1:
                break

    return n1, n2


##########################################################
def format_data(rows):
    '''format the list: use the maximum length of the column
    rows: a list of n rows and m columns
    returns a formated list of n rows and m columns.
    '''

    newrow = []

    fmt = []  # the writing format
    nrow, ncol = len(rows), len(rows[0])
    for i in range(ncol):  # put column format in a list
        tmp = []
        for j in range(nrow):
            # print i,j,rows[j][i]
            n = len(rows[j][i])
            tmp.append(n)
        fmt.append(max(tmp))

    for x in rows:  # re-write atom_site (left formated)
        rw = []
        for i, y in enumerate(x):
            ss = y.ljust(fmt[i] + 1)
            rw.append(ss)

        newrow.append(rw)

    return newrow


##########################################################
def sort_column(fp, k1, k2):
    ''' Sort list (fp=[[], [] ..]); k1,k2 column position of [])
    '''

    d1 = [[x[:k1], x[k1:k2], x[k2:]] for x in fp]  # separate 3 column
    d1.sort(key=lambda y: y[1])  # sort 2th column

    return d1


##########################################################
def program_other_bin(prog):

    if os.path.exists(config.PATH2 + prog):
        program = config.PATH2 + prog
    else:
        program = prog

    if 'MAPMAN_BINARY' in os.environ and prog == 'lx_mapman':
        program = os.environ['MAPMAN_BINARY']

    return program


##########################################################
def get_software_version(logfile, progid, dic):
    '''extract program versions and insert it to dic['software']
    '''

    if not check_file(10, logfile):
        return

    version = '?'
    fp = open(logfile, 'r').readlines()
    for x in fp:
        if 'MAPMAN' in progid.upper() and 'Created by MAPMAN V.' in x:
            version = x.split(' V.')[1].split()[0]
            ss = ['MAPMAN', version, config.RSR]
            if ss not in dic['software']:
                dic['software'].append(ss)
            break

        elif (('PDBREP' in progid.upper() or 'PDBCAL' in progid.upper())
              and 'REMARK   3   PROGRAM     : ' in x):

            if 'PDBREP' in progid.upper():
                info = config.VALID_REP
            elif 'PDBCAL' in progid.upper():
                info = config.VALID_CAL

            linesplit = x[26:].split()
            prog = '?'
            if len(linesplit) > 0:
                prog = linesplit[0].upper()

            if 'PHENIX' in x and len(linesplit) > 2:
                version = linesplit[2].replace(')', '')
            else:
                if len(linesplit) == 1 or len(linesplit) == 0:
                    version = '?'
                else:
                    version = linesplit[1]

            ss = [prog, version, info]
            if 'software' in dic and ss not in dic['software']:
                dic['software'].append(ss)
            break

        elif ('CNS' in progid.upper() and 'Version:' in x):
            version = x.split('Version:')[1].split()[0].strip()
            ss = ['CNS', version, config.VALID_CAL]
            if ss not in dic['software']:
                dic['software'].append(ss)
            break

        elif ('SHELX' in progid.upper() and 'Sheldrick' in x and 'Release' in x):
            version = x.split('Release')[1].split()[0]
            ss = ['SHELX', version, config.VALID_CAL]
            if ss not in dic['software']:
                dic['software'].append(ss)
            break

        elif ('SF_CONVERT' in progid.upper() and 'VERSION' in x):
            version = x.split('VERSION')[1].split('"')[1].split(':')[0]
            ss = ['SF_CONVERT', version, config.SF_FMT]
            if ss not in dic['software']:
                dic['software'].append(ss)
            break

        elif ('SFCHECK' in progid.upper() and 'SFCHECK' in x and 'Vers ' in x):
            version = x.split('Vers ')[1].split(';')[0]
            ss = ['SFCHECK', version, config.VALID_CAL]
            if ss not in dic['software']:
                dic['software'].append(ss)
            break

        elif ('EDSTAT' in progid.upper() and ' EDSTATS v' in x):
            version = x.split(' EDSTATS v')[1].split()[0]
            ss = ['EDSTAT', version, config.RSR]
            if ss not in dic['software']:
                dic['software'].append(ss)
            break

        elif ('POINTLESS' in progid.upper() and 'version ' in x):
            version = x.split('version ')[1].split(':')[0]
            ss = ['POINTLESS', version, config.DATA]
            if ss not in dic['software']:
                dic['software'].append(ss)
            break


##########################################################
def get_phenix_version_from_dic(dic):
    '''returns phenix version for dic['software']
       returns 1000*maj + min
       returns None if not found
    '''
    ret = 0
    soft = dic['software']

    for s in soft:
        if s[0] == 'PHENIX' and 'model vs diffraction data validation' in s[2]:
            vers = s[1].split('-')[0]
            maj = int(vers.split('.')[0])
            min = int(vers.split('.')[1])
            return maj * 1000 + min

    return ret
