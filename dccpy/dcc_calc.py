#!/usr/bin/env python
# #!/usr/bin/env /apps/python-2.6.1/bin/python  #(for RCSB)

import os
import sys
import shutil
from time import time
import config
import tls
import ligand
import sdsc_ligmap
import prog
import parse
import matt
import util
import ncs
import write_cif
import cifparse as cif
# import pychecker.checker


##########################################################
def process(*args):
    """ Calculate R factors and density correlation and real_space R  using:
    SFCHECK (modified version 7.02.4)
    REFMAC5
    Phenix  (model_vs_data  and xtriage)
    CNS/XPLOR (model_stat)
    LX_MAPMAN
    Other utility programs are:
    CCP4 (sub programs: cif2mtz, fft, cad, mapmask, tlsanl)
    SF_CONVERT
    """

    usage = """

===========================  sf-valid  version %s  ===========================
Simple Usage:   dcc  -pdb  xyzfile  -sf  sffile

  Options below may be added after the above arguments

  -o,       followed by an output file name to hold the calculated statistics.
            otherwise, the default output (xyzfile + _rcc_sum.cif) is assigned.
  -diags,   followed by a logfile name to hold error/warning messages.

  -verb,       add it to keep some of the intermediate files.
  -sfcheck,    add it to validate xray data by SFCHECK.
  -refmac,     add it to validate xray data by REFMAC.
  -phenix_x,   add it to validate xray data by PHENIX.
  -phenix_n,   add it to validate neutron data by PHENIX.
  -phenix_xn,  add it to validate neutron & xray hybrid data by PHENIX.
  -cns,        add it to validate xray data by CNS.
  -all,        add it to calculate R/Rfree by all programs(sfcheck,refmac,cns,phenix)
  -auto,       add it to validate xray data by the program in coordinate.
               If the refinement program fails, then try refmac or others.

  -rsr_all,    add it to calculate rsr rscc for  (overall, side/main chain, phosphate)
  -edstat,     add it to calculate rsr ZDIFF for (overall, side/main chain, phosphate)
  -map,        add it to calculate maps (mFo-DFc, 2mFo-DFc, ccp4 format).
  -mapsize,     add it to reduce map size (2.1 + 0.15*resolution**(3/2.), if>3, assign 3.0)
  -map_dsn6,   add it to calculate maps (mFo-DFc, 2mFo-DFc, dsn6 format).
  -ligmap,     add it to generate density (table,html,jmol) around all ligands.
  -omitmap,    add it to calculate RsR and density corr. after omitting all ligands.
  -fem,        add it to calculate RsR and density corr from the feature-enhanced-map.
  -mtzmap,     calculate map from mtz containing phase (dcc -mtzmap  mtzfile -pdb pdbfile).
  -mapcut,     cut a big map to small pieces. (dcc -mapcut mapfile coordfile id)
               Here, id = model_compound_chainID_resnumber_alter_insertion
  -omit,       omit a given ligand (e.g. dcc -pdb pdbfile -sf sffile -omit A_3:4 ).

  -bfull,       convert residual to full B factors: (dcc -bfull coor_cif).
  -no_xtriage,  add it to not include the xtriage calculation.
  -lib ,  followed by ligand restraint file needed by phenix
  -------------- Some examples are given below -----------------

  Runing REFMAC/MAPMAM for RSR/Zscore for each residue :
       dcc   -pdb pdbfile -sf sffile  -o output

  Runing REFMAC/MAPMAM for RSR/Zscore for main/side chain, residue :
       dcc -rsr_all  -pdb pdbfile -sf sffile  -o output

  Runing REFMAC/EDSTAT for RSR/Zscore for main/side chain, residue :
       dcc -edstat  -pdb pdbfile -sf sffile  -o output

  Using sfcheck only :
       dcc  -sfcheck -pdb pdbfile -sf sffile  -o output

  Using phenix (for xray) only :
       dcc  -phenix_x -pdb pdbfile -sf sffile  -o output

  Using phenix (for neutron) only :
       dcc  -phenix_n -pdb pdbfile -sf sffile  -o output

  Using phenix (for xray & neutron) only :
       dcc  -phenix_xn -pdb pdbfile -sf sffile  -o output

  Using CNS (version 1.3) only :
       dcc  -cns  -pdb pdbfile -sf sffile  -o output

  Running multiple program together(sfcheck, refmac, cns, phenix_x):
       dcc  -pdb pdbfile -sf sffile -all  -o output

    """ % config.VERSION

    '''
  items not shown for public.  (internal use)
  -fofc,   followed by a given file name of Fo-FC map.
  -2fofc,  followed by a given file name of 2Fo-FC map.
  -path,   followed by a full path (e.g. ?/sf-valid/) to locate the executables.
           (if the DCCPY environment is set, do not give path).

  -noeds,       add it (do not calculate EDS).
  -scale,       add it to use 'BULK LSSC ANISO EXPE', otherwise use defult.
  -refmac_tls,  add it to use refmac to get full B factors, then do validation.
  -twin,        add it to force refmac to calculate R/Rfree in twin mode.
  -refine,      add it to do refinement (4 cycles) by refmac.
  -sdsc_map,    generate ligmap required by sdsc.
  -ligmapcif,   add it to generate density (table,html,jmol) around all ligands.(D&A)
  -ligtable,    generate table,html,jmol by inputing dcc_file, dir, mapfile, pdbfile

  -cc,          add it to calculate cc/RsR using perfect FC (Biso=0).
  -lldf,        add it to calculate LLDF (testing).
  -rsrd,        add it to calculate zscores using database. (testing).
  -cif2pdb     convert cif to pdb (dcc -cif2pdb ciffile)
  -cif2cif     reformat cif to standard (dcc -cif2cif ciffile)
    '''

    # dic['software'].append(['DCC', config.VERSION.split()[0],config.DCC])

    print('\n============== sf-valid  version %s  ==============\n' % config.VERSION)

    start_time = time()
    if len(args[0]) < 2:
        print(usage)
        sys.exit()

    pdbfile, sffile = '', ''
    dic = initialize_dic()
    get_argument(args[0], dic)  # all args are stored in dic.
    check_env(dic)  # check all the needed environments
    pdbfile = check_xyz(dic)  # check xyz & generate a new xyz file & update dic

    if (dic['ligmap'] and not dic['ligmapcif']):
        if ((not dic['sdsc_map'] and not check_ligand(pdbfile, 0))
                or (dic['sdsc_map'] and not check_ligand(pdbfile, 1))):
            print('Warning: program exit, No ligand in file=%s' % pdbfile)
            sys.exit()

    if dic['mtzmap']:  # run fft with the mtz containing phase
        mapout1 = '%s_2FO-FC.map' % dic['outfile']
        mapout2 = '%s_FO-FC.map' % dic['outfile']
        prog.mtz2map(dic['mtzmap'], pdbfile, mapout1, '2FO_FC', dic)
        prog.mtz2map(dic['mtzmap'], pdbfile, mapout2, 'FO_FC', dic)
        if dic['display']:
            command = "coot --pdb %s --mtz %s & " % (pdbfile, dic['mtzmap'])
            os.system(command)
        if (dic['verb'] == 0):
            util.delete_file('get_mtz2map.csh')
        print("Time elapsed : %.2f(sec)" % (time() - start_time))

        sys.exit()

    sffile = check_sf(dic)  # get new sf (also for neutron) & update dic
    if util.check_file(10, pdbfile, sffile) == 0:
        return

    if dic['omit'] == 1 and 'omit_residue' in dic:
        pdbfile = omit_residue(dic, pdbfile)

    if (dic['auto'] and (dic['exp'] == 'x' or dic['exp'] == 'e')):
        dic_all = do_validation_auto(dic, pdbfile, sffile)
    else:
        dic_all = do_validation(dic, pdbfile, sffile)

    write_all(dic, dic_all)  # write final results

    map1 = '%s.PDB_sig_2fofc.map' % dic['xyzfile_orig']
    map2 = '%s.PDB_sig_fofc.map' % dic['xyzfile_orig']
    if dic['map'] and util.check_file(40, map1):
        shutil.move(map1, '%s_map-2fofc_P1.map' % dic['xyzfile_orig'])
    if dic['map'] and util.check_file(40, map2):
        shutil.move(map2, '%s_map-fofc_P1.map' % dic['xyzfile_orig'])

    if (dic['verb'] == 0):  # clean files
        util.delete_file('sf_information.cif', 'NATIVE.refmac', dic['cifrefm'])
        if dic['pdb_tls']:
            util.delete_file(dic['pdb_tls'])

        if dic['rsr_all'] or dic['edstat']:
            util.delete_file(dic['cifrefm1'], dic['cifrefm2'], dic['cifrefm3'])
        util.delete_file(dic['sf_xray'], dic['sffile'], dic['pdbfile'], dic['pdbfile_sig'])
        if (dic['rsrd']):
            util.delete_file('TMP_PROMOTIF_FILE.???')
        if dic['xyz_type'] == 'cif':
            util.delete_file(dic['pdbfile_orig'])
        if dic['ligmapcif']:
            util.delete_file(dic['xyzfile_orig'] + '.PDB_new', dic['xyzfile_orig'] + '.PDB_new_new')
            if dic['omitmap']:
                util.delete_file(dic['xyzfile_orig'] + '.PDB_new_NOLIG_new')

    print("Time elapsed : %.2f(sec)" % (time() - start_time))


##########################################################
def get_argument(arg, dic):
    '''get all the input arguments
    '''

    narg = len(arg)
    for k in range(narg):
        if (arg[k].upper() == "-PDB" or arg[k].upper() == "-CIF"):
            dic['xyzfile_orig'] = arg[k + 1]

        elif (arg[k].upper() == "-SF"):
            dic['sffile_orig'] = arg[k + 1]

        elif (arg[k].upper() == "-PDBID"):
            dic['xyzfile_orig'], dic['sffile_orig'] = util.get_file_by_pdbid(arg[k + 1], 'pdbid')

        elif (arg[k].upper() == "-CIFID"):
            dic['xyzfile_orig'], dic['sffile_orig'] = util.get_file_by_pdbid(arg[k + 1], 'cifid')

        elif (arg[k].upper() == "-O"):
            dic['outfile'] = arg[k + 1]

        elif (arg[k].upper() == "-MAPIN"):
            dic['mapfile'] = arg[k + 1]
        elif (arg[k].upper() == "-FOFC"):
            dic['fofc'] = arg[k + 1]
        elif (arg[k].upper() == "-2FOFC"):
            dic['2fofc'] = arg[k + 1]
        elif (arg[k].upper() == "-OMITMAP"):
            dic['omitmap'] = 1
        elif (arg[k].upper() == "-OMIT"):
            dic['omit'] = 1
            dic['omit_residue'] = arg[k + 1]   # such as A2:2 (A2) or A2:3 (A2 or A3)
        elif (arg[k].upper() == "-NOEDS"):
            dic['noeds'] = 1
        elif (arg[k].upper() == "-DIAGS"):
            dic['diags'] = arg[k + 1]
        elif (arg[k].upper() == "-DISPLAY"):
            dic['display'] = 1
        elif (arg[k].upper() == "-MAP_DSN6"):
            dic['map_dsn6'] = 1
        elif (arg[k].upper() == "-MAP" or arg[k].upper() == "-MAP_CCP4"):  # CCP4 type
            dic['map'] = 1
        elif (arg[k].upper() == "-MTZMAP"):  # generate a map from a mtz
            dic['mtzmap'] = arg[k + 1]

        elif (arg[k].upper() == "-SDSC_MAP"):  # generate all ligand  maps for SDSC
            dic['sdsc_map'] = 1
        elif (arg[k].upper() == "-LIGMAP" or arg[k].upper() == "-LIGMAPCIF"):
            # generate all the ligand related stuff
            dic['ligmap'] = 1
            dic['refmac'] = 1
            if arg[k].upper() == "-LIGMAPCIF":
                dic['ligmapcif'] = 1
        elif (arg[k].upper() == "-DIR"):  # for CGI, a dir to hold the data
            dic['dir'] = arg[k + 1]

        elif (arg[k].upper() == "-VERB"):
            dic['verb'] = 1
        elif (arg[k].upper() == "-AUTO"):
            dic['auto'] = 1
        elif (arg[k].upper() == "-XYZLIM"):
            dic['xyzlim'] = 1

        elif (arg[k].upper() == "-ALL"):  # All xray
            dic['sfcheck'], dic['refmac'], dic['cns'], dic['phenix_x'] = 1, 1, 1, 1
        elif (arg[k].upper() == "-SFCHECK"):
            dic['sfcheck'] = 1
        elif (arg[k].upper() == "-REFMAC"):
            dic['refmac'] = 1
        elif (arg[k].upper() == "-CNS"):
            dic['cns'] = 1
        elif (arg[k].upper() == "-SHELX"):
            dic['shelx'] = 1
        elif (arg[k].upper() == "-PHENIX_X"):
            dic['phenix_x'] = 1
        elif (arg[k].upper() == "-PHENIX_N"):
            dic['phenix_n'] = 1
        elif (arg[k].upper() == "-PHENIX_XN"):
            dic['phenix_xn'] = 1
        elif (arg[k].upper() == "-PHENIX_XTRIAGE"):
            dic['phenix_xtriage'] = 1
        elif (arg[k].upper() == "-NO_XTRIAGE"):  # not do xtriage for fast calc.
            dic['no_xtriage'] = 1

        elif (arg[k].upper() == "-SCALE"):  # use 'type BULK LSSC  ANISO  EXPE'
            dic['scale'] = 1
        elif (arg[k].upper() == "-RESO"):  # use reported resolution
            dic['reso'] = 1
        elif (arg[k].upper() == "-REST"):  # do restraint
            dic['rest'] = 1
        elif (arg[k].upper() == "-WAVE"):  # consider wavelength
            dic['wave'] = 1
        elif (arg[k].upper() == "-ONE"):  # only use the default. no more trials
            dic['one'] = 1

        elif (arg[k].upper() == "-RSR_ALL"):  # do all rsr (main/side chain ..)
            dic['rsr_all'] = 1

        elif (arg[k].upper() == "-ADD_ASA"):  # add ASA to rsr_all
            dic['add_asa'] = 1

        elif (arg[k].upper() == "-ADD_MOTIF"):  # add motif to rsr_all
            dic['add_motif'] = 1

        elif (arg[k].upper() == "-EDSTAT"):  # do edstat from Ian Tickle (diff map)
            dic['edstat'] = 1

        elif (arg[k].upper() == "-EDDIR"):  # direct it to a folder for edstat calculation
            dic['eddir'] = arg[k + 1]

        elif (arg[k].upper() == "-LLDF"):  # same as the wwPDF report
            dic['lldf'] = 1
        elif (arg[k].upper() == "-RSRD"):  # add the zscores calculated by database
            dic['rsrd'] = 1
        elif (arg[k].upper() == "-CC"):  # add to use perfect FC
            dic['cc'] = 1
        elif (arg[k].upper() == "-FEM"):  # feature enhanced map
            dic['fem'] = 1
        elif (arg[k].upper() == "-REFMAC_TLS"):  # use tls correction by refmac
            dic['refmac_tls'] = 1
        elif (arg[k].upper() == "-TWIN"):
            dic['refmactwin'] = 1

        elif (arg[k].upper() == "-LIB"):
            dic['lib'] = arg[k + 1]

        elif (arg[k].upper() == "-MAPSIZE"):
            dic['mapsize'] = 1
        elif (arg[k].upper() == "-PATH"):  # user assigned bin containing the executables
            config.PATH1 = arg[k + 1] + '/bin/'
            config.PATH2 = arg[k + 1] + '/other_bin/'

        elif (arg[k].upper() == "-BFULL"):  # update full B factors
            dic['bfull'] = arg[k + 1]
            dic['sffile_orig'] = arg[k + 2]
            update_bfull(dic)
            sys.exit()

        elif (arg[k].upper() == "-LIGTABLE"):  # for CGI
            # input density_file, pid, mapfile, pdbfile
            dic['dir'] = arg[k + 2]
            ligand.cut_map_around_ligand_peptide(arg[k + 1], dic, arg[k + 3], arg[k + 4])
            sys.exit()

        elif (arg[k].upper() == "-MAPCUT"):  # input mapfile, coordinate, pid
            # pid=(model_compound_chainID_resnumber_alter_insertion) for one compound
            # repeat pid by adding : Ex. '1_ALA_A_23_._.:1_PRO_A_24_._.'
            ligand.map_around_compound(arg[k + 1], arg[k + 2], arg[k + 3])
            # util.delete_file('mapmask_TMP.csh')
            sys.exit()

        elif (arg[k].upper() == "-ASA"):
            _asa = prog.calc_asa_areaimol(arg[k + 1], 1)  # noqa: F841
            sys.exit()

        elif (arg[k].upper() == "-MOTIF"):
            _asa = prog.run_motif(arg[k + 1])  # a pdbfile  # noqa: F841
            sys.exit()

        elif (arg[k].upper() == "-MATT"):
            mat, sol = matt.matthew_coeff(arg[k + 1])
            print('The matt,solv= %.3f  %.3f ' % (mat, sol))
            sys.exit()

        elif (arg[k].upper() == "-NCS"):
            ncs.check_ncs(arg[k + 1])
            sys.exit()

        elif (arg[k].upper() == "-CIF2PDB"):
            pdbfile = cif.cif2pdb(arg[k + 1])
            print('The output file = %s' % pdbfile)
            sys.exit()

        elif (arg[k].upper() == "-CIF2CIF"):  # parse cif and rewrite in new cif
            newfile = cif.cif2cif(arg[k + 1])
            print('The output file = %s' % newfile)
            sys.exit()

        if arg[k][0] == '-' and options(arg[k]) == 0:
            util.perror('Error: The option (%s) is not correct. Please try other option.' % arg[k])
            sys.exit()

    if (not (dic['sfcheck'] or dic['phenix_x'] or dic['phenix_n'] or dic['refmac']
             or dic['phenix_xn'] or dic['cns'] or dic['shelx'] or dic['auto'] or dic['fem'])):
        dic['refmac'] = 1  # default

    if not dic['outfile']:
        dic['outfile'] = "%s_rcc_sum.cif" % os.path.basename(dic['xyzfile_orig'])
# -------------------done with input -----------------


##########################################################
def check_env(dic):
    '''check the environment used for different programs
    '''

    if not config.PATH1:
        util.perror('Error: DCCPY environment is not set. Stopped!')
        print('"setenv DCCPY  ?/sf_valid/" (C shell) OR "export DCCPY=?/sf-valid/" (B shell).')
        sys.exit()

    if (dic['phenix_xtriage'] or dic['phenix_x'] or dic['phenix_n']
            or dic['exp'] == 'n' or dic['exp'] == 'xn' or dic['phenix_xn'] or dic['fem']):
        if 'PHENIX' not in os.environ:
            print("\nWarning: PHENIX failed! To run phenix program, source 'phenix_env'.\n")
        else:
            if 'PHENIX_VERSION' in os.environ and 'PHENIX_RELEASE_TAG' not in os.environ:
                version = '%s' % (os.environ['PHENIX_VERSION'])  # newer version
            else:
                version = '%s-%s' % (os.environ['PHENIX_VERSION'], os.environ['PHENIX_RELEASE_TAG'])

            ss = ['PHENIX', version, config.VALID_CAL]
            if ss not in dic['software']:
                dic['software'].append(ss)
            if dic['no_xtriage'] == 0:
                ss = ['Xtriage', version, config.DATA]
                if ss not in dic['software']:
                    dic['software'].append(ss)

    if (dic['cns'] and 'CNS_SOLVE' not in os.environ):
        print("\nWarning: CNS failed! To run CNS program, source 'cns_solve_env'.\n")

    if (dic['refmac'] and 'CCP4' not in os.environ):
        print("\nWarning: REFMAC failed! To run CCP4 program, source 'ccp4_env'.\n")


##########################################################
def options(x):
    ''' contains all the options:
    '''

    n = 0
    opt = ('-PDB', '-CIF', '-SF', '-PDBID', '-CIFID', '-ALL', '-SFCHECK', '-REFMAC', '-CNS',
           '-PHENIX_X', '-PHENIX_N', '-PHENIX_XN', '-PHENIX_XTRIAGE', '-AUTO',
           '-NO_XTRIAGE', '-RESO', '-REFMAC_TLS', '-TEST_R', '-TWIN', '-PATH',
           '-MAPIN', '-O', '-MAP_DSN6', '-SCALE', '-MAP', '-MAP_CCP4', '-LIGMAP', '-MAPSIZE',
           '-SDSC_MAP', '-LIGMAPCIF', '-ASA', '-MOTIF', '-ADD_ASA', '-ADD_MOTIF', '-RSRD',
           '-REFINE', '-MAPCUT', '-LIGTABLE', '-VERB', '-SHELX' , '-DIR', '-LIB',
           '-MTZMAP', '-DISPLAY', '-REST', '-WAVE', '-ONE', '-BFULL', '-DIAGS',
           '-XYZLIM', '-FOFC', '-2FOFC', '-NOEDS', '-OMITMAP', '-RSR', '-LLDF',
           '-RSR_ALL', '-EDSTAT', '-EDDIR', '-FEM', '-CC', '-OMIT', '-OMIT_RESIDUE', '-BUSTER')

    if x.upper() in opt:
        n = 1
    return n


##########################################################
def check_xyz(dic):
    '''dic contains all the input args. It grows as needed.

    '''

    xyzfile = dic['xyzfile_orig']

    if util.is_cif(xyzfile):
        print('Input coordinate format=mmCIF')
        dic['xyz_type'] = 'cif'
        # check_cif(xyzfile, dic) #check cif and update dic
        tmpfile = cif.cif2pdb(xyzfile)
        dic['pdbfile_orig'] = tmpfile  #
        dic['pdbfile'] = check_pdb(tmpfile, dic)

    else:
        print('Input coordinate format=PDB')
        dic['pdbfile_orig'] = xyzfile
        dic['pdbfile'] = check_pdb(xyzfile, dic)  # add coord if NCS

    if dic['omitmap']:
        dic['pdbfile'] = ligand.remove_ligand(dic['pdbfile'])
    pdbfile = dic['pdbfile']

    tls.check_tls(dic['pdbfile_orig'], 0)

    tlsok, pdb_tls = 0, ''  # run  TLSANL to convert B_partial to B_full
    if dic['tls_attempt'] == 'Y':  # from check_pdb
        tlsok, pdb_tls = pdb_bfull(dic['pdbfile_orig'])
        if tlsok:
            dic['tls_ok'] = 'Y'
            dic['pdb_tls'] = pdb_tls
            if dic['omit'] == 1 and 'omit_residue' in dic:
                dic['pdb_tls'] = omit_residue(dic, pdb_tls)

        else:
            util.perror('Warning: TLS correction failed.')

    check_error(dic['pdbfile_orig'], dic['sffile_orig'], dic)
#    print 'all coord=', dic['pdbfile_orig'], dic['xyzfile_orig'],  dic['pdbfile']
    print('Experiment type = %s' % dic['exp'])

    dic['pdbfile_sig'] = dic['pdbfile_orig'] + '_sig'
    fp = open(dic['pdbfile_orig'], 'r')
    fw = open(dic['pdbfile_sig'], 'w')
    for x in fp:
        if 'ENDMDL' in x[:6]:
            break

        if ((('UNL' in x[17:20] or 'UNK' in x[17:20] or 'UNX' in x[17:20])
             and (' U' in x[76:78] or '  ' in x[76:78])) or ' X' in x[76:78]):
            x1 = x[:12] + " C   " + x[17:76] + " C" + x[78:]
            fw.write('%s \n' % x1.strip())
            continue

        fw.write(x)
    fp.close()
    fw.close()
    return pdbfile


##########################################################
def omit_residue(dic, pdbfile):
    '''omit residues
    '''

    alist = open(pdbfile, 'r').readlines()

    fw = open(pdbfile + '_new', 'w')

    lig = pdbfile + '_lig.pdb'
    fwlig = open(lig, 'w')

    resid = dic['omit_residue'].split('_')
    ch = resid[0]
    res1 = int(resid[1].split(':')[0])
    res2 = int(resid[1].split(':')[1])
    print('omit residue = %s_%d:%d' % (ch, res1, res2))

    for x in alist:
        if x[:6] == 'CRYST1':
            fwlig.write(x)
        if (x[:4] == 'ATOM' or x[:6] == 'HETATM' or 'ANISOU' in x[:6]):
            ch1, res = x[20:22].strip(), int(x[22:26])
            if ch == ch1 and res1 <= res <= res2:
                fwlig.write(x)
                continue

            fw.write(x)
        else:
            fw.write(x)

    fw.close()
    fwlig.close()

    return pdbfile + '_new'


##########################################################
def do_validation_auto(dic, pdb_new, sf_new):
    '''validate the entry by the program for refinement.
    If failed, pick a different one (refmac ->phenix_x ->cns ->sfcheck)
    '''

    print('\nRuning DCC in the auto-mode...')
    rcut = 0.02
    if dic['model'] == 'Y':
        rcut = 0.05

    dic['refmac'] = 0  # refmac is by default, must be removed.
    prog = dic['prog'].upper()
    pp = ''
    if 'PHENIX' in prog:
        dic['phenix_x'] = 1
        pp = 'phenix_x'
    elif 'PLOR' in prog or 'CNS' in prog:
        dic['cns'] = 1
        pp = 'cns'
    else:
        dic['refmac'] = 1
        pp = 'refmac'

    dd_all = do_validation(dic, pdb_new, sf_new)  # by refinement prog
    if 'n' in dic['exp'] or dic['phenix_n'] or dic['phenix_xn']:
        return dd_all  # neutron diffr

    nn = best_of_solutions(dd_all)
    dd = dd_all[nn]

    print('Rwork/Rfree (reported=%s/%s by %s  : calc=%s/%s by %s)\n'
          % (dic['rfact'], dic['rfree'], prog, dd['rfact'], dd['rfree'], pp))

    cond_cut = (util.is_number(dic['rfact']) and util.is_number(dd['rfact'])
                and (float(dd['rfact']) - float(dic['rfact']) <= rcut))
    if cond_cut:
        return dd_all

    # prog_test=['refmac', 'phenix_x', 'cns', 'sfcheck']
    prog_test = ['phenix_x', 'refmac']

    dic[pp] = 0  # exclude the refinement program
    dd_all_new = []
    if '?' not in dd['rfact']:
        dd_all_new.append([dd, dd_all])
    for x in prog_test:
        if x == pp:
            continue
        dic[x] = 1
        print('\nTesting program  %s ...' % x)

        dd_all_in = do_validation(dic, pdb_new, sf_new)
        nn = best_of_solutions(dd_all_in)
        dd_in = dd_all_in[nn]
        dic[x] = 0  # remove after calculations

        cond_cut = (util.is_number(dic['rfact']) and util.is_number(dd_in['rfact'])
                    and (float(dd_in['rfact']) - float(dic['rfact']) <= rcut))
        if cond_cut or '?' in dic['rfact']:
            return dd_all_in
        elif '?' not in dd_in['rfact']:
            dd_all_new.append([dd_in, dd_all_in])

    if len(dd_all_new) == 0:
        return dd_all  # more test failed

    rfact = []
    for x in dd_all_new:  # pick the best from the test if non passed criteria
        # print x[0]['rfact']
        rfact.append(float(x[0]['rfact']))
    rmin = min(rfact)
    n = 0
    for i, x in enumerate(dd_all_new):
        if float(x[0]['rfact']) == rmin:
            n = i
            break

    return dd_all_new[n][1]


##########################################################
def do_validation(dic, pdb_new, sf_new):
    ''' Do various validations and write all the results in a cif file
        dic: contains all input informatin, pdb_new/sf_new is working PDB/SF format;
    '''

    pdbfile = dic['pdbfile_orig']  # orginal
    pdb_tls = dic['pdb_tls']  # TLS corrected file (if partial)

    tlsok = 0
    if dic['tls_ok'] == 'Y':
        tlsok = 1

    #  run sfcheck
    dic1 = initialize_dic1()
    dic2 = initialize_dic1()
    cifsfch1, cifsfch2 = '', ''
    if (dic['sfcheck']):
        (dic1, cifsfch1) = prog.run_sfcheck(pdb_new, sf_new, dic, 'NOTLS')
        util.get_software_version('sfcheck.log', 'SFCHECK', dic)
        if tlsok:
            (dic2, cifsfch2) = prog.run_sfcheck(pdb_tls, sf_new, dic, 'TLS')
        nn = best_of_solutions([dic1, dic2])
        if tlsok:
            tell_btype(nn, dic)
        dic['cifsfch'] = cifsfch1
        if nn == 1:
            dic['cifsfch'] = cifsfch2

    #  run refmac : multiple tries (non-large file) if needed.
    dic3 = initialize_dic1()
    dic4 = initialize_dic1()
    dic5 = initialize_dic1()

    reso_in = '  '
    mtzout = dic['outfile'] + '.mtz'
    if (dic['refmac']):
        mtz = prog.sf_convertor(sf_new, pdb_new, 'mtz')
        prog.run_pointless(mtz, dic)  # L test

        dic3['detail'] = 'without TLS correction'
        mtz3, mtz4, mtz5, pdb3, pdb4, pdb5, log3, log4, log5 = '', '', '', '', '', '', '', '', ''
        mtz3, pdb3, log3, scr = prog.run_refmac(pdb_new, mtz, dic, 1, reso_in)
        parse.refmac_log(log3, dic3)  # update  dic3 from log3

        print('Rwork/Rfree: reported=%s/%s, calculated=%s/%s, fom=%s '
              % (dic['rfact'], dic['rfree'], dic3['rfact'], dic3['rfree'], dic3['fom']))

        if tlsok:  # tls exist
            dic4['detail'] = 'with    TLS correction'
            mtz4, pdb4, log4, scr = prog.run_refmac(pdb_tls, mtz, dic, 2, reso_in)
            parse.refmac_log(log4, dic4)

        nn = best_of_solutions([dic3, dic4])

        if tlsok:
            tell_btype(nn, dic)
        elif 'PARTIAL' in dic['isob'] :
            util.perror('Warning: Input file shows PARTIAL B. But TLS failed, double check B type.')

        if not (dic['mtr'] > 1 or dic['model'] == 'Y' or len(dic['split1']) > 1 or len(dic['split2']) > 1 or dic['one'] == 1):
            if nn == 0:  # try scale, reso, wave, twin options
                mtz5, pdb5, log5, scr, dic5 = prog.run_refmac_more(pdb_new, mtz, dic, dic3)
            else:
                mtz5, pdb5, log5, scr, dic5 = prog.run_refmac_more(pdb_tls, mtz, dic, dic4)

            parse.refmac_log(log5, dic5)

        dics = [dic3, dic4, dic5]
        mtzs = [mtz3, mtz4, mtz5]
        pdbs = [pdb3, pdb4, pdb5]
        logs = [log3, log4, log5]

        nn = best_of_solutions(dics)
        mtzo = mtzs[nn]

        edsout = ''
        if dic['rsr_all'] > 0:  # all, main/side chain, update cifrefm in dic
            prog.run_edsall(mtzo, dic['pdbfile_sig'], dic, dics[nn])
        elif dic['edstat'] > 0:  # Tickle for diff map
            prog.run_edstat(mtzo, dic['pdbfile_sig'], dic, dics[nn])
            map_2fofc = dic['outfile'] + '-2FoFc.map'
            if dic['map']:
                prog.mtz2map(mtzo, pdbfile, map_2fofc, '2FO_FC', dic)
                if not dic['verb']:
                    util.delete_file('get_mtz2map.csh')
        else:
            if not dic['cc']:
                dic['cifrefm'], edsout = prog.run_eds(mtzo, dic['pdbfile_sig'], dic, dics[nn])  # use new
            else:
                mtzn = merge_fc2fwt(mtzo, pdb_new)
                dic['cifrefm'], edsout = prog.run_eds(mtzn, dic['pdbfile_sig'], dic, dics[nn])  # use new

        shutil.move(mtzo, mtzout)

        dic['k_sol'] = dics[nn]['k_sol']
        dic['b_sol'] = dics[nn]['b_sol']
        dic['twin_op'] = dics[nn]['twin_op']
        dic['twin_fr'] = dics[nn]['twin_fr']
        dic['fom'] = dics[nn]['fom']
        dic['dpi_fr'] = dics[nn]['dpi_fr']
        dic['dpi_xyz'] = dics[nn]['dpi_xyz']
        dic['realr'] = dics[nn]['realr']
        dic['dcorr'] = dics[nn]['dcorr']

        for i in range(len(mtzs)):
            if i == nn and dic['verb']:
                continue
            util.delete_file(mtzs[i], pdbs[i], logs[i])
        if not dic['verb']:
            util.delete_file(edsout)

    # run ligand map cut (always position it after runing refmac!)
    map2fofc = pdbfile + '.map'

    if dic['ligmap'] or dic['ligmapcif']:
        prog.mtz2map(mtzout, pdbfile, map2fofc, '2FO_FC', dic)
        ligand.cut_map_around_ligand_peptide(dic['cifrefm'], dic, map2fofc, pdbfile)
        if not dic['verb']:
            util.delete_file(map2fofc)

    elif dic['sdsc_map']:
        prog.mtz2map(mtzout, pdbfile, map2fofc, '2FO_FC', dic)
        sdsc_ligmap.sdsc_ligand(dic, map2fofc, pdbfile)

    #  run phenix (model_vs_data/xtriage)

    dic6 = initialize_dic1()
    dic7 = initialize_dic1()
    dic8 = initialize_dic1()
    dic9 = initialize_dic1()

    if not dic['no_xtriage'] and not dic['fem'] and 'PHENIX' in os.environ:  # always calc
        _phlog1 = prog.run_phenix(dic, dic, 'xtriage')  # noqa: F841

    if (dic['phenix_x']):
        _phlog2 = prog.run_phenix(dic, dic6, 'xray')  # noqa: F841
        dic6['detail'] = 'Calculated Xray'
    if (dic['phenix_n']):
        _phlog3 = prog.run_phenix(dic, dic7, 'neut')  # noqa: F841
        dic7['detail'] = 'Calculated Neut'
    if (dic['phenix_xn'] or dic['exp'] == 'xn' or dic['exp'] == 'nx'):
        _phlog4 = prog.run_phenix(dic, dic8, 'xray')  # noqa: F841
        dic8['detail'] = 'Calculated Xray'
        _phlog5 = prog.run_phenix(dic, dic9, 'neut')  # noqa: F841
        dic9['detail'] = 'Calculated Neut'

    # run cns program (model_stat)

    dic10 = initialize_dic1()
    dic11 = initialize_dic1()

    if dic['cns'] and 'CNS_SOLVE' in os.environ:
        cns_sf = prog.sf_convertor(sf_new, pdb_new, 'cns')

        cnsout1, cnslog1 = prog.run_cns(dic, pdb_new, cns_sf)
        parse.model_stat(cnsout1, cnslog1, dic10)
        util.get_software_version(cnslog1, 'CNS', dic)
        dic10['detail'] = 'without TLS correction'

        if tlsok:
            (cnsout2, cnslog2) = prog.run_cns(dic, pdb_tls, cns_sf)
            parse.model_stat(cnsout2, cnslog2, dic11)
            dic11['detail'] = 'with    TLS correction'

        nn = best_of_solutions([dic10, dic11])
        if tlsok:
            tell_btype(nn, dic)
        if (dic['verb'] == 0):
            util.delete_file(' model_stats.list model_stats.log')

    # run shelx program (LS refinement)

    dic12 = initialize_dic1()
    if dic['shelx']:  # use original SF !!
        shelx_sf = prog.sf_convertor(sf_new, pdb_new, 'shelx')
        shelx_log = prog.run_shelx(pdb_new, shelx_sf)
        util.get_software_version(shelx_log, 'SHELX', dic)
        parse.log_from_shelx(shelx_log, dic12)
        dic12['detail'] = 'without TLS correction'

    dic13 = initialize_dic1()
    if dic['fem']:  # calc RsR/CC by the FEM
        fem_mtz = prog.calc_fem(pdb_new, sf_new)
        fc_mtz = prog.calc_fc(pdb_new, sf_new)
        sf_mtz = 'fem_fc_all.mtz'
        prog.merge_mtz(fem_mtz, fc_mtz, sf_mtz)
        dic['cifrefm'], edsout = prog.run_eds_fem(sf_mtz, pdb_new, dic, dic13)  # use new
        mapout = pdbfile + '_fem_2fofc.map'
        prog.mtz2map(fem_mtz, pdbfile, mapout, 'FEM', dic)

    all_dic = [dic1, dic2, dic3, dic4, dic5, dic6, dic7, dic8, dic9, dic10, dic11, dic12, dic13]
    return all_dic


##########################################################
def write_all(dic, all_dic):
    # ---------  write the final data

    outfile = dic['outfile']

    fw = open(dic['outfile'], 'w')
    fw.write('data_' + dic['pdbid'] + '\n#\n')
    write_cif.dcc_residual(fw, dic['cifsfch'])

    if util.check_file(300, dic['cifrefm']):  # write RSR
        if dic['rsr_all']:
            write_cif.rsr_edsall(fw, dic)
        elif dic['edstat']:
            write_cif.rsr_edstat(fw, dic)
        else:
            write_cif.dcc_residual(fw, dic['cifrefm'])

    nn = best_of_solutions(all_dic)
    check_result(dic, all_dic[nn])

    write_cif.summary1(fw, dic)  # non-looped summary
    if dic['phenix_xn'] or dic['phenix_n'] or dic['phenix_x']:
        vers = util.get_phenix_version_from_dic(dic)
        if vers < 1013:
            write_cif.summary2(fw, dic)  # non-looped summary (geometry)

    write_cif.write_software(fw, dic)
    write_cif.write_dcc_sf(fw, dic)

    vers = dic['software']
    write_cif.final_cif_item(fw, {}, 0, vers)  # looped summary
    np = 0

    if '?' in dic['prog']:
        dic['prog'] = 'AUTH'  # the report program
    if not dic['phenix_xn'] and dic['exp'] != 'xn' and dic['exp'] != 'nx':
        write_cif.final_cif_item(fw, dic, 1, vers)
        np = 1
    else:
        version = '?'
        for x in vers:
            if 'PHENIX' in x[0]:
                version = x[1].replace(' ', '')
                break

        fw.write("1 %8s %6s %6s %6s %6s %6s %6s %6s %6s ?  ?  ?  %8s 'Reported Xray'\n"
                 % (dic['prog'], dic['resh'], dic['resl'], dic['rall'],
                    dic['rfact'], dic['rfree'], dic['nref'], dic['comp'],
                    dic['nfree'], version))

        fw.write("2 %8s %6s %6s %6s %6s %6s %6s %6s %6s ?  ?  ?  %8s 'Reported Neut'\n"
                 % (dic['prog'], dic['resh_n'], dic['resl_n'], dic['rall_n'],
                    dic['rfact_n'], dic['rfree_n'], dic['nref_n'], dic['comp_n'],
                    dic['nfree_n'], version))
        np = 2

    for i, dd in enumerate(all_dic):
        if '?' in dd['rfact']:
            continue
        np = np + 1
        write_cif.final_cif_item(fw, dd, np, vers)

    all_dic[nn]['detail'] = 'Best solution'
    write_cif.final_cif_item(fw, all_dic[nn], np + 1, vers)

    fw.close()

    if dic['rsrd'] == 0:
        add_modelid_to_outfile(dic, outfile)
    if dic['lldf']:
        add_lldf_to_outfile(dic, outfile)

    newfile = cif.cif2cif(outfile)
    shutil.move(newfile, outfile)
    print("\nThe final mmCIF file = %s\n" % outfile)

    if dic['edstat'] > 0:
        util.delete_file('NEW_sugar.pdb* NEW_side.pdb* NEW_phos.pdb*')

    if (dic['verb'] == 0):  # clean files
        util.delete_file('SF_4_validate.cif sf_information.txt  sf_format_guess.text ')
        util.delete_file('sfcheck.xml  SF_4_validate.cif_mtz ')
        util.delete_file('extract_map_sfcheck.mmcif_NOTLS  extract_map_sfcheck.mmcif_TLS ')
        util.delete_file('mapcut_TMP.csh  mapmask_TMP.csh get_mtzmap.csh  mtz2eds.csh run_refmac.csh ')
        if dic['omitmap']:
            util.delete_file(dic['pdb_nolig'], 'EDS_OMIT.OUT.cif')


##########################################################
def best_of_solutions(dd):
    ''' pick best of the solution from the list of dictionary dd
    return position of the best solution

    '''

    n = 0  #

    if (len(dd) == 2):
        if (util.is_number(dd[0]['rfact']) and util.is_number(dd[1]['rfact'])
                and float(dd[0]['rfact']) > float(dd[1]['rfact']) + 0.001):
            n = 1
        elif '?' in dd[0]['rfact'] and util.is_number(dd[1]['rfact']):
            n = 1

        return n

    d1 = []
    for x in dd:
        if not util.is_number(x['rfact']):
            continue
        d1.append(float(x['rfact']))

    if len(d1) == 0:
        return n
    min_r = min(d1)
    # max_r = max(d1)

    for i, x in enumerate(dd):
        if not util.is_number(x['rfact']):
            continue
        if min_r == float(x['rfact']):
            n = i
            break

    return n


##########################################################
def add_data2_dic(dic6, dic):
    '''copy some data from dic6 to dic for the summary
    '''

    dic['outlier'] = dic6['outlier']
    dic['allowed'] = dic6['allowed']
    dic['favored'] = dic6['favored']
    dic['rot_outlier'] = dic6['rot_outlier']
    dic['cbeta'] = dic6['cbeta']
    dic['clashscore'] = dic6['clashscore']
    dic['oscore'] = dic6['oscore']

    dic['bond'] = dic6['bond']
    dic['angle'] = dic6['angle']
    dic['dihedral'] = dic6['dihedral']
    dic['chirality'] = dic6['chirality']
    dic['planarity'] = dic6['planarity']
    dic['non_bond'] = dic6['non_bond']
    # dic['fofc_3sig'] = dic6['fofc_3sig']

    return


##########################################################
def merge_fc2fwt(mtz_in, pdbfile):
    '''get new mtz file that contains the perfect FC and FWT from refmac
    '''

    pdbnew = add_const_to_biso(pdbfile, 0, id=0)  # make Biso=0
    mtz_fc = prog.calc_fc(pdbnew, mtz_in)

    mtz1_lab = ['FC', 'PHIC', 'FP', 'SIGFP']  # from pdbfile
    mtz2_lab = ['FWT', 'PHWT', 'DELFWT', 'PHDELWT']  # from mtz_in

    mtzo = '%s_FC.mtz' % mtz_in

    prog.cat_mtz(mtz_fc, mtz_in, mtz1_lab, mtz2_lab, mtzo)

    return mtzo


##########################################################
def tell_btype(nn, dic):
    '''
    '''
    if nn == 1:
        dic['isob'] = 'PARTIAL'
    else:
        if 'PARTIAL' in dic['isob']:
            util.perror('Warning: PDB header shows PARTIAL B, the calculated is FULL B, double check.')
        dic['isob'] = 'FULL'


##########################################################
def get_segid_cif(xyzf):
    '''xyzf is a cif file
    '''

    flist = open(xyzf, 'r').readlines()

    items, values = cif.cifparse(flist, '_atom_site.')  # a loop
    asym = cif.parse_values(items, values, "_atom_site.auth_asym_id")
    comp = cif.parse_values(items, values , "_atom_site.auth_comp_id")
    seq = cif.parse_values(items, values, "_atom_site.auth_seq_id")
    ins = cif.parse_values(items, values, "_atom_site.pdbx_PDB_ins_code")
    if not ins:
        ins = cif.parse_values(items, values, "_atom_site.ndb_ins_code")
    alt = cif.parse_values(items, values, "_atom_site.label_alt_id")

    nlen = len(asym)
    id0 = ' '
    segid = []
    for i in range(nlen):
        inst = ins[i]
        alter = alt[i]
        if ins[i] == '?':
            inst = '.'
        if alt[i] == '?':
            alter = '.'

        id = '_'.join([asym[i][0], comp[i], seq[i], inst, alter])
        if id != id0:
            segid.append([asym[i], id])
            id0 = id

    return segid


##########################################################
def get_segid_pdb(xyzf):
    '''xyzf is a pdb file, must have segid at column 72-73
    '''

    fp = open(xyzf, 'r')

    id0 = ' '
    segid = []
    for x in fp:
        seg = x[72:74].strip()
        if 'ATOM' in x[:4] or 'HETA' in x[:4] and seg:
            comp = x[17:20].strip()
            ch = x[20:22].strip()
            seq = x[22:26].strip()
            ins = '.'
            tmp = x[26:27].strip()
            if tmp:
                ins = tmp
            alt = '.'
            tmp = x[16:17].strip()
            if tmp:
                alt = tmp

            id = '_'.join([ch, comp, seq, ins, alt])
            if id != id0:
                segid.append([seg, id])
                id0 = id

    fp.close()
    return segid


##########################################################
def add_lldf_to_outfile(dic, ofile):
    '''add local ligand density fit (LLDF) to the final cif file.
    '''

    print('Doing LLDF ..')
    flist = open(ofile, 'r').readlines()
    dcc = parse.parse_dcc(flist, '_pdbx_rscc_mapman.')
    # dcc is a list = nseq,chid,comp,alt,cc,rsr,biso,occ,pdbid,zrsr,ins
    #                 0    1    2    3   4  5  6    7    8    9    10

    fw = open(ofile, 'w')

    xyzfile = dic['pdbfile']
    ch_pep, chr_pep, ch_lig, chr_lig, ch_wat, chr_wat = tls.chain_res_range(xyzfile)
    #    print 'ch_lig,chr_lig=', ch_lig,chr_lig

    all_resid = []
    all_lldf = []
    for x in dcc:

        nseq, chid = int(x[0]), x[1]

        if chid in ch_lig.keys() and nseq in ch_lig[chid]:  # found ligand
            ligid = [x[0], x[1], x[2]]
            logfile = prog.get_contact(chid, nseq, xyzfile, 5)
            cont_resid = get_contact_resid(logfile)
            all_resid.append(ligid)
            all_resid.extend(cont_resid)
            lldf = ligand_lldf(dcc, all_resid)
            x.append('%.2f' % lldf)
            all_lldf.append(x)
            # print 'all_resid=', lldf,  all_resid
            all_resid = []

    for i, x in enumerate(flist):  # writing
        fw.write(x)
        if ('#' in x[:10] and 'Overall properties from' in flist[i + 1]
                and len(all_lldf) > 0):
            cifhead2 = write_cif.dcc_ciftoken('pdbx_non_poly', dic)
            fw.write(cifhead2)

            for y in all_lldf:
                t = "%4s %s %s %4s %s %s %s %s %s %s %s\n" % (y[8], y[1], y[2], y[0], y[3], y[4],
                                                              y[5], y[9], y[11], y[6], y[7])
                fw.write(t)
            fw.write('#\n')

    fw.close()


##########################################################
def ligand_lldf(dcc, all_resid):
    '''get the LLDF
    '''
    res = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'CYS', 'MET', 'PRO', 'PHE', 'TRP', 'TYR',
           'HIS', 'LYS', 'ARG', 'ASP', 'GLU', 'ASN', 'GLN', 'THR', 'SER', 'MSE',
           'A', 'G', 'C', 'T', 'DA', 'DG', 'DC', 'DT', 'DI', 'U', 'I']

    rsr = []
    # lignam = '%s_%s%s' % (all_resid[0][2], all_resid[0][1], all_resid[0][0])
    for x in dcc:  # the first one
        if x[0] == all_resid[0][0] and x[1] == all_resid[0][1]:
            rsr.append([float(x[5])])
            break

    for y in all_resid[1:]:

        if y[2] not in res:
            continue
        for x in dcc:
            if x[0] == y[0] and x[1] == y[1]:
                rsr.append([float(x[5])])
                # print y,x
                break

    avg, dev, mini, maxi = util.mean_dev(rsr[1:], 0)  # remove ligand
    # ncount = len(rsr[1:])
    if dev > 0.0:
        lldf = (rsr[0][0] - avg) / dev
    else:
        lldf = 999
        util.perror('Warning: sigma=0 around the ligand. No LLDF calculated!')

        # print '%6s:: Ncount=%d, avg=%.3f, dev=%.3f, mini=%.3f, maxi=%.3f, LLDF=%.2f' \
        #     %(lignam, ncount, avg, dev, mini, maxi, lldf)

    return lldf


##########################################################
def get_contact_resid(logfile):
    '''get residues which contact with ligand
    '''

    data = []
    if not util.check_file(100, logfile):
        return data

    fp = open(logfile, 'r')
    for x in fp:
        if len(x) > 56 and '.' in x and '[' in x and ']' in x and ':' in x:
            y = x[27:44]
            if y not in data:
                data.append(y)

    ndata = []
    for x in data:
        t = x.split('/')
        chid = t[2].strip()
        nres = t[3].split('(')[0].strip()
        comp = t[3].split('(')[1].split(')')[0].strip()
        if 'HOH' in comp:
            continue
        ndata.append([nres, chid, comp])

    fp.close()
    return ndata


##########################################################
def add_modelid_to_outfile(dic, ofile):
    '''add model id to the final cif file. Skip if only one model
    '''

    # if dic['model'] != 'Y' : return

    fp = open(ofile, 'r').readlines()

    firstl = ''
    for i, x in enumerate(fp):
        if (i > 0 and (('rscc_mapman.' in fp[i - 1] and '_mapman.' not in x)
                       or ('_sfcheck.' in fp[i - 1] and '_sfcheck.' not in x))):
            firstl = x
            break
    if not firstl:
        return

    t = firstl.split()
    id_mapman = [t[1], t[2], t[3], t[4], t[5]]
    id_sfcheck = [t[0], t[1], t[2], t[3]]
    fw = open(ofile, 'w')  # re-write the output with same name!
    n1, m, k1, k2, nmod = 0, 0, 0, 0, 0
    for x in fp:
        n1 = n1 + 1
        t = x.split()
        ss = x
        if k1 == 0 and ('_pdbx_rscc_mapman.' in x):
            fw.write(' _pdbx_rscc_mapman.model_id\n')
            k1 = k1 + 1
        elif (k2 == 0 and '_pdbx_rscc_sfcheck.' in x):
            fw.write(' _pdbx_rscc_sfcheck.model_id\n')
            k2 = k2 + 1
        elif len(t) > 7 and (k1 > 0 or k2 > 0):
            if k1 > 0:
                id = [t[1], t[2], t[3], t[4], t[5]]
                if id == id_mapman:
                    nmod = nmod + 1
            elif k2 > 0:
                id = [t[0], t[1], t[2], t[3]]
                if id == id_sfcheck:
                    nmod = nmod + 1
            else:
                nmod = 1

            ss = '%-d ' % nmod + x
            fw.write(ss)
            m = m + 1
            continue
        elif m > 1 and '#' in x:
            break

        fw.write(ss)

    for ln in fp[n1 - 1:]:
        fw.write(ln)

    fw.close()


##########################################################
def add_segid_to_outfile(xyzf, ofile):
    '''add pdbid for the split entries (pdb/cif input)
    sfcheck/refmac output file
    '''

    if util.is_cif(xyzf):
        segid = get_segid_cif(xyzf)
    else:
        segid = get_segid_pdb(xyzf)

    # for x in segid: print(x)

    fp = open(ofile, 'r').readlines()

    dcc, sfc, n = [], 0, 0
    for x in fp:
        t = x.split()
        nlen = len(t)
        if '_pdbx_rscc_mapman.' in x:
            n = n + 1
        elif '_pdbx_rscc_sfcheck.' in x:
            n = n + 1
            sfc = 1
        elif n > 0 and ('#' in x or nlen < 6):
            break
        elif nlen > 6 and n > 0:
            dcc.append(x)

    if sfc:  # rid of waters
        t = [x for x in segid if 'HOH' not in x[1]]
        segid = t

    len_seq, len_dcc = len(segid), len(dcc)
    if not (len_seq == len_dcc and len_seq > 0):
        util.perror('Warning: may have problem to map chainID to each residue(pdb=%d; dcc=%d).' % (len_seq, len_dcc))

    fw = open(ofile, 'w')  # re-write the output with same name!
    n1, m, k1, k2 = 0, 0, 0, 0
    for x in fp:
        n1 = n1 + 1
        t = x.split()
        ss = x
        if k1 == 0 and ('_pdbx_rscc_mapman.' in x):
            fw.write('_pdbx_rscc_mapman.segid\n')
            k1 = k1 + 1
        elif (k2 == 0 and '_pdbx_rscc_sfcheck.' in x):
            fw.write('_pdbx_rscc_sfcheck.segid\n')
            k2 = k2 + 1
        elif len(t) > 7 and (k1 > 0 or k2 > 0):
            s1 = ' '.join(t)
            ss = '  '.join([segid[m][0], s1]) + '\n'
            fw.write(ss)
            m = m + 1
            continue
        elif m > 10 and '#' in x:
            break

        fw.write(ss)

    for ln in fp[n1 - 1:]:
        fw.write(ln)

    fw.close()


##########################################################
def add_pdbid_to_outfile(dic, file):
    '''add pdbid for the split entries (input is multiple PDB files).
    only for input/output pdbfile
    '''

    if not util.check_file(file):
        return

    pdbid, resid = [], {}
    for x in dic['split1']:
        pdbid.append(x[0])
        resid[x[0]] = [x[1], x[2]]
    # print('pdbid=', pdbid,resid)
    fp = open(file, 'r')
    lis = fp.readlines()
    fo = open(file, 'w')

    n1, j, n, m = 0, 0, 0, 0
    end = resid[pdbid[j]][1]
    for ln in lis:
        n1 = n1 + 1
        tmp = ln.split()
        if len(tmp) == 0:
            continue
        s = ln
        if n == 0 and '_pdbx_rscc_mapman.' in tmp[0]:
            fo.write(' _pdbx_rscc_mapman.pdbid\n')
            n = n + 1
        elif len(tmp) > 7 and n > 0:
            m = m + 1
            s = pdbid[j] + ' ' + ln

            if (tmp[4] == end[2].strip() and tmp[2] == end[1].strip()
                    and tmp[3] == end[0].strip()):
                j = j + 1
                if j >= len(pdbid):
                    j = len(pdbid) - 1
                end = resid[pdbid[j]][1]
        elif m > 10 and '#' in ln:
            break

        fo.write(s)
    for ln in lis[n1 - 1:]:
        fo.write(ln)

    fo.close()
    fp.close()


##########################################################
def check_sf(dic):
    ''' check sf & update dic!
    In addition to the original input SF file dic['sffile_orig'], three more
    SF files are generated here.
    1. sf_valid (add F for refmac, if missing)
    2. sf_xray : first data block (could be the same as sf_valid if F exist.
    3. sf_neut : second block for neutron.

    '''
    sff = '%s/sf-convert/include/sf_convert.h' % os.environ['DCCPY']
    util.get_software_version(sff, 'SF_CONVERT', dic)

    pdbfile = dic['pdbfile_orig']
    sffile = dic['sffile_orig']

    sf_valid = 'SF_4_validate.cif'
    sf_stat = 'sf_information.cif'
    if os.path.exists(sf_valid):
        util.delete_file(sf_valid)
    if os.path.exists(sf_stat):
        util.delete_file(sf_stat)
    sf_new = prog.sf_convertor(sffile, pdbfile, 'mmcif')
    if not util.check_file(500, sf_valid):
        util.perror('Error: sf conversion is not successful. Check SF file (%s)' % (dic['sffile_orig']))
        sys.exit()

    dic['sf_xray'] = sffile + '_xray'
    dic['sf_neut'] = sffile + '_neut'

    fp = open(sf_new, 'r')

    fw1 = open(dic['sf_xray'], 'w')  # get first block
    n1 = 0
    for line in fp:
        if (' x ' in line or ' h ' in line or ' l ' in line
                or ' < ' in line or ' - ' in line):
            continue
        if '_refln.pdbx_r_free_flag' in line:
            line = '_refln.pdbx_r_free_flagx\n'
        if '_refln.pdbx_F_calc_with_solvent' in line:
            line = '_refln.fake2\n'
        if '_refln.pdbx_phase_calc_with_solvent' in line:
            line = '_refln.fake3\n'
        n1 = n1 + 1
        if n1 > 200 and "data_" in line:
            break
        fw1.write(line)
    fw1.close()

    if dic['exp'] == 'xn':  # X-ray & Neutron diffr.
        fw2 = open(dic['sf_neut'], 'w')
        fw2.write("data_neut\n")
        n2 = 0
        for line in fp:
            if (' x ' in line or ' h ' in line or ' l ' in line
                    or ' < ' in line or ' - ' in line):
                continue
            if '_refln.pdbx_r_free_flag' in line:
                line = '_refln.pdbx_r_free_flagx\n'
            if '_refln.pdbx_F_calc_with_solvent' in line:
                line = '_refln.fake2\n'
            if '_refln.pdbx_phase_calc_with_solvent' in line:
                line = '_refln.fake3\n'
            n2 = n2 + 1
            fw2.write(line)
        fw2.close()

        if (n2 < 100):
            util.perror('Error: Only one data set in SF file.\n \
      Please put Xray data set first and Neutron data set second.\n')

    get_sf_stat(dic, sf_stat)

    dic['sffile'] = sf_valid

    fp = open(dic['sf_xray'], 'r')
    n = 0
    for x in fp:
        n = n + 1
        if '_refln.F_meas_au' in x:
            dic['sf_f'] = 1
        elif '_refln.F_meas_sigma_au' in x:
            dic['sf_sigf'] = 1
        elif '_refln.intensity_meas' in x:
            dic['sf_i'] = 1
        elif '_refln.intensity_sigma' in x:
            dic['sf_sigi'] = 1
        elif '_refln.pdbx_F_plus' in x:
            dic['anom_f'] = 1
        elif '_refln.pdbx_I_plus' in x:
            dic['anom_i'] = 1
        elif n > 400:
            break

    fp.close()

    if not dic['verb']:
        util.delete_file(sf_new)
    return sf_valid


##########################################################
def get_sf_stat(dic, sf_stat):
    '''parse sf statistics
    '''
    #    print 'Extracting SF statistics from %s' %sf_stat

    dic['status_used'] = 'N'
    if not util.check_file(100, sf_stat):
        return ''

    fp = open(sf_stat, 'r').readlines()
    for x in fp:
        if '_sf_convert.error' in x:
            continue
        if 'warn' in x.lower() or 'error' in x.lower():
            if 'Warning: No wavelength value was found in SF file.' in x:
                util.perror('Warning: No wavelength value was found in SF file.')
            else:
                util.perror(x)

    for x in fp:
        if 'Use data with resolution' in x:
            break
        if ('Number of observed' in x or "reflections (status='o')" in x):
            nref = x.split('=')[-1].strip()
            dic['sf_nwork'] = nref

        elif ('Number for free set ' in x or "reflections (status='f')" in x):
            nfree = x.split('=')[-1].strip()
            dic['sf_nfree'] = nfree

            if util.is_number(nfree) and int(nfree) > 1:
                dic['status_used'] = 'Y'
                dic['status'] = 'Y'

        elif 'Total number of observed reflections = ' in x:
            dic['sf_nobs'] = x.split('=')[-1].strip()

        elif 'Percentage for free set ' in x:
            dic['sf_pfree'] = x.split('=')[-1].strip()
        elif 'Total number of Friedel pairs (F+/F-)' in x :
            dic['sf_nfpair'] = x.split('=')[-1].strip()
        elif 'Total number of observed F+' in x:
            dic['sf_nfpair_p'] = x.split('=')[-1].strip()
        elif 'Total number of observed F-' in x:
            dic['sf_nfpair_m'] = x.split('=')[-1].strip()
        elif 'Sum of observed F+ and F-' in x:
            dic['sf_nfpair_pm'] = x.split('=')[-1].strip()
        elif 'Total number of Friedel pairs (I+/I-)' in x :
            dic['sf_nipair'] = x.split('=')[-1].strip()
        elif 'Total number of observed I+' in x:
            dic['sf_nipair_p'] = x.split('=')[-1].strip()
        elif 'Total number of observed I-' in x:
            dic['sf_nipair_m'] = x.split('=')[-1].strip()
        elif 'Sum of observed I+ and I-' in x:
            dic['sf_nipair_pm'] = x.split('=')[-1].strip()

        elif 'Lowest resolution' in x:
            dic['resl_sf'] = x.split('=')[1].split(';')[0].strip()
        elif 'Highest resolution' in x:
            dic['resh_sf'] = x.split('=')[1].split(';')[0].strip()
        elif 'Number of reflections for validation' in x:
            break


##########################################################
def separate_sf(sffile):

    fp = open(sffile, 'r')
    sf_xray = sffile + '_xray'
    sf_neut = sffile + '_neut'

    fw1 = open(sf_xray, 'w')
    fw2 = open(sf_neut, 'w')

    n1 = 0
    for ln in fp:
        n1 = n1 + 1
        if n1 > 200 and "data_" in ln:
            break
        fw1.write(ln)

    fw2.write("data_sf2\n")
    n2 = 0
    for ln in fp:
        n2 = n2 + 1
        fw2.write(ln)

    fw1.close()
    fw2.close()
    fp.close()

    if (n1 > 100 and n2 < 100):
        util.perror("Error: Only one data set in SF file.\n \
      Please put Xray data set first and Neutron data set second.\n")

    return (sf_xray, sf_neut)


##########################################################
def change_cif_tag(line):
    ''' Change the cif tag after TLS correction '''

    # for sfcheck :  extract_map_sfcheck.mmcif
    if '_pdbx_rscc_sfcheck.' in line:
        tmp = line.split('sfcheck.')
        line = "%ssfcheck_TLS_correct.%s\n" % (tmp[0], tmp[1].strip())

    elif '_pdbx_rscc_sfcheck_overall.' in line:
        tmp = line.split('overall.')
        line = "%soverall_TLS_correct.%s\n" % (tmp[0], tmp[1].strip())

    elif '_pdbx_rscc_sfcheck_disorder.' in line:
        tmp = line.split('disorder.')
        line = "%sdisorder_TLS_correct.%s\n" % (tmp[0], tmp[1].strip())

    # for refmac:PDBFILE_new_eds_out.cif
    elif '_pdbx_rscc_mapman.' in line:
        tmp = line.split('mapman.')
        line = "%smapman_TLS_correct.%s\n" % (tmp[0], tmp[1].strip())

    elif line.lstrip()[0:6] == "_pdbx_" and "overall." in line:
        tmp = line.split('overall.')
        line = "%soverall_TLS_correct.%s\n" % (tmp[0], tmp[1].strip())

    return line


##########################################################
def check_cif(xyzfile, dic):
    """ Parse items from the cif file and return dic containning
    resolution, Rfactor, Rfree,correlation of Fc os Fo,molecular
    type, Mtrix, twin, TLS, NCS, Anisou ... .
    """

    # print 'Checking cif file=%s' %xyzfile

    if not util.check_file(100, xyzfile):
        return
    dic['detail'] = 'PDB reported'

    flist = open(xyzfile, 'r').readlines()

    items, values = cif.cifparse(flist, '_refine.')
    pdbid = cif.parse_values(items, values, '_refine.entry_id')
    if (pdbid):
        dic['pdbid'] = pdbid[0]
    nref = cif.parse_values(items, values, '_refine.ls_number_reflns_obs')
    if (nref):
        dic['nref'] = nref[0]
    if (nref) and len(nref) > 1:
        dic['nref_n'] = nref[1]

    resh = cif.parse_values(items, values, '_refine.ls_d_res_high')
    if (resh):
        dic['resh'] = resh[0]
    if (resh) and len(resh) > 1:
        dic['resh_n'] = resh[1]

    resl = cif.parse_values(items, values, '_refine.ls_d_res_low')
    if (resl):
        dic['resl'] = resl[0]
    if (resl) and len(resl) > 1:
        dic['resl_n'] = resl[1]

    comp = cif.parse_values(items, values, '_refine.ls_percent_reflns_obs')
    if (comp):
        dic['comp'] = comp[0]
    if (comp) and len(comp) > 1:
        dic['comp_n'] = comp[1]

    rall = cif.parse_values(items, values, '_refine.ls_R_factor_obs')
    if (rall):
        dic['rall'] = rall[0]
    if (rall) and len(rall) > 1:
        dic['rall_n'] = rall[1]

    rfact = cif.parse_values(items, values, '_refine.ls_R_factor_R_work')
    if (rfact):
        dic['rfact'] = rfact[0]
    if (rfact) and len(rfact) > 1:
        dic['rfact_n'] = rfact[1]

    rfree = cif.parse_values(items, values, '_refine.ls_R_factor_R_free')
    if (rfree):
        dic['rfree'] = rfree[0]
    if (rfree) and len(rfree) > 1:
        dic['rfree_n'] = rfree[1]

    pfree = cif.parse_values(items, values, '_refine.ls_percent_reflns_R_free')
    if (pfree):
        dic['pfree'] = pfree[0]
    if (pfree) and len(pfree) > 1:
        dic['pfree_n'] = pfree[1]

    nfree = cif.parse_values(items, values, '_refine.ls_number_reflns_R_free')
    if (nfree):
        dic['nfree'] = nfree[0]
    if (nfree) and len(nfree) > 1:
        dic['nfree_n'] = nfree[1]

    fcc = cif.parse_values(items, values, '_refine.correlation_coeff_Fo_to_Fc')
    if (fcc):
        dic['fcc'] = fcc[0]
    biso = cif.parse_values(items, values, '_refine.B_iso_mean')
    if (biso):
        dic['biso'] = biso[0]
    bmin = cif.parse_values(items, values, '_refine.B_iso_min')
    if (bmin):
        dic['bmin'] = bmin[0]
    bmax = cif.parse_values(items, values, '_refine.B_iso_max')
    if (bmax):
        dic['bmax'] = bmax[0]
    occ_min = cif.parse_values(items, values, '_refine.occupancy_min')
    if (occ_min):
        dic['occ_min'] = occ_min[0]
    occ_max = cif.parse_values(items, values, '_refine.occupancy_max')
    if (occ_max):
        dic['occ_max'] = occ_max[0]

    # ksol = cif.parse_values(items, values, '_refine.solvent_model_param_ksol')
    # bsol = cif.parse_values(items, values, '_refine.solvent_model_param_bsol')

    exp = cif.parse_values(items, values, '_refine.pdbx_refine_id')
    if len(exp) > 1:
        if ('X-RAY' in exp[0] and 'NEUTRON' in exp[1]):
            dic['exp'] = 'xn'
        elif ('X-RAY' in exp[1] and 'NEUTRON' in exp[0]):
            dic['exp'] = 'nx'
    elif len(exp) == 1:
        dic['exp'] = 'x'
        if 'NEUTRON' in exp[0]:
            dic['exp'] = 'n'

    ###
    items, values = cif.cifparse(flist, '_computing.')
    prog = cif.parse_values(items, values, '_computing.structure_refinement')
    if prog:
        dic['prog'] = prog[0].split()[0]

    ###
    c = cif.get_cell(flist)
    dic['cell'] = '%.3f %.3f %.3f %.2f %.2f %.2f' % (c[0], c[1], c[2], c[3], c[4], c[5])
    dic['spg'] = cif.get_symm(flist)

    ###
    items, values = cif.cifparse(flist, '_pdbx_database_related.')
    dbid = cif.parse_values(items, values, '_pdbx_database_related.db_id')
    dbtype = cif.parse_values(items, values, '_pdbx_database_related.content_type')
    if len(dbid) > 1 and len(dbtype) > 1:
        dic['split2'] = ''
        for i, x in enumerate(dbid):
            if dbtype[i] == 'split':
                dic['split2'] = dic['split2'] + dbid[i] + ' '

    ###
    items, values = cif.cifparse(flist, '_atom_site.')
    model = cif.parse_values(items, values, '_atom_site.pdbx_PDB_model_num')
    if not model:
        model = cif.parse_values(items, values, '_atom_site.ndb_model')
    if model:
        if int(model[-1]) > 1:
            dic['model'] = 'Y'
        print('Number of models in the cif file =  %s' % model[-1])
    ###
    #    items,values = cif.cifparse(flist, '_database_PDB_remark.')
    #    text=cif.parse_values(items,values, 'text')
    #    print text


###############################################################################
def check_pdb(pdbfile, dic):
    """ parse items from pdbfile and write a new one: return dic containning
    resolution, Rfactor, Rfree,correlation of Fc os Fo,... .

    """
    if not util.check_file(50, pdbfile):
        return 'XXXX'

    util.get_software_version(pdbfile, 'PDBREP', dic)

    fr = open(pdbfile, "r")
    pdb_new = os.path.basename(pdbfile) + '_new'
    fw = open(pdb_new, "w")

    prob_b, prob_occ = 0, 0
    biso_r, bmat, nmat, xatom, natom, = '', '', 0, 0, 0
    btrix, mtrix, m285, atom_record, biso, occ = [], [], [], [], [], []
    dic['detail'] = 'PDB reported'

    neut = 0
    for line in fr:

        if 'REMARK   3  NEUTRON DATA.' in line:
            neut = 1
        elif 'REMARK   3  X-RAY DATA.' in line:
            neut = 0

        if line[0:6] == 'HEADER':
            dic['pdbhead'] = line[9:49].strip()
            dic['pdbid'] = line[62:66].strip()
            if not len(dic['pdbid']) or len(dic['pdbid']) != 4:
                dic['pdbid'] = 'xxxx'

        elif line[0:6] == 'SPLIT ':
            dic['split2'] = line[7:].strip()

        elif line[0:6] == 'EXPDTA':
            if (('X-RAY' in line[6:] or 'XRAY' in line[6:] or 'X_RAY' in line[6:])
                    and 'NEUTRON' not in line[6:].upper()):
                dic['exp'] = 'x'
            elif ('ELECTRON CRYSTALLOGRAPHY' in line[6:].upper()):
                dic['exp'] = 'e'
            elif ('NEUTRON' in line[6:].upper() and 'X-RAY' not in line[6:].upper()):
                dic['exp'] = 'n'
                dic['phenix_n'] = 1
            elif ('NEUTRON' in line[6:].upper() and 'X-RAY' in line[6:].upper()):
                dic['exp'] = 'xn'
                dic['phenix_xn'] = 1
                continue  # phenix auto use neutron scattering even if it is xray
            else:
                util.perror("Error: The entry (%s) is %s" % (pdbfile, line[6:].strip()))
                print("       Only Xray or Neutron structure can be validated!\n")
                fw.close()
                util.delete_file(pdb_new)
                sys.exit()

        elif "REMARK   3   PROGRAM     : " in line:
            if line[26:].strip():
                dic['prog'] = line[26:].split()[0].upper()

        elif "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) :" in line:
            if not neut:
                dic['resh'] = util.get_value_after_id(line, ":")
            if neut:
                dic['resh_n'] = util.get_value_after_id(line, ":")

        elif "REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) :" in line:
            if not neut:
                dic['resl'] = util.get_value_after_id(line, ":")
            if neut:
                dic['resl_n'] = util.get_value_after_id(line, ":")

        elif "REMARK   3   COMPLETENESS " in line:
            if not neut:
                dic['comp'] = util.get_value_after_id(line, ":")
            if neut:
                dic['comp_n'] = util.get_value_after_id(line, ":")

        elif ("REMARK   3   NUMBER OF REFLECTIONS             " in line
              or 'NUMBER OF REFLECTIONS   (NO CUTOFF)' in line):
            if not neut:
                dic['nref'] = util.get_value_after_id(line, ":")
            if neut:
                dic['nref_n'] = util.get_value_after_id(line, ":")

        elif (" MEAN B VALUE      (OVERALL, A**2) : " in line):
            biso_r = util.get_value_after_id(line, ":")

        elif (('REMARK   3   R VALUE' in line[:21] and '(WORKING + TEST SET)' in line[22:47])
              or ('R VALUE   (WORKING + TEST SET, NO CUTOFF) :' in line)):
            if not neut:
                dic['rall'] = util.get_value_after_id(line, ':')
            if neut:
                dic['rall_n'] = util.get_value_after_id(line, ':')

        elif (('REMARK   3   R VALUE' in line[:21] and ' (WORKING SET)' in line[:47])
              or ('R VALUE          (WORKING SET, NO CUTOFF) : ' in line)):
            if not neut:
                dic['rfact'] = util.get_value_after_id(line, ':')
            if neut:
                dic['rfact_n'] = util.get_value_after_id(line, ':')

        elif ('REMARK   3   FREE R VALUE                    ' in line
              or 'FREE R VALUE                  (NO CUTOFF' in line):
            if not neut:
                dic['rfree'] = util.get_value_after_id(line, ':')
            if neut:
                dic['rfree_n'] = util.get_value_after_id(line, ':')

        elif 'REMARK   3   FREE R VALUE TEST SET SIZE  ' in line:
            if not neut:
                dic['pfree'] = util.get_value_after_id(line, ':')
            if neut:
                dic['pfree_n'] = util.get_value_after_id(line, ':')

        elif ('REMARK   3   FREE R VALUE TEST SET COUNT    ' in line
              or 'FREE R VALUE TEST SET COUNT   (NO CUTOFF)' in line):
            if not neut:
                dic['nfree'] = util.get_value_after_id(line, ':')
            if neut:
                dic['nfree_n'] = util.get_value_after_id(line, ':')

        elif "CORRELATION COEFFICIENT FO-FC      :" in line:
            dic['fcc'] = util.get_value_after_id(line, ":")

        elif "NUMBER OF TLS GROUPS" in line:
            tlsi = util.get_value_after_id(line, ":")
            if tlsi.isdigit() and int(tlsi) > 0:
                dic['tls'] = int(tlsi)
                dic['tls_report'] = 'Y'

        elif ('ATOM RECORD CONTAINS RESIDUAL B FACTORS ONLY' in line):
            dic['isob'] = 'PARTIAL'  # possibly partial

        elif ('REMARK   3  U VALUES      : RESIDUAL ONLY' in line):
            util.perror('Warning: PDB header shows B factor "RESIDUAL ONLY". Check needed!')

        elif ' NCS GROUPS :' in line:
            ncs = util.get_value_after_id(line, ":")
            if ncs.isdigit() and int(ncs) > 0:
                dic['ncs'] = int(ncs)

        elif ('REMARK   3' in line
              and ('  FRACTION:' in line or 'TWIN FRACTION :' in line)):
            fract = line.split(':')[1]
            if util.is_number(fract) and float(fract) > 0:
                dic['twin'] = 'Y'

        elif 'REMARK sg=' in line and 'gamma' in line:  # possible CNS format
            dic['prog'] = 'CNS'
            sg = util.str_between_id(line, 'sg=', 'a=')
            a = util.str_between_id(line, 'a=', 'b=')
            b = util.str_between_id(line, 'b=', 'c=')
            c = util.str_between_id(line, 'c=', 'alpha=')
            alpha = util.str_between_id(line, 'alpha=', 'beta=')
            beta = util.str_between_id(line, 'beta=', 'gamma=')
            gamma = util.str_after_id(line, 'gamma=')
            sg = sg.replace('(', '').replace(')', '')
            c = [float(x) for x in (a, b, c, alpha, beta, gamma)]
            line = ("CRYST1" + '%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %s\n'
                    % (c[0], c[1], c[2], c[3], c[4], c[5], sg))

        elif 'REMARK refinement resolution:' in line[:29] in line:  # CNS
            dic['resl'] = util.value_between_string(line, ':', ' - ')
            dic['resh'] = util.value_between_string(line, ' - ', 'A')

        elif 'REMARK final    r=' in line[:19] and 'free_r=' in line:  # CNS
            dic['rfact'] = util.get_value_after_id(line[:26], 'r=')
            dic['rfree'] = util.get_value_after_id(line[26:], 'r=')

        elif 'REMARK bulk solvent: density level' in line[:36]:  # CNS
            dic['k_sol'] = util.value_between_string(line, 'level=', 'e')
            dic['b_sol'] = util.value_between_string(line, 'factor=', 'A')

        elif 'REMARK total number of reflections used:' in line[:42]:  # CNS
            dic['nref'] = util.value_between_string(line, ':', '(')
            dic['comp'] = util.value_between_string(line, '(', '%')

        elif 'REMARK number of reflections in test set:' in line[:42]:  # CNS
            dic['nfree'] = util.value_between_string(line, ':', '(')
            dic['pfree'] = util.value_between_string(line, '(', '%')

        elif 'REMARK 200  WAVELENGTH OR RANGE        (A) :' in line[:45]:
            t = line[44:54].replace(' ', '').split(',')[0]
            if not util.is_number(t):
                t = '?'
            dic['wavelength'] = t

        elif 'REMARK 285 X0  1' in line[:16]:
            mxo1 = [float(m) for m in line[16:62].split()]
        elif 'REMARK 285 X0  2' in line[:16]:
            mxo2 = [float(m) for m in line[16:62].split()]
        elif 'REMARK 285 X0  3' in line[:16]:
            mxo3 = [float(m) for m in line[16:62].split()]
            m285 = [mxo1, mxo2, mxo3]
        elif 'REMARK 285 ' in line[:12] and '(BIOMT ' in line:
            bmat = util.value_between_string(line, '(BIOMT', ')').strip()

        elif 'REMARK 350   BIOMT1' in line[:19]:
            bt1 = [float(m) for m in line[23:68].split()]
        elif 'REMARK 350   BIOMT2' in line[:19]:
            bt2 = [float(m) for m in line[23:68].split()]
        elif 'REMARK 350   BIOMT3' in line[:19]:
            bt3 = [float(m) for m in line[23:68].split()]
            tmp = [bt1, bt2, bt3]
            btrix.append(tmp)

        elif line[:6] == 'MTRIX1' and '1' not in line[55:].strip():
            mt1 = [float(m) for m in line[10:55].split()]
            continue
        elif line[:6] == 'MTRIX2' and '1' not in line[55:].strip():
            mt2 = [float(m) for m in line[10:55].split()]
            continue
        elif line[:6] == 'MTRIX3' and '1' not in line[55:].strip():
            mt3 = [float(m) for m in line[10:55].split()]
            tmp = [mt1, mt2, mt3]
            mtrix.append(tmp)
            nmat = nmat + 1
            continue

        elif line[:6] == "CRYST1":
            cell = dic['cell'] = line[6:54].strip()
            spg = dic['spg'] = line[54:66].strip()
            if (len(cell) > 30 and len(spg) > 1):
                check_cell_spg(cell, spg)

        elif line[:6] == 'ENDMDL' or line[:6] == 'MODEL ':
            dic['model'] = 'Y'
            continue

        elif line[:3] == 'END' and nmat != 0:
            continue

        elif (('ATOM' in line[:4] or 'HETA' in line[:4] or 'ANISOU' in line[:6])
              and len(line.strip()) > 64):

            if '**' in line:
                util.perror('Error: ** in %s' % line)
                continue

            elif line[20:22] == '  ':
                line1 = line[:20] + " ." + line[22:]
                fw.write(line1)
                continue

            # util.perror('Error in pdb! %s' %line)

            elif 'ANISOU' in line[:6]:
                dic['anis'] = 'Y'
                if ' 0 ' in line[27:36]:
                    t = line[27:].split()
                    if (t[0] == '0' and t[1] == '0' and t[2] == '0' and t[3] == '0' and t[4] == '0' and t[5] == '0'):
                        util.perror('Error: Uij has all 0 values (%s).' % (line[:26]))

            else:  # atom/hetatm
                natom = natom + 1
                if (len(m285) == 3 and len(bmat) == 0):
                    line = new_xyz(m285, line)

                occt, bisot = float(line[54:60]), float(line[60:66])

                if (0 < occt <= 1 and (bisot < 0.001 or bisot > 900)):
                    prob_b = prob_b + 1
                    if (prob_b < 5):
                        util.perror('Warning: B factor problems (B=%6s :%s)' % (line[60:66], line[:26]))
                elif occt <= 0 or occt > 1:
                    prob_occ = prob_occ + 1
                    if (prob_occ < 5):
                        util.perror('Warning: Occupancy problems (OCC=%6s :%s)' % (line[54:60], line[:26]))

                occ.append(occt)
                biso.append(bisot)
                atom_record.append(line)

            if ((('UNL' in line[17:20] or 'UNK' in line[17:20] or 'UNX' in line[17:20])
                 and (' U' in line[76:78] or '  ' in line[76:78])) or ' X' in line[76:78]):
                util.perror('Warning: Unknown ATOM (%s) for (%s). C is assigned for validation!' % (line[76:78], line[0:26]))

                line1 = line[:12] + " C   " + line[17:76] + " C" + line[78:]
                xatom = xatom + 1
                fw.write(line1)
                continue

            elif (' D' in line[76:78] and 'exp' in dic
                  and ((dic['exp'] == 'x' or dic['exp'] == 'e') and dic['phenix_n'] == 0 and dic['phenix_xn'] == 0)):
                line1 = line[:76] + " H" + line[78:]  # for refmac & sfcheck
                fw.write(line1)
                continue

        fw.write(line)
    if (prob_b > 0):
        util.perror('Warning: Number of Bfactor problems =%d' % prob_b)
    if (prob_occ > 0):
        util.perror('Warning: Number of Occupancy problems =%d' % prob_occ)

    if dic['model'] == 'Y':
        print("Note: The PDB file %s has multiple models." % pdbfile)

    if xatom > 0:
        util.perror('Warning: %d X atoms were assigned to C atom for validation!' % xatom)

    if len(m285) == 3:
        print("Note: Matrix from REMARK 285 is applied to coordinates!")

    if (dic['tls'] > 0 and dic['anis'] != 'Y' and 'REFMAC' in dic['prog']):
        dic['tls_attempt'] = 'Y'

    if dic['anis'] == 'Y' and dic['isob'] == 'PARTIAL':
        util.perror('Warning: File header shows partial B, but ANISOU records exist. Check needed!')

    dic['mtr'] = nmat

    if dic['rfact'] == '?' and dic['rall'] != '?':
        dic['rfact'] = dic['rall']

    matt_coeff = 0
    # matt_coeff ,solvent = prog.get_matthew_coeff(pdbfile)  #using a CCP4 module
    matt_coeff, solvent = matt.matthew_coeff(pdbfile)  # use residue/atom/seq
    dic['matt'], dic['solvent'] = '%.2f' % matt_coeff, '%.2f' % solvent
    print('matt_coeff = %s ; solvent = %s ' % (dic['matt'], dic['solvent']))

    if (nmat > 0 and len(bmat) == 0):
        print("Note! New copy of coordinates is generated in ASU..\n")
        for i, mt in enumerate(mtrix):
            fw.write("REMARK  generated with matrix %d\n" % (i + 1))
            if (mt[0][0] == 1.0 and mt[1][1] == 1.0 and mt[2][2] == 1.0
                    and mt[0][3] == 0.0):
                continue
            for line in atom_record:
                natom = natom + 1
                line_new = new_xyz(mt, line)
                fw.write(line_new)

    fr.close()
    fw.close()
    if (len(bmat) > 0):
        natom_bio, pdb_new = coord_by_biomt(pdb_new, m285, btrix, bmat)
        natom = natom + natom_bio

    if len(biso) > 0:
        dic['bmin'] = str(min(biso))
        dic['bmax'] = str(max(biso))
        b = round(sum(biso) / len(biso), 3)
        dic['biso'] = str(b)
        if util.is_number(biso_r):
            if abs(float(biso_r) - b) > 10 or abs(float(biso_r) - b) / b > 0.10:
                util.perror('Warning: Large difference of mean B factor. reported (%.2f) and calculated(%.2f).' % (float(biso_r), b))

    if len(occ) > 0:
        dic['occ_min'] = str(min(occ))
        dic['occ_max'] = str(max(occ))
        dic['occ_avg'] = str(round(float(sum(occ)) / len(occ), 3))

    if natom > 99999:
        dic['one'] = 1
    print('Number of atoms = %d ' % natom)

    dic['natom'] = natom
    return pdb_new


##########################################################
def check_cell_spg(cell, spg):
    ''' check consistence of the crystal system and cell parameters.
    '''

    info = []

    cel = [float(x) for x in cell.split()]
    if len(cel) != 6:
        info.append('Cell dimension problem (%s)' % cell)
        return

    a, b, c, alpha, beta, gamma = cel[0], cel[1], cel[2], cel[3], cel[4], cel[5]
    cryst = util.space_group_crystal_match(spg)
    if cryst == 2:  # MONOCLINIC unique axis on b. alpha, gamma =90
        if alpha != 90.0 or gamma != 90.0:
            print("alpha=%s gamma=%s" % (alpha, gamma))
            info.append('Error:(b setting) alpha or gamma is not 90 for (%s).\n(cell=%s)' % (spg, cell))

    elif cryst == 20:  # MONOCLINIC unique axis on c. alpha,beta =90
        if alpha != 90.0 or beta != 90.0:
            info.append('Error:(c setting) alpha or beta is not 90 for (%s).\n(cell=%s)' % (spg, cell))

    elif cryst == 3 :  # ORTHORHOMBIC
        if not (alpha == beta == gamma == 90.0):
            info.append('Error: cell angle is not 90 for (%s). \n(cell=%s)' % (spg, cell))

    elif cryst == 4:  # TETRAGONAL
        if not (alpha == beta == gamma == 90.0):
            info.append('Error: cell angle is not 90 for (%s).\n(cell=%s)' % (spg, cell))
        if a != b:
            info.append('Error: cell size problem (a != b) for (%s). \n(cell=%s)' % (spg, cell))

    elif cryst == 5:  # TRIGONAL (H setting)
        if not (a == b):
            info.append('Error: problem with cell size for (%s).\n(cell=%s)' % (spg, cell))
        if not (alpha == beta == 90.0 and gamma == 120.0):
            info.append('Error: problem with cell angles for (%s).\n(cell=%s)' % (spg, cell))

    elif cryst == 50 :  # TRIGONAL (rhombohedral axis, a=b=c & alpha=beta=gamma)
        if not (a == b == c):
            info.append('Error: problem with cell size for (%s).\n(cell=%s)' % (spg, cell))
        if not (alpha == beta == gamma):
            info.append('Error: problem with cell angles for (%s).\n(cell=%s)' % (spg, cell))

    elif cryst == 6 :  # HEXAGONAL
        if not (a == b):
            info.append('Error: problem with cell size for (%s).\n(cell=%s)' % (spg, cell))
        if not (alpha == beta == 90.0 and gamma == 120.0):
            info.append('Error: problem with cell angles for (%s).\n(cell=%s)' % (spg, cell))

    elif cryst == 7:  # CUBIC
        if not (alpha == beta == gamma):
            info.append('Error: problem with cell angles for (%s).\n(cell=%s)' % (spg, cell))

        if not (a == b == c):
            info.append('Error: problem with cell size for (%s).\n(cell=%s)' % (spg, cell))
    if len(info) > 0:
        for x in info:
            print(x)
        config.ERRLOG.extend(info)


##########################################################
def coord_by_biomt(pdbfile, m285, btrix, bmat):
    ''' get full xyz using the CRYSTAL AU = (X0) * (BIOMT 1-30, ?,?) * CHAINS
    '''

    print('Transforming PDB using CRYSTAL AU = (X0) * (BIOMT ?).')

    fp = open(pdbfile, 'r').readlines()
    fw = open(pdbfile, 'w')
    biom = bmat.split(',')
    pdbhead, atom = [], []
    natom = 0
    for x in fp:
        if (x[:4] == 'ATOM' or x[:6] == 'HETATM'):
            atom.append(x)
        elif ('CONECT' not in x[:6] and 'MASTER' not in x[:6]
              and 'END' not in x[:3]):
            fw.write(x)
            pdbhead.append(x)

    for x in biom:
        nr = [int(m) for m in x.split('-')]
        for y in range(nr[0] - 1, nr[1]):
            fw.write("REMARK  generated with matrix (X0) * (BIOMT  %d)\n" % (y))
            for ln in atom:
                natom = natom + 1
                ln_bio = new_xyz(btrix[y], ln)
                ln_xo = new_xyz(m285, ln_bio)
                fw.write(ln_xo)
            fw.write('TER     \n')
    fw.close()

    return natom, pdbfile


##########################################################
def new_xyz(mt, line):
    '''Applying matrix (mt) to xyz and generate a new line
    '''

    if line[:4] != 'ATOM' and line[:6] != 'HETATM' or len(line) < 54:
        return line

    x = float(line[29:38])
    y = float(line[38:46])
    z = float(line[46:54])

    xp = mt[0][0] * x + mt[0][1] * y + mt[0][2] * z + mt[0][3]
    yp = mt[1][0] * x + mt[1][1] * y + mt[1][2] * z + mt[1][3]
    zp = mt[2][0] * x + mt[2][1] * y + mt[2][2] * z + mt[2][3]

    xyz = '%8.3f%8.3f%8.3f' % (xp, yp, zp)
    line_new = line[:30] + xyz + line[54:]

    return line_new


##########################################################
def initialize_dic():
    ''' give initial values for the dictionary (reporting)
    '''

    dic = {'map': 0, 'map_dsn6': 0, 'verb': 0, 'refmac': 0, 'refmac_tls': 0,
           'refmactwin': 0, 'scale': 0, 'outfile': '', 'no_xtriage': 0,
           'reso': 0, 'sfcheck': 0, 'wave': 0, 'phenix_x': 0, 'phenix_n': 0,
           'phenix_xn': 0, 'phenix_xtriage': 0, 'refine': 0, 'cns': 0,
           'ligmap': 0, 'sdsc_map': 0, 'ligmapcif': 0, 'shelx': 0,
           'display': 0,
           'mtzmap': '', 'rest': 0, 'cif': 0, 'add_asa': 0, 'add_motif': 0, 'mapsize': 0,
           'one': 0, 'auto': 0, 'bfull': '', 'mapfile': '', 'lib': '',
           'tls_report': 'N', 'tls_attempt': 'N', 'anis': 'N',
           'split1': '', 'split2': '', 'tls_ok': 'N', 'status': 'N',
           'status_used': 'N', 'isob': 'FULL', 'mtr': 0, 'exp': 'x',
           'error': '', 'twin_op': '', 'twin_fr': '', 'ncs': 0, 'tls': 0,
           'xyzlim': 0, 'noeds': 0, 'omitmap': 0,
           'rsrd': 0, 'rsr_all': 0, 'edstat': 0, 'cc': 0, 'anom_f': 0,
           'anom_i': 0, 'sf_i': 0, 'sf_sigi': 0, 'sf_f': 0, 'sf_sigf': 0,
           'fem': 0, 'omit': 0, 'lldf': 0, 'software': [], 'natom': 0}  # options

    dd = ['resh_r', 'resl_r', 'bwilson_r', 'comp_r', 'twin_r', 'nref_r',
          'nfree_r', 'rall_r', 'rfact_r', 'rfree_r', 'prog_r', 'fcc_r',
          'tls_r', 'realr_r', 'dcorr_r', 'detail_r', 'rsrz_p',
          'resh_c', 'resl_c', 'bwilson_c', 'comp_c', 'twin_c', 'nref_c',
          'nfree_c', 'rall_c', 'rfact_c', 'rfree_c', 'prog_c', 'fcc_c',
          'tls_c', 'realr_c', 'dcorr_c', 'detail_c',

          'resh_n', 'resl_n', 'comp_n', 'nref_n', 'rall_n', 'rfact_n',
          'rfree_n', 'nfree_n', 'pfree_n',

          'prog', 'pdbid', 'detail', 'cell', 'spg', 'resh', 'resl', 'bwilson',
          'bwilson_s', 'nref', 'comp', 'nfree', 'pfree', 'rall', 'rfact', 'rfree',
          'fcc', 'rfcc', 'realr', 'realr_d', 'dcorr', 'dcorr_d', 'k_sol', 'b_sol',
          'I2', 'F2', 'E2', 'L', 'L2', 'Z', 'twin', 'twin_ph', 'twin_op_ph',
          'twin_fr_ph', 'twin_type', 'twin_r', 'L2_pt', 'spg_pt',

          'bmin', 'bmax', 'resh_sf', 'resl_sf',
          'biso', 'occ_avg', 'occ_min', 'occ_max', 'dpi_xyz', 'dpi_fr', 'fom',
          'wavelength', 'matt', 'solvent', 'i_sigi_mean', 'i_sigi_resh',
          'i_sigi_resl', 'i_sigi_slop', 'dpi_xyz', 'dpi_fr',
          'pvalue', 'stran', 'ice_ring', 'aniso', 'zscore', 'pdbfile_sig', 'pdbfile_orig',
          'pdbhead', 'pdbfile', 'pdb_tls', 'sffile', 'split', 'dir', 'xyz_type', 'model',
          'xyzfile_orig', 'sffile_orig', 'sf_xray', 'sf_neut', 'fofc', '2fofc', 'pdb_nolig',
          'cifsfch', 'cifrefm', 'cifrefm1', 'cifrefm2', 'cifrefm3',

          'rf_rw', 'bm_bw', 'outlier', 'allowed', 'favored', 'rot_outlier',
          'outlier_num', 'allowed_num', 'favored_num', 'rot_outlier_num',
          'cbeta', 'clashscore', 'oscore', 'bond', 'angle', 'bond_lig',
          'angle_lig', 'bond_lig_max', 'angle_lig_max',
          'dihedral', 'bond_max', 'angle_max', 'dihedral_max', 'chirality',
          'chirality_max', 'planarity', 'planarity_max', 'non_bond', 'fofc_3sig_p',
          'fofc_3sig_n', 'fofc_6sig_p', 'fofc_6sig_n', 'eddir'
          ]  # data

    for x in dd:
        dic[x] = '?'
    dic['pdbid'] = 'xxxx'
    dd = ['sf_nobs', 'sf_nwork', 'sf_nfree', 'sf_pfree', 'sf_nfpair', 'sf_nfpair_p',
          'sf_nfpair_m', 'sf_nfpair_pm', 'sf_nipair', 'sf_nipair_p', 'sf_nipair_m',
          'sf_nipair_pm']
    for x in dd:
        dic[x] = ''
    return dic


##########################################################
def initialize_dic1():
    '''A dictionary for calculations
    '''

    dic1 = {
        'tls_report': 'N', 'tls_attempt': 'N', 'anis': 'N', 'split1': '', 'split2': '',
        'tls_ok': 'N', 'status': 'N', 'status_used': 'N', 'isob': 'FULL',
        'exp': 'x', 'error': '', 'twin_op': '', 'twin_fr': '', 'ncs': 0, 'mtr': 0,
        'tls': 0, 'phenix_n': 0, 'phenix_xn': 0, 'rsrd': 0, 'rsr_all': 0,
        'edstat': 0, 'cc': 0, 'fem': 0
    }

    dd = ['prog', 'pdbid', 'detail', 'cell', 'spg', 'resh', 'resl', 'bwilson',
          'bwilson_s', 'nref', 'comp', 'nfree', 'pfree', 'rall', 'rfact', 'rfree',
          'fcc', 'rfcc', 'realr', 'realr_d', 'dcorr', 'dcorr_d', 'k_sol', 'b_sol',
          'I2', 'F2', 'E2', 'L', 'L2', 'Z', 'twin', 'twin_ph', 'bmin', 'bmax',
          'biso', 'occ_avg', 'occ_min', 'occ_max', 'dpi_xyz', 'dpi_fr', 'fom',
          'wavelength', 'split', 'model', 'pdbfile', 'pdb_new',
          'rf_rw', 'bm_bw', 'outlier', 'allowed', 'favored', 'rot_outlier',
          'cbeta', 'clashscore', 'oscore', 'bond', 'angle', 'dihedral', 'bond_max',
          'angle_max', 'dihedral_max', 'chirality', 'chirality_max', 'planarity',
          'planarity_max', 'non_bond', 'fofc_3sig_p', 'fofc_3sig_n', 'fofc_6sig_p',
          'fofc_6sig_n'
          ]

    for x in dd:
        dic1[x] = '?'
    dic1['pdbid'] = 'xxxx'

    return dic1


##########################################################
def update_bfull(dic):
    '''Generate full B factors in the cif/pdb file
    '''

    file, outfile = dic['bfull'], dic['outfile']
    if not outfile:
        outfile = file + '_bfull'
    pdbfile = dic['bfull'] + '_bfull'
    shutil.copy(dic['bfull'], pdbfile)
    iscif = 0
    pdbfile = dic['bfull']
    if util.is_cif(file):
        pdbfile = cif.cif2pdb(file)
        iscif = 1

    pdb_new = check_pdb(pdbfile, dic)  # tls, also update dic
    tls.check_tls(pdbfile, 0)

    if dic['tls'] == 0 or 'refmac' not in dic['prog'].lower():
        print('Warning: There is no TLS groups or not refined by refmac.')
        return

    if int(dic['tls']) > 0 and 'Y' not in dic['anis']:  # regular update
        wk = tls.proc_tlsanl(pdb_new, pdbfile)
        if wk <= 0:
            print('Error: Failed to generate B full. Check the TLS groups.')
            return
        if iscif == 0:
            print('Note: The new output pdb file = %s\n' % pdbfile)
            return
        else:
            update_cif_biso_uij(file, pdbfile, outfile, 0)

    elif int(dic['tls']) > 0 and 'Y' in dic['anis']:  # complicated update
        print('Warning: ANISOU records exist (Is it a case of TLS + anisoutropic refinement ?).')
        print('If yes, use command  "dcc -bfull  coordinate_file  sf_file" to update Bfull.')

        if not util.check_file(100, dic['sffile']):
            print('Error: Please provide both coordinate and SF file and try again.')
            return

        sf_new = check_sf(dic)  # update dic
        mtz = prog.sf_convertor(sf_new, pdbfile, 'mtz')
        mtz, pdb, log, scr = prog.run_refmac(pdbfile, mtz, dic, 1, ' -tls -addu -all0 ')

        if not util.check_file(100, mtz):
            return

        if iscif == 0 :  # update the pdbfile
            print('Use a cif file.')
            return
            # update_pdb_biso_uij(file, pdb, outfile)
        else:
            update_cif_biso_uij(file, pdb, outfile, 1)


##########################################################
def update_cif_biso_uij(file, pdb, outfile, idd):
    '''idd=0: add anisou;  idd== 1: update anisou
    file: original file;  pdb: the new pdbfile with new Uij, Biso;
    '''

    bnew, u11, u22, u33, u12, u13, u23 = [], [], [], [], [], [], []
    atomid, atom, alt, comp, chain, nres, ins, symbol = [], [], [], [], [], [], [], []

    anis = '''
#
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.type_symbol
_atom_site_anisotrop.pdbx_auth_atom_id
_atom_site_anisotrop.pdbx_label_alt_id
_atom_site_anisotrop.pdbx_auth_comp_id
_atom_site_anisotrop.pdbx_auth_asym_id
_atom_site_anisotrop.pdbx_auth_seq_id
_atom_site_anisotrop.pdbx_PDB_ins_code
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
_atom_site_anisotrop.pdbx_PDB_model_num
'''

    fp = open(pdb, 'r')
    for x in fp:  # get new B
        if 'ATOM' in x[:4] or 'HETA' in x[:4]:
            bnew.append(x[60:66])
        if 'ANISOU' in x[:6]:
            u11.append(float(x[28:35]) * 0.0001)
            u22.append(float(x[35:42]) * 0.0001)
            u33.append(float(x[42:49]) * 0.0001)
            u12.append(float(x[49:56]) * 0.0001)
            u13.append(float(x[56:63]) * 0.0001)
            u23.append(float(x[63:70]) * 0.0001)
            atomid.append(x[6:11].strip())
            atom.append(x[12:16].strip())
            alt.append(x[16:17].strip())
            comp.append(x[17:20].strip())
            chain.append(x[20:22].strip())
            nres.append(x[22:26].strip())
            ins.append(x[26:27].strip())
            symbol.append(x[76:78].strip())
    fp.close()

    fw = open(outfile, 'w')  # write new file
    flist = open(file, 'r').readlines()

    if idd == 0:  # add aniso
        new_flist = update_biso(flist, bnew)
        for x in new_flist:
            fw.write(x)

        fw.write(anis)
        nlen = len(u11)
        for i in range(nlen):

            if not symbol[i]:
                symbol[i] = '.'
            if not alt[i]:
                alt[i] = '.'
            if not chain[i]:
                chain[i] = '.'
            if not ins[i]:
                ins[i] = '.'

            tmp = '%6s %2s %4s %s %3s %2s %4s %s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f 1\n' \
                % (atomid[i], symbol[i], atom[i], alt[i], comp[i], chain[i], nres[i],
                   ins[i], u11[i], u22[i], u33[i], u12[i], u13[i], u23[i])

            fw.write(tmp)

    else:  # update uij

        items, values = cif.cifparse(flist, '_atom_site_anisotrop.')  # a loop
        rows = cif.get_rows(items, values)
        if '_atom_site_anisotrop.U[1][1]' not in items:
            print('Warning: _atom_site_anisotrop.U[1][1] is not found in cif, no correction is applied.')
            return

        if len(rows) != len(u11):
            print('Warning: Number of rows in cif and new B not equal, no correction is applied.')
            return flist

        nu11 = items.index('_atom_site_anisotrop.U[1][1]')
        nu22 = items.index('_atom_site_anisotrop.U[2][2]')
        nu33 = items.index('_atom_site_anisotrop.U[3][3]')
        nu12 = items.index('_atom_site_anisotrop.U[1][2]')
        nu13 = items.index('_atom_site_anisotrop.U[1][3]')
        nu23 = items.index('_atom_site_anisotrop.U[2][3]')

        # print nu11,nu22,nu33,nu12,nu13,nu23

        nrow = len(u11)
        for i in range(nrow):  # upate full uij
            rows[i][nu11] = '%.4f' % u11[i]
            rows[i][nu22] = '%.4f' % u22[i]
            rows[i][nu33] = '%.4f' % u33[i]
            rows[i][nu12] = '%.4f' % u12[i]
            rows[i][nu13] = '%.4f' % u13[i]
            rows[i][nu23] = '%.4f' % u23[i]

        nst, nend = util.position_cif_table(flist, '_atom_site_anisotrop.')  # start - end position
        new_rows = util.format_data(rows)

        newlist = []  # put all lines to a new list
        newlist.extend(flist[:nst])
        newlist.append('#\nloop_\n')
        for x in items:
            newlist.append('%s \n' % x)
        for x in new_rows:
            s1 = ''.join(x) + '\n'
            newlist.append(s1)

        newlist.extend(flist[nend:])  # the last one

        new_flist = update_biso(newlist, bnew)
        for x in new_flist:
            fw.write(x)

    fw.close()

    print('Note: the new output file = %s\n' % outfile)


##########################################################
def update_uij(flist, bnew):
    '''update full B factors in the atom_site.
    flist is a list ['line1','line2',.. ]; bnew is list of new B factors
    '''

    items, values = cif.cifparse(flist, '_atom_site_anisotrop.')  # a loop
    rows = cif.get_rows(items, values)

    nb = -1
    b = '_atom_site.B_iso_or_equiv'
    if b in items:
        nb = items.index(b)

    # print 'rows=', len(rows), len(bnew)

    if nb < 0:
        print('Warning: _atom_site.B_iso_or_equiv is not found in cif, no correction is applied.')
        return flist
    if len(rows) != len(bnew):
        print('Warning: Number of rows in cif and new B not equal, no correction is applied.')
        return flist

    nrow = len(rows)
    for i in range(nrow):  # upate full B
        rows[i][nb] = bnew[i]

    nst, nend = util.position_cif_table(flist, '_atom_site.')  # start - end position

    new_rows = util.format_data(rows)

    newlist = []
    newlist.extend(flist[:nst])
    newlist.append('#\nloop_\n')
    for x in items:
        newlist.append('%s \n' % x)
    for x in new_rows:
        s1 = ''.join(x) + '\n'
        newlist.append(s1)

    newlist.extend(flist[nend:])

    return newlist


##########################################################
def update_biso(flist, bnew):
    '''update full B factors in the atom_site.
    flist is a list ['line1','line2',.. ]; bnew is list of new B factors
    '''

    items, values = cif.cifparse(flist, '_atom_site.')  # a loop
    rows = cif.get_rows(items, values)

    nb = -1
    b = '_atom_site.B_iso_or_equiv'
    if b in items:
        nb = items.index(b)

    # print 'rows=', len(rows), len(bnew)

    if nb < 0:
        print('Warning: _atom_site.B_iso_or_equiv is not found in cif, no correction is applied.')
        return flist
    if len(rows) != len(bnew):
        print('Warning: Number of rows in cif and new B not equal, no correction is applied.')
        return flist

    nrow = len(rows)
    for i in range(nrow):  # upate full B
        rows[i][nb] = bnew[i]

    nst, nend = util.position_cif_table(flist, '_atom_site.')  # start - end position

    new_rows = util.format_data(rows)

    newlist = []
    newlist.extend(flist[:nst])
    newlist.append('#\nloop_\n')
    for x in items:
        newlist.append('%s \n' % x)
    for x in new_rows:
        s1 = ''.join(x) + '\n'
        newlist.append(s1)

    newlist.extend(flist[nend:])

    return newlist


##########################################################
def pdb_bfull(pdb_new):
    ''' convert B to full if TLS records exist and ANISOU records missing
    '''

    pdb_tls = pdb_new + '__tls'
    wk = tls.proc_tlsanl(pdb_new, pdb_tls)
    bval = check_bfactor(pdb_tls)

    if wk == 1 and min(bval) <= 0:  # only for the older version of refmac5
        minb = 2 - min(bval)
        print('Correction: Add %.2f to B_full and re-run tlsanl' % minb)
        pdb_newb = add_const_to_biso(pdb_new, minb, 1)
        wk = tls.proc_tlsanl(pdb_newb, pdb_tls)  # do TLS again

    return wk, pdb_tls


##########################################################
def check_result(d1, d2):
    '''d1 is the reported,  d2 is the calculated.
    '''
    # d1=[d1['resh'], d1['resl'], d1['rfact'], d1['rfree'], d1['nref'],d1['nfree'],d1['comp'] ] reported
    # d2=[d2['resh'], d2['resl'], d2['rfact'], d2['rfree'], d2['nref'],d2['nfree'],d2['comp'] ] calculated

    rcut = 0.02  # if abs(Rcal-Rrep)>0.02, report warning messages.
    if d1['model'] == 'Y':
        rcut = 0.05

    info = []

    if (util.is_number(d1['nfree']) and util.is_number(d1['sf_nfree'])):
        if int(d1['nfree']) > int(d1['sf_nfree']) + 1:
            info.append('Warning: The reported number of free set in coordinates (%s) > the total number of free set in SF (%s).' % (d1['nfree'], d1['sf_nfree']))

    if util.is_number(d1['wavelength']):  # wavelength
        wave = float(d1['wavelength'])
        if wave > 2.0 or wave < 0.5:
            info.append('Warning: Wavelength (%s) is abnormal (double check).' % wave)

    if (util.is_number(d1['resh']) and util.is_number(d1['rfact'])):  # Rwork reported
        if float(d1['resh']) < 1.7:
            if float(d1['rfact']) > 0.23:
                info.append('Warning: The reported Rwork (%s) is too large at resolution (%s)' % (d1['rfact'], d1['resh']))
        elif float(d1['resh']) < 2.5:
            if float(d1['rfact']) > 0.28:
                info.append('Warning: The reported Rwork (%s) is too large at resolution (%s)' % (d1['rfact'], d1['resh']))
        elif float(d1['resh']) < 3.2:
            if float(d1['rfact']) > 0.33:
                info.append('Warning: The reported Rwork (%s) is too large at resolution (%s)' % (d1['rfact'], d1['resh']))
        elif float(d1['rfact']) > 0.36:
            info.append('Warning: The reported Rwork (%s) is large too at resolution (%s)' % (d1['rfact'], d1['resh']))

    if d1['tls_report'] == 'Y' and d1['anis'] == 'N':
        info.append('Warning: TLS reported but there is no ANISOU records.')

    if (util.is_number(d1['resh']) and float(d1['resh']) > 1.7
            and d1['tls_report'] == 'N' and d1['anis'] == 'Y' and d1['model'] == '?'):
        info.append('Warning: ANISOU records exists but no TLS records (resolution=%sA).' % d1['resh'])

    if (util.is_number(d1['resh']) and float(d1['resh']) < 1.5 and d1['tls_report'] == 'Y'
        and util.is_number(d1['rfact']) and util.is_number(d2['rfact'])
            and float(d2['rfact']) - float(d1['rfact']) > 0.015 and 'REFMAC' in d1['prog'].upper()):
        info.append('Warning: Possibly B_partial and mixed (TLS + individual_B) refinement.')

    if util.is_number(d1['resh']) and util.is_number(d2['resh']):  # resol
        if float(d1['resh']) - float(d2['resh']) < -0.01:
            info.append('Warning: Reported resolution (%s) higher than calculated (%s).' % (d1['resh'], d2['resh']))

        if abs(float(d1['resh']) - float(d2['resh'])) > 0.3:
            info.append('Warning: Large difference of resolution: reported (%s), calculated (%s).' % (d1['resh'], d2['resh']))
    if util.is_number(d1['resh']) and util.is_number(d1['resl']):  # resol
        if float(d1['resl']) - float(d1['resh']) < 4.0:
            info.append('Warning: Too small difference between high/low (%s , %s) resolution.' % (d1['resh'], d1['resl']))

    if util.is_number(d2['resh']) and util.is_number(d2['resl']):  # resol
        if float(d2['resl']) - float(d2['resh']) < 4.0:
            info.append('Warning: Too small difference between high/low (%s , %s) resolution.' % (d2['resh'], d2['resl']))

    if d1['mtr'] > 0 and 0 < float(d1['matt']) < 1.8:
        util.perror('Warning: Matrix should be flagged as 1. Double check the NCS MTRIX records.\n')

    if util.is_number(d1['rfact']) and util.is_number(d2['rfact']):  # R_work
        rwk_diff = float(d2['rfact']) - float(d1['rfact'])

        # print 'The new L2 & rdiff=', rwk_diff,  d1['L2_pt']

        if (abs(rwk_diff) > rcut and 'SFCHECK' not in d2['prog']) or (abs(rwk_diff) > 0.05 and 'SFCHECK' in d2['prog']):
            info.append('Warning: Large difference of R_work: reported (%s), calculated (%s).' % (d1['rfact'], d2['rfact']))

        if 'phenix' in d1['prog'].lower() and 'phenix' not in d2['prog'].lower() and rwk_diff > rcut:
            info.append('Note: It is suggested to do validation by phenix for this entry.')

        if (util.is_number(d1['L2_pt']) and float(d1['L2_pt']) < 0.16 and rwk_diff > 0.03):  # swapped F -> I
            info.append('Warning: possible swapped SF file F->I! (twin test L2=%s)' % d1['L2_pt'])

        if (util.is_number(d1['L2_pt']) and float(d1['L2_pt']) > 0.45 and rwk_diff > 0.20):  # swapped I -> F
            info.append('Warning: possible swapped SF file I->F! (twin test L2=%s)' % d1['L2_pt'])

        if (util.is_number(d1['matt'])):
            if float(d1['matt']) > 8.1 or float(d1['matt']) < 1.4:
                info.append('Warning: Matthew_coefficient(%s) is abnormal. May have crystal packing problem. ' % d1['matt'])
            if (float(d1['matt']) > 5.2 and rwk_diff > 0.15 and float(d2['rfact']) > 0.33):  # split
                info.append('Warning: Possible incomplete content of ASU (or a split entry). (Matt. coeff.=%s)' % d1['matt'])

        if d1['twin'] != 'Y' and '2:' in d1['twin_op'] and '2:' in d1['twin_fr'] and rwk_diff < -0.03 :  # twin by refmac
            info.append('Warning: twin was not reported, but is detected by refmac with the operators and fraction as below:')
            info.append(d1['twin_op'])
            info.append(d1['twin_fr'])

        if d1['twin'] != 'Y' and d1['twin_ph'] == 'Y' and rwk_diff < -0.03:  # twin by xtriage
            info.append('Warning: twin is detected by xtriage.')

        if d1['twin'] == 'Y' and rwk_diff > rcut:
            info.append('Warning: The deposited SF file might be detwinned. (double check!)')

        if 'use reported reso' in d2['detail'].lower() and abs(rwk_diff) <= 0.01:
            info.append('Warning: The reported and calculated resolution differs. Rfactor was calculated using the reported resolution range(%s - %s).' % (d1['resl'], d1['resh']))
            # info.append('         Data may have problem outside of the resolution range(%s - %s).' % (d1['resl'], d1['resh']))

    if util.is_number(d1['rfree']) and util.is_number(d2['rfree']):  # R_free
        rfwk_diff = abs(float(d1['rfree']) - float(d2['rfree']))
        if (rfwk_diff > rcut and 'SFCHECK' not in d2['prog']) or (rfwk_diff > 0.05 and 'SFCHECK' in d2['prog']) :
            info.append('Warning: Large difference of R_free: reported (%s) calculated (%s).' % (d1['rfree'], d2['rfree']))

    if util.is_number(d1['rfree']) and util.is_number(d1['rfact']):  # R_free & R_work (reported)
        val = float(d1['rfree']) - float(d1['rfact'])
        if val <= 0.005 and val >= 0:
            info.append('Warning: Too small difference between the reported R_free (%s) and R_work(%s).' % (d1['rfree'], d1['rfact']))
        if val < 0:
            info.append('Warning:  The reported R_free(%s) is smaller than R_work (%s).' % (d1['rfree'], d1['rfact']))

        if val > 0.10:
            info.append('Warning:  Large difference between the reported R_free (%s) and R_work(%s). Possible over/biased refinement!' % (d1['rfree'], d1['rfact']))

    if util.is_number(d2['rfree']) and util.is_number(d2['rfact']):  # R_free & R_work (calculated)
        val = float(d2['rfree']) - float(d2['rfact'])

        if val < 0:
            info.append('Warning: R_free(%s) is smaller than R_work(%s).' % (d2['rfree'], d2['rfact']))

        if val <= 0.01 and util.is_number(d2['resh']) and float(d2['resh']) > 1.5:
            info.append('Warning: There may be problem with free set. Please try free set 1 and test again.')

    if util.is_number(d1['rfree']) and not util.is_number(d2['rfree']):  # R_free
        info.append('Warning: R_free reported (%s), but not calculated.' % (d1['rfree']))

    if util.is_number(d1['comp']) and float(d1['comp']) < 65:  # comp
        info.append('Warning: Too small for reported data completeness (%s%%)' % (d1['comp']))
    if util.is_number(d1['comp']) and float(d1['comp']) > 100:  # comp
        info.append('Error: Too large for reported data completeness (%s%%)' % (d1['comp']))

    if 'comp' in d2 and util.is_number(d2['comp']) and float(d2['comp']) < 65:  # comp
        info.append('Warning: Too small for calculated data completeness (%s%%)' % (d2['comp']))

    if 'comp' in d2 and util.is_number(d1['comp']) and util.is_number(d2['comp']) and abs(float(d2['comp']) - float(d1['comp'])) > 20:
        info.append('Warning: Large difference of completeness reported(%s), calculated(%s).' % (d1['comp'], d2['comp']))

    # if util.is_number(d1['nref']) and util.is_number(d2['nref']):  #Nref
    #    if abs(float(d1['nref']) - float(d2['nref'])) >3000:
    #        info.append('Warning: Large difference of reflection: reported (%s), calculated (%s).' % (d1['nref'],d2['nref']))

    if (util.is_number(d1['nfree']) and util.is_number(d2['nfree'])
            and util.is_number(d1['resh']) and util.is_number(d2['resh'])):  # Nfree
        if abs(float(d1['resh']) - float(d2['resh'])) < 0.04 and abs(float(d1['nfree']) - float(d2['nfree'])) > 200:
            info.append('Warning: Large difference of free set: reported (%s), calculated (%s).' % (d1['nfree'], d2['nfree']))

    if (util.is_number(d2['k_sol']) and float(d2['k_sol']) < 0.1):  # k_solvent
        info.append('Warning: calculated k_sol (bulk solvent) (%s) is too small.' % (d1['k_sol']))
    if (util.is_number(d2['b_sol']) and float(d2['b_sol']) < 1.0):  # b_solvent
        info.append('Warning: calculated b_sol (bulk solvent) (%s) is too small.' % (d1['b_sol']))

    if d1['stran'] == 'Y':  # pseudo translational
        info.append('Note: Translational pseudo symmetry is detected by xtriage.')

    if (util.is_number(d1['biso']) and util.is_number(d1['bwilson'])):  # Biso with B_wilson
        br, bw = float(d1['biso']), float(d1['bwilson'])
        if br > 1.7 * bw or bw > 1.7 * br:
            info.append('Warning: Too much difference between Biso (%.2f) and B_Wilson(%.2f)' % (br, bw))

    if (util.is_number(d1['i_sigi_slop']) and float(d1['i_sigi_slop']) > 10):  # I/SigI (diff of last shell)
        info.append('Warning: I/SigI differs too much (%s) in the last two shells.' % d1['i_sigi_slop'])

    if (util.is_number(d1['L2']) and (float(d1['L2']) > 0.45 or float(d1['L2']) < 0.23)):  # twin
        info.append('Warning: Abnormal value for Padilla-Yeates_L2 test (%s). SF twinning or other problems.' % d1['L2'])
    elif (util.is_number(d1['L2_pt']) and (float(d1['L2_pt']) > 0.45 or float(d1['L2_pt']) < 0.23)):
        info.append('Warning: Abnormal value for Padilla-Yeates_L2 test (%s). SF twinning or other problems.' % d1['L2_pt'])

    if (util.is_number(d1['fom']) and float(d1['fom']) < 0.5):  # FOM
        info.append('Warning: Figure of merit is too small (FOM=%s). May have refinement problem.' % d1['fom'])

    if (util.is_number(d1['dpi_fr']) and float(d1['dpi_fr']) > 0.8):  # DPI
        info.append('Warning: Cruickshank_dpi (%s) is large. May have refinement problem.' % d1['dpi_fr'])

    if (util.is_number(d1['aniso']) and float(d1['aniso']) > 0.8):  # ANIS
        info.append('Warning: SF data is very anisotropic (ratio=%s).' % d1['aniso'])

    if '?' in d2['rfact']:
        info.append('Error: SF validation failed.')

    if len(info) > 0:
        for x in info:
            print(x)
        config.ERRLOG.extend(info)


##########################################################
def check_error(pdbfile, sffile, dic):
    '''check errors in sf/pdb file and the consistence
    print all error/warning message to dic['error']
    '''

    if not util.check_file(500, pdbfile, sffile):
        return

    cell_sf, spg_sf = ' ' , ' '
    flist = open(sffile, 'r').readlines()

    cell_sf = cif.get_cell(flist)
    spg_sf = cif.get_symm(flist)

    if '?' in dic['cell']:
        util.perror('Warning: input file (%s) has no cell parameters!' % pdbfile)
        cell = []
    else:
        cell = [float(x) for x in dic['cell'].split()]

    spg_pdb = dic['spg'].strip().upper()

    if len(cell) == len(cell_sf) == 6 and min(cell_sf) > 0.01:
        tmp1 = [abs(cell[i] - cell_sf[i]) for i in range(3)]
        tmp2 = [abs(cell[i] - cell_sf[i]) for i in range(3, 6)]
        tmp1_max, tmp2_max = max(tmp1), max(tmp2)
        if tmp1_max > 0.01 or max(tmp2) > 0.01:
            t1 = 'Error: Cell mismatch in PDB and SF (cell_max=%.2f, angle_max=%.2f).' % (tmp1_max, tmp2_max)
            util.perror(t1)

    if len(spg_sf) > 1 and len(spg_pdb) > 1:
        print('\nSpace group (%s) in coordinate;  (%s) in SF;   Resolution=%s' % (spg_pdb, spg_sf, dic['resh']))
        d1, d2 = spg_sf.split(), spg_pdb.split()
        if (len(d1) == 4 and d1[1] == d1[3] == '1'):
            spg_sf = d1[0] + d1[2]
        if (len(d2) == 4 and d2[1] == d2[3] == '1'):
            spg_pdb = d2[0] + d2[2]
        if (len(d1) == 3 and d1[2] == '1' and len(d1[0]) == 2 and d1[0][1] == '1'):
            spg_sf = d1[0][0] + d1[1]
        if (len(d2) == 3 and d2[2] == '1' and len(d2[0]) == 2 and d2[0][1] == '1'):
            spg_pdb = d2[0][0] + d2[1]

        spg_sf = spg_sf.replace(' ', '')
        spg_pdb = spg_pdb.replace(' ', '')

        if spg_pdb.upper() != spg_sf.upper():
            # update_sf_space_group(sffile, pdbfile, spg_pdb) #no longer active
            util.perror('Error: Space group mismatch in PDB (%s) and SF (%s).' % (spg_pdb, spg_sf))
            # util.perror('Warning: Space group in SF is corrected (%s)->(%s).' % (spg_sf, spg_pdb))

    if util.is_number(dic['bmin']) and float(dic['bmin']) < 0:
        util.perror('Error: B factors <0 in coordinate file. ')

    if util.is_number(dic['bmax']) and float(dic['bmax']) > 999:
        util.perror('Warning: B factors >999 in coordinate file.')

    if util.is_number(dic['occ_min']):
        if float(dic['occ_min']) == 0:
            util.perror('Warning: Atoms have 0 occupancy in coordinate file.')
        elif float(dic['occ_min']) < 0:
            util.perror('Error: Atoms have negative occupancy in coordinate file.')

    if util.is_number(dic['occ_max']) and float(dic['occ_max']) > 1:
        util.perror('Error: Occupancy >1 in the model file (occ_max=%s).' % dic['occ_max'])


##########################################################
def update_sf_space_group(sffile, pdbfile, spg_pdb):
    '''Replace space group and symmetry operators in PDB
    '''

    fr = open(sffile, 'r').readlines()
    fw = open(sffile, 'w')

    patt = 'REMARK .*290.*555'  # symmetry code in PDB
    item = prog.grep_items(pdbfile, patt)

    symm = '''#
loop_
_symmetry_equiv.id
_symmetry_equiv.pos_as_xyz
'''

    fr = remove_onecif_table(fr, '_symmetry_equiv.')
    n1, n2 = 0, len(fr)
    for i in range(n1, n2):
        if '_symmetry.space_group_name_H-M' in fr[i]:
            tmp = "_symmetry.space_group_name_H-M  '%s' \n" % spg_pdb
            fw.write(tmp)
            if len(item) > 1:  # add _symmetry_equiv from PDB
                fw.write(symm)
                m = 0
                for x in item:
                    m = m + 1
                    if len(x.strip()) < 4:
                        continue
                    y = x[22:].replace(' ', '')
                    tt = '%2d %s\n' % (m, y)
                    fw.write(tt)

            continue
        elif '_symmetry.Int_Tables_number ' in fr[i]:
            continue

        fw.write(fr[i])

    fw.close()


##########################################################
def remove_onecif_table(cif, table):
    '''a simplye way to remove a table in cif (not a safe way)
    '''

    n1, n2, m = 0, 0, len(cif)
    for i in range(0, m):
        if table in cif[i].strip()[:len(table) + 1]:
            n1 = i
            if 'loop_' in cif[i - 1]:
                n1 = i - 1
            for j in range(i, m):
                if ('#' in cif[j][:1] or 'loop_' in cif[j][:5]
                        or ('_' in cif[j][:1] and table not in cif[j])):
                    n2 = j
                    break
            if n2 > n1:
                break
    # print(len(table), len(cif),  n1, n2)
    if n2 < len(cif) and n2 > 0:
        del (cif[n1:n2])
    return cif


##########################################################
def check_ligand(pdbfile, id):
    '''Check if there is ligand in the pdbfile
    id=0, for ligand and peptide; id= 1, only for ligand;
    '''

    pdb = open(pdbfile, 'r').readlines()
    ch_pep, chr_pep, ch_lig, chr_lig, ch_wat, chr_wat = tls.chain_res_range(pdbfile)

    # print ch_lig,chr_pep
    if id == 0 or id == 1:
        for k, v in ch_lig.items():
            for y in v :
                natm = 0
                for x in pdb:
                    if (('ATOM' in x[:4] or 'HETATM' in x[:6])
                            and k == x[20:22].strip() and int(x[22:26]) == y):
                        natm = natm + 1
                        # print natm, k, v, x
                        if natm > 1:
                            return 1
    if id == 0:
        for k, v in ch_pep.items():
            if len(v) < 15:
                return 1

    return 0


##########################################################
def check_bfactor(pdbfile):
    ''' check B factors in the PDB. if B_full<0 add a constant
    '''

    bval, negb, bigb = [], [], 0
    if os.path.exists(pdbfile) and os.path.getsize(pdbfile) > 1000:
        fr = open(pdbfile, "r")
        for x in fr:
            if ((x[0:5] == "ATOM " or x[0:6] == 'HETATM') and '*' in x[60:66]):
                bigb = bigb + 1
                continue
            elif ((x[0:5] == "ATOM " or x[0:6] == 'HETATM') and float(x[60:66]) < 0):
                negb.append(float(x[60:66]))

            if x[0:4] == "ATOM" or x[0:4] == 'HETA':
                bval.append(float(x[60:66]))

        fr.close()
    n = len(negb)
    if bigb > 0:
        util.perror('Error: There are %d atom with B_full>999.' % (bigb))
    if n > 0:
        util.perror('Error: There are (%d) atom with B_full<0.' % (n))

    return bval


##########################################################
def add_const_to_biso(pdbfile, minb, id):
    ''' Add minb to Biso of each atom
    if id==0, make Biso all zero;  id==1, add minb
    '''

    fp = open(pdbfile, 'r').readlines()
    out = pdbfile + '_newb'
    fw = open(out, 'w')
    for x in fp:
        if id == 0 and 'ANISOU' in x[:6]:
            continue
        if 'ATOM' in x[:4] or 'HETA' in x[:4]:
            if id > 0:
                t1, t2, t3 = x[:60], float(x[60:66]) + minb, x[66:].rstrip()
                tmp = '%s%6.2f%s\n' % (t1, t2, t3)
            else:
                # t1, t2, t3, t4=x[:54], 1, 0, x[66:].rstrip() #both occ/Biso
                # tmp='%s%6.2f%6.2f%s\n' % (t1, t2, t3, t4)
                t1, t2, t3 = x[:60], 0, x[66:].rstrip()  # both occ/Biso
                tmp = '%s%6.2f%s\n' % (t1, t2, t3)
            fw.write(tmp)
        else:
            fw.write(x)
    fw.close()
    return out


##########################################################
if __name__ == '__main__':
    process(sys.argv)
    # get_rsr_database()
    # add_segid_to_outfile('test.cif', 'test.sum')
    # prog.calc_asa_areaimol('1fggu.pdb', 'asa.out')
    # check_pdb(sys.argv[1], dic={})
    # check_tls(sys.argv[1],0)
    sys.exit()
