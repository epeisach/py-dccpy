import os

# ===========================================================
# this module is only used as a global variable
# ===========================================================
PATH1 = ""  # for CGI, give full path for sfcheck, sf_convert
PATH2 = ""  # only for CGI, give full path for lx_mapman, shelx,
PATH3 = ""  # contains database of RsR for each residue and resolution
PATH4 = ""  # contains the src
if "DCCPY" in os.environ:
    PATH1 = os.environ["DCCPY"] + "/bin/"
    PATH2 = os.environ["DCCPY"] + "/other_bin/"
    PATH3 = os.environ["DCCPY"] + "/data/"
    PATH4 = os.environ["DCCPY"] + "/dccpy/"

ERRLOG = []  # contains all the error/warning messages

VERSION = "2.28  (2023-09-04))"

# ===================Some pre-defined messages =====================#
DCC = "A wrapper of crystallographic applications"
VALID_REP = "refinement program reported in coordinate file."
VALID_CAL = "model vs diffraction data validation (program used by DCC)."
RSR = "calculate density metrics such as RsR/RSCC/RSRZ"
DATA = "check structure factor quality"
SF_FMT = "convert structure factor format"

# ====================modifications below ==========================#

# ##############################################################################
#  started from version 1.0 (2009-05-20 created ), Current version same as name.
#  modifications:
#  1. If TLS correction is bad (e.g. Btot does not behaves well due to improper
#     range or refinement), the TLS correction will be skipped.
#  2. Add virus calculations for REFMAC.
#     (If MTRIX exist, generate full PDB for refmac)
#  3. Correct reported Rfactor extraction.
#  4. Add option to seperate sfcheck and REFMAC calculation
#  5. Add option '-pdbid ?' to replace '-pdb ? -sf ?'
#  6. Add -noed to skip density calculation (only for validation)
#  7. Add -refmac or -sfcheck for individual calculations
#  8. If multi-model exist in PDB, only the first model will be used.
#  9. add -test_r to test R, Rfree, Fcc by swaping I->F or F->I
# 10. corrected write_summary function to use only the resolution & Rfactor
#     as a indicator for writing refmac and sfcheck output.
# 11. improved sf_convert for mmcif to mmcif validation convertion.
# 12. If MTRIX record exist in coordinate file, full contents will be
#     generated for validation based on the value of Matthews Coefficient
#     if MC>5 and MTRIX exist, full contents in ASU are generated.else not.
# 13. Add -twin to force refmac to do twin validation (11/03/2009)
# 14. Add -o  to give user output name (11/03/2009)
# 15. Add -phenix_x for Xray validation using model_vs_data (11/30/2009)
# 16. Add -phenix_n for neutron validation using model_vs_data (11/30/2009)
# 17. Add -phenix_xn for Xray & neutron validation (11/30/2009)
# 18. Add -xtriage for SF Quality Assessment (11/30/2009)
# 19. Add -path (only for CGI, give full path containing sfcheck,sf_convert
#                ,lx_mapman (11/30/2009))
# 20. converted output of phenix into mmcif (12/04/2009)
# 21. Improve the TLS correction. If using defaut residue range (-? 9999),
#     the correct range will be put to the NEW pdb file (12/21/09)
# 22. Change unknow atom (X) to C for REFMAC & Phenix,. (12/22/09)
# 23. Changed threshold value Matt. Coeff=3.4 for NCS (02/02/2010)
# 24. TLS range correction also include HETATM (rare case) (12/03/2010)
# 25. If ANISOU records exist, refine type is given as 'bref ANIS' (03/15/2010)
# 26. parse data from xtriage and merge it to the summary cif table(04/06/2010)
# 27. parse matrix from REMARK 285 and do validation (04/15/2010)
# 28. Add a status check for NCS matrix(if 1 in 7th column, skip)(04/28/2010)
# 29. improved B factor status (Full or Partial) identification (05/11/2010)
# 30. TLS residue range include ligands (05/12/2010)
# 31. replace tlsextract which had problem with spacing(07/11/2010)
# 32. add cif item to hold all the error/warning message(08/15/2010)
# 33. add cif item for all phenix.model_vs_data(10/02/2010)
# 35. make REMARK 285 general. if 285 and BIOMT exist, new mtrix will
#     be calculated and used. If only 285 exist but not BIOMT, only
#     apply it to coordinates. (11/05/2010)
# 36. improved the error check of SG between SF and PDB.
#     The space group in SF is updated for correct one. (11/09/2010)
# 37. Added -refmac_tls to do tls correction by refmac (not tlsanl)(12/04/2010)
# 38. use new tls residue range correction from tls.py(12/21/2010)
# 39. Add checking TLS origin (the reported and calculated) (03/28/2011)
# 40. Add checking consistence between cystal system and cells(03/28/2011)
# 41. revised xtriage and neutron-xray hybrid (removed -i mmcif)(05/11/2011)
# 42. interpret the result (R/Rfree) and report errors in output(07/28/2011)
# 43. cut BIG map to small (around ligands)(08/31/2011)
# 44. generate small map around ligands (08/31/2011)
# 45. add validation for SF file itself. modify phenix (02/16/2012)
# 46. include mixed TLS (phenix_buster) selection (03/02/2012)
# 47. add -scale use 'BULK LSSC ANISO EXPE', else all by default refmac(03/15/2012)
# 48. modified tls parser for phenix (include not) (03/19/2012)
# 49. Add option -refine to do refinement by refmac (04/27/2012)
# 50. Modified ligand.py to handle jmol-13 (10/24/2012)
# 51. Modified CNS to parse both the log & stat.list (10/25/2012)
# 52. Add an option '-all' to do all programs (01/16/2013)
# 53. Add '-no_xtriage' to exclude the xtriage (02/25/2013)
# 54. Add '-one' for bigcif file to map 2 char chainID to residue.(04/10/2013)
#     If the cif file (-pdb ciffile) is bigcif, '-one' is automated.
# 55. TLS consistency check (between xyz and header of pdb).(04/15/2013, v2.02)
# 56. Add -auto for ebi (try refine program first, if failed, try refmac.(06/10/2013)
# 57. validation for ensemble refinement by refmac...(09/05/2013)
# 58. Add validation for NCS records (09/26/2013)
# 59. output more data from phenix (02/21/2014)
# 60. make occupancy in f6.3 format in -cif2pdb: make maximum model 180;(06/16/2014)
# 61. fix a bug for large entry (cifparse.cif2pdb) (06/17/2014)
# 62. improved calculation for many model xray (remove the restraint)  (06/17/2014)
# 63. add -omitmap to calculate RsR/dencorr. with Ligand omitted, Adjust Rcut=0.2(06/24/2014)
# 64. use all default of refmac to get 2FoFc map and FC by sfall. (calculate only once)
#     add -ligmapcif option for D&A (same style kept)  (07/16/2014)
# 65. if setenv MAPMAN_BINARY path/lx_mapman, specified binary is used (07/25/2014)
# 66. Corrected writing TLS in the buster format  (08/10/2014)
# 67. Add Zscore wrt to the whole database,similar resolu (08/19/2014)
# 68. Add option '-fem' to get FEM map and  CC/RsR  (09/11/2014)
# 69. separate geometry and pdbx_density, add more infor  (09/19/2014)
# 70. add new cif table for bad RsR and put cif tokens to write_cif (10/09/2014)
# 71. corrected cifparse.py to generate correct PDB (1/09/2015)
# 72. add -rsr_all option to get RSR in all,main/side chain, phosphate (2/1/2015)
# 73. add -edstat option to get RSR/ZDIFF in all,main/side chain, phosphate(3/10/2015)
# 74. fixed TLS problem for two letter chainID(4/25/2015)
# 75. sorted options to public for releasing DCC (9/10/2015)
# 76. modified tls.py to detect TLS truncation errors (1/21/2016)
# 77. change pdbx_DCC_software to pdbx_dcc_software (5/12/2016)
# 78. add table _pdbx_dcc_sf (12 items) (5/13/2016)
# 79. export error/warning messages for -auto as others (5/13/2016)
# 80. correct parsing of pointless log file (after 1.10.13)(5/20/2016)
# 81. fixed errors from Xray-Neutron validation and in tool options. (6/28/2016)
# 82. refactoring a lot of functions and export standard CIF (6/28/2016)
# 83. correct extraction of phenix version (7/9/2016)
# 84. export residues for one copy for EDS (map from all copies) (7/9/2016)
# 85. make faster for '-auto' option when R-report and R-calc diff a lot. (7/14/2016)
#     fix bug for sfcheck (connecting two large numbers in cif file.
# 86. fixed a new bug (can not find labels. stoped if only anomalous data) (8/5/2016)
# 87. Corrected the format of BUSTER TLS origin. (%10.4f) (8/8/2016)
# 88. If no chain-ID, given '.' and issue a warning message. (8/8/2016)
# 89. modified cifparse to handle B>1000 (cif->pdb). (9/1/2016)
# 90. fix omitmap in OneDep (dcc -cif ? -sf ? -ligmapcif -no_xtriage -omitmap) (9/12/2016)
# 91. Add error checking if reported NFREE>calculated NFREE in SF (10/20/2016)
# 92. modified -auto option. if phenix/refmac worked, no more continue regardless
#     of the Rfactor. (10/20/2016)
# 93. Add space group (P2 21 2). (11/2/2016)
# 94. correct program extraction (software first, then computing). (11/22/2016)
# 95. Correct write_cif.py, use asym chain only, no segid (11/22/2016)
# 96. modified (-auto) option to run phenix first for unsuported programs(11/29/2016)
# 97. modified cifparse.py to change '?' to 'NULL' when back to PDB format(1/11/2017)
# 98. revised cifparse.py, to handle '?' in Matrix (none stop) (2/1/2017)
# 99. revised tls.py, to handle weird buster TLS group (none stop) (2/1/2017)
# 100. revised cifparse.py to handle weird buster TLS group (none stop) (3/6/2017)
# 101. revised dcc-calc.py to handle multi-model with ENDMDL entrie  (3/9/2017)
# 102. revised dcc-calc.py to handle twin problem refined by phenix in auto mode. (4/25/2017)
# 103. corrected words for misleading message (compare free set in xyz & SF). (4/25/2017)
# 104  Support fore phenix 1.13 parsing of model_vs_data
# 105. Crash in neutron/X-ray hybrid with missing refinement data (4xpv) (06/26/2018)
# 106. When generating map for rsr/edstat calculation, fix retry fallback with larger map size
# 107. If EDSTATSBINDIR is set, run edstats.pl from there and pass as -exec argument.
# 108. Correct edstats log file parsing for newer version of edstats
# 109. Python 3 compatibility (12/31/2019)
# 110. Code cleanup. For Phenix, if PDBx/mmCIF provided, use it. (04/10/2023)
#
# Examples:
# (Note: '-pdbid id'  is the same as '-pdb pdbfile -sf sffile')
# ln -s dcc-calc-v?.py dcc  (link script to dcc)
#
# dcc  -pdbid id        (Use refmac. if B_partial, use TLSANL to get B_full)
# dcc  -pdbid id -refmac  (refmac, the same as above)
# dcc  -pdbid id -refmac_tls  (Use refmac. if B_partial, use REFMAC to get B_full)
# dcc  -pdbid id -sfcheck  (Use sfcheck, if B_partial, use TLSANL to get B_full)
# dcc  -pdbid id -phenix_x  (Use phenix for xray, auto handle TLS)
# dcc  -pdbid id -phenix_n  (Use phenix for nutron, auto handle TLS)
# dcc  -pdbid id -phenix_xn  (Use phenix for xray & nutron, auto handle TLS)
# dcc  -pdbid id -phenix_xtriage  (Use phenix to test quality of SF)
# dcc  -mapcut   bigmap  pdbfile -o small_map_name  (small box for ligand map)
# dcc  -maplig   -pdb pdbfile -sf sffile  -o small_map_name (not finished)
# dcc  -pdbid id -test_r   (test R,by swaping I->F or F->I)
#
# If add -auto: calculate R by refine program. If failed, try refmac, ...
# If add -noed: only do the regular validation (no real space R and CC)
# If add -map:  both ccp4 and dsn6 maps will be generated
# If add -verb: intermediate files will be kept.
# If add -twin: force refmac to do twin validation
# If add -o:  (-o output_name)
#
#
# ##############################################################################
