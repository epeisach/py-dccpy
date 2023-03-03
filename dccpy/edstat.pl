#!/usr/bin/perl -w
if(!@ARGV)
{ use File::Basename;
  print "
Usage:  ${\basename($0)}  <required arguments>  <options>

Required arguments:

-hklin   filename.mtz  - Input MTZ file (output from a refinement).
-xyzin   filename.pdb  - Input PDB file (output by the _same_ refinement
                         run that produced the MTZ file).

Optional arguments:

-flabel  <string>      - Column label for Fobs (by default tries F, FP,
                         FOSC, F_xxx, FP_xxx where xxx = any string).
-grid    <integer>     - Grid sampling (default = 4.5; must be >= 4).
-byatom                - Use by-atom averages (default is by-residue).
-select  <string>      - Residue selection (see program documentation
                         for details, except enclosing quotes are not
                         needed here unless the string contains spaces).
-bfac    <string>      - Statistic to place in B factor column of output
                         PDB file: R, RG, SRG, CCS, CCP, ZCC, ZD-, ZD+,
                         ZD, ZO or NONE (default = NONE).
-occu    <string>      - Statistic to place in occupancy column of
                         output PDB file (options as for -bfac except
                         default = ZD-).
-log     <filename>    - Log file (default is standard output).
-output  <filename>    - Output data file (default = 'edstats.out').
-hisout  <filename>    - Output density histogram plot (default = none).
-ppdout  <filename>    - Output P-P difference plot (default = none).
-qqdout  <filename>    - Output Q-Q difference plot (default = none).
-mapout1 <filename>    - Output rescaled rho_obs map (default = none).
-mapout2 <filename>    - Output rescaled rho_diff map (default = none).
-xyzout  <filename>    - Output PDB file with atom Z scores in occupancy
                         column (default = none).
-tmpdir  <dir>         - Directory for temporary files
                         (default = '/tmp/\$USER/').
-exec    <path>        - Path for executables (default = none).
-test    <integer>     - Test level (default = none).
-noerror               - Use input MTZ file after error in MTZFIX.
-debug                 - Debug (don't delete temp dir).

The columns and column labels in the output file are:

 1. RT:   Residue type (3-letter code).
 2. CI:   Chain ID (reset to underscore if blank).
 3. RN:   Residue number (including insertion code if present).

Then for the main-chain atoms in the residue (including CB):
 4. BAm:  Weighted average Biso.
 5. NPm:  No of statistically independent grid points covered by atoms.
 6. Rm:   Real-space R factor (RSR).
 7. RGm:  Real-space RG factor (RSRG).
 8. SRGm: Standard uncertainty of RSRG.
 9. CCSm: Real-space sample correlation coefficient (RSCC).
10. CCPm: Real-space 'population' correlation coefficient (RSPCC).
11. ZCCm: Z-score of real-space sample correlation coefficient.
12. ZOm:  Real-space Zobs score (RSZO).
13. ZDm:  Real-space Zdiff score (RSZD) i.e. max(-RSZD-,RSZD+).
14. ZD-m: Real-space Zdiff score for negative differences (RSZD-).
15. ZD+m: Real-space Zdiff score for positive differences (RSZD+).

Columns 16-27 (with labels 'BAs' etc.) contain the same information as
columns 4-15 above (i.e. add 12) but for the side-chain atoms
(excluding CB) if present.  Columns 28-39 (with labels 'BAa' etc.
contain the corresponding all-atom (main + side-chain) statistics.

";
  exit;
}

# Default settings.
$opt_grid = 4.5;
$opt_output = 'edstats.out';
$opt_exec = $opt_noerror = '';

# GetOptions zaps the argument list, so make a copy.
@args = @ARGV;

use Getopt::Long;

GetOptions('hklin=s','xyzin=s','flabel=s','grid=i','byatom','select=s',
  'bfac=s','occu=s','log=s','output=s','hisout:s','ppdout:s','qqdout:s',
  'mapout1:s','mapout2:s','xyzout:s','tmpdir=s','exec=s','test:i',
  'noerror','debug') or
  die 'ERROR in processing command line arguments,';

# If debugging echo command line for reassurance.
if($opt_debug)
{ $command = $0;
  $command .= /^$|[\s\#]/?" '$_'":" $_" for @args;
  print "\nCommand line: $command\n\n";
}

# Check we have CCP4 and all input files.
$ENV{CCP4} or die 'ERROR: you must set up the CCP4 environment,';

$opt_hklin && $opt_xyzin or
  die 'ERROR: MTZ and/or PDB file not specified,';

if(!-e $opt_hklin)
{ if(-e "$opt_hklin.mtz") {$opt_hklin .= '.mtz'}
  else {die "ERROR: input MTZ file $opt_hklin not found,"}
}

if(!-e $opt_xyzin)
{ if(-e "$opt_xyzin.pdb") {$opt_xyzin .= '.pdb'}
  else {die "ERROR: input PDB file $opt_xyzin not found,"}
}

# $opt_exec will be the full path to the EDSTATS executable (filename =
# 'edstats'); $path is the directory part for MTZINFO & MTZFIX only.
if(!($path = $opt_exec)) {$opt_exec = 'edstats'}
elsif(-d $path) {$opt_exec =~ s=/?$=/edstats=}
else {$path =~ s=[^/]+$==}

# Make sure executables can be found somewhere.
use Env '@PATH';

$mtzinfo='mtzinfo';
grep(-x "$_/mtzinfo",@PATH) || -x ($mtzinfo="${path}mtzinfo") or
  die "ERROR: mtzinfo executable '${path}mtzinfo' not found in path,";

$mtzfix='mtzfix';
grep(-x "$_/mtzfix",@PATH) || -x ($mtzfix="${path}mtzfix") or
  die "ERROR: mtzfix executable '${path}mtzfix' not found in path,";

grep(-x "$_/$opt_exec",@PATH) || -x $opt_exec or
  die "ERROR: edstats executable '$opt_exec' not found in path,";

# Options for MTZFIX.
$opt_flabel = $opt_flabel?" FLABEL $opt_flabel":'';

# Options for FFT.
$opt_grid=~/^\d+(\.\d*)?$/ && $opt_grid>=4 or
  die 'ERROR: grid sampling must be >= 4,';

# Options for EDSTATS.
$opt_byatom = $opt_byatom?',MAIN=atom,SIDE=atom':'';
$opt_select = $opt_select?",SELE=\\'$opt_select\\'":'';
$opt_bfac = $opt_bfac?",BFAC=$opt_bfac":'';
$opt_occu = $opt_occu?",OCCU=$opt_occu":'';

$opt_hisout = 'edstats.his' if defined($opt_hisout) && !$opt_hisout;
$opt_hisout = $opt_hisout?" HISOUT $opt_hisout":'';

$opt_ppdout = 'edstats.ppd' if defined($opt_ppdout) && !$opt_ppdout;
$opt_ppdout = $opt_ppdout?" PPDOUT $opt_ppdout":'';

$opt_qqdout = 'edstats.qqd' if defined($opt_qqdout) && !$opt_qqdout;
$opt_qqdout = $opt_qqdout?" QQDOUT $opt_qqdout":'';

$opt_mapout1 = 'edstats-fo.map' if defined($opt_mapout1) &&
  !$opt_mapout1;
$opt_mapout1 = $opt_mapout1?" MAPOUT1 $opt_mapout1":'';

$opt_mapout2 = 'edstats-df.map' if defined($opt_mapout2) &&
  !$opt_mapout2;
$opt_mapout2 = $opt_mapout2?" MAPOUT2 $opt_mapout2":'';

$opt_xyzout = 'edstats.pdb' if defined($opt_xyzout) && !$opt_xyzout;
$opt_xyzout = $opt_xyzout?" XYZOUT $opt_xyzout":'';

$opt_log = $opt_log?">$opt_log":'';

# Create temporary files directory if it doesn't exist.
if($opt_tmpdir)
{ if(-d $opt_tmpdir) {$keep = 1}
  else {mkdir $opt_tmpdir,0755 or die "ERROR in mkdir $opt_tmpdir,"}
}
else
{ $opt_tmpdir = "/tmp/$ENV{USER}";
  mkdir $opt_tmpdir,0755 or die "ERROR in mkdir $opt_tmpdir,"
    if !-e $opt_tmpdir;
  use File::Temp qw(tempdir);
  $opt_tmpdir = tempdir('edstats-XXXXXX',DIR=>$opt_tmpdir);
}

# Get column labels & resolution limits from the MTZ header.
open IN,"$mtzinfo $opt_hklin |";

while(<IN>)
{ print;
  if(/^LABELS/i) {@a = split}
  elsif(/^XDATA/i) {($rl,$rh) = (split)[7,8]}
}
close IN;

# Set up column labels for FFT.
if(grep(/2FOFCWT/,@a) && grep(/PH2FOFCWT/,@a) && grep(/FOFCWT/,@a) &&
  grep(/PHFOFCWT/,@a))
  {($fo,$po,$fd,$pd) = ('2FOFCWT','PH2FOFCWT','FOFCWT','PHFOFCWT')}

elsif(grep(/FWT/,@a) && grep(/PHWT/,@a) && grep(/DELFWT/,@a) &&
  grep(/PHDELWT/,@a))
  {($fo,$po,$fd,$pd) = ('FWT','PHWT','DELFWT','PHDELWT')}

else {die 'ERROR - standard map coefficient labels not found,'}

# Fix up the map coefficients: FLABEL specifies the label for Fobs &
# sigma(Fobs) (defaults are F/SIGF or FOSC/SIGFOSC or FP/SIGFP or
# F_xxx/SIGF_xxx or FP_xxx/SIGFP_xxx).  Here, 'in.mtz' is the output
# reflection file from the refinement program in MTZ format.
!system "rm -f $opt_tmpdir/mtzfix.mtz; $mtzfix$opt_flabel".
  " HKLIN $opt_hklin HKLOUT $opt_tmpdir/mtzfix.mtz".
  " >$opt_tmpdir/mtzfix.log" or $opt_noerror or die 'ERROR in mtzfix,';
  
# If no fix-up was needed, or error in MTZFIX, use the original file.
if(!-e "$opt_tmpdir/mtzfix.mtz")
{ if($opt_hklin!~/^\//)
  { use Cwd;
    $opt_hklin = &cwd.'/'.$opt_hklin;
  }
  symlink $opt_hklin,"$opt_tmpdir/mtzfix.mtz";
}

# Compute the 2mFo-DFc map (mFo for centrics).  Note that EDSTATS needs
# only 1 asymmetric unit (but will also work with more).  Grid sampling
# must be at least 4.
!system "echo 'LABI F1=$fo PHI=$po\nXYZL asu\nGRID samp $opt_grid'".
  " |fft HKLIN $opt_tmpdir/mtzfix.mtz MAPOUT $opt_tmpdir/fo.map".
  " >$opt_tmpdir/fft-fo.log" or die 'ERROR in fft(1),';

# Compute the 2(mFo-DFc) map (mFo-DFc for centrics).
!system "echo 'LABI F1=$fd PHI=$pd\nXYZL asu\nGRID samp $opt_grid'".
  " |fft HKLIN $opt_tmpdir/mtzfix.mtz MAPOUT $opt_tmpdir/df.map".
  " >$opt_tmpdir/fft-df.log" or die 'ERROR in fft(2),';

# Test option for EDSTATS.
$opt_test = -1 if defined($opt_test) && !$opt_test;
$opt_test = $opt_test?",TEST=$opt_test":'';

# Main- & side-chain residue statistics.
!system "echo RESL=$rl,RESH=$rh$opt_byatom$opt_select$opt_bfac".
  "$opt_occu$opt_test |$opt_exec XYZIN $opt_xyzin".
  " MAPIN1 $opt_tmpdir/fo MAPIN2 $opt_tmpdir/df OUT $opt_output".
  "$opt_hisout$opt_ppdout$opt_qqdout$opt_mapout1$opt_mapout2".
  "$opt_xyzout$opt_log" or die "ERROR in $opt_exec,";

# Delete temp directory if not required.
$opt_debug || $keep or system "rm -fr $opt_tmpdir";
