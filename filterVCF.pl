#!/usr/bin/perl


# 19/09/2018
# NTM

# UPDATE 15/10/2018: making a single unified script that can do 
# everything that the previous 5 or 6 scripts in Old/ did,
# using GetOptions.
# See @good* comments below, and examples of usage in the README.
#
# UPDATE 18/10/2018: allow multi-allelic sites.
#
# UPDATE 31/10/2018: add --onesNotRef option.
#
# Parse on stdin a VCF file, with FORMAT and data columns as produced
# by gvcf2vcf.pl.
# Find good candidates variants (either for parent strain 1 or 2).
#
# - all @goodZeroes samples must be x/x with some x (or ./. if $allowNocalls)
# - all @goodOnes samples must be y/y with y>=1, y!=x (or ./. if $allowNocalls)
# - all @goodLow samples must be NOCALL, x/x, or HET x/? with freq(x) >= 1-$goodLowFreq
#   (this is a relaxed version of goodZeroes)
# - all @goodMedium samples must be NOCALL or have freq(x) in
#   [$goodMinMedium,100-$goodMinMedium] , ie usually HET
#
# Strains are as follows:
# A : WT
# B: delta cox12
# UPDATE 12/10/18 C and D are actually clones from the backcross experiments,
# ie C is one if the clones in pool E and D is one of the clones in pool G.
# This doesn't change the way we look for good candidates though.
# C: evolved parent strain 1 (clonal)
# D: evolved parent strain 2 (clonal)
# E: back-crossed strain 1, growing (pool of 15 haploid samples)
# F: back-crossed strain 1, not growing (pool of 8 haploid samples)
# G: back-crossed strain 2, growing (pool of 15 haploid samples)
# H: back-crossed strain 2, not growing (pool of 15 haploid samples)
#
# So in strain 1 we look for variants with:
# 0% in A, B, D
# 100% in C
# 0% in F G H
# 100% in E.
# Similarly in strain 2 we want:
# 0% in A B C E F H
# 100% in D G
# This is all expressed with:
# @goodZeroes = (0,1,3,5,6,7) and @goodOnes = (2,4) for strain 1;
# @goodZeroes = (0,1,2,4,5,7) and @goodOnes = (3,6) for strain 2.
# UPDATE 19/10/2018: args for samples can now be A-H or a-h

use strict;
use warnings;
use Getopt::Long;

# Parental strain we are interested in: 1 or 2, leave blank if we don't care
my $goodStrain = '';
# relax stringency: if true we allow samples in @goodZeroes and @goodOnes to
# be NOCALLs, but in any case there must be at least one true homo call of each class
my $allowNocalls = '';
# by default "ones" can be REF if some other strain is non-REF; if $onesNotRef
# is true, impose "ones" to be non-ref.
my $onesNotRef = '';
# Each @good* array contains a list of samples == data column indexes (starting at 0)
my @goodZeroes = ();
my @goodOnes = ();
my @goodLow = ();
my $goodLowFreq = 25;
my @goodMedium = ();
my $goodMinMedium = 20;

# temp strings to produce @good* arrays, for GetOptions
my ($zero,$one,$low,$medium) = ("","","","");
GetOptions ("strain=i" => \$goodStrain,
	    "nocalls" => \$allowNocalls,
	    "onesNotRef" => \$onesNotRef,
	    "zero=s" => \$zero,
	    "one=s" => \$one,
	    "low=s" => \$low,
	    "medium=s" => \$medium)
    or die("Error in command line arguments\n");
# strain must be 1 or 2 if defined
($goodStrain eq '') || ($goodStrain == 1) || ($goodStrain == 2) || 
    die "if you specifiy --strain it must be 1 or 2\n";
# translate sample IDs A-H or a-h into indexes
$zero =~ tr/A-H/0-7/; $zero =~ tr/a-h/0-7/;
$one =~ tr/A-H/0-7/; $one =~ tr/a-h/0-7/;
$low =~ tr/A-H/0-7/; $low =~ tr/a-h/0-7/;
$medium =~ tr/A-H/0-7/; $medium =~ tr/a-h/0-7/;
# fill arrays
@goodZeroes = split(/,/,$zero);
@goodOnes = split(/,/,$one);
@goodLow = split(/,/,$low);
@goodMedium = split(/,/,$medium);
{
    # sanity : make sure we don't have the same sample twice
    my @seen = (0) x 8;
    foreach my $thisCol (@goodZeroes, @goodOnes, @goodLow, @goodMedium) {
	($thisCol =~ /^\d$/) || die "argument $thisCol is not a sample ID A-H or a-h, try again\n";
	(($thisCol >= 0) && ($thisCol <= 7)) || die "argument $thisCol is not in [0;7], impossible?\n";
	($seen[$thisCol]) && die "Sample $thisCol present in several arguments, try again!\n";
	$seen[$thisCol] = 1;
    }
}

# skip header
while(my $line = <STDIN>) {
    if ($line =~ /^##/) {
    }
    elsif ($line =~ /^#CHROM/) {
	chomp($line);
	# add info with full command line run
	my $com = qx/ps -o args= $$/;
	chomp($com);
	$com .= " < ".`readlink -f /proc/$$/fd/0` ;
	chomp($com);
	$com .= " > ".`readlink -f /proc/$$/fd/1` ;
	chomp($com);
	# I don't want stderr here
	#$com .= " 2> ".`readlink -f /proc/$$/fd/2` ;
	#chomp($com);
	print "$line\tcommandLine=\"$com\"\n";
	last;
    }
    else {
	die "E: parsing header, found bad line:\n$line";
    }
}

# now parse data lines
 LINE: while(my $line = <STDIN>) {
     chomp($line);
     my @data = split(/\t/, $line);
     (@data == 17) || die "not 8 sample columns in line?\n$line\n";

     # keep this line if it is a good candidate for $goodStrain
     my $good = 1;
     # allele numbers x and y: usually x should be 0 and y 1, but for the
     # multiallelic cases we need to grab x and y in the data
     # default values must be different and negative, so it works if
     # no goodZeroes or goodOnes are requested or found
     my $zeroAllele = -2;
     my $oneAllele = -1;
     # number of zeroes found in the @goodZeroes (others may have been ./.)
     my $zeroesFound = 0;
     # same for number of ones
     my $onesFound = 0;
     foreach my $i (@goodZeroes)  {
	 my $thisData = $data[9+$i];
	 # accept nocalls but don't increment counter
	 if ($thisData eq './.') {
	     next;
	 }
	 elsif ($thisData =~ m~^(\d)/\1:~ ) {
	     # we are HOMO
	     my $thisAllele = $1;
	     # if this is the first non-noCall among goodZeroes, this defines allele x
	     ($zeroAllele >= 0) || ($zeroAllele = $thisAllele);
	     if ($thisAllele == $zeroAllele) {
		 $zeroesFound++;
	     }
	     else {
		 # HOMOZYG but not with the zeroAllele
		 $good=0;
		 last;
	     }
	 }
	 else{
	     # HET can't be goodZero
	     $good=0;
	     last;
	 }
     }

     foreach my $i (@goodOnes)  {
	 my $thisData = $data[9+$i];
	 # accept nocalls but don't increment counter
	 if ($thisData eq './.') {
	     next;
	 }
	 elsif ($thisData =~ m~^(\d)/\1:~ ) {
	     # we are HOMO
	     my $thisAllele = $1;
	     # if this is the first non-noCall among goodOnes it defines allele y
	     ($oneAllele >= 0) || ($oneAllele = $thisAllele);
	     # we can't have x==y
	     if ($thisAllele == $zeroAllele) {
		 $good=0;
		 last;
	     }
	     elsif ($onesNotRef && ($thisAllele == 0)) {
		 # onesNotRef mode: don't allow "one" allele to be REF
		 $good=0;
		 last;
	     }
	     elsif ($thisAllele == $oneAllele) {
		 $onesFound++;
	     }
	     else {
		 # HOMOZYG but not with the oneAllele
		 $good=0;
		 last;
	     }
	 }
	 else{
	     # HET can't be goodOne
	     $good=0;
	     last;
	 }
     }

     # for goodLow and goodMedium we will consider that x is 0 if we didn't
     # set $zeroAllele above, except if $oneAllele is 0 (in which case
     # we arbitrarily set $zeroAllele=1)
     if ($zeroAllele < 0) {
	 $zeroAllele = 0;
	 ($oneAllele==0) && ($zeroAllele=1);
     }

     foreach my $i (@goodLow) {
	 my $thisData = $data[9+$i];
	 if ($thisData eq './.') {
	     # accept nocalls
	     next;
	 }
	 elsif ($thisData =~ m~^$zeroAllele/$zeroAllele:~) {
	     # HOMO for zeroAllele, accept
	     next;
	 }
	 elsif ($thisData =~ m~^$zeroAllele/\d:\d+:(\d+),\d+$~) {
	     # HET with zeroAllele major, examine its AF
	     if ($1 >= 100 - $goodLowFreq) {
		 next;
	     }
	     else {
		 # freq(x) not high enough
		 $good = 0;
		 last;
	     }
	 }
	 else {
	     # other cases can't be good
	     $good = 0;
	     last;
	 }
     }

     foreach my $i (@goodMedium) {
	 my $thisData = $data[9+$i];
	 if ($thisData eq './.') {
	     # accept nocalls
	     next;
	 }
	 elsif ( ($thisData =~ m~^$zeroAllele/\d:\d+:(\d+),\d+$~) ||
		 ($thisData =~ m~^\d/$zeroAllele:\d+:\d+,(\d+)$~) ) {
	     # HET with zeroAllele, in both cases $1 is its AF
	     if (($1 >= $goodMinMedium) && ($1 <= 100 - $goodMinMedium)) {
		 next;
	     }
	     else {
		 # freq(x) not medium enough
		 $good = 0;
		 last;
	     }
	 }
	 else {
	     # other cases can't be good
	     $good = 0;
	     last;
	 }
     }

     if ($good) {
	 # might be good...
	 if ($allowNocalls) {
	     # at least one non-nocall if any zeroes or ones were requested
	     (! @goodZeroes) || ($zeroesFound) || ($good=0);
	     (! @goodOnes) || ($onesFound) || ($good=0);
	 }
	 else {
	     # no tolerance, everything requested must have been found
	     (@goodZeroes == $zeroesFound) || ($good=0);
	     (@goodOnes == $onesFound) || ($good=0);
	 }
     }
     # are we still goood?
     if ($good) {
	 # found a good candidate! if $goodStrain add at start of INFO
	 if ($goodStrain) {
	     my $oldInfo = $data[7] ;
	     ($oldInfo =~ /GoodForStrain/) &&
		 die "WTF, current variant is good for both strains? impossible!\n$line\n";
	     my $newInfo = "GoodForStrain=$goodStrain";
	     if ($oldInfo ne '.') {
		 $newInfo .= ";".$oldInfo ;
	     }
	     $data[7] = $newInfo;
	 }
	 print join("\t",@data)."\n";
     }
}
