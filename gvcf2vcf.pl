#!/usr/bin/perl

# NTM
# 30/08/2018

# Parses on stdin a GVCF file with one or more sample data columns;
# prints to stdout a VCF file where we:
# - skip lines where the only ALT allele is <NON_REF>
# - replace FORMAT with GT:DP:AF
# - data columns will have for GT:
#   * ./. if DP < $minDP or if we would have called NON_REF, this is a NOCALL;
#   * x/x if a single allele (whether REF or ALT, but not NON_REF) has 
#     read frequency >= $minFracToCall, in which case AF has this single value;
#   * x/y where x has the highest AF (first number in AF) and y has the 
#     second highest (second AF value) and it's >= $minFracToCall.
# - data columns will have zero (geno is ./.) or one int for DP, 
#   and zero (geno is ./.), one (geno is x/x) or two (geno is x/y) 
#   ints for AF, these are the allele freqs as percentages (ignoring decimals)
#
# Any ALT allele that doesn't appear in a GT field is removed from
# the ALT column.
# <NON_REF> is also removed from ALTs, but if it was called this will be NOCALL.
# If the ALT column is then empty (ie every sample is REF or NOCALL),
# the line is skipped (but counted).
#
# Since data columns A-D are clonal haploid strains, they should be 
# HOMOREF or HOMOVAR for one ALT. So, any line where one of samples
# A-D is not x/x is suspicious.

# Takes as args:
# --minDP int (must have DP >= $minDP to make call, eg 10),
# --strandDisc [NOT IMPLEMENTED YET], 
# --minFracToCall int (eg 15, this is a percentage)
#


use strict;
use warnings;
use Getopt::Long;

# default values for args:
# min DP to make a call
my $minDP = 10;
# min read freq to call an allele (as percentage)
my $minFracToCall = 15;

GetOptions ("minDP=i" => \$minDP,
	    "minFracToCall=i" => \$minFracToCall)
    or die("Error in command line arguments\n");

# parse header, just copy it except for FORMAT lines, which we fix
while(my $line = <STDIN>) {
    if ($line =~ /^##FORMAT/) {
	# FORMAT: copy GT and DP
	if (($line =~ /^##FORMAT=<ID=GT,/) || ($line =~ /^##FORMAT=<ID=DP,/)) {
	    print $line;
	}
	# skip other FORMAT lines: NOOP
    }
    elsif ($line =~ /^##/) {
	print $line;
    }
    elsif ($line =~ /^#CHROM/) {
	# add FORMAT definition for new AF field
	print "##FORMAT=<ID=AF,Number=.,Type=Integer,Description=\"Fraction of reads supporting the allele(s) present in GT, in the same order\">\n";
	# add info with full command line run
	my $com = qx/ps -o args= $$/;
	chomp($com);
	$com .= " < ".`readlink -f /proc/$$/fd/0` ;
	chomp($com);
	$com .= " > ".`readlink -f /proc/$$/fd/1` ;
	chomp($com);
	$com .= " 2> ".`readlink -f /proc/$$/fd/2` ;
	chomp($com);
	print "##gvcf2vcf=<commandLine=\"$com\">\n";
	print $line;
	last;
    }
    else {
	die "E: parsing header, found bad line:\n$line";
    }
}

# counters for log
# lines skipped because no ALT allele
my $skippedNoAlt = 0;
# lines skipped because all samples have same HOMOZYG genotype
my $skippedNoVar = 0;
# lines with a sample where no allele has sufficient readFrac to call
my $warnedBadFrac = 0;
# samples where NONREF is the most frequent allele, it's a nocall
my $nonRefFirst = 0;
# samples where NONREF is the second highest frac and would have been called
# (we call HOMO for the first allele, but log and count)
my $warnNonRefSecond = 0;

# now parse data lines
 LINE: while(my $line = <STDIN>) {
     chomp($line);
     my @data = split(/\t/, $line);
     (@data == 17) || die "not 8 sample columns in line?\n$line\n";
     my @toPrint;
     # first 4 fields are just copied
     foreach my $i (1..4) {
	 push(@toPrint, shift(@data));
     }
     # col 5 is ALT
     my $alts = shift(@data);
     # skip line if only ALT is NON_REF
     if ($alts eq "<NON_REF>") {
	 $skippedNoAlt++;
	 next LINE;
     }
     #otherwise save alts, we will only print the ones that are called
     my @alts = split(/,/,$alts);
     # make sure NON_REF is last, we use this later (look for $#ADs)
     ($alts[$#alts] eq "<NON_REF>") || die "line where NON_REF isnt the last ALT, WTF?\n$line\n";
     # for now push empty ALT value
     push(@toPrint, "");
     # fields 6-8 are copied
     foreach my $i (6..8) {
	 push(@toPrint, shift(@data));
     }
     # col 9 is FORMAT: replace with ours
     my $format = shift(@data);
     push(@toPrint, "GT:DP:AF");
     # %format: key is a FORMAT key (eg GQ), value is the index of that key in $format
     my %format;
     { 
	 my @format = split(/:/, $format);
	 foreach my $i (0..$#format) {
	     $format{$format[$i]} = $i ;
	 }
     }
     # sanity: make sure the fields we need are there
     ((defined $format{"GT"}) && ($format{"GT"}==0)) || 
	 die "no GT key or GT not first in FORMAT string for line:\n$line\n";
     (defined $format{"DP"}) || die "no DP key in FORMAT string for line:\n$line\n";
     (defined $format{"AD"}) || die "no AD key in FORMAT string for line:\n$line\n";

     # @allelesCalled: for each allele in (REF,ALTs), value 1 iff allele is called
     # (so, the indexes of @allelesCalled correspond to the numbers in GT)
     my @allelesCalled = (0) x (1+scalar(@alts));

     # now deal with actual data fields
     foreach my $sample (1..8) {
	 my $data = shift(@data);
	 # examine content and apply filters
	 my @thisData = split(/:/, $data) ;

	 if ( (! defined $thisData[$format{"DP"}]) || ($thisData[$format{"DP"}] eq '.') || 
	      ($thisData[$format{"DP"}] < $minDP)  || ($thisData[$format{"DP"}] == 0) ) {
	     # DP undefined or too low (or zero when $minDP==0 must still be skipped), 
	     # change to NOCALL
	     push(@toPrint, "./.");
	     next;
	 }

	 else {
	     # coverage is good, make call ie find the 2 biggest ADs (and their 
	     # indexes in @ADs)
	     my @ADs = split(/,/, $thisData[$format{"AD"}]);
	     my $maxAdIndex = 0;
	     my $maxAd = $ADs[0];
	     my $secondAdIndex;
	     my $secondAd = -1;
	     foreach my $i (reverse (1..$#ADs)) { # reverse because we want to start with NON_REF
		 # it's super important to start with NON_REF, so if ADs are equal it sticks
		 if ($ADs[$i] > $maxAd) {
		     $secondAd = $maxAd;
		     $secondAdIndex = $maxAdIndex;
		     $maxAd = $ADs[$i];
		     $maxAdIndex = $i;
		 }
		 elsif ($ADs[$i] > $secondAd) {
		     $secondAd = $ADs[$i];
		     $secondAdIndex = $i;
		 }
	     }
	     my $fracReadsMax = int(100 * $maxAd / $thisData[$format{"DP"}]) ;
	     my $fracReadsSecond = int(100 * $secondAd / $thisData[$format{"DP"}]) ;
	     if ($maxAdIndex == $#ADs) {
		 # most frequent allele is NON_REF, this is a NOCALL
		 push(@toPrint, "./.");
		 $nonRefFirst++;
	     }
	     elsif (($fracReadsSecond >=  $minFracToCall) && ($secondAdIndex < $#ADs)) {
		 # can call 2 "real" (not NONREF) alleles here
		 push(@toPrint, "$maxAdIndex/$secondAdIndex:".$thisData[$format{"DP"}].
		      ":$fracReadsMax,$fracReadsSecond");
		 # save called alleles
		 $allelesCalled[$maxAdIndex] = 1;
		 $allelesCalled[$secondAdIndex] = 1;
		 next;
	     }
	     elsif ($fracReadsSecond >=  $minFracToCall) {
		 # second most frequent allele is above cutoff but it's NON_REF
		 # call it HOMO for the first allele but log
		 push(@toPrint, "$maxAdIndex/$maxAdIndex:".$thisData[$format{"DP"}].":$fracReadsMax");
		 $allelesCalled[$maxAdIndex] = 1;
		 # too many lines and seems legit, don't log the lines just count them
		 # warn "CALLING HOMO BUT NONREF IS SECOND AND ABOVE CUTOFF (sample $sample, fracReadsSecond $fracReadsSecond): $line\n";
		 $warnNonRefSecond++;
	     }
	     elsif ($fracReadsMax >= $minFracToCall) {
		 # single allele (HOMOZYG) called here
		 push(@toPrint, "$maxAdIndex/$maxAdIndex:".$thisData[$format{"DP"}].":$fracReadsMax");
		 # save called allele
		 $allelesCalled[$maxAdIndex] = 1;
		 next;
	     }
	     else {
		 # BAD sample, DP sufficient but no single allele is >= $minFrac!
		 # NOCALL this sample but log as BAD
		 push(@toPrint, "./.");
		 warn "BAD SAMPLE (sample $sample , fracReadsMax $fracReadsMax): $line\n";
		 $warnedBadFrac++;
		 next;
	     }
	 }
     }


     # do we have at least 2 different called alleles?
     my $numCalledAlleles = 0;
     foreach my $i (0..$#allelesCalled) {
	 ($allelesCalled[$i]) && ($numCalledAlleles++);
     }
     if ($numCalledAlleles < 2) {
	 # all samples are noCall or HOMO for a single allele (possibly REF),
	 # skip line but count
	 $skippedNoVar++;
	 next LINE;
     }

     # build list of ALTs to keep
     my @newAlts = ();
     # $allelesCalled[$i] becomes the new number to use in GT for that allele
     # (ie 0 for REF and 1+index in @newAlts for ALTs)
     $allelesCalled[0] = 0;
     foreach my $i (0..$#alts) {
	 if ($allelesCalled[1+$i]) {
	     push(@newAlts, $alts[$i]);
	     $allelesCalled[1+$i] = 1 + $#newAlts;
	 }
     }

     # fix ALT column (index 4)
     $toPrint[4] = join(',', @newAlts);

     # adjust all GTs (cols indexed 9-16)
     foreach my $i (9..16) {
	 if ($toPrint[$i] eq './.') {
	     # NOCALL, nothing to do
	     next;
	 }
	 elsif ($toPrint[$i] =~ /^(\d+)\/(\d+)(:.*)$/) {
	     # called, fix geno
	     my ($g1,$g2,$rest) = ($1,$2,$3);
	     $toPrint[$i] = $allelesCalled[$g1]."/".$allelesCalled[$g2].$rest;
	 }
	 else {
	     die "cannot decompose data in col $i: ".$toPrint[$i]."\nLINE IS:$line\n";
	 }
     }

     # OK all done, print line
     print join("\t", @toPrint)."\n";
}

# log counts of skipped lines
warn "\nLines skipped because no ALT allele: $skippedNoAlt\n";
warn "\nLines skipped because all samples are NOCALL or HOMO for a single allele: $skippedNoVar\n";
warn "\nLines with a BAD sample (good cov but all readFracs low, NOCALL'ed but logged): $warnedBadFrac\n";
warn "\nCalls where NONREF is first, we NOCALL'd but counted them: $nonRefFirst\n";
warn "\nCalls where NONREF is second and would have been called: $warnNonRefSecond\n";


