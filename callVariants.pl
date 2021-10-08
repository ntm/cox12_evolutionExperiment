#!/usr/bin/perl

# 19/07/2018
# PKS
# Launching GATK Haplotypecaller command.
# Generates g.vcf file containing raw variants from ready to analyze-bam file.
# resulting files are produced in current directory.


use strict;
use warnings;


# key == string present in BAM files (will also be present in the created GVCFs),
# value == ref genome file (with PATH)
my %refGenomes = ("S288C" => "../RefGenome/GCF_000146045.2_R64_genomic.fna",
		  "CENPK" => "../RefGenome/CEN.PK2-1Ca_SGD_2015_JRIV01000000.fasta",
		  "VEP" => "../RefGenome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa");

# full path to GATK jar file
my $gatk = "/home/nthierry/Software/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar";

# read ready to analyze-bam file from $inDir
my $inDir = "../AlignFastqsAndSortBAMs/";
(-d $inDir) || die "indir $inDir is not a directory\n";

opendir(INDIR, "$inDir") || die "cannot opendir $inDir";


foreach my $inFile (sort readdir(INDIR)) {
    # silently skip bai files and files starting with .
    (($inFile =~ /sorted\.bam\.bai$/) || ($inFile =~ /^\./)) && next;
    ($inFile !~ /sorted\.bam$/) && (warn "we only want bam files, skipping infile $inFile\n") && next;    
    ($inFile =~ /^([A-H])\.(\w+)\.sorted\.bam$/) ||
	die "cannot grab sample id or strain name from filename $inFile\n";
    my ($sample,$strain) = ($1,$2);

    (defined $refGenomes{$strain}) || die "don't have a refGenome entry for strain $strain\n";
    my $refGenome = $refGenomes{$strain};
    my $outFile = "$sample.$strain.g.vcf";

    if (-e $outFile) {
	warn "outFile $outFile already exists, skipping inFile $inFile\n";
	next;
    }

    my $command = "java -Xmx16G -jar $gatk -T HaplotypeCaller -S LENIENT -R $refGenome -I $inDir/$inFile -ERC BP_RESOLUTION -o $outFile -dt NONE -nct 8 -A StrandAlleleCountsBySample";

    warn "will run: $command\n";
    system($command);
}
