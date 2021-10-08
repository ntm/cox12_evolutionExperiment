#!/usr/bin/perl

# 13/07/2018
# NTM

# Align paired-end trimmed FASTQs from $inDir, then sort the BAMs,
# then index them.
# We are aligning against several different reference genomes,
# see %refGenomes.
# We also add RG tags ("read groups", see $RG below) to the BAM, GATK 
# needs them downstream.
# Resulting files are produced in current dir.


use strict;
use warnings;


# key == ref genome file (with PATH), value == string added to output files
#my %refGenomes = ("../RefGenome/GCF_000146045.2_R64_genomic.fna.gz" => "S288C",
#		  "../RefGenome/CEN.PK2-1Ca_SGD_2015_JRIV01000000.fsa.gz" => "CENPK");
# Running again with the VEP/Ensembl version of S288, calling it VEP:
my %refGenomes = ("../RefGenome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz" => "VEP");

# read fastqs from $inDir
my $inDir = "../Fastq_PerSample_Trimmed/";
(-d $inDir) || die "indir $inDir is not a directory\n";


# Number of threads, both for bwa and for samtools sort
my $numThreads = 12;
# max memory per thread for samtools sort
my $memPerThread = "5G";


#############################################

opendir(INDIR, "$inDir") || die "cannot opendir $inDir";



foreach my $inFile (sort readdir(INDIR)) {
    ($inFile !~ /R1\.trimmed\.fastq\.gz$/) && (warn "we only want R1 fastqs, skipping infile $inFile\n") && next;
    ($inFile =~ /^([A-H])\.R1\.trimmed\.fastq\.gz$/) ||
	die "cannot grab sample id from filename $inFile\n";
    my $fileStart = $1;

    foreach my $refGenome (sort keys %refGenomes) {
	# this script aligns AND THEN sorts the BAM
	my $outFile = "$fileStart.".$refGenomes{$refGenome}.".sorted.bam";

	my $inFileMate = $inFile;
	($inFileMate =~ s/\.R1\./.R2./) || die "cannot make R2 file from infile $inFile\n";

	# read group info to add to BAM headers
	my $RG = '@RG\tID:cox12\tSM:'.$fileStart.'\tPL:Illumina';

	my $command = "bwa mem -t $numThreads -R \'$RG\' $refGenome $inDir/$inFile $inDir/$inFileMate | samtools sort -m $memPerThread -\@".$numThreads." -o $outFile -";
	# also build samtools index command
	my $indexCom = "samtools index $outFile";

	# skip BAM files that already exist, but still index if needed
	if (-e $outFile) {
	    warn "outfile $outFile already exists, skipping\n";
	    (-e "$outFile.bai") || ((warn "index missing, still indexing with: $indexCom\n") && system($indexCom));
	    next;
	}


	warn "\n###############################################\nStarting to run: $command\n";
	system($command);
	warn "\n####################\Indexing with: $indexCom\n";
	system($indexCom);
    }
}

