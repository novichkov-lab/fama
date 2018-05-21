#! /usr/bin/perl
use strict;
use warnings;

########################################################
###  Parse Spades output and filter list of contigs  ###
########################################################

my $fasta_file = "/mnt/data2/SEED/test_dir/test_comparative/assembly/contigs.fasta";
my $sam_file = "/mnt/data2/FEBA/nitrogen_v6/bowtie2_out.txt";
my $outfile = "/mnt/data2/SEED/test_dir/test_comparative/contigs_data.txt";
my $fasta_outfile = "/mnt/data2/SEED/test_dir/test_comparative/contigs_filtered.fasta";

# Parameters
my $coverage_cutoff = 3;
my $read_count_cutoff = 7;
my $length_cutoff = 300;

if (@ARGV == 4) {
	$fasta_file = $ARGV[0];
	$sam_file = $ARGV[1];
	$outfile = $ARGV[2];
	$fasta_outfile = $ARGV[3];
} else {
	print "Usage: perl filter_spades_contiglist.pl <fasta file name> <bowtie output> <output file> <output fasta file>\n";
#	exit(0);
};

if (!(-e $fasta_file)) {
	print "File $fasta_file not found!\n";
	exit(1);
}

if (!(-e $sam_file)) {
	print "File $sam_file not found!\n";
	exit(1);
}

my %sample_ids = ();
my %contigs = ();
my %coverage = ();
my %read_counts = ();

#read list of contigs
read_contigs($fasta_file);

#read Bowtie output file
read_sam_file($sam_file);

open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile");
print OUTFILE "#Contig\tNumber\tLength\tCoverage_spades\tCoverage_effective\tFiltered\t";
print OUTFILE join("\t", sort keys %sample_ids), "\tTotal\n";
foreach my $contig (sort keys %contigs){
	my $coverage_eff = 0;
	if (defined $coverage{$contig}) {
		$coverage_eff = $coverage{$contig}/$contigs{$contig}{"length"};
	}
	if ($coverage_eff < $coverage_cutoff){
		$contigs{$contig}{"filtered"} = "true";
	}
	if ((!exists $read_counts{$contig}{"Total"})||($read_counts{$contig}{"Total"} < $read_count_cutoff)){
		$contigs{$contig}{"filtered"} = "true";
	}
	if ($contigs{$contig}{"length"} < $length_cutoff){
		$contigs{$contig}{"filtered"} = "true";
	}
	print OUTFILE $contig . "\t" . $contigs{$contig}{"number"} . "\t" . 
					$contigs{$contig}{"length"} . "\t" . 
					$contigs{$contig}{"coverage"} . "\t" . 
					$coverage_eff . "\t" . 
					$contigs{$contig}{"filtered"} . "\t";
	foreach my $sample_id (sort keys %sample_ids){
		if (exists $read_counts{$contig}{$sample_id}){
			print OUTFILE $read_counts{$contig}{$sample_id} . "\t";
		} else {
			print OUTFILE "0\t";
		}
	}
	if (exists $read_counts{$contig}{"Total"}){
		print OUTFILE $read_counts{$contig}{"Total"};
	} else {
		print OUTFILE "0";
	}
	print OUTFILE "\n";
}
print OUTFILE "Unmapped reads\tN/A\tN/A\tN/A\t0\tfalse\t";
foreach my $sample_id (sort keys %sample_ids){
	if (exists $read_counts{"*"}{$sample_id}){
		print OUTFILE $read_counts{"*"}{$sample_id} . "\t";
	} else {
		print OUTFILE "0\t";
	}
}
if (exists $read_counts{"*"}{"Total"}){
	print OUTFILE $read_counts{"*"}{"Total"};
} else {
	print OUTFILE "0";
}
print OUTFILE "\n";
close OUTFILE;

my $flag = 0;
open (INFILE, $fasta_file) or die ("Unable to open file $fasta_file");
open (OUTFILE, ">$fasta_outfile") or die ("Unable to open file $fasta_outfile");
while (my $line = <INFILE>){
	if ($line =~ /^>/) {
		$flag = 0;
		my $contig = $line;
		chomp $contig;
		$contig =~ s/^>//;
		unless (exists $contigs{$contig}{"filtered"}){
			print "ERROR: Contig $contig not found.\n";
			exit(1);
		}
		unless ($contigs{$contig}{"filtered"} eq "true"){
			$flag = 1;
			print OUTFILE $line;
		}
	} elsif ($flag) {
		print OUTFILE $line;
	}
	
}


#######################
###   SUBROUTINES   ###
#######################

sub read_contigs{
	my ($infile) = @_;
	open (INFILE, $infile) or die ("Unable to open file $infile");
	while (my $line = <INFILE>){
		chomp $line;
		unless ($line) {
			next;
		}
		if ($line =~ /^>/){
			my $contig = $line;
			$contig =~ s/^>//;
			my (undef, $number, undef, $length, undef, $coverage) = split(/\_/,$contig);
			$contigs{$contig}{"number"} = $number;
			$contigs{$contig}{"length"} = $length;
			$contigs{$contig}{"coverage"} = $coverage;
			$contigs{$contig}{"filtered"} = "false";
		}
	}
	close INFILE;
}

sub read_sam_file{
	my ($infile) = @_;
	open (INFILE, $infile) or die ("Unable to open file $infile");
	while (my $line = <INFILE>){
		chomp $line;
		unless ($line) {
			next;
		}
		if ($line =~ /^@/){
			next;
		}
		my ($read, undef, $contig, undef, undef, undef, undef, undef, undef, $segment_seq) = split(/\t/, $line);
		if (exists $coverage{$contig}){
			$coverage{$contig} += length($segment_seq);
		} else {
			$coverage{$contig} = length($segment_seq);
		}
		my ($sample) = split(/_/, $line);
		$sample_ids{$sample} = 1;
		if (exists $read_counts{$contig}{"Total"}) {
			$read_counts{$contig}{"Total"} += 1;
		} else {
			$read_counts{$contig}{"Total"} = 1;
		}
		if (exists $read_counts{$contig}{$sample}){
			$read_counts{$contig}{$sample} += 1;
		} else {
			$read_counts{$contig}{$sample} = 1;
		}
	}
}

