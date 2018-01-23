#! /usr/bin/perl
use strict;
use warnings;

########################################################
###  Parse Spades output and filter list of contigs  ###
########################################################

my $fasta_file = "/mnt/data2/SEED/test_dir/test_comparative/assembly/contigs.fasta";
my $outfile = "/mnt/data2/SEED/test_dir/test_comparative/contigs_filtered1.fasta";

# Parameters
my $length_cutoff = 300;

if (@ARGV == 2) {
	$fasta_file = $ARGV[0];
	$outfile = $ARGV[1];
} else {
	print "Usage: perl filter_contigs_by_length.pl <fasta file name> <output fasta file>\n";
#	exit(0);
};

if (!(-e $fasta_file)) {
	print "File $fasta_file not found!\n";
	exit(1);
}

my @seq = ();
my $id = "";
my $contig_length = 0;

open (INFILE, $fasta_file) or die ("Unable to open file $fasta_file");
open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile");
while (my $line = <INFILE>){
	if ($line =~ /^>/) {
		if ($contig_length > $length_cutoff){
			print OUTFILE $id;
			foreach my $seq_line(@seq){
				print OUTFILE $seq_line."\n";
			}
		}
		$id = $line;
		@seq = ();
		$contig_length = 0;
	} else {
		chomp $line;
		push @seq, $line;
		$contig_length += length ($line);
	}
	
}
if ($contig_length > $length_cutoff){
	print OUTFILE $id;
	foreach my $seq_line(@seq){
		print OUTFILE $seq_line."\n";
	}
}
close OUTFILE;
close INFILE;
