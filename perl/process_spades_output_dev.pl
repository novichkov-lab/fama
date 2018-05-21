#! /usr/bin/perl
use strict;
use warnings;

####################################################################################################
###  Parse MetaSPAdes output and create a list of contigs with counts of reads from each sample  ###
####################################################################################################

# Paths
my $contigs_datafile = "/mnt/data2/FEBA/nitrogen_v6/11016/spades/contigs_data.txt";
my $diamond_output = "/mnt/data2/FEBA/nitrogen_v6/11016/contigs_uniref100_besthit_annotated.txt";
my $out_file = "/mnt/data2/FEBA/nitrogen_v6/11016/11016.contigs.data.txt";

# Parameters
my $coverage_cutoff = 3;
my $read_count_cutoff = 7;

if (@ARGV == 3) {
	$contigs_datafile = $ARGV[0];
	$diamond_output = $ARGV[1];
	$out_file = $ARGV[2];
} else {
#	print "Usage: perl process_spades_output.pl <contigs table> <annotated DIAMOND output file> <output file name>\n";
#	print "This script parses DIAMOND output and creates a list of contigs with hits from UniRef database.\n";
#	exit(0);
};

unless (-e $contigs_datafile) {
	print "File $contigs_datafile not found!\n";
	exit(1);
}

unless (-e $diamond_output) {
	print "File $diamond_output not found!\n";
	exit(1);
}

#Data structures
my %sample_ids = ();
my %contigs=();
my %hits = ();
my %read_stats = ();
my %read_counts = ();
my %coverage = ();
my $header = "";

#read list of contigs
read_contigs_table($contigs_datafile);
$header .= "\t";
#read DIAMOND output file
read_diamond_output($diamond_output);


#write output
open (OUTFILE, ">$out_file") or die;
print OUTFILE $header;
foreach my $contig (sort keys %contigs){
	print OUTFILE $contigs{$contig}{"data"}, "\t";
	if ($contigs{$contig}{"filtered"} eq "true"){
		print OUTFILE "\tHits not expected\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n";
	} elsif (exists $contigs{$contig}{"hit"}) {
		print OUTFILE $contigs{$contig}{"hit"}, "\n";
	} else {
		print OUTFILE "\tHits not found\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n";
	}
}
close OUTFILE;

#######################
###   SUBROUTINES   ###
#######################

sub read_contigs_table{
	my ($infile) = @_;
	open (INFILE, $infile) or die ("Unable to open file $infile");
	$header = <INFILE>;
	chomp $header;
	my @tokens = split(/\t/, $header);
	while(){
		my $token = shift @tokens;
		if ($token eq "Filtered"){
			last;
		}
		unless (defined $token){
			print "ERROR: Filtered not found\n";
			exit(1);
		}
	}
	foreach my $token (@tokens) {
		$sample_ids{$token} = 1;
	}
	while (my $line = <INFILE>){
		chomp $line;
		unless ($line) {
			next;
		}
		my ($contig, undef, undef, undef, undef, $filtered) = split(/\t/, $line);
		$contigs{$contig}{"filtered"} = $filtered;
		$contigs{$contig}{"data"} = $line;
	}
	close INFILE;
}

sub read_diamond_output{
	my ($infile) = @_;
	open (INFILE, $infile) or die ("Unable to open file $infile");
	$header .= <INFILE>;
	while (my $line = <INFILE>){
		chomp $line;
		unless ($line) {
			next;
		}
		if ($line =~ /^#/) {
			next;
		}
		my ($contig) = split(/\t/, $line);
		if (exists $contigs{$contig}{"data"}){
			$contigs{$contig}{"hit"} = $line;
		} else {
			print "ERROR: Contig $contig from DIAMOND output is missing from FASTA file\n";
			exit(1);
		}
	}
	close INFILE;
}

