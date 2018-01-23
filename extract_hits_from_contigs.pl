#! /usr/bin/perl
use strict;
use warnings;

###################################################################
###  Parse DIAMOND output and report protein sequences of hits  ###
###################################################################

my $diamond_output_file = "/mnt/data2/FEBA/nitrogen_cycle_roles/1728_assembly/AmoB1_assembly/AmoB1_d_results/AmoB1_out.unpadded.fasta_uniref100_besthit_annotated.txt";
my $fasta_file = "/mnt/data2/FEBA/nitrogen_cycle_roles/1728_assembly/AmoB1_assembly/AmoB1_d_results/AmoB1_out.unpadded.fasta";
my $outfile = "/mnt/data2/FEBA/nitrogen_cycle_roles/1728_assembly/AmoB1_assembly/AmoB1_d_results/AmoB1_out.unpadded.proteins.fasta";

my $identity_cutoff = 40;
my $alignment_length_cutoff = 15;

if (@ARGV == 3) {
	$diamond_output_file = $ARGV[0];
	$fasta_file = $ARGV[1];
	$outfile = $ARGV[2];
} else {
	print "Usage: perl extract_hits_from_contigs.pl <diamond output file name> <fasta file name> <output file name>\n";
#	exit(0);
};

if (!(-e $diamond_output_file)) {
	print "File $diamond_output_file not found!\n";
	exit(1);
}

if (!(-e $fasta_file)) {
	print "File $fasta_file not found!\n";
	exit(1);
}

my %contig_data = ();

open (INFILE, $diamond_output_file) or die ("Unable to open file $diamond_output_file");
while (my $line = <INFILE>){
	chomp $line;
	unless ($line) {
		next;
	}
		
	my ($contig, undef, $identity, $alignment_length, undef, undef, $start, $end) = split(/\t/, $line);
	no warnings 'numeric';
	$identity = int($identity);
	unless ($identity > $identity_cutoff) {
		next;
	}
	unless ($alignment_length > $alignment_length_cutoff) {
		next;
	}
	my @arr = ($start, $end);
	$contig_data{$contig} = \@arr;

}
close INFILE;

my $seq = "";
my $current_contig = "";
open (INFILE, $fasta_file) or die ("Unable to open file $fasta_file");
open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile");
while (my $line = <INFILE>){
	chomp $line;
	if ($line =~ /^>/) {
		if (exists $contig_data{$current_contig}){
			$seq =~ s/\s//g;
			my $protein = get_protein(${$contig_data{$current_contig}}[0], ${$contig_data{$current_contig}}[1], $seq);
			if ($protein) {
				#print "Contig $current_contig\n";
				print OUTFILE ">" . $current_contig . "_" . ${$contig_data{$current_contig}}[0] . "_" . ${$contig_data{$current_contig}}[1] . "\n" . $protein . "\n";
			} else {
				print "Hit for contig $current_contig not found\n";
			}
		} else {
			print "Contig $current_contig not found\n";
		}
		$current_contig = $line;
		$current_contig =~ s/^>//;
		$seq = "";
	} else {
		$seq .= $line;
	}
	
}

if (exists $contig_data{$current_contig}){
	$seq =~ s/\s//g;
	my $protein = get_protein(${$contig_data{$current_contig}}[0], ${$contig_data{$current_contig}}[1], $seq);
	if ($protein) {
		#print "Contig $current_contig\n";
		print OUTFILE ">" . $current_contig . "_" . ${$contig_data{$current_contig}}[0] . "_" . ${$contig_data{$current_contig}}[1] . "\n" . $protein . "\n";
	} else {
		print "Protein for contig $current_contig not found\n";
	}
}

close OUTFILE;
close INFILE;

print "Run cdhit -i " . $outfile . " -o " . $outfile . ".clusters -c 0.9 -d 0\n";


sub get_protein {
	my ($start, $end, $seq) = @_;
#	print "Start $start End $end \n";
	$start += 0;
	$end += 0;
	if ($start < $end) {
		my $offset = $start - 1;
		my $length = $end - $start + 1;
		$seq = substr $seq, $offset, $length;
	} else {
		my $offset = $end - 1;
		my $length = $start - $end + 1;
#		print $seq . "\n\n";
		$seq = substr $seq, $offset, $length;
#		print $seq . "\n\n";
		$seq = reverse_complement($seq);
#		print $seq . "\n\n";
	}		

	my $protein='';
	my $codon;
	for(my $i=0;$i<(length($seq)-2);$i+=3){
		$codon=substr($seq,$i,3);
		$protein.=&codon2aa($codon);
	}

	return $protein;
}

sub codon2aa{
	my ($codon) = @_;
	my $undetermined_aa = "X";
	$codon= uc $codon;
#	my (%g) = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');

	# for CD-HIT, stops should be masked with 'X'
	my (%g) = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'X','TAG'=>'X','TGC'=>'C','TGT'=>'C','TGA'=>'X','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
	if (exists $g{$codon}){
		return $g{$codon};
	} else {
#		print STDERR "Bad codon \"$codon\"!!\n";
		return $undetermined_aa;
#		exit(1);
	}
}

sub reverse_complement {
	my ($seq) = @_;
	$seq = lc $seq;
	my @nucleotides = split(//, $seq);
	@nucleotides = reverse @nucleotides;
	$seq = "";
	foreach my $nucleotide (@nucleotides){
		if ($nucleotide eq "a") {
			$seq .= "t";
		} elsif ($nucleotide eq "c") {
			$seq .= "g";
		} elsif ($nucleotide eq "g") {
			$seq .= "c";
		} elsif ($nucleotide eq "t") {
			$seq .= "a";
		} else {
			$seq .= "n";
		}
	}
	return $seq;	
}
