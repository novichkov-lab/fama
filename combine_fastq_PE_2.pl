#! /usr/bin/perl
use strict;
use warnings;

#########################################################################################
###  Combines reads mapped to the same SEED functional role into a single FASTQ file  ###
#########################################################################################

# Paths
my $role = "";
my $outdir = "";
my %samples = ();

if (@ARGV > 3) {
	$role = shift;
	$outdir = shift;
	if (scalar(@ARGV)%2 == 0){
		for (my $i = 0; $i < @ARGV; $i+=2){
			unless (-e $ARGV[$i+1]) {
				print "ERROR: Directory $ARGV[$i+1] not found!\n";
				exit(1);
			}
			if (exists $samples{$ARGV[$i]}) {
				push @{$samples{$ARGV[$i]}}, $ARGV[$i+1];
			} else {
				my @arr = ($ARGV[$i+1]);
				$samples{$ARGV[$i]} = \@arr;
			}
		}
	} else {
		print "ERROR: number of sample IDs maust be equal to number of directories";
		exit(1);
	}
} else {
	print "Usage: perl combine_fastq_PE_2.pl <output directory> <list of sample ids and  directories with paired-end fastq files>\n";
	print "In FASTQ identifier, the last symbol must designate the member of a pair, i.e. either 1 or 2.\n";
	print "This script reads list of best hits, extract reads and create separate fasta for each functional role.\n";
#	exit(0);
};

# Constants
my $file1_extension = ".1.reads.fastq";
my $file2_extension = ".2.reads.fastq";
my $info_file = "sample_info.txt";

# Data structures
my %reads = ();

unless (-e $outdir) {
	print "Directory $outdir not found!\n";
	exit(1);
}
my $outfile1 = $outdir . "/" . $role . $file1_extension;
my $outfile2 = $outdir . "/" . $role . $file2_extension;
$info_file = $outdir . "/" . $info_file;


open (INFOFILE, ">>$info_file") or die ("Unable to create file $info_file \n");

print INFOFILE "Paramters:\nFunctional role $role\nDirectories:\n";
my $sample = 0;

foreach my $sample (sort keys %samples){
	foreach my $directory (@{$samples{$sample}}){
		print INFOFILE $sample . "\t" . $directory;
		my $file1 = $directory . "/" . $role . $file1_extension;
		my $file2 = $directory . "/" . $role . $file2_extension;
		unless (-e $file1) {
			print "File $file1 not found! Skipping file $file2\n";
			next;
		}
		unless (-e $file2) {
			print "File $file2 not found! Skipping file $file1\n";
			next;
		}
		print INFOFILE "\t" . read_fastq_file($file1, $sample);
		print INFOFILE "\t" . read_fastq_file($file2, $sample) . "\n";
	}
}

open (OUTFILE1, ">$outfile1") or die ("Unable to create file $outfile1 \n");
open (OUTFILE2, ">$outfile2") or die ("Unable to create file $outfile2 \n");
my $read_count = 0;
foreach my $sample (sort keys %samples){
	foreach my $read_id (sort keys %{$reads{$sample}}){
		unless (exists $reads{$sample}{$read_id}{"1"}){
			print "ERROR: First end for read $read_id not found.\n";
			exit(1);
		}
		unless (exists $reads{$sample}{$read_id}{"2"}){
			print "ERROR: Second end for read $read_id not found.\n";
			exit(1);
		}
		print OUTFILE1 $reads{$sample}{$read_id}{"1"};
		print OUTFILE2 $reads{$sample}{$read_id}{"2"};
		$read_count++;
	}
}
print INFOFILE "Total number of reads written:\t$read_count\t$read_count\n\n";
close OUTFILE1;
close OUTFILE2;
close INFOFILE;

print "File $outfile1 created.\nFile $outfile2 created.\n";


sub read_fastq_file{
	my ($infile, $sample) = @_;
	open (INFILE, $infile) or die ("Unable to open file $infile");
	my $current_read = "";
	my $flag = "";
	my $read_count = 0;
	my $line_count = 0;
	while (my $line = <INFILE>){
		chomp $line;
		if ($line_count == 0){
			$line_count++;
			$line =~ s/^@//;
			$current_read = "@" . $sample . "_" . $line;
			if (($current_read =~ /\/1$/)||($current_read =~ /\/2$/)){
				$flag = chop $current_read;
			} elsif (($current_read =~ /\s1/)||($current_read =~ /\s2/)) {
				my $tmp_id = $current_read;
				$tmp_id =~ s/.*\s//;
				($flag) = split("",$tmp_id);
			} else {
				die("ERROR: Unable to find end identifier in string \"$line\".\n");
			}
			$current_read =~ s/\s.*//;
			if (exists $reads{$sample}{$current_read}{$flag}){
				$current_read = "";
				$flag = "";
			} else {
				$reads{$sample}{$current_read}{$flag} = "@" . $sample . "_" . $line . "\n";
				$read_count++;
			}
		} else {
			$line_count++;
			if ($line_count == 4){
				$line_count = 0;
			}
			if ($flag ne "") {
				$reads{$sample}{$current_read}{$flag} .= $line . "\n";
			}
		}
	}
	close INFILE;
	return $read_count;
}
