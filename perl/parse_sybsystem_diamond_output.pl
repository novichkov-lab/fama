#! /usr/bin/perl
use strict;
use warnings;
use IO::Zlib qw(:gzip_external 1);
use FileHandle;

#############################################################
###  Parse DIAMOND output of a metagenome vs. subsystems  ###
#############################################################

my $diamond_outfile = "";
my $outfile = "";
my $nucleotide_outfile = "";
my $metagenome_file = "";
my $identity_cutoff = 50;
my $alignment_length_cutoff = 15;
my $hits_overlap = 10;
my %blast_results = ();
my $line = "";

if (@ARGV == 4) {
	$diamond_outfile = $ARGV[0];
	$metagenome_file = $ARGV[1];
	$outfile = $ARGV[2];
	$nucleotide_outfile = $ARGV[3];
} else {
	print "Usage: perl parse_subsystem_diamond_fastq.pl <diamond output file name> <metagenome fastq file name> <output file name> <hits fastq file name>\n";
	exit(0);
};

unless (-e $metagenome_file) {
	print "File $metagenome_file not found!\n";
	exit(1);
}

open (INFILE, $diamond_outfile) or die ("Unable to open file $diamond_outfile");
while ($line = <INFILE>){
	chomp $line;
	if ($line =~ /^\# /) {
	# skip comments
	} elsif ($line eq "") {
	# skip empty lines
	} else {
		
		my @entry = split(/\t/, $line);
		my $q_start = $entry[6] + 0;
		my $q_end = $entry[7] + 0;
		# hits my overlap by $hits_overlap base pairs
		if ($q_start < $q_end) {
			$q_start = $q_start + $hits_overlap;
			$q_end = $q_end - $hits_overlap;
		} else {
			$q_start = $q_start - $hits_overlap;
			$q_end = $q_end + $hits_overlap;		
		}
		my $identity = $entry[2] + 0;
		unless ($identity > $identity_cutoff) {
			next;
		}
		my $alignment_length = $entry[3] + 0;
		unless ($alignment_length >= $alignment_length_cutoff) {
			next;
		}
		my $bitscore = $entry[11] + 0;
#		print "Next line start=$q_start end=$q_end identity=$identity bit-score=$bitscore\n";
		if (exists $blast_results{$entry[0]})  {
			
			my $flag = 1; #if the flag is set to 0, the line will not be added to list of hits
			my @existing_hits = @{$blast_results{$entry[0]}};
			my $i = 0;
			foreach my $existing_hit (@existing_hits) {

#				print "Entry: $ex_entry \n";
				my @existing_entry = split (/\t/, $existing_hit);
				my $exist_q_start = $existing_entry[6] + 0;
				my $exist_q_end = $existing_entry[7] + 0;
				my $exist_identity = $existing_entry[2] + 0;
				my $exist_bitscore = $existing_entry[11] + 0;
				# print "Existing line in hash $existing_entry[1]: start=$exist_q_start end=$exist_q_end identity=$exist_identity bit-score=$exist_bitscore\n";

				if (($q_start < $q_end) && ($exist_q_start > $exist_q_end)) {
					#do nothing if hits on different strands
#					print "hits on different strands \n";
				} elsif (($q_start > $q_end) && ($exist_q_start < $exist_q_end)) {
					#do nothing if hits on different strands
#					print "hits on different strands \n";
				} elsif ($q_start < $q_end) {
					#direct strand
					if ( ($exist_q_end < $q_start) || ($exist_q_start > $q_end) ) { 
						#direct strand, non-overlapping hits
#						print "direct strand, non-overlapping hits \n";
					} elsif ($bitscore <= $exist_bitscore) {
						#direct strand, overlapping hits, existing hit is better than new one
#						print "direct strand, overlapping hits, low score \n";
						$flag = 0;
					} else {
						#direct strand, overlapping hits, new hit is better than existing one
						$blast_results{$entry[0]}[$i] = $line;
#						print "direct strand, better overlapping hit found!\n";
						$flag = 0; #new line replaced exising entry, no need to add it
					}
				} elsif ($q_start > $q_end) {
					#reverse strand
					if ( ($exist_q_end > $q_start) || ($exist_q_start < $q_end) ) { 
						#reverse strand, non-overlapping hits
#						print "reverse strand, non-overlapping hits \n";
					} elsif ($bitscore <= $exist_bitscore) {
						#reverse strand, overlapping hits, existing hit is better than new one
#						print "reverse strand, overlapping hits, low score \n";
						$flag = 0;
					} else {
						#reverse strand, overlapping hits, new hit is better than existing one
						$blast_results{$entry[0]}[$i] = $line;
#						print "reverse strand, better overlapping hit found!\n";
						$flag = 0; #new line replaced exising entry, no need to add it
					}
				} else {
					print "WARNING: something wrong with line \n$line\n";
				}
				$i++;
			}
			if ($flag) {
				my @arr = ($line);
#				print "Putting line to hash: " .$line."\n";
				push @{$blast_results{$entry[0]}}, @arr;
#				print "Line in hash:" . $blast_results{$entry[0]}[1] . "\n";
			}
		} else {
			my @arr = ($line);
#			print "First hit in hash\n";
			$blast_results{$entry[0]} = \@arr;
			#print "Line in hash:" . $blast_results{$entry[0]}[0] . "\n";
		}
	}
}
close INFILE;

open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile ");
foreach my $entries (keys %blast_results) {
	foreach my $entry (@{$blast_results{$entries}}) {
		print OUTFILE $entry."\n";
	}
}
close OUTFILE;

#exit(0);

my %sequences = ();
my $fh;
if ($metagenome_file =~ /\.gz$/) {
	$fh = IO::Zlib->new($metagenome_file, "rb");
} else {
	$fh = FileHandle->new($metagenome_file, "r");
}
my $sequence_id = "";
my $flag = 0;
if (defined $fh) {
	while ($line = $fh->getline()){
		chomp $line;
		if ($line =~ /^@/) {
			$flag = 0;
			$sequence_id = $line;
			$sequence_id =~ s/@//;
			my @arr = split(/ /, $sequence_id);
			$sequence_id = $arr[0];
			if (exists $blast_results{$sequence_id}) {
				my $read_identifier =$sequence_id;
				$read_identifier =~ s/\/1|\/2//;
				$sequences{$sequence_id}{'line1'} = $read_identifier . "/" . get_end_number($line);
				$flag = 1;
			} else {
				$sequence_id = "";
			}
		} elsif ($flag == 1) {
			$sequences{$sequence_id}{'line2'} = $line;
			$flag = 2;
		} elsif ($flag == 2) {
			$sequences{$sequence_id}{'line3'} = $line;
			$flag = 3;
		} elsif ($flag == 3) {
			$sequences{$sequence_id}{'line4'} = $line;
			$flag = 0;
		}
	}
	$fh->close;
}

open (OUTFILE, ">$nucleotide_outfile") or die ("Unable to open file $nucleotide_outfile");

foreach my $sequence_id (keys %blast_results) {
	foreach my $entries (@{$blast_results{$sequence_id}}) {
		my @entry = split (/\t/, $entries);
		if (exists $sequences{$sequence_id}) {
			print OUTFILE "@" . $sequences{$sequence_id}{'line1'} . "|" . $entry[6] . "|" . $entry[7] . "\n";
			my $start = $entry[6] + 0;
			my $end = $entry[7] + 0;		
			my $hit_sequence = $sequences{$sequence_id}{'line2'};
			my $hit_quality = $sequences{$sequence_id}{'line4'};
			if ($start < $end) {
				my $offset = $start - 1;
				my $length = $end - $start;
				$hit_sequence = substr $hit_sequence, $offset, $length;
				$hit_quality = substr $hit_quality, $offset, $length;
			} else {
				my $offset = $end - 1;
				my $length = $start - $end;
				$hit_sequence = substr $hit_sequence, $offset, $length;
				$hit_quality = substr $hit_quality, $offset, $length;
			}
			print OUTFILE $hit_sequence."\n";
			print OUTFILE $sequences{$sequence_id}{'line3'} . "\n";
			print OUTFILE $hit_quality."\n";
		} else {
			print "WARNING: Sequence $sequence_id not found in the metagenome\n";
		}
	}
}

close OUTFILE;
exit(0);

#######################
###   SUBROUTINES   ###
#######################

sub get_end_number{
	my ($read_id) = @_;
	my $ret_val = 0;
	if ($read_id =~ /\/1$|\/2$/) { # old Illumina format
		$ret_val = chop $read_id; 
	} else{
		my (undef, $id) = split(/ /, $read_id);
		if ($id =~ /^1:|2:/){
			$ret_val = substr($id, 0, 1);
		} else {
			return "";
		}
	}
	return $ret_val;
}
