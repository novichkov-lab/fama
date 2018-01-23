#! /usr/bin/perl
use strict;
use warnings;

###############################################################
###  Parse BLAST output of search 2, generate output report  ###
###############################################################

my $search1_hits_file = "";
my $search2_blast_output = "";
my $outfile = "";
my $rolesfile = "/mnt/data2/SEED/seed_roles/nr_roles_list.txt";
my $functional_proteins_collection = "nitrate_roles_db/nitrate_genes_roles_filtered_nr95.txt";
my $hits_fasta_file = "";
my $fasta_dir = "";
my $focus_sh="runfocus.sh";
my $identity_cutoff = 50; #
my $length_cutoff = 15; #
my $reads_total = 0; #Total number of sequence reads in fasta file
my $narrow_cutoff = 10;
my $uncertain_role = "*_*";
my $fasta_file_name = "all_rbh.fna";
my $carma_file_suffix = ".tax";
my $out_roles_file = "roles_stat.txt";
my $out_hist_file = "hits_identity_hist.txt";

my %genes_roles = ();
my %search1_results_list = ();
my %search2_results_list = ();
my %selected_genes = ();
my %skipped_genes = ();
my %different_hits_same_role_genes = ();
my %relevant_function_genes = ();
my %taxonomy_HoA = ();
my @missing_genes = ();
my %roles = ();
my %roles_dict = ();
my %roles_scores = (); 
my %roles_counts = (); 
my %sequences = ();
my %identity_distribution = ();
my $identity_distribution_step = 10;
$roles_dict{$uncertain_role} = "Functional role uncertain: this role assigned to proteins that have multiple role mappings";


if (@ARGV == 7) {
	$search1_hits_file = $ARGV[0];
	$search2_blast_output = $ARGV[1];
	$functional_proteins_collection = $ARGV[2]; 
	$hits_fasta_file = $ARGV[3];
	$outfile = $ARGV[4];
	$fasta_dir = $ARGV[5];
	$reads_total = $ARGV[6];
} else {
	print "Usage: perl parse_proteomes_blastx.pl <file with list of search1 hits> <file with blastx tabular outputof search2> <text file with list of functional roles> <search2 fasta query file name> <out file name> <directory for fasta files> <number of sequence reads>\n";
	print "This script compares results of two sequence similarity searches. Search 1 should be done against a database of functional proteins \n
(proteins mapped to SEED functional roles of interest). Search 2 should be done against a database of complete bacterial proteomes. \n
Script outputs a text file with list of best hits and separate fasta files for individual functional roles.\n";
	exit(0);
};

$out_roles_file = $fasta_dir . "/" . $out_roles_file;
$out_hist_file = $fasta_dir . "/" . $out_hist_file;

#initialize distribution hash
for (my $i = 0; $i < 10; $i++){
	my $bin = ($i+1)*$identity_distribution_step;
	$identity_distribution{$bin} = 0;
}


#check if required files and directories exist
if (!(-e $fasta_dir)) {
	print "Directory $fasta_dir not found!\n";
	exit(1);
}

if (!(-e $rolesfile)) {
	print "File $rolesfile not found!\n";
	exit(1);
}

if (!(-e $functional_proteins_collection)) {
	print "File $functional_proteins_collection not found!\n";
	exit(1);
}

#read list of genes in the database of functional genes (search 1 db))
open (INFILE, $functional_proteins_collection) or die ("Unable to open file $functional_proteins_collection");
my $line ="";
while ($line = <INFILE>){
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ s/^>//;
		my ($roles, $gene) = split(/\_/, $line);
		$genes_roles{$gene} = $roles;
	}
}
close INFILE;


#read results of blastx search 1 (vs. database of functional genes)
open (INFILE, $search1_hits_file) or die ("Unable to open file $search1_hits_file");
while ($line = <INFILE>){
	chomp $line;
	my @entry = split(/\t/, $line);
	if (exists $search1_results_list{$entry[0]}) {
		my @arr = ($line);
		push @{$search1_results_list{$entry[0]}}, @arr;
	} else {
		my @arr = ($line);
		$search1_results_list{$entry[0]} = \@arr;
	}
}
close INFILE;

#read results of blastx search 2
open (INFILE, $search2_blast_output) or die ("Unable to open file $search2_blast_output");
while ($line = <INFILE>){
	chomp $line;
	if ($line =~ /^\# /) {
	# skip comments
	} elsif ($line eq "") {
	# skip empty lines
	} else {
		my @entry = split(/\t/, $line);
		unless ($entry[2] >= $identity_cutoff) {
			next;
		}
		unless ($entry[3] >= $length_cutoff) {
			next;
		}
		if (exists $search2_results_list{$entry[0]}) {
			my @arr = ($line);
			push @{$search2_results_list{$entry[0]}}, @arr;
		} else {
			my @arr = ($line);
			$search2_results_list{$entry[0]} = \@arr;
		}
	}
}
close INFILE;


open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile ");
print OUTFILE "True best hits:\n";

for my $seq_id (sort keys %search1_results_list) {
	my @search1_hits_per_protein = @{$search1_results_list{$seq_id}};
	foreach my $search1_hit (@search1_hits_per_protein) {
		my @search1_hit_fields = split(/\t/, $search1_hit);
		#my ($subsystem_id, $role_id, $search1_hit_subject) = split(/_/, $search1_hit_fields[1]);
		my ($role_id, $search1_hit_subject) = split(/_/, $search1_hit_fields[1]);
		my $bitscore_of_besthit1_in_search2 = 0;
		my $search2_query_id = $search1_hit_fields[0] . "|" . $search1_hit_fields[6] . "|" . $search1_hit_fields[7];
		if (exists $search2_results_list{$search2_query_id}) {
			my $search2_best_hit = "";
			my $search2_best_hit_subject = "";
			my $search2_best_hit_identity = "";
			my $search2_best_hit_bitscore = 0;
			my $search2_best_hit_evalue = "";
			foreach my $search2_hit_line (@{$search2_results_list{$search2_query_id}}) {
				my @search2_hit_fields = split(/\t/, $search2_hit_line);
				my $search2_hit_bitscore = $search2_hit_fields[11] + 0;
				if (($search2_hit_fields[1] eq $search1_hit_subject) && ($bitscore_of_besthit1_in_search2 == 0)) {
					$bitscore_of_besthit1_in_search2 = $search2_hit_bitscore;
				}
				if ($search2_hit_bitscore > $search2_best_hit_bitscore) {
					$search2_best_hit = $search2_hit_line;
					$search2_best_hit_subject = $search2_hit_fields[1];
					$search2_best_hit_identity = $search2_hit_fields[2];
					$search2_best_hit_evalue = $search2_hit_fields[10];
					$search2_best_hit_bitscore = $search2_hit_bitscore;
				}
			}
			
			if ($bitscore_of_besthit1_in_search2 == 0) {
				if (exists $genes_roles{$search2_best_hit_subject}) {
					my $rpkm_score = calculate_rpkm($search2_best_hit);
					$relevant_function_genes{$search2_query_id} = $search2_best_hit . "\t" . $rpkm_score . "\t" . $search1_hit;

					#collect data for taxonomic profile
					my $tax_id = get_taxid_id($search2_best_hit_subject);
					my $carma_line = "$search2_query_id\t$search2_best_hit_identity\t$rpkm_score\t$tax_id\tundef\t$search2_best_hit_evalue\n";
					if (exists $taxonomy_HoA{$genes_roles{$search1_hit_subject}}) {
						push @{$taxonomy_HoA{$genes_roles{$search2_best_hit_subject}}}, $carma_line;
					} else {
						my @arr = ($carma_line);
						$taxonomy_HoA{$genes_roles{$search2_best_hit_subject}} = \@arr;
					}

					
				} else {
					#best hit from search1 was not found in search2 results
#					print "Search 1 best hit $search1_hit_subject for query $search2_query_id not found in search 2: removed from list\n";
					$skipped_genes{$search2_query_id} = $search2_best_hit . "\tSEARCH 1 BEST HIT NOT FOUND IN SEARCH 2\t" . $search1_hit;
				}
			} elsif ($search2_best_hit_bitscore <= $bitscore_of_besthit1_in_search2) {
				#Best hit in search 1 is the best hit in search 2. This is really THE BEST hit
				my $rpkm_score = calculate_rpkm($search2_best_hit);
				build_identity_distribution($search1_hit_fields[2]);
				my $output_entry = $search2_best_hit . "\t" . $rpkm_score . "\t" . $role_id;
				print OUTFILE $output_entry . "\n";
				$selected_genes{$search2_query_id} = $output_entry;
				# put the gene to list of roles
				if (exists $roles{$role_id}) {
					$roles{$role_id} = $roles{$role_id} . "\t" . $search2_query_id;
				} else {
					$roles{$role_id} = $search2_query_id;
				}

				#collect data for taxonomic profile
				my $tax_id = get_taxid_id($search1_hit_subject);
				my $carma_line = "$search2_query_id\t$search2_best_hit_identity\t$rpkm_score\t$tax_id\tundef\t$search2_best_hit_evalue\n";
				if (exists $taxonomy_HoA{$role_id}) {
					push @{$taxonomy_HoA{$role_id}}, $carma_line;
				} else {
					my @arr = ($carma_line);
					$taxonomy_HoA{$role_id} = \@arr;
				}
			} elsif ($search2_best_hit_bitscore < $bitscore_of_besthit1_in_search2 + $narrow_cutoff) {
				if (exists $genes_roles{$search2_best_hit_subject}) {
					#Best hit in search 1 and best hit in search 2 are different, but both hits have functional roles of interest
					my $rpkm_score = calculate_rpkm($search2_best_hit);
					if ($genes_roles{$search2_best_hit_subject} eq $genes_roles{$search1_hit_subject}) {
						$different_hits_same_role_genes{$search2_query_id} = $search2_best_hit . "\t" . $rpkm_score . "\t" . $search1_hit; 
					} else {
						$relevant_function_genes{$search2_query_id} = $search2_best_hit . "\t" . $rpkm_score . "\t" . $search1_hit;
					}
					
					#collect data for taxonomic profile
					my $carma_line = "$search2_query_id\t$search2_best_hit_identity\t$rpkm_score\t0\tunknown\t$search2_best_hit_evalue\n";
					if (exists $taxonomy_HoA{$genes_roles{$search2_best_hit_subject}}) {
						push @{$taxonomy_HoA{$genes_roles{$search2_best_hit_subject}}}, $carma_line;
					} else {
						my @arr = ($carma_line);
						$taxonomy_HoA{$genes_roles{$search2_best_hit_subject}} = \@arr;
					}

				} else {

					#Best hit in search 1 was found in search 2, but better hit found in search 2 with close bitscore
					$skipped_genes{$search2_query_id} = $search2_best_hit . "\tBEST HITS IN SEARCHES 1 AND 2 HAVE CLOSE SCORES\t" . $search1_hit;
				}
			
			} elsif ($search2_best_hit_bitscore > $bitscore_of_besthit1_in_search2) {
				if (exists $genes_roles{$search2_best_hit_subject}) {
					my $rpkm_score = calculate_rpkm($search2_best_hit);
					$relevant_function_genes{$search2_query_id} = $search2_best_hit . "\t" . $rpkm_score . "\t" . $search1_hit; 
					
					#collect data for taxonomic profile
					my $carma_line = "$search2_query_id\t$search2_best_hit_identity\t$rpkm_score\t0\tunknown\t$search2_best_hit_evalue\n";
					if (exists $taxonomy_HoA{$genes_roles{$search2_best_hit_subject}}) {
						push @{$taxonomy_HoA{$genes_roles{$search2_best_hit_subject}}}, $carma_line;
					} else {
						my @arr = ($carma_line);
						$taxonomy_HoA{$genes_roles{$search2_best_hit_subject}} = \@arr;
					}

				} else {
					#Best hit in search 1 was found in search 2, but better hit found in search 2
					$skipped_genes{$search2_query_id} = $search2_best_hit . "\tBETTER HIT FOUND IN SEARCH 2\t" . $search1_hit;
				}
			} else {
				print "ERROR: Hit $search1_hit\n";
			}
		} else {
			#Query missing from search 2 results
			#print "$search2_query_id not found in search 2 blastx results\n";
			push @missing_genes, $search2_query_id . "\t" . $search1_hit;
		}
	}
}


print OUTFILE "\nGenes with different best hits in two searches having identical functional roles:\n";
for my $gene (sort keys %different_hits_same_role_genes) {
	print OUTFILE $different_hits_same_role_genes{$gene} . "\n";
}

print OUTFILE "\nGenes with different best hits in two searches having different functional roles:\n";
for my $gene (sort keys %relevant_function_genes) {
	print OUTFILE $relevant_function_genes{$gene} . "\n";
}


print OUTFILE "\nGenes excluded from RBH list:\n";
for my $gene (sort keys %skipped_genes) {
	print OUTFILE $skipped_genes{$gene} . "\n";
}

print OUTFILE "\nMissing genes:\n";
for my $gene (sort @missing_genes) {
	print OUTFILE $gene . "\n";
}

print OUTFILE "\nList of genes for each role:\n";
#for my $role (sort {$a <=> $b} keys %roles) {
for my $role (sort keys %roles) {
	print OUTFILE $role . "\t" . $roles{$role} . "\n";
}

open (INFILE, $rolesfile) or die ("Unable to open file $rolesfile");
while ($line = <INFILE>){
	chomp $line;
	my @entry = split(/\t/, $line);
	$roles_dict{$entry[0]} = $entry[1];
}
close INFILE;
close OUTFILE;

open (OUTFILE, ">$out_roles_file") or die ("Unable to open file $out_roles_file ");
print OUTFILE "\nFunctional roles statistics:\n";
print OUTFILE "\nRole\tRead count\tRPKM\tName\n";
#for my $role (sort {$a<=>$b} keys %roles_scores) {
for my $role (sort keys %roles_scores) {
	print OUTFILE $role . "\t" . $roles_counts{$role} . "\t" . $roles_scores{$role};
	if (exists $roles_dict{$role}) {
	print OUTFILE "\t" . $roles_dict{$role};	
	}
	print OUTFILE "\n";	
}
close OUTFILE;

open (OUTFILE, ">$out_hist_file") or die ("Unable to open file $out_hist_file");

print OUTFILE "\nDistribution of hits identity:\n";
my $hits_count = 0;

for my $bin (keys %identity_distribution){
	$hits_count = $hits_count + $identity_distribution{$bin};
}

for my $bin (sort {$a <=> $b} keys %identity_distribution){
	my $bin_fraction = int($identity_distribution{$bin}*10/$hits_count);
	print OUTFILE "\|";
	for (my $i = 0; $i < $bin_fraction; $i++) {
		print OUTFILE "*";
	}
	for (my $i = $bin_fraction; $i < 10; $i++) {
		print OUTFILE " ";
	}
	print OUTFILE "\| " . $bin . "\%  \(" . $identity_distribution{$bin} . " hits\) \n";
}



#for my $role (sort keys %roles) {
#	my @ids = split(/_/, $role);
#	my @gene_number = split(/\t/, $roles{$role});
#	print OUTFILE $ids[0] . "\t" . $ids[1] . "\t" . scalar(@gene_number) . "\t" . $roles_dict{$role} . "\n";
#}



close OUTFILE;

open (INFILE, $hits_fasta_file) or die ("Unable to open file $hits_fasta_file");
#open (OUTFILE, ">$rbh_fasta_file") or die ("Unable to open file $rbh_fasta_file");
my $sequence_id = "";
my $role_id = "";
while ($line = <INFILE>){
	chomp $line;
	if ($line =~ /^>/) {
		$sequence_id = $line;
		$sequence_id =~ s/>//;
		if (exists $selected_genes{$sequence_id}) {
#			print OUTFILE $line."\n";
			my @entry = split(/\t/, $selected_genes{$sequence_id});
			$role_id = $entry[13];
			if (exists $sequences{$role_id}) {
				my @arr = ($line);
				push @{$sequences{$role_id}}, @arr;
			} else {
				my @arr = ($line);
				$sequences{$role_id} = \@arr;
			}
		} else {
			$sequence_id = "";
			$role_id = "";
		}
	} elsif ($sequence_id) {
		if ($role_id) {
				my @arr = ($line);
				push @{$sequences{$role_id}}, @arr;	
		}
#		print OUTFILE $line."\n";
	}
}
close INFILE;

$fasta_file_name = $fasta_dir . "/" . $fasta_file_name;
open (OUTFILE, ">$fasta_file_name") or die ("Unable to open file $fasta_file_name");
$focus_sh = $fasta_dir . "/" . $focus_sh;
open (SHELLFILE, ">$focus_sh") or die ("Unable to open file $focus_sh");
print SHELLFILE "python focus.py -q $fasta_file_name -m 0.1\n";
my $role_fasta_file = "";
foreach $role_id (keys %sequences) {
	$role_fasta_file = $fasta_dir . "/" . $role_id . ".fna";
	$role_fasta_file =~ s/\|/_/g;
	print SHELLFILE "python focus.py -q \"$role_fasta_file\" -m 0.1\n";
	open (ROLEFILE, ">$role_fasta_file") or die ("Unable to open file $role_fasta_file");
	foreach my $line (@{$sequences{$role_id}}) {
		print OUTFILE $line."\n";
		print ROLEFILE $line."\n";
	}
	close ROLEFILE;
}

close OUTFILE;
close SHELLFILE;

create_carma_inputfiles($fasta_dir);
exit (0);

#####################
###  SUBROUTINES  ###
#####################
sub get_taxid_id {
	my ($tax_id) = @_;
	$tax_id =~ s/^fig\|//;
	$tax_id =~ s/.peg.*//;
	$tax_id =~ s/\..*//;
	return $tax_id;
}


sub calculate_rpkm {
	my ($hit) = @_;
	my @entry = split (/\t/, $hit);
	my $score = 1000000000 / (($entry[5] - $entry[3]) * 3 * $reads_total);
	if ($score < 0) {
		$score = 0 - $score
	}
	store_score ($entry[1], $score);
	store_count ($entry[1]);
	return $score;
}

sub store_score {
	my ($gene, $score) = @_;
	if (exists $genes_roles{$gene}){
		my @roles_list = split(/\|/, $genes_roles{$gene});
		foreach my $role (@roles_list){
			if (!($role)) {
				print "Undefined role " . $role ."\n";
			}
			if (exists $roles_scores{$role}){
				$roles_scores{$role} = $roles_scores{$role} + ($score/scalar(@roles_list));
			} else {
				$roles_scores{$role} = $score/scalar(@roles_list);
			}
		}
	}
}

sub store_count {
	my ($gene) = @_;
	if (exists $genes_roles{$gene}){
		my @roles_list = split(/\|/, $genes_roles{$gene});
		foreach my $role (@roles_list){
			if (!($role)) {
				print "Undefined role " . $role ."\n";
			}
			if (exists $roles_counts{$role}){
				$roles_counts{$role} = $roles_counts{$role} + (1/scalar(@roles_list));
			} else {
				$roles_counts{$role} = 1/scalar(@roles_list);
			}
		}
	}
}

sub build_identity_distribution {
	my ($identity) = @_;
	
	foreach my $bin (keys %identity_distribution){
		if ((($identity + 0) >= $bin) &&(($identity + 0) < ($bin + $identity_distribution_step))){
			$identity_distribution{$bin} = $identity_distribution{$bin} + 1;
			last;
		}
	}
}

sub create_carma_inputfiles {
	my ($carma_dir) = @_;
	my @carma_files = ();
	foreach my $role (keys %taxonomy_HoA) {
		my $carma_file_name = $carma_dir . "/" . $role . $carma_file_suffix;
		$carma_file_name =~ s/\|/_/g;
		open (OUTFILE, ">$carma_file_name") or die ("Unable to open file $carma_file_name");
		foreach my $line (@{$taxonomy_HoA{$role}}){
			print OUTFILE $line;
		}
		close OUTFILE;
		my $taxprofile_file = "$carma_dir/$role.taxprofile.tsv";
		$taxprofile_file =~ s/\|/_/g;		
		push @carma_files, "perl /mnt/data2/SEED/scripts_v2/taxonomy/getTaxonomicProfileRPKM_idstats.pl -o $taxprofile_file -e 1 -l 10000 $carma_file_name\n";
	}
	my $carma_script = $carma_dir . "/run_carma.sh";
	open (OUTFILE, ">$carma_script") or die ("Unable to open file run_carma.sh");
	print OUTFILE "#! /bin/sh\n";
	foreach my $file (@carma_files) {
		print OUTFILE $file;	
	}
	close OUTFILE;
	chmod(0775, $carma_script);
}
