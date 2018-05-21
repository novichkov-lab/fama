#! /usr/bin/perl
use strict;
use warnings;
#BEGIN
#{
	#use File::Basename;
	#use Cwd 'abs_path';
	#use lib dirname(abs_path($0));
	#use KronaTools;
#}

my $infile = "/mnt/data2/SEED/test_dir/test_comparative/contigs_uniref100_besthit.txt";
my $uniref_file = "/mnt/data2/Databases/Uniref100_20170719/uniref100.tsv";
my $outfile = "/mnt/data2/SEED/test_dir/test_comparative/contigs_uniref100_besthit_annotated.txt";

my $taxonomy_file = "/mnt/data2/Databases/Krona/taxonomy.tab";
my $taxa_exceptions_file = "/mnt/data2/SEED/scripts_v2/taxonomy/taxa_exceptions.txt";

#Filtering parameters
my $alignment_length_cutoff = 15;
my $identity_cutoff = 40;
my $bitscore_fraction_cutoff = 0.95;

my %proteins = ();
my %hits = ();
my %annotations = ();
my %taxonomy_exceptions = ();

#################
# Lookup tables #
#################

my @taxDepths;
my @taxParents;
my @taxRanks;
my @taxNames;
my %missingTaxIDs;

my %taxInfoByID;
my %taxIDByAcc;
my %invalidAccs;
my %missingAccs;


if (@ARGV == 4) {
	$infile = shift;
	$uniref_file = shift;
	$taxonomy_file = shift;
	$outfile = shift;
} else {
	print "Usage: perl annotate_uniref_hits_dev.pl <besthit file> <UniRef fasta file> <Krona taxonomy.tab file> <output file>\n";
	print "This scripts assigns NCBI taxonomy ID to best hits and produces separate files for functional roles of interest, for subsequent construction of taxonomical profiles.\n";
	exit(0);
};

my $tempfile = "/mnt/data2/Databases/Uniref100_20170719/tempfile";

# my @taxNames = loadTaxonomy();
loadTaxonomy();
load_taxonomy_exceptions();

my $count = 0;
open (INFILE, "$infile") or die ("File $infile not found!");
while (my $line = <INFILE>) {
	chomp $line;
	unless ($line) {
		next;
	}
	my ($contig, $protein, $identity, $length) = split(/\t/, $line);
	if ($identity < $identity_cutoff){
		next;
	}
	if ($length < $alignment_length_cutoff){
		next;
	}
	$proteins{$protein} = 1;
	if (exists $hits{$contig}){
		push @{$hits{$contig}}, $line;
		$count++;
	} else {
		my @arr = ($line);
		$hits{$contig} = \@arr;
		$count++;
	}
}
close INFILE;

print scalar(keys %hits) . " contigs, " . $count . " hits found.\n";

open (INFILE, $uniref_file) or die ("Unable to open input file $uniref_file");
my $skip_line = <INFILE>; #skip first line
while (my $line = <INFILE>){
	chomp $line;
	unless (%proteins) {
		last;
	}
	my ($prot_id, $name, $members_number, $tax, $taxid, $repid) = split (/\t/, $line);
	if (exists $proteins{$prot_id}) {
		$annotations{$prot_id}{"name"} = $name;
		$annotations{$prot_id}{"n"} = $members_number;
		$annotations{$prot_id}{"tax"} = $tax;
		$annotations{$prot_id}{"taxid"} = $taxid;
		$annotations{$prot_id}{"repid"} = $repid;
		delete $proteins{$prot_id};
	}
}
close INFILE;


print "Protein annotations: " . scalar(keys %annotations) . " entries found.\n";


open (OUTFILE, ">$outfile") or die ("Unable to open file $outfile ");
print OUTFILE "#Contig\tBest hit\tIdentity\tLength\tMismatch\tSubj_length" 
	. "\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tEvalue" 
	. "\tBitscore\tName\tMembers\tTax\tTaxID\tRepID\tPhylum\tLCA_hits\tLCA_taxid\tLCA_name\n";
foreach my $contig (sort keys %hits){
	print OUTFILE get_besthit($contig), "\n";
}
close OUTFILE;


sub load_taxonomy_exceptions{
	open (INFILE, "$taxa_exceptions_file") or die ("File $taxa_exceptions_file not found!");
	while (my $line = <INFILE>) {
		chomp $line;
		unless ($line) {
			next;
		}
		$taxonomy_exceptions{$line} = 1;
	}
	close INFILE;
}

sub get_besthit{
	my ($contig) = @_;
	my $first_hit = ${$hits{$contig}}[0];
#	my ($contig, $best_protein, $best_identity, $best_length, $best_mismatch, $best_slen, $best_qstart, $best_qend, $best_sstart, $best_send, $best_evalue, $best_bitscore) = split ("\t", $first_hit);
	my @ret_val = split ("\t", $first_hit);
	my $best_bitscore = $ret_val[11];
	my $best_protein = $ret_val[1];
	if (exists $annotations{$best_protein}{"name"}){
		push @ret_val, $annotations{$best_protein}{"name"};
		push @ret_val, $annotations{$best_protein}{"n"};
		push @ret_val, $annotations{$best_protein}{"tax"};
		push @ret_val, $annotations{$best_protein}{"taxid"};
		push @ret_val, $annotations{$best_protein}{"repid"};
		if ($annotations{$best_protein}{"taxid"}) {
			push @ret_val, get_phylum_name($annotations{$best_protein}{"taxid"});
		} else {
			push @ret_val, "";
		}
	}
	my %candidate_taxids = ();
	foreach my $next_hit (@{$hits{$contig}}){
		my (undef, $protein, $identity, $length, undef, undef, undef, undef, undef, undef, undef, $bitscore) = split ("\t", $next_hit);
		if ($bitscore/$best_bitscore >= $bitscore_fraction_cutoff){
			if (exists $annotations{$protein}{"taxid"}) {
				unless (($annotations{$protein}{"taxid"} eq "0")||($annotations{$protein}{"taxid"} eq "1")){
					my $next_hit_taxid = $annotations{$protein}{"taxid"};
					if (exists $taxNames[$next_hit_taxid]) {
						$candidate_taxids{$next_hit_taxid} = 1;
					} else {
						print "WARNING: TaxID $next_hit_taxid not found in Krona database.\n";
					}
				}
			} else {
				print "ERROR: TaxID for protein $protein not found.\n";
			}
		}
	}
	my $lca_taxid = "1";
	push @ret_val, scalar(keys %candidate_taxids); # Number of hits for LCA
	if (scalar(keys %candidate_taxids) > 1) {
		$lca_taxid = get_lca(keys %candidate_taxids);
	} elsif (scalar(keys %candidate_taxids) == 1) {
		$lca_taxid = $annotations{$best_protein}{"taxid"};
	}
	if ($lca_taxid eq "1") {
		push @ret_val, "LCA not found";
		push @ret_val, "LCA not found";
	} else {
		push @ret_val, $lca_taxid;
		push @ret_val, get_lca_name($lca_taxid);
	}
	my $return_value = join("\t", @ret_val);
	return $return_value;
}

sub get_lca{
	return taxLowestCommonAncestor(@_);
}

sub get_lca_name{
	my ($taxid) = @_;
	return $taxNames[$taxid];
}

sub get_phylum_name{
	my ($taxid) = @_;
	my $depth = "unset";
	my $phylum_name = "Unknown";
	until ($taxid == 1){
		unless (defined $taxParents[$taxid]){
			last;
		}
		if (exists $taxonomy_exceptions{$taxParents[$taxid]}){
			$phylum_name = $taxNames[$taxid];
			last;
		} elsif ($taxRanks[$taxid] eq "phylum"){
			$phylum_name = $taxNames[$taxid];
			last;
		} else {
			$taxid = $taxParents[$taxid];
		}
	}
	return $phylum_name;
}

######################################
#### KronaTools code starts here  ####
######################################

sub loadTaxonomy
{
	open INFO, $taxonomy_file or die
		"Taxonomy file $taxonomy_file not found";
	
	while ( my $line = <INFO> )
	{
		chomp $line;
		my ($id, $depth, $parent, $rank, $name) = split /\t/, $line;
		
		$taxParents[$id] = $parent;
		$taxDepths[$id] = $depth;
		$taxRanks[$id] = $rank;
		$taxNames[$id] = $name;
	}
	
	close INFO;
}


sub taxLowestCommonAncestor
{
	my @nodes = @_;
	
	# walk the nodes up to an equal depth
	#
	my $minDepth;
	#
	foreach my $node ( @nodes )
	{
		if ( ! taxIDExists($node) )
		{
			$missingTaxIDs{$node} = 1;
			$node = 1;
		}
		
		if ( ! defined $minDepth || getTaxDepth($node) < $minDepth )
		{
			$minDepth = getTaxDepth($node);
		}
	}
	#
	foreach my $node ( @nodes )
	{
		while ( getTaxDepth($node) > $minDepth )
		{
			$node = getTaxParent($node);
		}
	}
	
	my $done = 0;
	
	while ( ! $done )
	{
		$done = 1;
		
		my $prevNode;
		
		foreach my $node ( @nodes )
		{
			if ( defined $prevNode && $prevNode != $node )
			{
				$done = 0;
				last;
			}
			
			$prevNode = $node;
		}
		
		if ( ! $done )
		{
			for ( my $i = 0; $i < @nodes; $i++ )
			{
				if ( ! defined getTaxParent($nodes[$i]) )
				{
					ktDie("Undefined parent for taxID $nodes[$i]");
					return;
				}
				
				$nodes[$i] = getTaxParent($nodes[$i]);
			}
		}
	}
	
	return $nodes[0];
}

sub taxIDExists
{
	my ($taxID) = @_;
	
	return defined getTaxParent($taxID);
}

sub getTaxDepth
{
	my ($taxID) = @_;
	
	if ( @taxDepths )
	{
		return $taxDepths[$taxID];
	}
	else
	{
		return (getTaxInfo($taxID))[1];
	}
}

sub getTaxParent
{
	my ($taxID) = @_;
	
	if ( @taxParents )
	{
		return $taxParents[$taxID];
	}
	else
	{
		print "Invoked getTaxInfo($taxID)\n";
		return (getTaxInfo($taxID))[2];
	}
}
