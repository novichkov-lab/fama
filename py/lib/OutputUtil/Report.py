#!/usr/bin/python
import os
from collections import defaultdict,Counter,OrderedDict
from lib.DiamondParser.DiamondParser import cleanup_protein_id
from lib.ReferenceLibrary.TaxonomyData import TaxonomyData


def generate_report(parser):
    outfile = os.path.join(parser.project.get_project_dir(parser.sample), parser.project.get_output_subdir(parser.sample),parser.sample + '_' + parser.end + '_'+ parser.project.get_report_name())
    with open(outfile, 'w') as of:
        # Write general info
        of.write('\nRun info\n\n')
        of.write('Sample ID:\t' + parser.project.get_sample_id(parser.sample) + '\n')
        of.write('Paired end:\t' + parser.end + '\n')
        of.write('FASTQ file:\t' + parser.project.get_fastq_path(parser.sample, parser.end) + '\n')
        of.write('Total number of reads:\t' + str(parser.project.get_fastq1_readcount(parser.sample)) + '\n')
        of.write('*****************************************\n\n')

        # Write read statistics
        of.write('\nRead statistics\n\n')
        read_stats = Counter()
        for read in sorted(parser.reads.keys()):
            read_stats[parser.reads[read].get_status()] += 1
        for status in OrderedDict(read_stats.most_common()):
            if status == 'unaccounted':
                of.write('Reads missing from background DB search result\t' + str(read_stats[status]) + '\n')
            elif status == 'nofunction':
                of.write('Reads not mapped to any function\t' + str(read_stats[status]) + '\n')
            elif status == 'function':
                of.write('Reads mapped to a function of interest\t' + str(read_stats[status]) + '\n')
            elif status == 'function,besthit':
                of.write('Reads mapped to a function of interest and a taxon\t' + str(read_stats[status]) + '\n')
            else:
                of.write(status + '\t' + str(read_stats[status]) + '\n')

        of.write('*****************************************\n\n')

        # Write function scores
        of.write('\nFunction statistics\n')
        func_stats = defaultdict(float)
        func_counts = Counter()
        func_identity = defaultdict(float)
        func_hit_counts = Counter()
        for read in parser.reads.keys():
            if parser.reads[read].get_status() == 'function,besthit' or parser.reads[read].get_status() == 'function':
                functions = parser.reads[read].get_functions()
                for function in functions:
                    func_stats[function] += functions[function]
                    func_counts[function] += 1/len(functions)
                for hit in parser.reads[read].get_hit_list().get_hits():
                    for function in hit.get_functions():
                        func_identity[function] += hit.get_identity()
                        func_hit_counts[function] += 1
        for function in func_identity:
            func_identity[function] = func_identity[function]/func_hit_counts[function]
        of.write('\nFunction\tDefinition\tRPKM score\tRead count\tAvg. identity\n')
        for function in sorted(func_stats.keys()):
            of.write(function + '\t' 
                    + parser.ref_data.lookup_function_name(function) + '\t' 
                    + str(func_stats[function]) + '\t'
                    + str(func_counts[function]) + '\t'
                    + str(func_identity[function]) + '\n')
        of.write('*****************************************\n\n')

        # Write group scores
        of.write('\nFunction statistics by category\n')
        func_stats = defaultdict(float)
        func_counts = Counter()
        func_identity = defaultdict(float)
        func_hit_counts = Counter()
        for read in parser.reads.keys():
            if parser.reads[read].get_status() == 'function,besthit' or parser.reads[read].get_status() == 'function':
                functions = parser.reads[read].get_functions()
                for function in functions:
                    func_stats[parser.ref_data.lookup_function_group(function)] += functions[function]
                    func_counts[parser.ref_data.lookup_function_group(function)] += 1/len(functions)
                for hit in parser.reads[read].get_hit_list().get_hits():
                    for function in hit.get_functions():
                        func_identity[parser.ref_data.lookup_function_group(function)] += hit.get_identity()
                        func_hit_counts[parser.ref_data.lookup_function_group(function)] += 1
        for function in func_identity:
            func_identity[function] = func_identity[function]/func_hit_counts[function]
        of.write('\nCategory\tRPKM score\tRead count\tAvg. identity\n')
        for function in sorted(func_stats.keys()):
            of.write(function + '\t' 
                    + str(func_stats[function]) + '\t'
                    + str(func_counts[function]) + '\t'
                    + str(func_identity[function]) + '\n')

        of.write('*****************************************\n\n')
        # Write taxonomy stats
        of.write('\nTaxonomy statistics for best hits\n\n')
        tax_stats = Counter()
        identity_stats = defaultdict(float)
        rpkm_stats = defaultdict(float)
        for read in parser.reads.keys():
#            print (read, parser.reads[read].get_status())
#            if parser.reads[read].get_status() == 'function,besthit' or parser.reads[read].get_status() == 'function':
            if parser.reads[read].get_status() == 'function,besthit':
#                print ('\t',parser.reads[read].get_status())
                hits = parser.reads[read].get_hit_list().get_hits()
                for hit in hits:
#                    print (hit.get_subject_id())
#                    print (parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id())))
#                    print (hit.get_identity())
                    protein_taxid = parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                    tax_stats[protein_taxid] += 1
                    identity_stats[protein_taxid] += hit.get_identity()
                if len(hits) == 1:
                    read_functions = parser.reads[read].get_functions()
                    for function in read_functions:
                        rpkm_stats[parser.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].get_subject_id()))] += read_functions[function]
                else:
                    read_functions = parser.reads[read].get_functions()
                    protein_taxids = {}
                    for hit in hits:
                        hit_taxid = parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                        hit_functions = hit.get_functions()
                        for hit_function in hit_functions:
                            protein_taxids[hit_taxid] = hit_function
                    for taxid in protein_taxids:
                        if protein_taxids[taxid] in read_functions:
                            rpkm_stats[taxid] += read_functions[protein_taxids[taxid]]
                        
        tax_data = TaxonomyData(parser.config)
        tax_data.load_taxdata(parser.config)
        counts_per_rank, identity_per_rank, rpkm_per_rank = tax_data.get_taxonomy_profile(counts=tax_stats, identity=identity_stats, scores = rpkm_stats)

        ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        for rank in ranks:
            of.write('Taxonomy report for rank ' + rank + '\n\n')
            of.write('Taxon\tRead count\tRPKM score\tAverage identity\n')
            
            for tax in OrderedDict(Counter(counts_per_rank[rank]).most_common()):
                #print (tax + '\t' + str(counts_per_rank[rank][tax]) + '\t' + str(identity_per_rank[rank][tax]))
                of.write(rank + '\t' + tax + '\t' 
                        + str(counts_per_rank[rank][tax]) + '\t' 
                        + str(rpkm_per_rank[rank][tax]) + '\t' 
                        + str(identity_per_rank[rank][tax]) + '\n')
            of.write('*****************************************\n\n')

                    
        of.write('\nList of reads\n')
        # Write list of reads
        for read in sorted(parser.reads.keys()):
            of.write(read + '\t' + parser.reads[read].get_status() + '\n')
            of.write('\t' + str(parser.reads[read].get_functions()) + '\n')
            for hit in parser.reads[read].get_hit_list().get_hits():
                of.write('\t' + str(hit) + '\n')
                
        of.write('\n\n*** End of report ***\n')
        of.closed

