import os
from collections import defaultdict,Counter,OrderedDict

from Fama.utils import autovivify,cleanup_protein_id,sanitize_file_name
from Fama.DiamondParser.hit_utils import get_efpk_score
from Fama.TaxonomyProfile import TaxonomyProfile
from Fama.OutputUtil.XlsxUtil import generate_function_sample_xlsx, generate_function_taxonomy_sample_xlsx, generate_sample_taxonomy_function_xlsx
from Fama.OutputUtil.KronaXMLWriter import generate_functional_taxonomy_chart,generate_taxonomy_series_chart

#ENDS = ['pe1','pe2']
def generate_fastq_report(parser):
    outfile = os.path.join(parser.options.get_project_dir(parser.sample.sample_id), parser.options.get_output_subdir(parser.sample.sample_id),parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_report_name())
    with open(outfile, 'w') as of:
        # Write general info
        of.write('\nRun info\n\n')
        of.write('Sample ID:\t' + parser.options.get_sample_name(parser.sample.sample_id) + '\n')
        of.write('Paired end:\t' + parser.end + '\n')
        of.write('FASTQ file:\t' + parser.options.get_fastq_path(parser.sample.sample_id, parser.end) + '\n')
        of.write('Total number of reads:\t' + str(parser.options.get_fastq1_readcount(parser.sample.sample_id)) + '\n')
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
#            elif status == 'function,besthit':
#                of.write('Reads mapped to a function of interest and a taxon\t' + str(read_stats[status]) + '\n')
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
        of.write('\nFunction\tDefinition\tRPK score\tRead count\tAvg. identity\n')
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
        of.write('\nCategory\tRPK score\tRead count\tAvg. identity\n')
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
            if parser.reads[read].get_status() == 'function':
#            if parser.reads[read].get_status() == 'function':
#                print ('\t',parser.reads[read].get_status())
                taxonomy = parser.reads[read].taxonomy
                if taxonomy is None:
                    print ('No taxonomy ID assigned to ', read)
                    continue
                tax_stats[taxonomy] += 1
                read_functions = parser.reads[read].get_functions()
                for f in read_functions:
                    rpkm_stats[taxonomy] += read_functions[f]
                identity_stats[taxonomy] += sum(list(hit.get_identity() for hit in parser.reads[read].get_hit_list().get_hits()))
                    
                #~ hits = parser.reads[read].get_hit_list().get_hits()
                #~ for hit in hits:
#~ #                    print (hit.get_subject_id())
#~ #                    print (parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id())))
#~ #                    print (hit.get_identity())
                    #~ protein_taxid = parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                    #~ tax_stats[protein_taxid] += 1
                    #~ identity_stats[protein_taxid] += hit.get_identity()
                #~ if len(hits) == 1:
                    #~ read_functions = parser.reads[read].get_functions()
                    #~ for function in read_functions:
                        #~ rpkm_stats[parser.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].get_subject_id()))] += read_functions[function]
                #~ else:
                    #~ read_functions = parser.reads[read].get_functions()
                    #~ protein_taxids = {}
                    #~ for hit in hits:
                        #~ hit_taxid = parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.get_subject_id()))
                        #~ hit_functions = hit.get_functions()
                        #~ for hit_function in hit_functions:
                            #~ protein_taxids[hit_taxid] = hit_function
                    #~ for taxid in protein_taxids:
                        #~ if protein_taxids[taxid] in read_functions:
                            #~ rpkm_stats[taxid] += read_functions[protein_taxids[taxid]]
                        
        tax_data = parser.taxonomy_data
        counts_per_rank, identity_per_rank, rpkm_per_rank = tax_data.get_taxonomy_profile(counts=tax_stats, identity=identity_stats, scores = rpkm_stats)

        ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus']
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
#            if parser.reads[read].get_status() != 'nofunction':
#                of.write(read + '\t' + parser.reads[read].get_status() + '\n')
#                of.write('\t' + str(parser.reads[read].get_functions()) + '\n')
            of.write(read + ': ' + parser.reads[read].get_status() + ': ' + ','.join(sorted(parser.reads[read].get_functions().keys())))
            if not parser.reads[read].taxonomy is None:
                of.write(' Taxonomy :' + parser.reads[read].taxonomy + '\n')
            else:
                of.write(' No taxonomy\n')
            for hit in parser.reads[read].get_hit_list().get_hits():
                of.write('\t' + str(hit) + '\n')
                
        of.write('\n\n*** End of report ***\n')
        of.closed

def generate_fasta_report(parser):
    outfile = os.path.join(parser.options.get_project_dir(parser.sample.sample_id), parser.options.get_output_subdir(parser.sample.sample_id),parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_report_name())
    with open(outfile, 'w') as of:
        # Write general info
        of.write('\nRun info\n\n')
        of.write('Sample ID:\t' + parser.options.get_sample_name(parser.sample.sample_id) + '\n')
        of.write('Sequence file:\t' + parser.options.get_fastq_path(parser.sample.sample_id, parser.end) + '\n')
        of.write('Total number of reads:\t' + str(parser.options.get_fastq1_readcount(parser.sample.sample_id)) + '\n')
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
            if parser.reads[read].get_status() == 'function':
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
        of.write('\nFunction\tDefinition\tRPK score\tRead count\tAvg. identity\n')
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
            if parser.reads[read].get_status() == 'function':
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
        of.write('\nCategory\tRPK score\tRead count\tAvg. identity\n')
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
            if parser.reads[read].get_status() == 'function':
#                print ('\t',parser.reads[read].get_status())
                taxonomy = parser.reads[read].taxonomy
                if taxonomy is None:
                    print ('No taxonomy ID assigned to ', read)
                    continue
                tax_stats[taxonomy] += 1
                read_functions = parser.reads[read].get_functions()
                for f in read_functions:
                    rpkm_stats[taxonomy] += read_functions[f]
                identity_stats[taxonomy] += sum(list(hit.get_identity() for hit in parser.reads[read].get_hit_list().get_hits()))
                    
        tax_data = parser.taxonomy_data
        counts_per_rank, identity_per_rank, rpkm_per_rank = tax_data.get_taxonomy_profile(counts=tax_stats, identity=identity_stats, scores = rpkm_stats)

        ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus']
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

                    
        of.write('\nList of sequences and hits\n')
        # Write list of reads
        for read in sorted(parser.reads.keys()):
#            if parser.reads[read].get_status() != 'nofunction':
#                of.write(read + '\t' + parser.reads[read].get_status() + '\n')
#                of.write('\t' + str(parser.reads[read].get_functions()) + '\n')
            of.write(read + ': ' + parser.reads[read].get_status() + ': ' + ','.join(sorted(parser.reads[read].get_functions().keys())))
            if not parser.reads[read].taxonomy is None:
                of.write(' Taxonomy :' + parser.reads[read].taxonomy + '\n')
            else:
                of.write(' No taxonomy\n')
            for hit in parser.reads[read].get_hit_list().get_hits():
                of.write('\t' + str(hit) + '\n')
                
        of.write('\n\n*** End of report ***\n')
        of.closed

def get_function_scores(project, sample_id = None, metrics=None):
    # This function actually returns read counts for readcount metrics,
    # RPK for rpkm and rpkg
    # or FPK for fpkm and fpkg
    # it also returns read count as 'count' metrics, sum of identity % of 
    # all best hits as 'identity' and number of best hits as 'hit_count' 
    # for calculation of average identity %
    # resulting data structure is ret_val[function_id][sample_id][metrics/'count'/'identity'/'hit_count']
    
    ret_val = autovivify(3,float)
    for s in project.list_samples():
        if sample_id is not None and s != sample_id:
            continue
        
        length_cutoff = project.config.get_length_cutoff(project.options.get_collection(s))
        average_read_length = project.samples[s].get_avg_read_length('pe1')
        # Check if reads were processed or imported for this sample
        if project.samples[s].reads is None or 'pe1' not in project.samples[s].reads or len(project.samples[s].reads['pe1']) == 0:
            raise ValueError ('No reads data loaded for sample',s,'end pe1')
        if project.samples[s].is_paired_end:
            if 'pe2' not in project.samples[s].reads or len(project.samples[s].reads['pe2']) == 0:
                raise ValueError ('No reads data loaded for sample',s,'end pe2')

        norm_factor = 0.0
        if metrics == 'readcount':
            norm_factor = 1.0
        elif metrics in ['rpkm','fpkm','erpkm','efpkm']:
            norm_factor = project.samples[s].rpkm_scaling_factor
        elif metrics in ['rpkg','fpkg','erpkg','efpkg']:
            norm_factor = project.samples[s].rpkg_scaling_factor
        
        if norm_factor == 0.0:
            raise ValueError('Cannot get normalization factor')

        # Calculate scores
        if metrics in ['readcount','rpkg','rpkm','erpkg','erpkm']:
            if project.samples[s].is_paired_end:
                raise ValueError('Read count, RPKG and RPKM metrics provided only for single-end sequences')
            
            for read_id, read in project.samples[s].reads['pe1'].items():
                if read.get_status() != 'function': # Filter unmapped reads
                    continue
                read_rpk_scores = read.get_functions()
                for function in read_rpk_scores:
                    if metrics == 'readcount':
                        ret_val[function][s][metrics] += 1.0
                    else:
                        ret_val[function][s]['count'] += 1.0
                        ret_val[function][s][metrics] += norm_factor * read_rpk_scores[function]
                        
                function_maxbitscores = defaultdict(dict)
                hits = [hit for hit in read.get_hit_list().get_hits()]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [function for function in hit.get_functions() if function in read_rpk_scores]:
                        if hit_function in function_maxbitscores:
                            if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                        else:
                            function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                            function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                # Count hits and sum identity
                for function in read_rpk_scores:
                    if function in function_maxbitscores:
                            ret_val[function][s]['hit_count'] += 1.0 
                            ret_val[function][s]['identity'] += function_maxbitscores[function]['identity']
                    else:
                        print('Function',function,'not found in hits of read', read_id)
                        
        elif metrics in ['fpkg','fpkm','efpkg','efpkm']:
            if not project.samples[s].is_paired_end:
                raise ValueError('FPKG and FPKM metrics require paired-end sequences')
            fragment_length = project.get_fragment_length(project.samples[s])
            reads_processed = set()
            
            for read_id in project.samples[s].reads['pe1'].keys():
                read_pe1 = project.samples[s].reads['pe1'][read_id]
                if read_pe1.get_status() != 'function':
                    continue
                reads_processed.add(read_id)
                                        
                if 'pe2' in project.samples[s].reads and read_id in project.samples[s].reads['pe2'].keys(): 
                    read_pe2 = project.samples[s].reads['pe2'][read_id]
                    if read_pe2.get_status() == 'function':
                        # Both ends are mapped
                        
                        fragment_functions = set() # List of functions assigned to the current read

                        read1_functions = read_pe1.get_functions()
                        read2_functions = read_pe2.get_functions()

                        fragment_functions.update(read1_functions.keys())
                        fragment_functions.update(read2_functions.keys())

                        # We filled list of functions. Let's calculate FPKM score
                            
                        #~ for function in fragment_functions:
                            #~ ret_val[function][s]['count'] += 1
                            #~ if function in read1_functions and function in read2_functions: # Take higher score
                                #~ ret_val[function][s][metrics] += norm_factor * (max (read1_functions[function], read2_functions[function]))
                            #~ elif function in read1_functions:
                                #~ ret_val[function][s][metrics] += norm_factor * read1_functions[function]
                            #~ elif function in read2_functions:
                                #~ ret_val[function][s][metrics] += norm_factor * read2_functions[function]

                        function_maxbitscores = defaultdict(dict)
                        hits1 = [hit for hit in read_pe1.get_hit_list().get_hits()]
                        hits2 = [hit for hit in read_pe2.get_hit_list().get_hits()]
                        # Find hit with max. bitscore for each function
                        for hit in hits1:
                            for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                        function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                        function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                        function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                    function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                    function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                        for hit in hits2:
                            for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                        function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                        function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                        function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                    function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                    function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                        # Count hits and sum identity
                        for function in fragment_functions:
                            ret_val[function][s]['count'] += 1
                            if function in function_maxbitscores:
                                    ret_val[function][s]['hit_count'] += 1.0 
                                    ret_val[function][s]['identity'] += function_maxbitscores[function]['identity']
                                    ret_val[function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, fragment_length=fragment_length)
                            else:
                                print('Function',function,'not found in hits of read', read_id)
                                    
                else: #Only end1 is mapped
                    read_functions = read_pe1.get_functions()
                    fragment_functions = set(read_functions.keys())
                    # Count FPKM
                    #~ for function in fragment_functions:
                        #~ ret_val[function][s]['count'] += 1
                        #~ ret_val[function][s][metrics] += norm_factor * read_functions[function]

                    function_maxbitscores = defaultdict(dict)
                    hits1 = [hit for hit in read_pe1.get_hit_list().get_hits()]

                    # Find max. bitscore for each function
                    for hit in hits1:
                        for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                            if hit_function in function_maxbitscores:
                                if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                    function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                    function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                            else:
                                function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                    # Count hits and sum identity
                    for function in fragment_functions:
                        ret_val[function][s]['count'] += 1
                        if function in function_maxbitscores:
                            ret_val[function][s]['hit_count'] += 1.0 
                            ret_val[function][s]['identity'] += function_maxbitscores[function]['identity']
                            ret_val[function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, fragment_length=fragment_length)
                        else:
                            print('Function',function,'not found in hits of read', read_id)

            for read_id in project.samples[s].reads['pe2'].keys():
                if read_id in reads_processed: 
                    continue #Skip read if it was already counted
                    
                read_pe2 = project.samples[s].reads['pe2'][read_id]
                if read_pe2.get_status() != 'function':
                    continue

                fragment_functions = set()
                read_functions = read_pe2.get_functions()
                fragment_functions.update(read_functions.keys())

                # Count FPKM
                #~ for function in fragment_functions:
                    #~ ret_val[function][s]['count'] += 1
                    #~ ret_val[function][s][metrics] += norm_factor * read_functions[function]

                function_maxbitscores = defaultdict(dict)
                hits = [hit for hit in read_pe2.get_hit_list().get_hits()]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                        if hit_function in function_maxbitscores:
                            if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                        else:
                            function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                            function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                            function_maxbitscores[hit_function]['length'] = hit.get_subject_length()

                # Count hits and sum identity
                for function in fragment_functions:
                    ret_val[function][s]['count'] += 1
                    if function in function_maxbitscores:
                        ret_val[function][s]['hit_count'] += 1.0 
                        ret_val[function][s]['identity'] += function_maxbitscores[function]['identity']
                        ret_val[function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, fragment_length=fragment_length)
                    else:
                        print('Function',function,'not found in hits of read', read_id)

    return ret_val

def get_function_taxonomy_scores(project, sample_id = None, metrics=None):
    # This function actually returns read counts for readcount metrics,
    # RPK for rpkm and rpkg
    # or FPK for fpkm and fpkg
    # it also returns read count as 'count' metrics, sum of identity % of 
    # all best hits as 'identity' and number of best hits as 'hit_count' 
    # for calculation of average identity %
    # resulting data structure is ret_val[taxonomy_id][function_id][sample_id][metrics/'count'/'identity'/'hit_count']
    
    ret_val = autovivify(4,float)
    for s in project.list_samples():
        if sample_id is not None and s != sample_id:
            continue

        # Check if reads were processed or imported for this sample
        if project.samples[s].reads is None or 'pe1' not in project.samples[s].reads or len(project.samples[s].reads['pe1']) == 0:
            raise ValueError ('No reads data loaded for sample',s,'end pe1')
        if project.samples[s].is_paired_end:
            if 'pe2' not in project.samples[s].reads or len(project.samples[s].reads['pe2']) == 0:
                raise ValueError ('No reads data loaded for sample',s,'end pe2')

        norm_factor = 0.0
        if metrics == 'readcount':
            norm_factor = 1.0
        elif metrics in ['rpkm','fpkm','erpkm','efpkm']:
            norm_factor = project.samples[s].rpkm_scaling_factor
        elif metrics in ['rpkg','fpkg','erpkg','efpkg']:
            norm_factor = project.samples[s].rpkg_scaling_factor
        
        if norm_factor == 0.0:
            raise ValueError('Cannot get normalization factor')


        # Calculate scores
        if metrics in ['readcount','rpkg','rpkm','erpkg','erpkm']:
            if project.samples[s].is_paired_end:
                raise ValueError('Read count, RPKG and RPKM metrics provided only for single-end sequences')
            
            for read_id, read in project.samples[s].reads['pe1'].items():
                if read.get_status() != 'function': # Filter unmapped reads
                    continue
                read_rpk_scores = read.get_functions()
                for function in read_rpk_scores:
                    ret_val[read.taxonomy][function][s]['count'] += 1.0
                    ret_val[read.taxonomy][function][s][metrics] += norm_factor * read_rpk_scores[function]
                        
                function_maxbitscores = {}
                hits = [hit for hit in read.get_hit_list().get_hits()]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [function for function in hit.get_functions() if function in read.functions]:
                        if hit_function in function_maxbitscores:
                            if hit.get_bitscore() > function_maxbitscores[hit_function]:
                                function_maxbitscores[hit_function] = hit.get_bitscore()
                        else:
                            function_maxbitscores[hit_function] = hit.get_bitscore()
                # Count hits and sum identity
                for hit in hits:
                    for hit_function in [function for function in hit.get_functions() if function in read.functions]:
                        if hit_function in function_maxbitscores and hit.get_bitscore() == function_maxbitscores[hit_function]:
                            ret_val[read.taxonomy][function][s]['hit_count'] += 1.0 
                            ret_val[read.taxonomy][function][s]['identity'] += hit.get_identity()
                            del function_maxbitscores[hit_function]
        
        elif metrics in ['fpkg','fpkm','efpkg','efpkm']:
            if not project.samples[s].is_paired_end:
                raise ValueError('FPKG and FPKM metrics require paired-end sequences')
            fragment_length = project.get_fragment_length(project.samples[s])
            reads_processed = set()
            length_cutoff = project.config.get_length_cutoff(project.options.get_collection(s))
            average_read_length = project.samples[s].get_avg_read_length('pe1')
            
            for read_id in project.samples[s].reads['pe1'].keys():
                read_pe1 = project.samples[s].reads['pe1'][read_id]
                if read_pe1.get_status() != 'function':
                    continue
                reads_processed.add(read_id)
                                        
                if 'pe2' in project.samples[s].reads and read_id in project.samples[s].reads['pe2'].keys(): 
                    read_pe2 = project.samples[s].reads['pe2'][read_id]
                    if read_pe2.get_status() == 'function':
                        # Both ends are mapped
                        fragment_taxonomy = project.taxonomy_data.get_lca([read_pe1.taxonomy, read_pe2.taxonomy])
                        
                        fragment_functions = set() # List of functions assigned to the current read

                        read1_functions = read_pe1.get_functions()
                        read2_functions = read_pe2.get_functions()

                        fragment_functions.update(read1_functions.keys())
                        fragment_functions.update(read2_functions.keys())

                        # We filled list of functions. Let's calculate FPKM score
                            
                        #~ for function in fragment_functions:
                            #~ ret_val[fragment_taxonomy][function][s]['count'] += 1
                            #~ if function in read1_functions and function in read2_functions: # Take higher score
                                #~ ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * (max (read1_functions[function], read2_functions[function]))
                            #~ elif function in read1_functions:
                                #~ ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * read1_functions[function]
                            #~ elif function in read2_functions:
                                #~ ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * read2_functions[function]

                        # Let's get hits data for calculation of identity% 
                        
                        function_maxbitscores = defaultdict(dict)
                        hits1 = [hit for hit in read_pe1.get_hit_list().get_hits()]
                        hits2 = [hit for hit in read_pe2.get_hit_list().get_hits()]
                        # Find max. bitscore for each function
                        for hit in hits1:
                            for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                        function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                        function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                        function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                    function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                    function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                        for hit in hits2:
                            for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                        function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                        function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                        function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                    function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                    function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                        # Count hits and sum identity for best hits for each function
                        for function in fragment_functions:
                            if function in function_maxbitscores:
                                ret_val[fragment_taxonomy][function][s]['count'] += 1
                                ret_val[fragment_taxonomy][function][s]['hit_count'] += 1.0 
                                ret_val[fragment_taxonomy][function][s]['identity'] += function_maxbitscores[function]['identity']
                                ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, fragment_length=fragment_length)
                            else:
                                print('Function',function,'not found in hits of read', read_id)


                else: #Only end1 is mapped
                    fragment_taxonomy = read_pe1.taxonomy
                    read_functions = read_pe1.get_functions()
                    fragment_functions = set(read_functions.keys())
                    # Count FPKM
                    #~ for function in fragment_functions:
                        #~ ret_val[fragment_taxonomy][function][s]['count'] += 1
                        #~ ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * read_functions[function]
                    function_maxbitscores = defaultdict(dict)
                    hits1 = [hit for hit in read_pe1.get_hit_list().get_hits()]
                    # Find max. bitscore for each function
                    for hit in hits1:
                        for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                            if hit_function in function_maxbitscores:
                                if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                    function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                    function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                            else:
                                function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                    # Count hits and sum identity for best hits for each function
                    for function in fragment_functions:
                        if function in function_maxbitscores:
                            ret_val[fragment_taxonomy][function][s]['count'] += 1
                            ret_val[fragment_taxonomy][function][s]['hit_count'] += 1.0 
                            ret_val[fragment_taxonomy][function][s]['identity'] += function_maxbitscores[function]['identity']
                            ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, fragment_length=fragment_length)
                        else:
                            print('Function',function,'not found in hits of read', read_id)

            for read_id in project.samples[s].reads['pe2'].keys():
                if read_id in reads_processed: 
                    continue #Skip read if it was already counted
                    
                read_pe2 = project.samples[s].reads['pe2'][read_id]
                if read_pe2.get_status() != 'function':
                    continue

                fragment_functions = set()
                fragment_taxonomy = read_pe2.taxonomy
                read_functions = read_pe2.get_functions()
                fragment_functions.update(read_functions.keys())

                # Count FPKM
                #~ for function in fragment_functions:
                    #~ ret_val[fragment_taxonomy][function][s]['count'] += 1
                    #~ ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * read_functions[function]
                # Build FPKM-based functional taxonomic profile
                function_maxbitscores = defaultdict(dict)
                hits = [hit for hit in read_pe2.get_hit_list().get_hits()]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [function for function in hit.get_functions() if function in fragment_functions]:
                        if hit_function in function_maxbitscores:
                            if hit.get_bitscore() > function_maxbitscores[hit_function]['bitscore']:
                                function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                                function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                                function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                        else:
                            function_maxbitscores[hit_function]['bitscore'] = hit.get_bitscore()
                            function_maxbitscores[hit_function]['identity'] = hit.get_identity()
                            function_maxbitscores[hit_function]['length'] = hit.get_subject_length()
                # Count hits and sum identity for best hits for each function
                for function in fragment_functions:
                    if function in function_maxbitscores:
                        ret_val[fragment_taxonomy][function][s]['count'] += 1
                        ret_val[fragment_taxonomy][function][s]['hit_count'] += 1.0 
                        ret_val[fragment_taxonomy][function][s]['identity'] += function_maxbitscores[function]['identity']
                        ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, fragment_length=fragment_length)
                    else:
                        print('Function',function,'not found in hits of read', read_id)
    return ret_val

def generate_sample_report(project, sample_id):
    #outfile = os.path.join(project.samples[sample_id].work_directory, sample_id + '_report.txt')
    outfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), sample_id + '_report.txt'))
    with open(outfile, 'w') as of:
        of.write('Report for ' + sample_id + '\n\n')
        of.write('Sample ID:\t' + project.options.get_sample_name(sample_id) + '\n')
        of.write('Replicate:\t' + str(project.samples[sample_id].replicate) + '\n\n')
        of.write('Project name:\t' + project.options.get_name() + '\n\n')
        of.write('Reference data set:\t' + project.options.get_collection(sample_id) + '\n\n')
        
        of.write('Source files\n')
        of.write('FASTQ file 1:\t' + project.options.get_fastq_path(sample_id, 'pe1') + '\n')
        of.write('    Number of reads:\t' + str(project.samples[sample_id].fastq_fwd_readcount) + '\n')
        of.write('    Number of bases:\t' + str(project.samples[sample_id].fastq_fwd_basecount) + '\n')
        
        if project.samples[sample_id].is_paired_end:
            of.write('FASTQ file 2:\t' + project.options.get_fastq_path(sample_id, 'pe2') + '\n')
            of.write('    Number of reads:\t' + str(project.samples[sample_id].fastq_rev_readcount) + '\n')
            of.write('    Number of bases:\t' + str(project.samples[sample_id].fastq_rev_basecount) + '\n')
        
        of.write('\nNormalization\n')
        if not project.samples[sample_id].rpkm_scaling_factor is None:
            of.write('RPKM normalization factor:\t' + str(project.samples[sample_id].rpkm_scaling_factor) + '\n')
        if not project.samples[sample_id].rpkg_scaling_factor is None:
            of.write('Average genome size:\t' + format(project.samples[sample_id].rpkg_scaling_factor * project.samples[sample_id].fastq_fwd_basecount, "0.0f") + '\n')
            of.write('RPKG normalization factor:\t' + str(project.samples[sample_id].rpkg_scaling_factor) + '\n')

        of.write('\nNumber of mapped reads\n')
        of.write('FASTQ file 1:\t' + str(len(project.samples[sample_id].reads['pe1'])) + '\n')
        if project.samples[sample_id].is_paired_end:
            of.write('FASTQ file 2:\t' + str(len(project.samples[sample_id].reads['pe2'])) + '\n')

        of.closed
    
    if project.samples[sample_id].is_paired_end:
        if not project.samples[sample_id].rpkg_scaling_factor is None:
            metrics = 'efpkg'
        elif not project.samples[sample_id].rpkm_scaling_factor is None:
            metrics = 'efpkm'
        else:
            metrics = 'readcount'
    else:
        if not project.samples[sample_id].rpkg_scaling_factor is None:
            metrics = 'erpkg'
        elif not project.samples[sample_id].rpkm_scaling_factor is None:
            metrics = 'erpkm'
        else:
            metrics = 'readcount'

    scores_function = get_function_scores(project, sample_id, metrics=metrics)
    scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id, metrics=metrics)
    generate_function_sample_xlsx(project, 
                                scores_function, 
                                metrics=metrics, 
                                sample_id = sample_id)
    generate_function_taxonomy_sample_xlsx(project, scores_function_taxonomy, metrics=metrics, sample_id = sample_id)
    sample_scores_taxonomy = autovivify(3, float)
    for tax in scores_function_taxonomy.keys():
        for f in scores_function_taxonomy[tax].keys():
            if sample_id in scores_function_taxonomy[tax][f]:
                for k,v in scores_function_taxonomy[tax][f][sample_id].items():
                    sample_scores_taxonomy[tax][f][k] = v

    tax_profile = TaxonomyProfile()
    outfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), sample_id + '_' + metrics + '_functional_taxonomy_profile.xml'))
    tax_profile.build_functional_taxonomy_profile(project.taxonomy_data, sample_scores_taxonomy)
    function_list = sorted(project.ref_data.functions_dict.keys())
    generate_taxonomy_series_chart(tax_profile, function_list, outfile, score=metrics)

    
def generate_project_report(project):
    is_paired_end = None
    for sample_id in project.list_samples():
        if project.samples[sample_id].is_paired_end:
            if is_paired_end is None:
                is_paired_end = True
            elif not is_paired_end:
                print ('Project contains both single-end and paired-end sequences. No comparative tables will be generated.')
                return
        else:
            if is_paired_end is None:
                is_paired_end = False
            elif is_paired_end:
                print ('Project contains both single-end and paired-end sequences. No comparative tables will be generated.')
                return

        # If there are no read data, try to load them from JSON
        if len(project.samples[sample_id].reads) == 0:
            project.import_reads_json(sample_id, project.ENDS)
    print ('Generating spreadsheets and interactive diagrams for the project...')
    if is_paired_end:
        metrics = 'efpkm'
        scores = get_function_scores(project, sample_id=None, metrics=metrics)
        generate_function_sample_xlsx(project, scores, metrics=metrics)
        scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id=None, metrics=metrics)
        generate_function_taxonomy_sample_xlsx(project, scores_function_taxonomy, metrics=metrics, sample_id=None)
        generate_sample_taxonomy_function_xlsx(project, scores_function_taxonomy, metrics=metrics, function_id=None)
        metrics = 'efpkg'
        scores = get_function_scores(project, sample_id=None, metrics=metrics)
        generate_function_sample_xlsx(project, scores, metrics=metrics)
        scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id=None, metrics=metrics)
        generate_function_taxonomy_sample_xlsx(project, scores_function_taxonomy, metrics=metrics, sample_id=None)
        generate_sample_taxonomy_function_xlsx(project, scores_function_taxonomy, metrics=metrics, function_id=None)
    else:
        metrics = 'erpkm'
        scores = get_function_scores(project, sample_id=None, metrics=metrics)
        generate_function_sample_xlsx(project, scores, metrics=metrics)
        scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id=None, metrics=metrics)
        generate_function_taxonomy_sample_xlsx(project, scores_function_taxonomy, metrics=metrics, sample_id=None)
        generate_sample_taxonomy_function_xlsx(project, scores_function_taxonomy, metrics=metrics, function_id=None)
        metrics = 'erpkg'
        scores = get_function_scores(project, sample_id=None, metrics=metrics)
        generate_function_sample_xlsx(project, scores, metrics=metrics)
        scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id=None, metrics=metrics)
        generate_function_taxonomy_sample_xlsx(project, scores_function_taxonomy, metrics=metrics, sample_id=None)
        generate_sample_taxonomy_function_xlsx(project, scores_function_taxonomy, metrics=metrics, function_id=None)
        metrics = 'readcount'
        scores = get_function_scores(project, sample_id=None, metrics=metrics)
        generate_function_sample_xlsx(project, scores, metrics=metrics)
        scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id=None, metrics=metrics)
        generate_function_taxonomy_sample_xlsx(project, scores_function_taxonomy, metrics=metrics, sample_id=None)
        generate_sample_taxonomy_function_xlsx(project, scores_function_taxonomy, metrics=metrics, function_id=None)

def generate_protein_project_report(project):
    for sample_id in project.list_samples():
        # If there are no read data, try to load them from JSON
        if len(project.samples[sample_id].reads) == 0:
            project.import_reads_json(sample_id, project.ENDS)
    print ('Generating spreadsheets and interactive diagrams for the project...')
    outfile = os.path.join(project.options.get_work_dir(), 'project_report.txt')
    with open (outfile, 'w') as of:
        of.write(project.options.get_name() + '\n\n')
        for sample_id in project.list_samples():
            of.write(sample_id + ':\t' + project.samples[sample_id].sample_name + '\tproteins mapped: ' + str(len(project.samples[sample_id].reads['pe1'])))
            of.write('\n')
        of.closed

    metrics = 'readcount'
    scores = get_function_scores(project, sample_id=None, metrics=metrics)
    generate_function_sample_xlsx(project, scores, metrics=metrics)
    scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id=None, metrics=metrics)
    generate_function_taxonomy_sample_xlsx(project, scores_function_taxonomy, metrics=metrics, sample_id=None)
    generate_sample_taxonomy_function_xlsx(project, scores_function_taxonomy, metrics=metrics, function_id=None)


def generate_functions_scores_table(project):
    functions_list = set()
    scores = defaultdict(lambda : defaultdict(float))
    read_counts = defaultdict(lambda : defaultdict(float))
    
    # fill tables of scores and read counts
    for sample in project.list_samples():
        for end in ENDS:
            project.import_reads_json(sample, [end,]) # Lazy load
            scaling_factor = 1.0
            if end == 'pe1':
                scaling_factor = project.options.get_fastq1_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            elif end == 'pe2':
                scaling_factor = project.options.get_fastq2_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            else:
                raise Exception('Unknown identifier of read end: ' + end)
            for read_id in project.samples[sample].reads[end]:
                read = project.samples[sample].reads[end][read_id]
                if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                    for function in read.functions:
                        functions_list.add(function)
                        scores[function][sample] += scaling_factor * read.functions[function]
                        read_counts[function][sample] += 1.0/len(read.functions)
            project.samples[sample].reads[end] = None

    # generate output
    samples_list = sorted(project.samples.keys())
    lines = ['Function\t' + '\t'.join(samples_list) + '\tDefinition',]
    for function in sorted(functions_list):
        line = function
        for sample in samples_list:
            if sample in scores[function]:
                line += '\t' + '{0:.2f}'.format(scores[function][sample])
            else:
                line += '\t0.00'
        line += '\t' + project.ref_data.lookup_function_name(function)
        lines.append(line)
    
    lines.append('\nFunction\t' + '\t'.join(samples_list) + '\tDefinition')
    for function in sorted(functions_list):
        line = function
        for sample in samples_list:
            if sample in read_counts[function]:
                line += '\t' + str(read_counts[function][sample])
            else:
                line += '\t0'
        line += '\t' + project.ref_data.lookup_function_name(function)
        lines.append(line)
    
    return '\n'.join(lines)
    
def generate_functions_scores_list(project):
    functions = defaultdict(dict)
    protein_counts = defaultdict(dict)
    # initialize list of functions
    for sample in project.list_samples():
        for end in project.samples[sample].reads:
            for read in project.samples[sample].reads[end]:
                if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                    for function in read.functions:
                        if sample in functions[function]:
                            functions[function][sample] += read.functions[function]
                            protein_counts[function][sample] += 1
                        else:
                            functions[function][sample] = read.functions[function]
                            protein_counts[function][sample] = 1
    
    samples_list = sorted(project.samples.keys())
    lines = ['Function\tSample\tScore\tDefinition',]
    for sample in samples_list:
        for function in sorted(functions.keys()):
            line = sample
            if sample in functions[function]:
                line +=  '\t' + function + '\t' + '{0:.2f}'.format(functions[function][sample])
            else:
                line +=  '\t' + function + '\t0.00'
            line += '\t' + project.ref_data.lookup_function_name(function)
            lines.append(line)
    
    lines.append('\nFunction\tSample\tCount\tDefinition')
    for sample in samples_list:
        for function in sorted(functions.keys()):
            line = sample
            if sample in functions[function]:
                line +=  '\t' + function + '\t' + str(protein_counts[function][sample])
            else:
                line +=  '\t' + function + '\t0'
            line += '\t' + project.ref_data.lookup_function_name(function)
            lines.append(line)
    
    return '\n'.join(lines)


def create_functions_markdown_document(project):

    functions_list = set()
    categories_list = set()
    scores = defaultdict(lambda : defaultdict(float))
    scores_cat = defaultdict(lambda : defaultdict(float))
    
    # calculate scores

    for sample in project.list_samples():
        for end in project.ENDS:
            project.samples[sample].reads[end] = None
            if end == 'pe2' and not project.samples[sample].is_paired_end:
                continue
            project.import_reads_json(sample, [end,])
            scaling_factor = 1.0
            if end == 'pe1':
                scaling_factor = project.options.get_fastq1_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            elif end == 'pe2':
                scaling_factor = project.options.get_fastq2_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            else:
                raise Exception('Unknown identifier of read end: ' + end)
            for read_id in project.samples[sample].reads[end]:
                read = project.samples[sample].reads[end][read_id]
                if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                    for function in read.functions:
                        functions_list.add(function)
                        scores[function][sample] += scaling_factor * read.functions[function]
                        category = project.ref_data.lookup_function_group(function)
                        categories_list.add(category)
                        scores_cat[category][sample] += scaling_factor * read.functions[function]
            project.samples[sample].reads[end] = None

    # write output
    samples_list = sorted(project.list_samples())
    outfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), 'index.md'))

    with open(outfile, 'w') as of:
        of.write('# ')
        of.write(project.options.get_name())
        of.write('\n\n')

        of.write('## Categories\n\nAbundance of functional categories (RPKM score)\n\n')

        of.write('Category|')
        of.write('|'.join(samples_list))
        of.write('|\n')
        of.write('|---'*(len(samples_list) + 1))
        of.write('|\n')
        
        for category in sorted(categories_list):
            of.write(category)
            for sample in samples_list:
                of.write('|')
                if sample in scores_cat[category]:
                    of.write('{0:.3f}'.format(scores_cat[category][sample]))
                else:
                    of.write('0.000')
            of.write('|\n')

        of.write('\n## Functions\n\nAbundance of individual functions (RPKM score)\n\n')


        of.write('Function|')
        of.write('|'.join(samples_list))
        of.write('|Definition|\n')
        of.write('|---'*(len(samples_list) + 2))
        of.write('|\n')
        
        for function in sorted(functions_list):
            of.write(function)
            for sample in samples_list:
                of.write('|')
                if sample in scores[function]:
                    of.write('{0:.3f}'.format(scores[function][sample]))
                else:
                    of.write('0.000')
            of.write('|')
            of.write(project.ref_data.lookup_function_name(function))
            of.write('|\n')

        targetfile = sanitize_file_name(os.path.join('data', project.options.get_name() + '_functions.xlsx'))

        of.write('\n<a href="')
        of.write(targetfile)
        of.write('" target="_blank">Download table of RPKM scores and read counts in XLSX format</a>\n\n')

        of.write('## Taxonomy profile for all reads mapped to nitrogen cycle genes\n\n')

        targetfile = sanitize_file_name(os.path.join('data', project.options.get_name() + '_taxonomy_profile.xml.html'))
        of.write('<div class="krona-wrapper">\n<iframe src="')
        of.write(targetfile)
        of.write('" height="800" width="100%">')
        of.write(project.options.get_name())
        of.write(' taxonomy profile</iframe>\n<br>\n<a href="')
        of.write(targetfile)
        of.write('" target="_blank">Open chart in a new window</a>\n</div>\n\n')
        
        of.write('## Taxonomy profiles for individual functions\n\n')

        targetfile = sanitize_file_name(os.path.join('data', project.options.get_name() + '_functions_taxonomy.xlsx'))
        of.write('<a href="')
        of.write(targetfile)
        of.write('" target="_blank">Download detailed taxonomic profile for all functions in all samples (XLSX format)</a>\n\n<div>\n')
        for sample in samples_list:
            targetfile = sanitize_file_name(os.path.join('data', sample + '_functional_taxonomy_profile.xml.html'))
            of.write('<a href="')
            of.write(targetfile)
            of.write('" target="_blank">Taxonomy profile for individual functions in sample ')
            of.write(sample)
            of.write(' (interactive chart)</a><br>\n')
        of.write('</div>\n')
        of.closed

    # write files for samples
    for sample in samples_list:
        outfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), sample + '.md'))
        with open(outfile, 'w') as of:
            of.write('# Sample ')
            of.write(sample)
            of.write('\n\n## Functional taxonomy profile for reads mapped to ')
            of.write(os.path.join(project.options.get_collection()))
            of.write(' dataset\n\n')
            
            targetfile = sanitize_file_name(os.path.join('..','data', sample + '_functional_taxonomy_profile.xml.html'))
            of.write('<div class="krona-wrapper">\n<iframe src="')
            of.write(targetfile)
            of.write('" height="800" width="100%">')
            of.write(project.options.get_name())
            of.write(' taxonomy profile</iframe>\n<br>\n<a href="')
            of.write(targetfile)
            of.write('" target="_blank">Open chart in a new window</a>\n</div>\n\n')
            
            of.write('## Reports for FASTQ files\n\n')
            if os.path.exists(os.path.join(project.options.get_project_dir(sample), project.options.get_output_subdir(sample),sample + '_pe1_'+ project.options.get_report_name() + '.pdf')):
                targetfile = os.path.join('..','data', sample + '_pe1_'+ project.options.get_report_name() + '.pdf')
                of.write('<a href="')
                of.write(targetfile)
                of.write('">Download report for read end 1 (PDF format)</a>\n<br>\n')
            if os.path.exists(os.path.join(project.options.get_project_dir(sample), project.options.get_output_subdir(sample),sample + '_pe2_'+ project.options.get_report_name() + '.pdf')):
                targetfile = os.path.join('..','data', sample + '_pe2_'+ project.options.get_report_name() + '.pdf')
                of.write('<a href="')
                of.write(targetfile)
                of.write('">Download report for read end 2 (PDF format)</a>\n<br>\n')
            of.closed
            


def generate_assembly_report(assembler):
    outfile = os.path.join(assembler.assembly_dir, 'out', 'assembly_report.txt')
    with open(outfile, 'w') as of:
        of.write('\nFunction statistics\n\n')
        for function in assembler.assembly.contigs:
            of.write(function + ':\n')
            for contig in assembler.assembly.contigs[function]:
                for gene_id in assembler.assembly.contigs[function][contig].genes:
                    gene = assembler.assembly.contigs[function][contig].genes[gene_id]
                    if gene.get_status() == 'function,besthit' or gene.get_status() == 'function':
                        of.write(function + '\t' + 
                                contig + '\t' + 
                                str(len(assembler.assembly.contigs[function][contig].sequence)) + '\t' + 
                                gene_id + '\t' + 
                                ','.join(gene.functions.keys()))
                        for sample in assembler.assembly.contigs[function][contig].read_count.keys():
                            of.write('\t' + sample + ':read count ' + str(assembler.assembly.contigs[function][contig].read_count[sample]) +
                                    ';coverage ' + str(assembler.assembly.contigs[function][contig].get_coverage(sample))  +
                                    ';rpkm ' + str(assembler.assembly.contigs[function][contig].get_rpkm(sample, assembler.project.options.get_fastq1_readcount(sample)))
                                    )
                        of.write('\n')
        
        of.write('\n\n')
        for function in assembler.assembly.reads:
            of.write(function + '\t' + str(len(assembler.assembly.reads[function])) + '\n')
                
        of.write('\n\n*** End of report ***\n')
        of.closed

