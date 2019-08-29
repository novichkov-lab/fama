import os
from collections import defaultdict,Counter,OrderedDict

from Fama import const
from Fama.utils import autovivify,cleanup_protein_id,sanitize_file_name
from Fama.DiamondParser.hit_utils import get_efpk_score,get_fpk_score
from Fama.TaxonomyProfile import TaxonomyProfile
from Fama.OutputUtil.XlsxUtil import generate_function_sample_xlsx, generate_function_taxonomy_sample_xlsx, generate_sample_taxonomy_function_xlsx
from Fama.OutputUtil.KronaXMLWriter import generate_taxonomy_series_chart

#ENDS = ['pe1','pe2']
def generate_fastq_report(parser):
    outfile = os.path.join(parser.options.get_project_dir(parser.sample.sample_id), parser.options.get_output_subdir(parser.sample.sample_id),parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.report_name)
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
            read_stats[parser.reads[read].status] += 1
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
            if parser.reads[read].status == 'function':
                functions = parser.reads[read].functions
                for function in functions:
                    func_stats[function] += functions[function]
                    func_counts[function] += 1/len(functions)
                for hit in parser.reads[read].hit_list.hits:
                    for function in hit.functions:
                        func_identity[function] += hit.identity
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
            if parser.reads[read].status == 'function':
                functions = parser.reads[read].functions
                for function in functions:
                    func_stats[parser.ref_data.lookup_function_group(function)] += functions[function]
                    func_counts[parser.ref_data.lookup_function_group(function)] += 1/len(functions)
                for hit in parser.reads[read].hit_list.hits:
                    for function in hit.functions:
                        func_identity[parser.ref_data.lookup_function_group(function)] += hit.identity
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
#            print (read, parser.reads[read].status)
            if parser.reads[read].status == 'function':
                taxonomy = parser.reads[read].taxonomy
                if taxonomy is None:
                    print ('No taxonomy ID assigned to ', read)
                    continue
                tax_stats[taxonomy] += 1
                read_functions = parser.reads[read].functions
                for f in read_functions:
                    rpkm_stats[taxonomy] += read_functions[f]
                identity_stats[taxonomy] += sum(list(hit.identity for hit in parser.reads[read].hit_list.hits))
                        
        tax_data = parser.taxonomy_data
        counts_per_rank, identity_per_rank, rpkm_per_rank = tax_data.get_taxonomy_profile(counts=tax_stats, identity=identity_stats, scores = rpkm_stats)

        ranks = const.RANKS[1:]
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
            of.write(read + ': ' + parser.reads[read].status + ': ' + ','.join(sorted(parser.reads[read].functions.keys())))
            if not parser.reads[read].taxonomy is None:
                of.write(' Taxonomy :' + parser.reads[read].taxonomy + '\n')
            else:
                of.write(' No taxonomy\n')
            for hit in parser.reads[read].hit_list.hits:
                of.write('\t' + str(hit) + '\n')
                
        of.write('\n\n*** End of report ***\n')
        of.closed

def generate_fasta_report(parser):
    outfile = os.path.join(parser.options.get_project_dir(parser.sample.sample_id), parser.options.get_output_subdir(parser.sample.sample_id),parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.report_name)
    with open(outfile, 'w') as of:
        # Write general info
        of.write('\nRun info\n\n')
        of.write('Sample ID:\t' + parser.options.get_sample_name(parser.sample.sample_id) + '\n')
        of.write('Sequence file:\t' + parser.options.get_fastq_path(parser.sample.sample_id, parser.end) + '\n')
        of.write('Total number of reads:\t' + str(parser.sample.fastq_fwd_readcount) + '\n')
        of.write('*****************************************\n\n')

        # Write read statistics
        of.write('\nRead statistics\n\n')
        read_stats = Counter()
        for read in sorted(parser.reads.keys()):
            read_stats[parser.reads[read].status] += 1
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
            if parser.reads[read].status == 'function':
                functions = parser.reads[read].functions
                for function in functions:
                    func_stats[function] += functions[function]
                    func_counts[function] += 1/len(functions)
                for hit in parser.reads[read].hit_list.hits:
                    for function in hit.functions:
                        func_identity[function] += hit.identity
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
            if parser.reads[read].status == 'function':
                functions = parser.reads[read].functions
                for function in functions:
                    func_stats[parser.ref_data.lookup_function_group(function)] += functions[function]
                    func_counts[parser.ref_data.lookup_function_group(function)] += 1/len(functions)
                for hit in parser.reads[read].hit_list.hits:
                    for function in hit.functions:
                        func_identity[parser.ref_data.lookup_function_group(function)] += hit.identity
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
            if parser.reads[read].status == 'function':
                taxonomy = parser.reads[read].taxonomy
                if taxonomy is None:
                    print ('No taxonomy ID assigned to ', read)
                    continue
                tax_stats[taxonomy] += 1
                read_functions = parser.reads[read].functions
                for f in read_functions:
                    rpkm_stats[taxonomy] += read_functions[f]
                identity_stats[taxonomy] += sum(list(hit.identity for hit in parser.reads[read].hit_list.hits))
                    
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
            of.write(read + ': ' + parser.reads[read].status + ': ' + ','.join(sorted(parser.reads[read].functions.keys())))
            if not parser.reads[read].taxonomy is None:
                of.write(' Taxonomy :' + parser.reads[read].taxonomy + '\n')
            else:
                of.write(' No taxonomy\n')
            for hit in parser.reads[read].hit_list.hits:
                of.write('\t' + str(hit) + '\n')
                
        of.write('\n\n*** End of report ***\n')
        of.closed

def get_function_scores(project, sample_id = None, metrics=None):
    # This function actually returns read counts for readcount metrics,
    # RPK for rpkm and rpkg
    # or FPK for fpkm and fpkg
    # it also always returns read count as 'count' metrics, sum of identity % of 
    # all best hits as 'identity' and number of best hits as 'hit_count' 
    # for calculation of average identity %
    # resulting data structure is ret_val[function_id][sample_id][metrics|'count'|'identity'|'hit_count']
    
    ret_val = autovivify(3,float)
    for s in project.list_samples():
        if sample_id is not None and s != sample_id:
            continue
        
        length_cutoff = project.config.get_length_cutoff(project.options.get_collection(s))
        average_read_length = project.samples[s].get_avg_read_length('pe1')
        # Check if reads were processed or imported for this sample
        if project.samples[s].reads is None or 'pe1' not in project.samples[s].reads:
            raise ValueError ('No reads data loaded for sample',s,'end pe1')
        if project.samples[s].is_paired_end:
            if 'pe2' not in project.samples[s].reads:
                raise ValueError ('No reads data loaded for sample',s,'end pe2')

        norm_factor = 0.0
        if metrics in ['readcount','erpk','fragmentcount','fpk','efpk']:
            norm_factor = 1.0
        elif metrics in ['fpkm','erpkm','efpkm']:
            norm_factor = project.samples[s].rpkm_scaling_factor
        elif metrics in ['fpkg','erpkg','efpkg']:
            norm_factor = project.samples[s].rpkg_scaling_factor
        else:
            raise ValueError('Unknown metrics:' + metrics)
        
        if norm_factor == 0.0:
            raise ValueError('Cannot get normalization factor')

        # Calculate scores
        if metrics in ['readcount','erpk','erpkg','erpkm']:
            if project.samples[s].is_paired_end:
                raise ValueError('Read count, RPKG and RPKM metrics provided only for single-end sequences')
            
            for read_id, read in project.samples[s].reads['pe1'].items():
                if read.status != 'function': # Filter unmapped reads
                    continue
                read_erpk_scores = read.functions
                for function in read_erpk_scores:
                    if metrics == 'readcount':
                        ret_val[function][s]['count'] += 1.0
                        ret_val[function][s][metrics] += 1.0
                    else:
                        ret_val[function][s]['count'] += 1.0
                        ret_val[function][s][metrics] += norm_factor * read_erpk_scores[function]
                        
                function_maxbitscores = defaultdict(dict)
                # Find max. bitscore for each function
                for hit in read.hit_list.hits:
                    for hit_function in [function for function in hit.functions if function in read_erpk_scores]:
                        if hit_function in function_maxbitscores:
                            if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                function_maxbitscores[hit_function]['identity'] = hit.identity
                        else:
                            function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                            function_maxbitscores[hit_function]['identity'] = hit.identity
                # Count hits and sum identity
                for function in read_erpk_scores:
                    if function in function_maxbitscores:
                            ret_val[function][s]['hit_count'] += 1.0 
                            ret_val[function][s]['identity'] += function_maxbitscores[function]['identity']
                    else:
                        print('Function',function,'not found in hits of read', read_id)
                        
        elif metrics in ['fragmentcount','fpk','efpk','fpkg','fpkm','efpkg','efpkm']:
            if not project.samples[s].is_paired_end:
                raise ValueError('Metrics based on fragment count require paired-end sequences')
            insert_size = project.get_insert_size(project.samples[s])
                
            reads_processed = set()
            
            for read_id in project.samples[s].reads['pe1'].keys():
                read_pe1 = project.samples[s].reads['pe1'][read_id]
                if read_pe1.status != 'function':
                    continue
                reads_processed.add(read_id)
                                        
                if 'pe2' in project.samples[s].reads and read_id in project.samples[s].reads['pe2'].keys(): 
                    read_pe2 = project.samples[s].reads['pe2'][read_id]
                    if read_pe2.status == 'function':
                        # Both ends are mapped
                        
                        fragment_functions = set() # List of functions assigned to the current read

                        read1_functions = read_pe1.functions
                        read2_functions = read_pe2.functions

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
                        hits1 = [hit for hit in read_pe1.hit_list.hits]
                        hits2 = [hit for hit in read_pe2.hit_list.hits]
                        # Find hit with max. bitscore for each function
                        for hit in hits1:
                            for hit_function in [function for function in hit.functions if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                        function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                        function_maxbitscores[hit_function]['identity'] = hit.identity
                                        function_maxbitscores[hit_function]['length'] = hit.s_len
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                        for hit in hits2:
                            for hit_function in [function for function in hit.functions if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                        function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                        function_maxbitscores[hit_function]['identity'] = hit.identity
                                        function_maxbitscores[hit_function]['length'] = hit.s_len
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                        # Count hits and sum identity
                        for function in fragment_functions:
                            ret_val[function][s]['count'] += 1
                            if function in function_maxbitscores:
                                    ret_val[function][s]['hit_count'] += 1.0 
                                    ret_val[function][s]['identity'] += function_maxbitscores[function]['identity']
                                    if metrics == 'fragmentcount':
                                        ret_val[function][s][metrics] += 1
                                    elif metrics in ['fpk','fpkg','fpkm']:
                                        ret_val[function][s][metrics] += norm_factor * get_fpk_score(function_maxbitscores[function]['length'])
                                    else:
                                        ret_val[function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, insert_size=insert_size)
                            else:
                                print('Function',function,'not found in hits of read', read_id)
                                    
                else: #Only end1 is mapped
                    read_functions = read_pe1.functions
                    fragment_functions = set(read_functions.keys())
                    # Count FPKM
                    #~ for function in fragment_functions:
                        #~ ret_val[function][s]['count'] += 1
                        #~ ret_val[function][s][metrics] += norm_factor * read_functions[function]

                    function_maxbitscores = defaultdict(dict)
                    hits1 = [hit for hit in read_pe1.hit_list.hits]

                    # Find max. bitscore for each function
                    for hit in hits1:
                        for hit_function in [function for function in hit.functions if function in fragment_functions]:
                            if hit_function in function_maxbitscores:
                                if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                            else:
                                function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                function_maxbitscores[hit_function]['identity'] = hit.identity
                                function_maxbitscores[hit_function]['length'] = hit.s_len
                    # Count hits and sum identity
                    for function in fragment_functions:
                        ret_val[function][s]['count'] += 1
                        if function in function_maxbitscores:
                            ret_val[function][s]['hit_count'] += 1.0 
                            ret_val[function][s]['identity'] += function_maxbitscores[function]['identity']
                            if metrics == 'fragmentcount':
                                ret_val[function][s][metrics] += 1
                            elif metrics in ['fpk','fpkg','fpkm']:
                                ret_val[function][s][metrics] += norm_factor * get_fpk_score(function_maxbitscores[function]['length'])
                            else:
                                ret_val[function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, insert_size=insert_size)
                        else:
                            print('Function',function,'not found in hits of read', read_id)

            for read_id in project.samples[s].reads['pe2'].keys():
                if read_id in reads_processed: 
                    continue #Skip read if it was already counted
                    
                read_pe2 = project.samples[s].reads['pe2'][read_id]
                if read_pe2.status != 'function':
                    continue

                fragment_functions = set()
                read_functions = read_pe2.functions
                fragment_functions.update(read_functions.keys())

                # Count FPKM
                #~ for function in fragment_functions:
                    #~ ret_val[function][s]['count'] += 1
                    #~ ret_val[function][s][metrics] += norm_factor * read_functions[function]

                function_maxbitscores = defaultdict(dict)
                hits = [hit for hit in read_pe2.hit_list.hits]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [function for function in hit.functions if function in fragment_functions]:
                        if hit_function in function_maxbitscores:
                            if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                function_maxbitscores[hit_function]['identity'] = hit.identity
                                function_maxbitscores[hit_function]['length'] = hit.s_len
                        else:
                            function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                            function_maxbitscores[hit_function]['identity'] = hit.identity
                            function_maxbitscores[hit_function]['length'] = hit.s_len

                # Count hits and sum identity
                for function in fragment_functions:
                    ret_val[function][s]['count'] += 1
                    if function in function_maxbitscores:
                        ret_val[function][s]['hit_count'] += 1.0 
                        ret_val[function][s]['identity'] += function_maxbitscores[function]['identity']
                        if metrics == 'fragmentcount':
                            ret_val[function][s][metrics] += 1
                        elif metrics in ['fpk','fpkg','fpkm']:
                            ret_val[function][s][metrics] += norm_factor * get_fpk_score(function_maxbitscores[function]['length'])
                        else:
                            ret_val[function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, insert_size=insert_size)
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
        if project.samples[s].reads is None or 'pe1' not in project.samples[s].reads:
            raise ValueError ('No reads data loaded for sample',s,'end pe1')
        if project.samples[s].is_paired_end:
            if 'pe2' not in project.samples[s].reads:
                raise ValueError ('No reads data loaded for sample',s,'end pe2')

        norm_factor = 0.0
        if metrics in ['readcount','erpk','fragmentcount','fpk','efpk']:
            norm_factor = 1.0
        elif metrics in ['rpkm','fpkm','erpkm','efpkm']:
            norm_factor = project.samples[s].rpkm_scaling_factor
        elif metrics in ['rpkg','fpkg','erpkg','efpkg']:
            norm_factor = project.samples[s].rpkg_scaling_factor
        
        if norm_factor == 0.0:
            raise ValueError('Cannot get normalization factor')


        # Calculate scores
        if metrics in ['readcount','rpkg','rpkm','erpk','erpkg','erpkm']:
            if project.samples[s].is_paired_end:
                raise ValueError('Read count, RPKG and RPKM metrics provided only for single-end sequences')
            
            for read_id, read in project.samples[s].reads['pe1'].items():
                if read.status!= 'function': # Filter unmapped reads
                    continue
                read_erpk_scores = read.functions
                for function in read_erpk_scores:
                    ret_val[read.taxonomy][function][s]['count'] += 1.0
                    if metrics == 'readcount':
                        ret_val[read.taxonomy][function][s][metrics] += 1
                    else:
                        ret_val[read.taxonomy][function][s][metrics] += norm_factor * read_erpk_scores[function]
                        
                function_maxbitscores = {}
                hits = [hit for hit in read.hit_list.hits]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [function for function in hit.functions if function in read.functions]:
                        if hit_function in function_maxbitscores:
                            if hit.bitscore > function_maxbitscores[hit_function]:
                                function_maxbitscores[hit_function] = hit.bitscore
                        else:
                            function_maxbitscores[hit_function] = hit.bitscore
                # Count hits and sum identity
                for hit in hits:
                    for hit_function in [function for function in hit.functions if function in read.functions]:
                        if hit_function in function_maxbitscores and hit.bitscore == function_maxbitscores[hit_function]:
                            ret_val[read.taxonomy][function][s]['hit_count'] += 1.0 
                            ret_val[read.taxonomy][function][s]['identity'] += hit.identity
                            del function_maxbitscores[hit_function]
        
        elif metrics in ['fragmentcount','fpk','efpk','fpkg','fpkm','efpkg','efpkm']:
            if not project.samples[s].is_paired_end:
                raise ValueError('FPKG and FPKM metrics require paired-end sequences')
            insert_size = project.get_insert_size(project.samples[s])
            reads_processed = set()
            length_cutoff = project.config.get_length_cutoff(project.options.get_collection(s))
            average_read_length = project.samples[s].get_avg_read_length('pe1')
            
            for read_id in project.samples[s].reads['pe1'].keys():
                read_pe1 = project.samples[s].reads['pe1'][read_id]
                if read_pe1.status != 'function':
                    continue
                reads_processed.add(read_id)
                                        
                if 'pe2' in project.samples[s].reads and read_id in project.samples[s].reads['pe2'].keys(): 
                    read_pe2 = project.samples[s].reads['pe2'][read_id]
                    if read_pe2.status == 'function':
                        # Both ends are mapped
                        fragment_taxonomy = project.taxonomy_data.get_lca([read_pe1.taxonomy, read_pe2.taxonomy])
                        
                        fragment_functions = set() # List of functions assigned to the current read

                        read1_functions = read_pe1.functions
                        read2_functions = read_pe2.functions

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
                        hits1 = [hit for hit in read_pe1.hit_list.hits]
                        hits2 = [hit for hit in read_pe2.hit_list.hits]
                        # Find max. bitscore for each function
                        for hit in hits1:
                            for hit_function in [function for function in hit.functions if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                        function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                        function_maxbitscores[hit_function]['identity'] = hit.identity
                                        function_maxbitscores[hit_function]['length'] = hit.s_len
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                        for hit in hits2:
                            for hit_function in [function for function in hit.functions if function in fragment_functions]:
                                if hit_function in function_maxbitscores:
                                    if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                        function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                        function_maxbitscores[hit_function]['identity'] = hit.identity
                                        function_maxbitscores[hit_function]['length'] = hit.s_len
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                        # Count hits and sum identity for best hits for each function
                        for function in fragment_functions:
                            if function in function_maxbitscores:
                                ret_val[fragment_taxonomy][function][s]['count'] += 1
                                ret_val[fragment_taxonomy][function][s]['hit_count'] += 1.0 
                                ret_val[fragment_taxonomy][function][s]['identity'] += function_maxbitscores[function]['identity']
                                if metrics == 'fragmentcount':
                                    ret_val[fragment_taxonomy][function][s][metrics] += 1
                                elif metrics in ['fpk','fpkg','fpkm']:
                                    ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * get_fpk_score(function_maxbitscores[function]['length'])
                                else:
                                    ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, insert_size=insert_size)
                            else:
                                print('Function',function,'not found in hits of read', read_id)


                else: #Only end1 is mapped
                    fragment_taxonomy = read_pe1.taxonomy
                    read_functions = read_pe1.functions
                    fragment_functions = set(read_functions.keys())
                    # Count FPKM
                    #~ for function in fragment_functions:
                        #~ ret_val[fragment_taxonomy][function][s]['count'] += 1
                        #~ ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * read_functions[function]
                    function_maxbitscores = defaultdict(dict)
                    hits1 = [hit for hit in read_pe1.hit_list.hits]
                    # Find max. bitscore for each function
                    for hit in hits1:
                        for hit_function in [function for function in hit.functions if function in fragment_functions]:
                            if hit_function in function_maxbitscores:
                                if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                            else:
                                function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                function_maxbitscores[hit_function]['identity'] = hit.identity
                                function_maxbitscores[hit_function]['length'] = hit.s_len
                    # Count hits and sum identity for best hits for each function
                    for function in fragment_functions:
                        if function in function_maxbitscores:
                            ret_val[fragment_taxonomy][function][s]['count'] += 1
                            ret_val[fragment_taxonomy][function][s]['hit_count'] += 1.0 
                            ret_val[fragment_taxonomy][function][s]['identity'] += function_maxbitscores[function]['identity']
                            if metrics == 'fragmentcount':
                                ret_val[fragment_taxonomy][function][s][metrics] += 1
                            elif metrics in ['fpk','fpkg','fpkm']:
                                ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * get_fpk_score(function_maxbitscores[function]['length'])
                            else:
                                ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, insert_size=insert_size)
                        else:
                            print('Function',function,'not found in hits of read', read_id)

            for read_id in project.samples[s].reads['pe2'].keys():
                if read_id in reads_processed: 
                    continue #Skip read if it was already counted
                    
                read_pe2 = project.samples[s].reads['pe2'][read_id]
                if read_pe2.status != 'function':
                    continue

                fragment_functions = set()
                fragment_taxonomy = read_pe2.taxonomy
                read_functions = read_pe2.functions
                fragment_functions.update(read_functions.keys())

                # Count FPKM
                #~ for function in fragment_functions:
                    #~ ret_val[fragment_taxonomy][function][s]['count'] += 1
                    #~ ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * read_functions[function]
                # Build FPKM-based functional taxonomic profile
                function_maxbitscores = defaultdict(dict)
                hits = [hit for hit in read_pe2.hit_list.hits]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [function for function in hit.functions if function in fragment_functions]:
                        if hit_function in function_maxbitscores:
                            if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                function_maxbitscores[hit_function]['identity'] = hit.identity
                                function_maxbitscores[hit_function]['length'] = hit.s_len
                        else:
                            function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                            function_maxbitscores[hit_function]['identity'] = hit.identity
                            function_maxbitscores[hit_function]['length'] = hit.s_len
                # Count hits and sum identity for best hits for each function
                for function in fragment_functions:
                    if function in function_maxbitscores:
                        ret_val[fragment_taxonomy][function][s]['count'] += 1
                        ret_val[fragment_taxonomy][function][s]['hit_count'] += 1.0 
                        ret_val[fragment_taxonomy][function][s]['identity'] += function_maxbitscores[function]['identity']
                        if metrics == 'fragmentcount':
                            ret_val[fragment_taxonomy][function][s][metrics] += 1
                        elif metrics in ['fpk','fpkg','fpkm']:
                            ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * get_fpk_score(function_maxbitscores[function]['length'])
                        else:
                            ret_val[fragment_taxonomy][function][s][metrics] += norm_factor * get_efpk_score(function_maxbitscores[function]['length'], average_read_length, length_cutoff, insert_size=insert_size)
                    else:
                        print('Function',function,'not found in hits of read', read_id)
    return ret_val

def generate_sample_report(project, sample_id, metrics = None):
    # This function creates output files only in project's working directory
    if metrics is None:
        if project.samples[sample_id].is_paired_end:
            metrics = 'efpkg'
        else:
            metrics = 'erpkg'
    #outfile = os.path.join(project.samples[sample_id].work_directory, sample_id + '_report.txt')
    outfile = sanitize_file_name(os.path.join(project.options.work_dir, sample_id + '_report.txt'))
    with open(outfile, 'w') as of:
        of.write('Report for ' + sample_id + '\n\n')
        of.write('Sample ID:\t' + project.options.get_sample_name(sample_id) + '\n')
        of.write('Replicate:\t' + str(project.samples[sample_id].replicate) + '\n\n')
        of.write('Project name:\t' + project.options.project_name + '\n\n')
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
        if project.samples[sample_id].is_paired_end:
            of.write('Average insert size:\t' + str(project.get_insert_size(project.samples[sample_id])) + '\n')

        of.write('\nNumber of mapped reads\n')
        of.write('FASTQ file 1:\t' + str(len(project.samples[sample_id].reads['pe1'])) + '\n')
        if project.samples[sample_id].is_paired_end:
            of.write('FASTQ file 2:\t' + str(len(project.samples[sample_id].reads['pe2'])) + '\n')

        of.closed

    if project.samples[sample_id].rpkg_scaling_factor is None and metrics in ['efpkg', 'fpkg', 'erpkg', 'rpkg']:
        raise ValueError('Not enough data to normalize by average genome size')
    if project.samples[sample_id].rpkm_scaling_factor is None and metrics in ['efpkm', 'fpkm', 'erpkm', 'rpkm']:
        raise ValueError('Not enough data to normalize by sample size')

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
    outfile = sanitize_file_name(os.path.join(project.options.work_dir, sample_id + '_' + metrics + '_functional_taxonomy_profile.xml'))
    tax_profile.build_functional_taxonomy_profile(project.taxonomy_data, sample_scores_taxonomy)
    function_list = sorted(project.ref_data.functions_dict.keys())
    generate_taxonomy_series_chart(tax_profile, function_list, outfile, project.config.krona_path, score=metrics)

def generate_protein_sample_report(project, sample_id, metrics = None):
    # This function creates output files only in project's working directory
    #outfile = os.path.join(project.samples[sample_id].work_directory, sample_id + '_report.txt')
    outfile = sanitize_file_name(os.path.join(project.options.work_dir, sample_id + '_report.txt'))
    with open(outfile, 'w') as of:
        of.write('Report for ' + sample_id + '\n\n')
        of.write('Sample ID:\t' + project.options.get_sample_name(sample_id) + '\n')
        of.write('Replicate:\t' + str(project.samples[sample_id].replicate) + '\n\n')
        of.write('Project name:\t' + project.options.project_name + '\n\n')
        of.write('Reference data set:\t' + project.options.get_collection(sample_id) + '\n\n')
        
        of.write('Source files\n')
        of.write('Input file:\t' + project.options.get_fastq_path(sample_id, 'pe1') + '\n')
        of.write('    Number of proteins:\t' + str(project.samples[sample_id].fastq_fwd_readcount) + '\n')
        of.write('    Number of amino acids:\t' + str(project.samples[sample_id].fastq_fwd_basecount) + '\n')
        
        of.write('\nNumber of mapped proteins:' + str(len(project.samples[sample_id].reads['pe1'])) + '\n')

        for read in sorted(project.samples[sample_id].reads['pe1'].keys()):
            protein = project.samples[sample_id].reads['pe1'][read]
            if protein.status == 'nofunction':
                continue
            of.write(read + ': ' + ','.join(sorted(protein.functions.keys())))
            tax_id = protein.taxonomy
            if tax_id is None:
                of.write(' No taxonomy\n')
            else:
                of.write(' Taxonomy :' + project.taxonomy_data.names[tax_id]['name'] + '(' + tax_id + ')\n')
            of.write(' Hits:\n')
            for hit in protein.hit_list.hits:
                of.write('\t' + str(hit) + '\n')
        
        of.closed

    if len(project.samples[sample_id].reads['pe1']) == 0:
        return

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
    outfile = sanitize_file_name(os.path.join(project.options.work_dir, sample_id + '_' + metrics + '_functional_taxonomy_profile.xml'))
    tax_profile.build_functional_taxonomy_profile(project.taxonomy_data, sample_scores_taxonomy)
    function_list = sorted(project.ref_data.functions_dict.keys())
    generate_taxonomy_series_chart(tax_profile, function_list, outfile, project.config.krona_path, score=metrics)

    
def generate_project_report(project, metrics = None):
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
    
    if metrics is None: 
        if is_paired_end:
            metrics = 'efpkg'
        else:
            metrics = 'erpkg'
    
    scores = get_function_scores(project, sample_id=None, metrics=metrics)
    generate_function_sample_xlsx(project, scores, metrics=metrics)
    scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id=None, metrics=metrics)
    generate_function_taxonomy_sample_xlsx(project, scores_function_taxonomy, metrics=metrics, sample_id=None)
    generate_sample_taxonomy_function_xlsx(project, scores_function_taxonomy, metrics=metrics, function_id=None)

    generate_project_markdown_document(project, scores, metrics=metrics)


def generate_protein_project_report(project):
    for sample_id in project.list_samples():
        # If there are no read data, try to load them from JSON
        if len(project.samples[sample_id].reads) == 0:
            project.import_reads_json(sample_id, project.ENDS)
    print ('Generating spreadsheets and interactive diagrams for the project...')
    outfile = os.path.join(project.options.work_dir, 'project_report.txt')
    with open (outfile, 'w') as of:
        of.write(project.options.project_name + '\n\n')
        for sample_id in project.list_samples():
            of.write(sample_id + ':\t' + project.samples[sample_id].sample_name + '\tproteins mapped: ' + str(len(project.samples[sample_id].reads['pe1'])))
            of.write('\n')
        of.write('\nList of mapped proteins\n')

        for sample_id in project.list_samples():
            for protein_id in sorted(project.samples[sample_id].reads['pe1'].keys()):
                protein = project.samples[sample_id].reads['pe1'][protein_id]
                if protein.status == 'nofunction':
                    continue
                for hit in protein.hit_list.hits:
                    of.write('\t'.join([sample_id, protein_id, project.taxonomy_data.names[protein.taxonomy]['name'], str(hit)]) + '\n')
        of.closed

    metrics = 'readcount'
    scores = get_function_scores(project, sample_id=None, metrics=metrics)
    generate_function_sample_xlsx(project, scores, metrics=metrics)
    scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id=None, metrics=metrics)
    generate_function_taxonomy_sample_xlsx(project, scores_function_taxonomy, metrics=metrics, sample_id=None)
    generate_sample_taxonomy_function_xlsx(project, scores_function_taxonomy, metrics=metrics, function_id=None)


def generate_project_markdown_document(project, scores, sample_id = None, metrics = None):

    functions_list = sorted(scores.keys())
    
    # write output
    if sample_id is None:
        samples_list = sorted(project.list_samples())
    else:
        samples_list = [sample_id]

    outfile = sanitize_file_name(os.path.join(project.options.work_dir, 'index.md'))

    with open(outfile, 'w') as of:
        of.write('# ')
        of.write(project.options.project_name)
        of.write('\n\n')

        #~ of.write('## Categories\n\nAbundance of functional categories (RPKM score)\n\n')

        #~ of.write('Category|')
        #~ of.write('|'.join(samples_list))
        #~ of.write('|\n')
        #~ of.write('|---'*(len(samples_list) + 1))
        #~ of.write('|\n')
        
        #~ for category in sorted(categories_list):
            #~ of.write(category)
            #~ for sample in samples_list:
                #~ of.write('|')
                #~ if sample in scores_cat[category]:
                    #~ of.write('{0:.3f}'.format(scores_cat[category][sample]))
                #~ else:
                    #~ of.write('0.000')
            #~ of.write('|\n')

        of.write('\n## Functions\n\n' + metrics + ' scores of individual functions\n\n')


        of.write('Function|')
        of.write('|'.join(samples_list))
        of.write('|Definition|\n')
        of.write('|---'*(len(samples_list) + 2))
        of.write('|\n')
        
        for function in sorted(functions_list):
            of.write(function)
            for sample in samples_list:
                of.write('|')
                if function in scores and sample in scores[function]:
                    of.write('{0:.3f}'.format(scores[function][sample][metrics]))
                else:
                    of.write('0.000')
            of.write('|')
            of.write(project.ref_data.lookup_function_name(function))
            of.write('|\n')

        targetfile = sanitize_file_name(os.path.join('data', project.options.project_name + '_functions.xlsx'))

        of.write('\n<a href="')
        of.write(targetfile)
        of.write('" target="_blank">Download table of RPKM scores and read counts in XLSX format</a>\n\n')

        of.write('## Taxonomy profile for all reads mapped to nitrogen cycle genes\n\n')

        targetfile = sanitize_file_name(os.path.join('data', project.options.project_name + '_taxonomy_profile.xml.html'))
        of.write('<div class="krona-wrapper">\n<iframe src="')
        of.write(targetfile)
        of.write('" height="800" width="100%">')
        of.write(project.options.project_name)
        of.write(' taxonomy profile</iframe>\n<br>\n<a href="')
        of.write(targetfile)
        of.write('" target="_blank">Open chart in a new window</a>\n</div>\n\n')
        
        of.write('## Taxonomy profiles for individual functions\n\n')

        targetfile = sanitize_file_name(os.path.join('data', project.options.project_name + '_functions_taxonomy.xlsx'))
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
        outfile = sanitize_file_name(os.path.join(project.options.work_dir, sample + '.md'))
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
            of.write(project.options.project_name)
            of.write(' taxonomy profile</iframe>\n<br>\n<a href="')
            of.write(targetfile)
            of.write('" target="_blank">Open chart in a new window</a>\n</div>\n\n')
            
            of.write('## Reports for FASTQ files\n\n')
            if os.path.exists(os.path.join(project.options.get_project_dir(sample), project.options.get_output_subdir(sample),sample + '_pe1_'+ project.options.report_name + '.pdf')):
                targetfile = os.path.join('..','data', sample + '_pe1_'+ project.options.report_name + '.pdf')
                of.write('<a href="')
                of.write(targetfile)
                of.write('">Download report for read end 1 (PDF format)</a>\n<br>\n')
            if os.path.exists(os.path.join(project.options.get_project_dir(sample), project.options.get_output_subdir(sample),sample + '_pe2_'+ project.options.report_name + '.pdf')):
                targetfile = os.path.join('..','data', sample + '_pe2_'+ project.options.report_name + '.pdf')
                of.write('<a href="')
                of.write(targetfile)
                of.write('">Download report for read end 2 (PDF format)</a>\n<br>\n')
            of.closed
            
def generate_functions_stamp_input(project, scores, metrics):
    scores_outfile = sanitize_file_name(os.path.join(project.options.work_dir, project.options.project_name + '_' + metrics + '_functions.stamp.tsv'))
    metadata_outfile = sanitize_file_name(os.path.join(project.options.work_dir, project.options.project_name + '_functions.metadata.stamp.tsv'))
    
    with open(scores_outfile,'w') as of:
        of.write('Category\tFunction\t' + '\t'.join(project.list_samples()) + '\n')
        for function in sorted(scores.keys()):
            of.write(project.ref_data.lookup_function_group(function) + '\t' + function)
            for sample_id in project.list_samples():
                if sample_id in scores[function]:
                    #of.write('\t' + str(1000 * scores[function][sample_id][metrics]))
                    of.write('\t' + str(scores[function][sample_id][metrics]))
                else:
                    of.write('\t0')
            of.write('\n')
        of.closed
                
    with open(metadata_outfile,'w') as of:
        of.write('SampleID\tSample_name\tReplicate\n')
        for sample_id in project.list_samples():
            of.write(sample_id + '\t' + project.samples[sample_id].sample_name + '\t' + project.samples[sample_id].replicate + '\n')
        of.closed


def print_stamp_dataseries_taxonomy(tax_profile, sample_list, metrics, taxonomy_id='1', prefix='', taxonomy_level=0):
    attribute_values = autovivify(2,float)
    
    if taxonomy_id not in tax_profile.tree.data:
        raise ValueError (taxonomy_id,'not found in the tree!!!')
    
    if taxonomy_id == '1': # top level 
        new_prefix = prefix
    elif prefix == '': # superkingdom level: do not start with tab
        new_prefix = tax_profile.tree.data[taxonomy_id].name
    else:
        new_prefix = prefix + '\t' + tax_profile.tree.data[taxonomy_id].name
    new_taxonomy_level = taxonomy_level + 1
    ret_val = ''
    if tax_profile.tree.data[taxonomy_id].children:
        for child_taxid in tax_profile.tree.data[taxonomy_id].children:
            child_node, child_values = print_stamp_dataseries_taxonomy(tax_profile, 
                                            sample_list, 
                                            metrics, 
                                            taxonomy_id=child_taxid, 
                                            prefix=new_prefix, 
                                            taxonomy_level=new_taxonomy_level)
            ret_val += child_node
            for datapoint in child_values.keys():
                for k,v in child_values[datapoint].items():
                    attribute_values[datapoint][k] += v

        unclassified_flag = False
#        print(tax_profile.tree.data[taxid].attributes)
        for s in sample_list:
            if s in tax_profile.tree.data[taxonomy_id].attributes:
                if attribute_values[s][metrics] < tax_profile.tree.data[taxonomy_id].attributes[s][metrics]:
                    unclassified_flag = True
                    break

        if unclassified_flag:
            if taxonomy_id == '1': # line sohuld not start with tab symbol
                ret_val += 'unclassified' + '\tunclassified' * (len(tax_profile.RANKS) - new_taxonomy_level - 1)
            else:
                ret_val += new_prefix + '\tUnclassified ' + tax_profile.tree.data[taxonomy_id].name + '\tunclassified' * (len(tax_profile.RANKS) - new_taxonomy_level - 1)
            for s in sample_list:
                if s in tax_profile.tree.data[taxonomy_id].attributes and attribute_values[s][metrics] < tax_profile.tree.data[taxonomy_id].attributes[s][metrics]:
                    ret_val += '\t' + str(tax_profile.tree.data[taxonomy_id].attributes[s][metrics] - attribute_values[s][metrics])
                else:
                    ret_val += '\t0'
            ret_val += '\n'

    else:
    #~ if new_taxonomy_level == len(tax_profile.RANKS):
        #~ ret_val += new_prefix 
    #~ else:
        ret_val += new_prefix + '\tunclassified' * (len(tax_profile.RANKS) - new_taxonomy_level)
        if tax_profile.tree.data[taxonomy_id].attributes:
            for sample_id in sample_list:
                if sample_id in tax_profile.tree.data[taxonomy_id].attributes and metrics in tax_profile.tree.data[taxonomy_id].attributes[sample_id]:
                    ret_val += '\t' + str(tax_profile.tree.data[taxonomy_id].attributes[sample_id][metrics])
                else:
                    ret_val += '\t0'
        else:
            ret_val += '\t0' * len(sample_list)
        ret_val += '\n'
        
    
    attribute_values = autovivify(1)
    for sample_id in sample_list:
        if sample_id in tax_profile.tree.data[taxonomy_id].attributes and metrics in tax_profile.tree.data[taxonomy_id].attributes[sample_id]:
            attribute_values[sample_id][metrics] = tax_profile.tree.data[taxonomy_id].attributes[sample_id][metrics]

    return ret_val, attribute_values


def generate_functions_taxonomy_stamp_input(project, scores, metrics):
    metadata_outfile = sanitize_file_name(os.path.join(project.options.work_dir, project.options.project_name + '_functions_taxonomy.metadata.stamp.tsv'))
    
    sample_list = project.list_samples()
    root_id = '1'

    for function in project.ref_data.list_functions():

        # Subsetting scores
        sample_scores = autovivify(3, float)
        for tax in scores.keys():
            if function in scores[tax].keys():
                for s in project.list_samples():
                    if s in scores[tax][function] and metrics in scores[tax][function][s]:
                        #sample_scores[tax][s][metrics] = 1000 * scores[tax][function][s][metrics]
                        sample_scores[tax][s][metrics] = scores[tax][function][s][metrics]
                    else:
                        sample_scores[tax][s][metrics] = 0.0
        
        if len(sample_scores) == 0:
            continue

        scores_outfile = sanitize_file_name(os.path.join(project.options.work_dir, project.options.project_name + '_' + metrics + '_' + function + '_taxonomy.stamp.tsv'))

        with open(scores_outfile,'w') as of:
            of.write('\t'.join(project.taxonomy_data.RANKS[1:-1]) + '\t' + '\t'.join(project.list_samples()) + '\n')
            tax_profile = TaxonomyProfile()
            tax_profile.build_functional_taxonomy_profile(project.taxonomy_data, sample_scores)
            child_nodes, attribute_values = print_stamp_dataseries_taxonomy(tax_profile, sample_list, metrics, taxonomy_id=root_id, prefix='', taxonomy_level=0)
            of.write(child_nodes)
        of.closed
                
    with open(metadata_outfile,'w') as of:
        of.write('SampleID\tSample_name\tReplicate\n')
        for sample_id in project.list_samples():
            of.write(sample_id + '\t' + project.samples[sample_id].sample_name + '\t' + project.samples[sample_id].replicate + '\n')
        of.closed
    
def generate_sample_html_report(project, sample_id, metrics = None):
    # This function creates output files only in project's working directory
    if metrics is None:
        if project.samples[sample_id].is_paired_end:
            metrics = 'efpkg'
        else:
            metrics = 'erpkg'
    #outfile = os.path.join(project.samples[sample_id].work_directory, sample_id + '_report.txt')
    outfile = sanitize_file_name(os.path.join(project.options.work_dir, sample_id + '_report.txt'))
    with open(outfile, 'w') as of:
        of.write('Report for ' + sample_id + '\n\n')
        of.write('Sample ID:\t' + project.options.get_sample_name(sample_id) + '\n')
        of.write('Replicate:\t' + str(project.samples[sample_id].replicate) + '\n\n')
        of.write('Project name:\t' + project.options.project_name + '\n\n')
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
        if project.samples[sample_id].is_paired_end:
            of.write('Average insert size:\t' + str(project.get_insert_size(project.samples[sample_id])) + '\n')

        of.write('\nNumber of mapped reads\n')
        of.write('FASTQ file 1:\t' + str(len(project.samples[sample_id].reads['pe1'])) + '\n')
        if project.samples[sample_id].is_paired_end:
            of.write('FASTQ file 2:\t' + str(len(project.samples[sample_id].reads['pe2'])) + '\n')

        of.closed

    if project.samples[sample_id].rpkg_scaling_factor is None and metrics in ['efpkg', 'fpkg', 'erpkg', 'rpkg']:
        raise ValueError('Not enough data to normalize by average genome size')
    if project.samples[sample_id].rpkm_scaling_factor is None and metrics in ['efpkm', 'fpkm', 'erpkm', 'rpkm']:
        raise ValueError('Not enough data to normalize by sample size')

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
    outfile = sanitize_file_name(os.path.join(project.options.work_dir, sample_id + '_' + metrics + '_functional_taxonomy_profile.xml'))
    tax_profile.build_functional_taxonomy_profile(project.taxonomy_data, sample_scores_taxonomy)
    function_list = sorted(project.ref_data.functions_dict.keys())
    generate_taxonomy_series_chart(tax_profile, function_list, outfile, project.config.krona_path, score=metrics)


def generate_assembly_report(assembler):
    outfile = os.path.join(assembler.assembly_dir, 'out', 'assembly_report.txt')
    with open(outfile, 'w') as of:
        of.write('\nFunction statistics\n\n')
        for function in assembler.assembly.contigs:
            of.write(function + ':\n')
            for contig in assembler.assembly.contigs[function]:
                for gene_id in assembler.assembly.contigs[function][contig].genes:
                    gene = assembler.assembly.contigs[function][contig].genes[gene_id]
                    if gene.status == 'function':
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

