import os
from collections import defaultdict,Counter,OrderedDict

from Fama.utils import autovivify,cleanup_protein_id,sanitize_file_name

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
            if parser.reads[read].get_status() == 'function,besthit' or parser.reads[read].get_status() == 'function':
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
            if parser.reads[read].get_status() != 'nofunction':
#                of.write(read + '\t' + parser.reads[read].get_status() + '\n')
#                of.write('\t' + str(parser.reads[read].get_functions()) + '\n')
                of.write(read + ': ' + parser.reads[read].get_status() + ': ' + ','.join(sorted(parser.reads[read].get_functions().keys())) + '\n')#' Taxonomy :' + parser.reads[read].taxonomy + '\n')
                for hit in parser.reads[read].get_hit_list().get_hits():
                    of.write('\t' + str(hit) + '\n')
                
        of.write('\n\n*** End of report ***\n')
        of.closed

def generate_fasta_report(project, sample_id):
    # TODO
    pass

def generate_sample_report(project, sample_id):
    outfile = os.path.join(project.samples[sample_id].work_directory, sample_id + '_report.txt')
    with open(outfile, 'w') as of:
        of.write('Report for ' + sample_id + '\n\n')
        of.write('Sample ID:\t' + project.options.get_sample_name(sample_id) + '\n')
        of.write('FASTQ file 1:\t' + project.options.get_fastq_path(sample_id, 'pe1') + '\n')
        if project.samples[sample_id].is_paired_end:
            of.write('FASTQ file 2:\t' + project.options.get_fastq_path(sample_id, 'pe2') + '\n')

        of.closed
    

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

