#!/usr/bin/python
import os
from collections import defaultdict,Counter,OrderedDict
import xlsxwriter

from Fama.DiamondParser.hit_utils import cleanup_protein_id,autovivify
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData

ENDS = ['pe1','pe2']
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

def generate_functions_scores_table(project):
    functions_list = set()
    scores = defaultdict(lambda : defaultdict(float))
    read_counts = defaultdict(lambda : defaultdict(float))
    
    # fill tables of scores and read counts
    for sample in project.list_samples():
        for end in ENDS:
            project.load_annotated_reads(sample, end) # Lazy load
            scaling_factor = 1.0
            if end == 'pe1':
                scaling_factor = project.options.get_fastq1_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            elif end == 'pe2':
                scaling_factor = project.options.get_fastq2_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            else:
                raise Exception('Unknown identifier of read end: ' + end)
            #reads = project.samples[sample][end]
            for read_id in project.samples[sample][end]:
                read = project.samples[sample][end][read_id]
                if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                    for function in read.functions:
                        functions_list.add(function)
                        scores[function][sample] += scaling_factor * read.functions[function]
                        read_counts[function][sample] += 1.0/len(read.functions)
            project.samples[sample][end] = None

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
    for sample in project.samples:
        for end in project.samples[sample]:
            reads = project.samples[sample][end]
            for read in reads:
                if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                    for function in reads[read].functions:
                        if sample in functions[function]:
                            functions[function][sample] += reads[read].functions[function]
                            protein_counts[function][sample] += 1
                        else:
                            functions[function][sample] = reads[read].functions[function]
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

def create_functions_xlsx(project):
    xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), project.options.get_name() + '_functions.xlsx'))    
    xlsxfile = xlsxfile.replace(' ', '_')
    xlsxfile = xlsxfile.replace("'", "")
    xlsxfile = xlsxfile.replace('"', '')
    workbook = xlsxwriter.Workbook(xlsxfile)
    bold = workbook.add_format({'bold': True})
    
    functions_list = set()
    categories_list = set()
    scores = defaultdict(lambda : defaultdict(float))
    read_counts = defaultdict(lambda : defaultdict(float))
    scores_cat = defaultdict(lambda : defaultdict(float))
    read_counts_cat = defaultdict(lambda : defaultdict(float))
    
    # calculate scores

    for sample in project.list_samples():
        for end in ENDS:
            if not project.samples[sample][end]:
                project.load_annotated_reads(sample, end) 
            scaling_factor = 1.0
            if end == 'pe1':
                scaling_factor = project.options.get_fastq1_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            elif end == 'pe2':
                scaling_factor = project.options.get_fastq2_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            else:
                raise Exception('Unknown identifier of read end: ' + end)
#            reads = project.samples[sample][end]
            for read_id in project.samples[sample][end]:
                read = project.samples[sample][end][read_id]
                if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                    for function in read.functions:
                        functions_list.add(function)
                        scores[function][sample] += scaling_factor * read.functions[function]
                        read_counts[function][sample] += 1.0/len(read.functions)
                        category = project.ref_data.lookup_function_group(function)
                        categories_list.add(category)
                        scores_cat[category][sample] += scaling_factor * read.functions[function]
                        read_counts_cat[category][sample] += 1.0/len(read.functions)
            project.samples[sample][end] = None
            
    # generate output
    samples_list = sorted(project.list_samples())
    
    # generate tables for functions
    scores_worksheet = workbook.add_worksheet('Functions RPKM')
    counts_worksheet = workbook.add_worksheet('Functions read count')
    
    row = 0
    col = 0
    
    scores_worksheet.write(row, col, 'Function', bold)
    counts_worksheet.write(row, col, 'Function', bold)
    
    for sample in samples_list:
        col += 1
        scores_worksheet.write(row, col, sample, bold)
        counts_worksheet.write(row, col, sample, bold)
    
    col += 1
    scores_worksheet.write(row, col, 'Definition', bold)
    counts_worksheet.write(row, col, 'Definition', bold)
    
    for function in sorted(functions_list):
        row += 1
        col = 0
        scores_worksheet.write(row, col, function, bold)
        counts_worksheet.write(row, col, function, bold)
        for sample in samples_list:
            col += 1
            if sample in scores[function]:
                #scores_worksheet.write(row, col, '{0:.3f}'.format(scores[function][sample]))
                #counts_worksheet.write(row, col, '{0:.0f}'.format(read_counts[function][sample]))
                scores_worksheet.write(row, col, scores[function][sample])
                counts_worksheet.write(row, col, read_counts[function][sample])
            else:
                scores_worksheet.write(row, col, 0.0)
                counts_worksheet.write(row, col, 0)
        col += 1
        scores_worksheet.write(row, col, project.ref_data.lookup_function_name(function))
        counts_worksheet.write(row, col, project.ref_data.lookup_function_name(function))

    # adjust column width
    scores_worksheet.set_column(0, 0, 10)
    counts_worksheet.set_column(0, 0, 10)    
    scores_worksheet.set_column(col, col, 50)
    counts_worksheet.set_column(col, col, 50)

    # generate tables for categories
    scores_cat_worksheet = workbook.add_worksheet('Categories RPKM')
    counts_cat_worksheet = workbook.add_worksheet('Categories read count')
    
    row = 0
    col = 0
    
    scores_cat_worksheet.write(row, col, 'Category', bold)
    counts_cat_worksheet.write(row, col, 'Category', bold)
    
    for sample in samples_list:
        col += 1
        scores_cat_worksheet.write(row, col, sample, bold)
        counts_cat_worksheet.write(row, col, sample, bold)
    
    col += 1
    
    for category in sorted(categories_list):
        row += 1
        col = 0
        scores_cat_worksheet.write(row, col, category, bold)
        counts_cat_worksheet.write(row, col, category, bold)
        for sample in samples_list:
            col += 1
            if sample in scores_cat[category]:
#                scores_cat_worksheet.write(row, col, '{0:.3f}'.format(scores_cat[category][sample]))
#                counts_cat_worksheet.write(row, col, '{0:.0f}'.format(read_counts_cat[category][sample]))
                scores_cat_worksheet.write(row, col, scores_cat[category][sample])
                counts_cat_worksheet.write(row, col, read_counts_cat[category][sample])
            else:
                scores_cat_worksheet.write(row, col, 0.0)
                counts_cat_worksheet.write(row, col, 0)

    # adjust column width
    scores_cat_worksheet.set_column(0, 0, 40)
    counts_cat_worksheet.set_column(0, 0, 40)    

    workbook.close()
    
def create_functions_rpkm_xlsx(project):
    xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), project.options.get_name() + '_functions.xlsx'))    
    xlsxfile = xlsxfile.replace(' ', '_')
    xlsxfile = xlsxfile.replace("'", "")
    xlsxfile = xlsxfile.replace('"', '')
    workbook = xlsxwriter.Workbook(xlsxfile)
    bold = workbook.add_format({'bold': True})
    
    functions_list = set()
    categories_list = set()
    scores = autovivify(3, float)#defaultdict(lambda : defaultdict(float))
    read_counts = autovivify(3, float)#defaultdict(lambda : defaultdict(float))
    scores_cat = autovivify(3, float)#defaultdict(lambda : defaultdict(float))
    read_counts_cat = autovivify(3, float)#defaultdict(lambda : defaultdict(float))
    
    # calculate scores

    for sample in project.list_samples():
        for end in ENDS:
            project.load_annotated_reads(sample, end) 
            scaling_factor = 1.0
            if end == 'pe1':
                scaling_factor = project.options.get_fastq1_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            elif end == 'pe2':
                scaling_factor = project.options.get_fastq2_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            else:
                raise Exception('Unknown identifier of read end: ' + end)
#            reads = project.samples[sample][end]
            for read_id in project.samples[sample][end]:
                read = project.samples[sample][end][read_id]
                if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                    for function in read.functions:
                        functions_list.add(function)
                        scores[function][sample][end] += scaling_factor * read.functions[function]
                        read_counts[function][sample][end] += 1.0/len(read.functions)
                        category = project.ref_data.lookup_function_group(function)
                        categories_list.add(category)
                        scores_cat[category][sample][end] += scaling_factor * read.functions[function]
                        read_counts_cat[category][sample][end] += 1.0/len(read.functions)
            project.samples[sample][end] = None
            
    # generate output
    samples_list = sorted(project.list_samples())
    
    # generate tables for functions
    scores_worksheet = workbook.add_worksheet('Functions RPKM')
    counts_worksheet = workbook.add_worksheet('Functions read count')
    
    row = 0
    col = 0
    
    scores_worksheet.write(row, col, 'Function', bold)
    counts_worksheet.write(row, col, 'Function', bold)
    
    for sample in samples_list:
        for end in ENDS:
            col += 1
            scores_worksheet.write(row, col, sample, bold)
            counts_worksheet.write(row, col, sample, bold)
    
    col += 1
    scores_worksheet.write(row, col, 'Definition', bold)
    counts_worksheet.write(row, col, 'Definition', bold)
    
    for function in sorted(functions_list):
        row += 1
        col = 0
        scores_worksheet.write(row, col, function, bold)
        counts_worksheet.write(row, col, function, bold)
        for sample in samples_list:
            for end in ENDS:
                col += 1
                if sample in scores[function]:
                    #scores_worksheet.write(row, col, '{0:.3f}'.format(scores[function][sample]))
                    #counts_worksheet.write(row, col, '{0:.0f}'.format(read_counts[function][sample]))
                    scores_worksheet.write(row, col, scores[function][sample][end])
                    counts_worksheet.write(row, col, read_counts[function][sample][end])
                else:
                    scores_worksheet.write(row, col, 0.0)
                    counts_worksheet.write(row, col, 0)
        col += 1
        scores_worksheet.write(row, col, project.ref_data.lookup_function_name(function))
        counts_worksheet.write(row, col, project.ref_data.lookup_function_name(function))

    # adjust column width
    scores_worksheet.set_column(0, 0, 10)
    counts_worksheet.set_column(0, 0, 10)    
    scores_worksheet.set_column(col, col, 50)
    counts_worksheet.set_column(col, col, 50)

    # generate tables for categories
    scores_cat_worksheet = workbook.add_worksheet('Categories RPKM')
    counts_cat_worksheet = workbook.add_worksheet('Categories read count')
    
    row = 0
    col = 0
    
    scores_cat_worksheet.write(row, col, 'Category', bold)
    counts_cat_worksheet.write(row, col, 'Category', bold)
    
    for sample in samples_list:
        for end in ENDS:
            col += 1
            scores_cat_worksheet.write(row, col, sample, bold)
            counts_cat_worksheet.write(row, col, sample, bold)
    
    col += 1
    
    for category in sorted(categories_list):
        row += 1
        col = 0
        scores_cat_worksheet.write(row, col, category, bold)
        counts_cat_worksheet.write(row, col, category, bold)
        for sample in samples_list:
            for end in ENDS:
                col += 1
                if sample in scores_cat[category]:
    #                scores_cat_worksheet.write(row, col, '{0:.3f}'.format(scores_cat[category][sample]))
    #                counts_cat_worksheet.write(row, col, '{0:.0f}'.format(read_counts_cat[category][sample]))
                    scores_cat_worksheet.write(row, col, scores_cat[category][sample][end])
                    counts_cat_worksheet.write(row, col, read_counts_cat[category][sample][end])
                else:
                    scores_cat_worksheet.write(row, col, 0.0)
                    counts_cat_worksheet.write(row, col, 0)

    # adjust column width
    scores_cat_worksheet.set_column(0, 0, 40)
    counts_cat_worksheet.set_column(0, 0, 40)    

    workbook.close()

def create_functions_markdown_document(project):

    functions_list = set()
    categories_list = set()
    scores = defaultdict(lambda : defaultdict(float))
    scores_cat = defaultdict(lambda : defaultdict(float))
    
    # calculate scores

    for sample in project.list_samples():
        for end in ENDS:
            if not project.samples[sample][end]:
                project.load_annotated_reads(sample, end) 
            scaling_factor = 1.0
            if end == 'pe1':
                scaling_factor = project.options.get_fastq1_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            elif end == 'pe2':
                scaling_factor = project.options.get_fastq2_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            else:
                raise Exception('Unknown identifier of read end: ' + end)
            reads = project.samples[sample][end]
            for read_id in project.samples[sample][end]:
                read = project.samples[sample][end][read_id]
                if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                    for function in reads[read].functions:
                        functions_list.add(function)
                        scores[function][sample] += scaling_factor * read.functions[function]
                        category = project.ref_data.lookup_function_group(function)
                        categories_list.add(category)
                        scores_cat[category][sample] += scaling_factor * read.functions[function]
            project.samples[sample][end] = None

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
            


def sanitize_file_name(filename):
    filename = filename.replace(' ', '_')
    filename = filename.replace("'", "")
    filename = filename.replace('"', '')
    return filename
    

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

def create_assembly_xlsx(assembler, taxonomy_data):
    xlsxfile = sanitize_file_name(os.path.join(assembler.project.options.get_work_dir(), 'assembly', assembler.project.options.get_name() + '_assembly.xlsx'))    
    xlsxfile = xlsxfile.replace(' ', '_')
    xlsxfile = xlsxfile.replace("'", "")
    xlsxfile = xlsxfile.replace('"', '')
    workbook = xlsxwriter.Workbook(xlsxfile)
    bold = workbook.add_format({'bold': True})
    cell_numformat1 = workbook.add_format()
    cell_numformat1.set_num_format('0.0')
    cell_numformat5 = workbook.add_format()
    cell_numformat5.set_num_format('0.00000')
    
    functions_list = set()
    samples_list = sorted(assembler.project.list_samples())
    function_read_counts = autovivify(2, float) # function_read_counts[function][sample]
    contig_scores = autovivify(3, float)#contig_scores[function][contig][sample]
    gene_rpkm = autovivify(3, float) # gene_rpkm[function][gene][sample], parameters are RPKM, coverage, identity 
    
    # count reads per function, per sample
    for function in assembler.assembly.reads:
        functions_list.add(function)
        for read in assembler.assembly.reads[function]:
            function_read_counts[function][assembler.assembly.reads[function][read]] += 1
    
    # collect RPKM scores for contigs per function, per sample (for contigs? for genes?)
    
    
    # calculate total read count
    total_read_count = 0 

    for sample in assembler.project.list_samples():
        total_read_count += assembler.project.options.get_fastq1_readcount(sample)
        #total_read_count += assembler.project.options.get_fastq2_readcount(sample)
            
            
    # generate output
    
    # make worksheet for read counts per function
    reads_worksheet = workbook.add_worksheet('Functions read count')
    
    row = 0
    col = 0
    reads_worksheet.write(row, col, 'Function', bold)

    for sample in samples_list:
            col += 1
            reads_worksheet.write(row, col, sample, bold)
    col += 1
    reads_worksheet.write(row, col, 'All samples', bold)
    col += 1
    reads_worksheet.write(row, col, 'Assembled reads', bold)
    col += 1
    reads_worksheet.write(row, col, 'Unassembled reads', bold)
    col += 1
    reads_worksheet.write(row, col, 'Definition', bold)

    for function in sorted(functions_list):
        row += 1
        col = 0
        reads_worksheet.write(row, col, function, bold)
        for sample in samples_list:
                col += 1
                if sample in function_read_counts[function]:
                    reads_worksheet.write(row, col, '{0:.0f}'.format(function_read_counts[function][sample]))
                else:
                    reads_worksheet.write(row, col, '0')
        col += 1
        all_reads = sum(function_read_counts[function].values())
        reads_worksheet.write(row, col, '{0:.0f}'.format(all_reads), bold)
        col += 1
        assembled_reads = 0
        if function in assembler.assembly.contigs:
            assembled_reads = sum([len(c.reads) for c in assembler.assembly.contigs[function].values()]) / 2
        reads_worksheet.write(row, col, '{0:.0f}'.format(assembled_reads), bold)
        col += 1
        reads_worksheet.write(row, col, '{0:.0f}'.format(all_reads - assembled_reads), bold)
        col += 1
        reads_worksheet.write(row, col, assembler.project.ref_data.lookup_function_name(function))

    # adjust column width
    reads_worksheet.set_column(0, 0, 10)
    reads_worksheet.set_column(col, col, 50)


    # make worksheet with contig data
    contigs_worksheet = workbook.add_worksheet('Contigs')
    
    row = 0
    col = 0
    contigs_worksheet.write(row, col, 'Contig', bold)
    col += 1
    contigs_worksheet.write(row, col, 'Function', bold)
    col += 1
    contigs_worksheet.write(row, col, 'Length', bold)
    col += 1
    contigs_worksheet.write(row, col, 'Read count', bold)
    col += 1
    contigs_worksheet.write(row, col, 'RPKM', bold)
    col += 1
    contigs_worksheet.write(row, col, 'Coverage', bold)
    col += 1
    contigs_worksheet.write(row, col, 'Number of genes', bold)

    for sample in samples_list:
        col += 1
        contigs_worksheet.write(row, col, sample, bold)
        col += 1
        contigs_worksheet.write(row, col, sample, bold)
        col += 1
        contigs_worksheet.write(row, col, sample, bold)

    col += 1
    contigs_worksheet.write(row, col, 'Definition', bold)

    row += 1
    col = 6
    for sample in samples_list:
        col += 1
        contigs_worksheet.write(row, col, 'Read count', bold)
        col += 1
        contigs_worksheet.write(row, col, 'RPKM', bold)
        col += 1
        contigs_worksheet.write(row, col, 'Coverage', bold)

    
    for function in sorted(functions_list):
        if function in assembler.assembly.contigs:
            for contig in sorted(assembler.assembly.contigs[function].keys()):
                row += 1
                col = 0
                contigs_worksheet.write(row, col, contig, bold)
                col += 1
                contigs_worksheet.write(row, col, function)
                col += 1
                contigs_worksheet.write(row, col, len(assembler.assembly.contigs[function][contig].sequence))
                col += 1
                contigs_worksheet.write(row, col, assembler.assembly.contigs[function][contig].get_read_count())
                col += 1
                contigs_worksheet.write(row, col, assembler.assembly.contigs[function][contig].get_rpkm(total_read_count), cell_numformat5)
                col += 1
                contigs_worksheet.write(row, col, assembler.assembly.contigs[function][contig].get_coverage(), cell_numformat1)
                col += 1
                contigs_worksheet.write(row, col, len(assembler.assembly.contigs[function][contig].genes))
                col += 1

                for sample in samples_list:
                    contigs_worksheet.write(row, col, assembler.assembly.contigs[function][contig].get_read_count(sample))
                    col += 1
                    contigs_worksheet.write(row, col, assembler.assembly.contigs[function][contig].get_rpkm(assembler.project.options.get_fastq1_readcount(sample),sample), cell_numformat5)
                    col += 1
                    contigs_worksheet.write(row, col, assembler.assembly.contigs[function][contig].get_coverage(sample), cell_numformat1)
                    col += 1
                contigs_worksheet.write(row, col, assembler.project.ref_data.lookup_function_name(function))

    # adjust column width
    contigs_worksheet.set_column(0, 1, 10)
    contigs_worksheet.set_column(col, col, 50)

    # make worksheet for genes
    genes_worksheet = workbook.add_worksheet('Genes')
    
    row = 0
    col = 0
    genes_worksheet.write(row, col, 'Gene', bold)
    col += 1
    genes_worksheet.write(row, col, 'Reads function', bold)
    col += 1
    genes_worksheet.write(row, col, 'Contig', bold)
    col += 1
    genes_worksheet.write(row, col, 'Gene length', bold)
    col += 1
    genes_worksheet.write(row, col, 'Read count', bold)
    col += 1
    genes_worksheet.write(row, col, 'RPKM', bold)
    col += 1
    genes_worksheet.write(row, col, 'Coverage', bold)
    col += 1
    genes_worksheet.write(row, col, 'Gene status', bold)
    col += 1
    genes_worksheet.write(row, col, 'Uniref best hit taxonomy ID', bold)
    col += 1
    genes_worksheet.write(row, col, 'Uniref best hit organism', bold)
    col += 1
    genes_worksheet.write(row, col, 'Function best hit taxonomy IDs', bold)
    col += 1
    genes_worksheet.write(row, col, 'Protein functions', bold)
    col += 1
    genes_worksheet.write(row, col, 'Identity', bold)
    col += 1
    genes_worksheet.write(row, col, 'CDS completeness', bold)

    for sample in samples_list:
        col += 1
        genes_worksheet.write(row, col, sample, bold)
        col += 1
        genes_worksheet.write(row, col, sample, bold)
        col += 1
        genes_worksheet.write(row, col, sample, bold)

    col += 1
    genes_worksheet.write(row, col, 'Definition', bold)

    row += 1
    col = 13
    for sample in samples_list:
        col += 1
        genes_worksheet.write(row, col, 'Read count', bold)
        col += 1
        genes_worksheet.write(row, col, 'RPKM', bold)
        col += 1
        genes_worksheet.write(row, col, 'Coverage', bold)

    for function in sorted(functions_list):
        if function in assembler.assembly.contigs:
            for contig in sorted(assembler.assembly.contigs[function].keys()):
                for gene_id in sorted(assembler.assembly.contigs[function][contig].genes.keys()):
                    gene = assembler.assembly.contigs[function][contig].genes[gene_id]
                    row += 1
                    col = 0
                    genes_worksheet.write(row, col, gene_id)
                    col += 1
                    genes_worksheet.write(row, col, function)
                    col += 1
                    genes_worksheet.write(row, col, contig)
                    col += 1
                    genes_worksheet.write(row, col, len(gene.protein_sequence) * 3)
                    col += 1
                    gene_read_count = assembler.assembly.contigs[function][contig].get_read_count() * len(gene.protein_sequence) * 3 / len(assembler.assembly.contigs[function][contig].sequence)
                    genes_worksheet.write(row, col, gene_read_count, cell_numformat1)
                    col += 1
                    gene_rpkm = assembler.assembly.contigs[function][contig].get_rpkm(total_read_count) * len(gene.protein_sequence) * 3 / len(assembler.assembly.contigs[function][contig].sequence)
                    genes_worksheet.write(row, col, gene_rpkm, cell_numformat5)
                    col += 1
                    genes_worksheet.write(row, col, assembler.assembly.contigs[function][contig].get_coverage(), cell_numformat1)
                    col += 1
                    genes_worksheet.write(row, col, gene.get_status())
                    taxonomy_id = gene.get_taxonomy_id()
                    col += 1
                    if taxonomy_id:
                        genes_worksheet.write(row, col, taxonomy_id)
                    else:
                        genes_worksheet.write(row, col, 'Not found')
                    col += 1
                    if taxonomy_id in taxonomy_data.names:
                        genes_worksheet.write(row, col, taxonomy_data.names[taxonomy_id]['name'])
                    else:
                        genes_worksheet.write(row, col, 'Not found')
                    col += 1

                    if gene.get_status() == 'function,besthit' or gene.get_status() == 'function':
                        gene_taxonomy = [assembler.project.ref_data.lookup_protein_tax(cleanup_protein_id(x.get_subject_id())) for x in gene.hit_list.get_hits()]
                        genes_worksheet.write(row, col, ','.join(gene_taxonomy))
                        col += 1
                        gene_functions = [y for x in gene.hit_list.get_hits() for y in x.functions]
                        genes_worksheet.write(row, col, ','.join(gene_functions))
                        col += 1
                        gene_identity = [x.get_identity() for x in gene.hit_list.get_hits()]
                        genes_worksheet.write(row, col, sum (gene_identity) / len (gene_identity), cell_numformat1)
                        col += 1
                        ref_length = [x.get_subject_length() for x in gene.hit_list.get_hits()]
                        genes_worksheet.write(row, col, len (gene.protein_sequence) * 100 / sum (ref_length), cell_numformat1)
                    else:
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')

                    for sample in samples_list:
                        col += 1
                        gene_read_count = assembler.assembly.contigs[function][contig].get_read_count(sample) * len(gene.protein_sequence) * 3 / len(assembler.assembly.contigs[function][contig].sequence)
                        genes_worksheet.write(row, col, gene_read_count, cell_numformat1)
                        col += 1
                        gene_rpkm = assembler.assembly.contigs[function][contig].get_rpkm(assembler.project.options.get_fastq1_readcount(sample),sample) * len(gene.protein_sequence) * 3 / len(assembler.assembly.contigs[function][contig].sequence)
                        genes_worksheet.write(row, col, gene_rpkm, cell_numformat5)
                        col += 1
                        genes_worksheet.write(row, col, assembler.assembly.contigs[function][contig].get_coverage(sample), cell_numformat1)
                    col += 1
                    genes_worksheet.write(row, col, assembler.project.ref_data.lookup_function_name(function))
                        


    # adjust column width
    genes_worksheet.set_column(0, 0, 20)
    genes_worksheet.set_column(1, 1, 10)
    genes_worksheet.set_column(7, 9, 15)
    genes_worksheet.set_column(col, col, 50)
    
    workbook.close()
