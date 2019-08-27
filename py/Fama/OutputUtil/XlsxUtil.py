import os
from collections import defaultdict,Counter,OrderedDict
import xlsxwriter
import pandas as pd

from Fama.utils import autovivify,cleanup_protein_id,sanitize_file_name
from Fama.TaxonomyProfile import TaxonomyProfile


def generate_function_sample_xlsx(project, scores, metrics, sample_id = None):
    # This function generates XLSX file for function scores for one or more samples.
    # Scores may be read counts, RPKM, RPKG, ERPKM, ERPKG (for single-end sequencing) or FPKM, FPKG, EFRKM, EFPKG for paired-end sequencing
    
    # scores must be dict of dicts, like scores[function_id][sample_id][metrics] = value
    # metrics must be 'readcount', 'rpkm', 'rpkg', fpkm', 'fpkg'
    if sample_id is None:
        xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), project.options.get_name() + '_' + metrics + '_functions.xlsx'))
    else:
        xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), sample_id + '_' + metrics + '_functions.xlsx'))
    
    print('Writing', xlsxfile)
    workbook = xlsxwriter.Workbook(xlsxfile)
    bold = workbook.add_format({'bold': True})
    
    functions_list = sorted(project.ref_data.functions_dict.keys())
    categories_list = sorted(list(set([project.ref_data.functions_dict[x]['group'] for x in project.ref_data.functions_dict.keys()])))
    scores_cat = autovivify(2,float)

    # generate tables for functions
    scores_worksheet = workbook.add_worksheet('Functions ' + metrics)
    
    row = 0
    col = 0
    
    scores_worksheet.write(row, col, 'Function', bold)
    
    for s in project.list_samples():
        if not sample_id is None and s != sample_id:
            continue
        col += 1
        scores_worksheet.write(row, col, s, bold)
    
    col += 1
    scores_worksheet.write(row, col, 'Definition', bold)
    
    for function in functions_list:
        category = project.ref_data.lookup_function_group(function)
        row += 1
        col = 0
        scores_worksheet.write(row, col, function, bold)
        for s in project.list_samples():
            if not sample_id is None and s != sample_id:
                continue
            col += 1
            if function in scores and s in scores[function]:
                scores_worksheet.write(row, col, scores[function][s][metrics])
                scores_cat[category][s] += scores[function][s][metrics]
            else:
                scores_worksheet.write(row, col, 0.0)

        col += 1
        scores_worksheet.write(row, col, project.ref_data.lookup_function_name(function))

    # adjust column width
    scores_worksheet.set_column(0, 0, 10)
    scores_worksheet.set_column(col, col, 50)

    # Write worksheet for categories
    scores_cat_worksheet = workbook.add_worksheet('Categories ' + metrics)
    row = 0
    col = 0
    scores_cat_worksheet.write(row, col, 'Categories', bold)
    
    for s in project.list_samples():
        if not sample_id is None and s != sample_id:
            continue
        col += 1
        scores_cat_worksheet.write(row, col, s, bold)
    
    for category in categories_list:
        row += 1
        col = 0
        scores_cat_worksheet.write(row, col, category, bold)
        for s in project.list_samples():
            if not sample_id is None and s != sample_id:
                continue
            col += 1
            if category in scores_cat and s in scores_cat[category]:
                scores_cat_worksheet.write(row, col, scores_cat[category][s])
            else:
                scores_cat_worksheet.write(row, col, 0.0)
    # adjust column width
    scores_cat_worksheet.set_column(0, 0, 50)

    workbook.close()

def generate_function_taxonomy_sample_xlsx(project, scores, metrics, sample_id = None, rank = None):
    # This function generates XLSX file for function and taxa scores for one or more samples.
    # Scores may be read counts, RPKM, RPKG (for single-end sequencing) or FPKM, FPKG for paired-end sequencing
    
    # scores must have the following structure: scores[taxonomy_id][function_id][metrics] = value
    # metrics must be 'readcount', 'rpkm', 'rpkg', fpkm', 'fpkg'
    
    # if rank parameter is not None, the resulting XLSX file will contain only entries for this rank
    if sample_id is None:
        if rank is None:
            xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), project.options.get_name() + '_' + metrics + '_functions_taxonomy.xlsx'))
        else:
            xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), project.options.get_name() + '_' + metrics + '_functions_' + rank + '_taxonomy.xlsx'))
    else:
        if rank is None:
            xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), sample_id + '_' + metrics + '_functions_taxonomy.xlsx'))
        else:
            xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), sample_id + '_' + metrics + '_functions_' + rank + '_taxonomy.xlsx'))

    print('Writing', xlsxfile)
    writer = pd.ExcelWriter(xlsxfile, engine='xlsxwriter')

    for s in project.list_samples():
        if not sample_id is None and s != sample_id:
            continue

        # Subsetting scores
        sample_scores = autovivify(3, float)
        for tax in scores.keys():
            for f in scores[tax].keys():
                if s in scores[tax][f]:
                    for k,v in scores[tax][f][s].items():
                        sample_scores[tax][f][k] = v
                

        tax_profile = TaxonomyProfile()
        tax_profile.build_functional_taxonomy_profile(project.taxonomy_data, sample_scores)

        if sample_id is None:
            df = tax_profile.convert_taxonomic_profile_into_score_df(score=metrics)
        else:
            df = tax_profile.convert_function_taxonomic_profile_into_df(score=metrics)
        
        if rank is None:
            df.to_excel(writer, sheet_name=s, merge_cells=False)
        else:
            df_filtered = df[df[('','Rank')] == rank]
            df_filtered.to_excel(writer, sheet_name=s, merge_cells=False)

        format_taxonomy_worksheet(writer, s)

    writer.save()

def format_taxonomy_worksheet(xlsx_writer, worksheet_label):
    workbook  = xlsx_writer.book
    worksheet = xlsx_writer.sheets[worksheet_label]
    superkingdom_format = workbook.add_format({'bg_color': '#FF6666'})
    phylum_format = workbook.add_format({'bg_color': '#FF9900'})
    class_format = workbook.add_format({'bg_color': '#FFCC99'})
    order_format = workbook.add_format({'bg_color': '#FFFFCC'})
    family_format = workbook.add_format({'bg_color': '#99FFCC'})
#    genus_format = workbook.add_format({'bg_color': '#99FFFF'})
    worksheet.conditional_format('B4:B1048560', {'type':     'text',
                           'criteria': 'containing',
                           'value':    'superkingdom',
                           'format':   superkingdom_format})
    worksheet.conditional_format('B4:B1048560', {'type':     'text',
                           'criteria': 'containing',
                           'value':    'phylum',
                           'format':   phylum_format})
    worksheet.conditional_format('B4:B1048560', {'type':     'text',
                           'criteria': 'containing',
                           'value':    'class',
                           'format':   class_format})
    worksheet.conditional_format('B4:B1048560', {'type':     'text',
                           'criteria': 'containing',
                           'value':    'order',
                           'format':   order_format})
    worksheet.conditional_format('B4:B1048560', {'type':     'text',
                           'criteria': 'containing',
                           'value':    'family',
                           'format':   family_format})
    #~ worksheet.conditional_format('B4:B1048560', {'type':     'text',
                           #~ 'criteria': 'containing',
                           #~ 'value':    'genus',
                           #~ 'format':   genus_format})
    worksheet.set_column(1, 1, 15)
    worksheet.set_column(2, 2, 35)
    

    
def generate_sample_taxonomy_function_xlsx(project, scores, metrics, function_id = None, rank = None):
    # This function generates XLSX file for taxa scores in all samples for one or more functions.
    # Scores may be read counts, RPKM, RPKG (for single-end sequencing) or FPKM, FPKG for paired-end sequencing
    
    # scores must have the following structure: scores[taxonomy_id][function_id][metrics] = value
    # metrics must be 'readcount', 'rpkm', 'rpkg', fpkm', 'fpkg'
    if function_id is None:
        if rank is None:
            xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), project.options.get_name() + '_' + metrics + '_samples_taxonomy.xlsx'))
        else:
            xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), project.options.get_name() + '_' + metrics + '_samples_' + rank + '_taxonomy.xlsx'))
            
    else:
        if rank is None:
            xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), function_id + '_' + metrics + '_samples_taxonomy.xlsx'))
        else:
            xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), function_id + '_' + metrics + '_samples_' + rank + '_taxonomy.xlsx'))

    print('Writing', xlsxfile)
    writer = pd.ExcelWriter(xlsxfile, engine='xlsxwriter')
    
    for function in sorted(project.ref_data.functions_dict.keys()):
        if not function_id is None and function != function_id:
            continue
    
        # Subsetting scores
        sample_scores = autovivify(3, float)
        for tax in scores.keys():
            if function in scores[tax].keys():
                for s in project.list_samples():
                    if s in scores[tax][function]:
                        for k,v in scores[tax][function][s].items():
                            sample_scores[tax][s][k] = v
                    else:
                        sample_scores[tax][s][metrics] = 0.0
                
        tax_profile = TaxonomyProfile()
        tax_profile.build_functional_taxonomy_profile(project.taxonomy_data, sample_scores)

        df = tax_profile.convert_taxonomic_profile_into_score_df(score=metrics)
        
        if rank is None:
            df.to_excel(writer, sheet_name=function, merge_cells=False)
        else:
            df_filtered = df[df[('','Rank')] == rank]
            df_filtered.to_excel(writer, sheet_name=function, merge_cells=False)
        
        format_taxonomy_worksheet(writer, function)
    
    # Make 'Average' sheet
    if function_id is None:
        sample_scores = autovivify(3, float)
        for tax in scores.keys():
            for function in sorted(project.ref_data.functions_dict.keys()):
                if function in scores[tax].keys():
                    for s in project.list_samples():
                        if s in scores[tax][function]:
                            for k,v in scores[tax][function][s].items():
                                sample_scores[tax][s][k] += v
                        else:
                            sample_scores[tax][s][metrics] += 0.0
        for tax in sample_scores:
            for s in sample_scores[tax]:
                sample_scores[tax][s][metrics] = sample_scores[tax][s][metrics]/len(project.ref_data.functions_dict.keys())

        tax_profile = TaxonomyProfile()
        tax_profile.build_functional_taxonomy_profile(project.taxonomy_data, sample_scores)

        df = tax_profile.convert_taxonomic_profile_into_score_df(score=metrics)
        
        if rank is None:
            df.to_excel(writer, sheet_name='Average', merge_cells=False)
        else:
            df_filtered = df[df[('','Rank')] == rank]
            df_filtered.to_excel(writer, sheet_name='Average', merge_cells=False)
        
        format_taxonomy_worksheet(writer, 'Average')
    
    writer.save()

#~ def create_functions_xlsx(project, metrics, sample_id_only = None):
    #~ xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), project.options.get_name() + '_' + metrics + '_functions.xlsx'))    
    #~ xlsxfile = xlsxfile.replace(' ', '_')
    #~ xlsxfile = xlsxfile.replace("'", "")
    #~ xlsxfile = xlsxfile.replace('"', '')
    #~ workbook = xlsxwriter.Workbook(xlsxfile)
    #~ bold = workbook.add_format({'bold': True})
    
    #~ functions_list = set()
    #~ categories_list = set()
    #~ scores = autovivify(3,float)
    #~ scores_cat = autovivify(3,float)
    #~ samples_list = []
    #~ # calculate scores

    #~ for sample_id in project.list_samples():
        #~ if not sample_id_only is None and sample_id_only != sample_id:
            #~ continue
        #~ samples_list.append(sample_id)
        #~ for end in project.ENDS:
            #~ # Skip reverse read for single-end samples
            #~ if end =='pe2' and not project.samples[sample_id].is_paired_end:
                #~ continue
            
            #~ # If there are no read data, import annotated reads from JSON
            #~ if end not in project.samples[sample_id].reads or project.samples[sample_id].reads[end] is None:
                #~ project.import_reads_json(sample_id, [end,])

        #~ norm_factor = 0.0
        #~ if metrics == 'Read count':
            
            #~ for end in project.samples[sample_id].reads.keys():

                #~ for read_id,read in project.samples[sample_id].reads[end].items():
                    #~ if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                        #~ for function in read.functions:
                            #~ functions_list.add(function)
                            #~ scores[function][sample_id][end] += 1.0/len(read.functions)
                            #~ category = project.ref_data.lookup_function_group(function)
                            #~ categories_list.add(category)
                            #~ scores_cat[category][sample_id][end] += 1.0/len(read.functions)
            
        #~ elif metrics in ['fpkm', 'efpkm', 'fpkg', 'efpkg']:
            #~ if not project.samples[sample_id].is_paired_end:
                #~ continue
            #~ if metrics in ['fpkm', 'efpkm']:
                #~ norm_factor = project.samples[sample_id].rpkm_scaling_factor
                #~ if norm_factor == 0.0:
                    #~ norm_factor = 1000000/project.options.get_fastq1_readcount(sample_id)
            #~ elif metrics in ['fpkg', 'efpkg']:
                #~ norm_factor = project.samples[sample_id].rpkg_scaling_factor
                #~ if norm_factor == 0.0:
                    #~ raise ValueError('FPKG scaling factor is missing')

            #~ reads_processed = set()
            #~ for read_id,read_fwd in project.samples[sample_id].reads['pe1'].items():
                #~ if read.get_status() == 'function':
                    #~ reads_processed.add(read)
                                            
                    #~ if read_id in project.samples[sample_id].reads['pe2']: 
                        #~ read_rev = project.samples[sample_id].reads['pe2'][read_id]
                        #~ if read_rev.get_status() == 'function':
                            #~ # Both ends are mapped
                            #~ fragment_scores = defaultdict(float) # Stores scores for the current read
                            #~ fragment_functions = set() # List of functions assigned to the current read
                            #~ if have_similar_functions(read_fwd, read_rev): # Ends have overlapping functions, count them together

                                #~ read1_functions = read_fwd.get_functions()
                                #~ read2_functions = read_rev.get_functions()

                                #~ fragment_functions.update(read1_functions.keys())
                                #~ fragment_functions.update(read2_functions.keys())
                                #~ # We filled list of functions. Let's calculate the score
                                
                                #~ for fragment_function in fragment_functions:
                                    #~ if fragment_function in read1_functions and fragment_function in read2_functions: # Calculate average score
                                        #~ fragment_scores[fragment_function] = max(read1_functions[fragment_function], read2_functions[fragment_function])
                                    #~ elif fragment_function in read1_functions:
                                        #~ fragment_scores[fragment_function] += read1_functions[fragment_function]
                                    #~ elif fragment_function in read2_functions:
                                        #~ fragment_scores[fragment_function] += read2_functions[fragment_function]
                            #~ else:
                                #~ fragment_functions = read_fwd.get_functions()
                                #~ for fragment_function in fragment_functions:
                                    #~ fragment_scores[fragment_function] += fragment_functions[fragment_function]
                                #~ fragment_functions = read_rev.get_functions()
                                #~ for fragment_function in fragment_functions:
                                    #~ fragment_scores[fragment_function] += fragment_functions[fragment_function]
                            
                            #~ # Add FPKM score of the read to FPKM score of the sample
                            #~ for function in fragment_scores.keys():
                                #~ functions_list.add(function)
                                #~ scores[function][sample_id][''] += norm_factor * fragment_scores[function]
                                #~ category = project.ref_data.lookup_function_group(function)
                                #~ categories_list.add(category)
                                #~ scores_cat[category][sample_id][''] += norm_factor * fragment_scores[function]
                            
                    #~ else: #Only end1 is mapped
                        #~ for function,rpk_score in read_fwd.get_functions().items():
                            #~ functions_list.add(function)
                            #~ scores[function][sample_id][''] += norm_factor * rpk_score
                            #~ category = project.ref_data.lookup_function_group(function)
                            #~ categories_list.add(category)
                            #~ scores_cat[category][sample_id][''] += norm_factor * rpk_score

            #~ for read_id,read_rev in project.samples[sample_id].reads['pe2'].items():
                #~ if read_id not in reads_processed: #Only end2 was mapped
                    #~ if read_rev.get_status() == 'function':
                        #~ read_functions = read_rev.get_functions()
                        #~ # Count FPKM
                        #~ for function,rpk_score in read_rev.get_functions().items():
                            #~ functions_list.add(function)
                            #~ scores[function][sample_id][''] += norm_factor * rpk_score
                            #~ category = project.ref_data.lookup_function_group(function)
                            #~ categories_list.add(category)
                            #~ scores_cat[category][sample_id][''] += norm_factor * rpk_score
                            
        #~ elif metrics in ['rpkm', 'erpkm', 'rpkg', 'erpkg']:
            #~ if metrics in ['rpkm', 'erpkm']:
                #~ norm_factor = project.samples[sample_id].rpkm_scaling_factor
                #~ if norm_factor == 0.0:
                    #~ norm_factor = 1000000/project.options.get_fastq1_readcount(sample_id)
            #~ elif metrics in ['rpkg', 'erpkg']:
                #~ norm_factor = project.samples[sample_id].rpkg_scaling_factor
                #~ if norm_factor == 0.0:
                    #~ raise ValueError('RPKG scaling factor is missing')

            #~ for end in project.samples[sample_id].reads.keys():
                #~ for read_id,read in project.samples[sample_id].reads[end].items():
                    #~ if read.get_status() == 'function':
                        #~ for function in read.functions:
                            #~ functions_list.add(function)
                            #~ scores[function][sample_id][end] += norm_factor * read.functions[function]
                            #~ category = project.ref_data.lookup_function_group(function)
                            #~ categories_list.add(category)
                            #~ scores_cat[category][sample_id][end] += norm_factor * read.functions[function]
        #~ else:
            #~ raise ValueError('Unknown metrics: ' + metrics)
        #~ # Free memory
        #~ project.samples[sample_id].reads = None
            
    #~ # generate output
    
    #~ # generate tables for functions
    #~ scores_worksheet = workbook.add_worksheet('Functions ' + metrics)
    
    #~ row = 0
    #~ col = 0
    
    #~ scores_worksheet.write(row, col, 'Function', bold)
    
    #~ for sample_id in samples_list:
        #~ if metrics in ['fpkm', 'efpkm', 'fpkg', 'efpkg']:
            #~ col += 1
            #~ scores_worksheet.write(row, col, sample_id, bold)
        #~ else:
            #~ for end in project.ENDS:
                #~ if end == 'pe2' and not project.samples[sample_id].is_paired_end:
                    #~ continue
                #~ col += 1
                #~ scores_worksheet.write(row, col, sample_id + '_' + end , bold)
    
    #~ col += 1
    #~ scores_worksheet.write(row, col, 'Definition', bold)
    
    #~ for function in sorted(functions_list):
        #~ row += 1
        #~ col = 0
        #~ scores_worksheet.write(row, col, function, bold)
        #~ for sample_id in samples_list:
            #~ if metrics in ['fpkm', 'efpkm', 'fpkg', 'efpkg']:
                #~ col += 1
                #~ if sample_id in scores[function]:
                    #~ scores_worksheet.write(row, col, scores[function][sample_id][''])
                #~ else:
                    #~ scores_worksheet.write(row, col, 0.0)
            #~ else:
                #~ for end in project.ENDS:
                    #~ if end == 'pe2' and not project.samples[sample_id].is_paired_end:
                        #~ continue
                    #~ col += 1
                    #~ if sample_id in scores[function] and end in scores[function][sample_id]:
                        #~ scores_worksheet.write(row, col, scores[function][sample_id][end])
                    #~ else:
                        #~ scores_worksheet.write(row, col, 0.0)
        #~ col += 1
        #~ scores_worksheet.write(row, col, project.ref_data.lookup_function_name(function))

    #~ # adjust column width
    #~ scores_worksheet.set_column(0, 0, 10)
    #~ scores_worksheet.set_column(col, col, 50)

    #~ # generate tables for categories
    #~ scores_cat_worksheet = workbook.add_worksheet('Categories RPKM')
    
    #~ row = 0
    #~ col = 0
    
    #~ scores_cat_worksheet.write(row, col, 'Category', bold)
    
    #~ for sample_id in samples_list:
        #~ if metrics in ['fpkm', 'efpkm', 'fpkg', 'efpkg']:
            #~ col += 1
            #~ scores_cat_worksheet.write(row, col, sample_id , bold)
        #~ else:
            #~ for end in project.ENDS:
                #~ if end == 'pe2' and not project.samples[sample_id].is_paired_end:
                    #~ continue
                #~ col += 1
                #~ scores_cat_worksheet.write(row, col, sample_id + '_' + end , bold)
    
#~ #    col += 1
    
    #~ for category in sorted(categories_list):
        #~ row += 1
        #~ col = 0
        #~ scores_cat_worksheet.write(row, col, category, bold)


        #~ for sample_id in samples_list:
            #~ if metrics in ['fpkm', 'efpkm', 'fpkg', 'efpkg']:
                #~ col += 1
                #~ if sample_id in scores_cat[category]:
                   #~ scores_cat_worksheet.write(row, col, scores_cat[category][sample_id][''])
                #~ else:
                    #~ scores_cat_worksheet.write(row, col, 0.0)
            
            #~ else:
                #~ for end in project.ENDS:
                    #~ if end == 'pe2' and not project.samples[sample_id].is_paired_end:
                        #~ continue
                    #~ col += 1
                    #~ if sample_id in scores_cat[category] and end in scores_cat[category][sample_id]:
                       #~ scores_cat_worksheet.write(row, col, scores_cat[category][sample_id][end])
                    #~ else:
                        #~ scores_cat_worksheet.write(row, col, 0.0)

    #~ # adjust column width
    #~ scores_cat_worksheet.set_column(0, 0, 40)

    #~ workbook.close()


#~ def create_functions_rpkm_xlsx(project):
    #~ xlsxfile = sanitize_file_name(os.path.join(project.options.get_work_dir(), project.options.get_name() + '_functions.xlsx'))    
    #~ xlsxfile = xlsxfile.replace(' ', '_')
    #~ xlsxfile = xlsxfile.replace("'", "")
    #~ xlsxfile = xlsxfile.replace('"', '')
    #~ workbook = xlsxwriter.Workbook(xlsxfile)
    #~ bold = workbook.add_format({'bold': True})
    
    #~ functions_list = set()
    #~ categories_list = set()
    #~ scores = autovivify(3, float)#defaultdict(lambda : defaultdict(float))
    #~ read_counts = autovivify(3, float)#defaultdict(lambda : defaultdict(float))
    #~ scores_cat = autovivify(3, float)#defaultdict(lambda : defaultdict(float))
    #~ read_counts_cat = autovivify(3, float)#defaultdict(lambda : defaultdict(float))
    
    #~ # calculate scores

    #~ for sample in project.list_samples():
        #~ for end in project.ENDS:
            #~ if end == 'pe2' and not project.samples[sample].is_paired_end:
                #~ continue
            #~ project.import_reads_json(sample, [end,])
            #~ scaling_factor = 1.0
            #~ if end == 'pe1':
                #~ scaling_factor = project.options.get_fastq1_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            #~ elif end == 'pe2':
                #~ scaling_factor = project.options.get_fastq2_readcount(sample)/(project.options.get_fastq1_readcount(sample) + project.options.get_fastq2_readcount(sample))
            #~ else:
                #~ raise Exception('Unknown identifier of read end: ' + end)
            #~ for read_id in project.samples[sample].reads[end]:
                #~ read = project.samples[sample][end].reads[read_id]
                #~ if read.get_status() == 'function,besthit' or read.get_status() == 'function':
                    #~ for function in read.functions:
                        #~ functions_list.add(function)
                        #~ scores[function][sample][end] += scaling_factor * read.functions[function]
                        #~ read_counts[function][sample][end] += 1.0/len(read.functions)
                        #~ category = project.ref_data.lookup_function_group(function)
                        #~ categories_list.add(category)
                        #~ scores_cat[category][sample][end] += scaling_factor * read.functions[function]
                        #~ read_counts_cat[category][sample][end] += 1.0/len(read.functions)
            #~ project.samples[sample].reads[end] = None
            
    #~ # generate output
    #~ samples_list = sorted(project.list_samples())
    
    #~ # generate tables for functions
    #~ scores_worksheet = workbook.add_worksheet('Functions RPKM')
    #~ counts_worksheet = workbook.add_worksheet('Functions read count')
    
    #~ row = 0
    #~ col = 0
    
    #~ scores_worksheet.write(row, col, 'Function', bold)
    #~ counts_worksheet.write(row, col, 'Function', bold)
    
    #~ for sample in samples_list:
        #~ for end in ENDS:
            #~ col += 1
            #~ scores_worksheet.write(row, col, sample, bold)
            #~ counts_worksheet.write(row, col, sample, bold)
    
    #~ col += 1
    #~ scores_worksheet.write(row, col, 'Definition', bold)
    #~ counts_worksheet.write(row, col, 'Definition', bold)
    
    #~ for function in sorted(functions_list):
        #~ row += 1
        #~ col = 0
        #~ scores_worksheet.write(row, col, function, bold)
        #~ counts_worksheet.write(row, col, function, bold)
        #~ for sample in samples_list:
            #~ for end in ENDS:
                #~ col += 1
                #~ if sample in scores[function]:
                    #~ #scores_worksheet.write(row, col, '{0:.3f}'.format(scores[function][sample]))
                    #~ #counts_worksheet.write(row, col, '{0:.0f}'.format(read_counts[function][sample]))
                    #~ scores_worksheet.write(row, col, scores[function][sample][end])
                    #~ counts_worksheet.write(row, col, read_counts[function][sample][end])
                #~ else:
                    #~ scores_worksheet.write(row, col, 0.0)
                    #~ counts_worksheet.write(row, col, 0)
        #~ col += 1
        #~ scores_worksheet.write(row, col, project.ref_data.lookup_function_name(function))
        #~ counts_worksheet.write(row, col, project.ref_data.lookup_function_name(function))

    #~ # adjust column width
    #~ scores_worksheet.set_column(0, 0, 10)
    #~ counts_worksheet.set_column(0, 0, 10)    
    #~ scores_worksheet.set_column(col, col, 50)
    #~ counts_worksheet.set_column(col, col, 50)

    #~ # generate tables for categories
    #~ scores_cat_worksheet = workbook.add_worksheet('Categories RPKM')
    #~ counts_cat_worksheet = workbook.add_worksheet('Categories read count')
    
    #~ row = 0
    #~ col = 0
    
    #~ scores_cat_worksheet.write(row, col, 'Category', bold)
    #~ counts_cat_worksheet.write(row, col, 'Category', bold)
    
    #~ for sample in samples_list:
        #~ for end in ENDS:
            #~ col += 1
            #~ scores_cat_worksheet.write(row, col, sample, bold)
            #~ counts_cat_worksheet.write(row, col, sample, bold)
    
    #~ col += 1
    
    #~ for category in sorted(categories_list):
        #~ row += 1
        #~ col = 0
        #~ scores_cat_worksheet.write(row, col, category, bold)
        #~ counts_cat_worksheet.write(row, col, category, bold)
        #~ for sample in samples_list:
            #~ for end in ENDS:
                #~ col += 1
                #~ if sample in scores_cat[category]:
    #~ #                scores_cat_worksheet.write(row, col, '{0:.3f}'.format(scores_cat[category][sample]))
    #~ #                counts_cat_worksheet.write(row, col, '{0:.0f}'.format(read_counts_cat[category][sample]))
                    #~ scores_cat_worksheet.write(row, col, scores_cat[category][sample][end])
                    #~ counts_cat_worksheet.write(row, col, read_counts_cat[category][sample][end])
                #~ else:
                    #~ scores_cat_worksheet.write(row, col, 0.0)
                    #~ counts_cat_worksheet.write(row, col, 0)

    #~ # adjust column width
    #~ scores_cat_worksheet.set_column(0, 0, 40)
    #~ counts_cat_worksheet.set_column(0, 0, 40)    

    #~ workbook.close()

def create_assembly_xlsx(assembler, taxonomy_data):
    xlsxfile = sanitize_file_name(os.path.join(assembler.project.options.get_assembly_dir(), 'out', assembler.project.options.get_name() + '_assembly.xlsx'))    
    xlsxfile = xlsxfile.replace(' ', '_')
    xlsxfile = xlsxfile.replace("'", "")
    xlsxfile = xlsxfile.replace('"', '')
    workbook = xlsxwriter.Workbook(xlsxfile)
    bold = workbook.add_format({'bold': True})
    cell_numformat0 = workbook.add_format()
    cell_numformat0.set_num_format('0')
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
                    reads_worksheet.write(row, col, function_read_counts[function][sample]*2, cell_numformat0)
                else:
                    reads_worksheet.write(row, col, 0, cell_numformat0)
        col += 1
        all_reads = sum(function_read_counts[function].values())*2
        reads_worksheet.write(row, col, all_reads, cell_numformat0)
        col += 1
        assembled_reads = 0
        if function in assembler.assembly.contigs:
            assembled_reads = sum([len(c.reads) for c in assembler.assembly.contigs[function].values()])
        reads_worksheet.write(row, col, assembled_reads, cell_numformat0)
        col += 1
        reads_worksheet.write(row, col, all_reads - assembled_reads, cell_numformat0)
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
    genes_worksheet.write(row, col, 'Gene start', bold)
    col += 1
    genes_worksheet.write(row, col, 'Gene end', bold)
    col += 1
    genes_worksheet.write(row, col, 'Gene length', bold)
    col += 1
    genes_worksheet.write(row, col, 'Gene strand', bold)
    col += 1
    genes_worksheet.write(row, col, 'Read count', bold)
    col += 1
    genes_worksheet.write(row, col, 'RPKM', bold)
    col += 1
    genes_worksheet.write(row, col, 'Coverage', bold)
    col += 1
    genes_worksheet.write(row, col, 'Fama gene status', bold)
    col += 1
    genes_worksheet.write(row, col, 'Fama function', bold)
    col += 1
    genes_worksheet.write(row, col, 'Fama identity', bold)
    col += 1
    genes_worksheet.write(row, col, 'CDS completeness', bold)
    col += 1
    genes_worksheet.write(row, col, 'Fama best hit', bold)
    col += 1
    genes_worksheet.write(row, col, 'Fama best hit taxonomy ID', bold)
    col += 1
    genes_worksheet.write(row, col, 'Fama best hit organism', bold)
    col += 1
    genes_worksheet.write(row, col, 'Fama best hit taxonomy', bold)
    col += 1
    genes_worksheet.write(row, col, 'Fama LCA taxonomy ID', bold)
    col += 1
    genes_worksheet.write(row, col, 'Fama LCA organism', bold)
    col += 1
    genes_worksheet.write(row, col, 'Fama LCA taxonomy', bold)

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
    col = 20
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
                    # Write Gene ID
                    genes_worksheet.write(row, col, gene_id)
                    col += 1
                    # Write Gene function from read mapping
                    genes_worksheet.write(row, col, function)
                    col += 1
                    # Write Contig ID
                    genes_worksheet.write(row, col, contig)
                    col += 1
                    # Write gene start
                    genes_worksheet.write(row, col, int(gene.start))
                    col += 1
                    # Write gene end
                    genes_worksheet.write(row, col, int(gene.end))
                    col += 1
                    # Write gene length
                    gene_length = int(gene.end) - int(gene.start) + 1
                    genes_worksheet.write(row, col, gene_length)
                    col += 1
                    # Write gene strand
                    genes_worksheet.write(row, col, gene.strand)
                    col += 1
                    # Write read count (calculated from read count of contig, adjusted by gene length)
                    gene_read_count = assembler.assembly.contigs[function][contig].get_read_count() * gene_length / len(assembler.assembly.contigs[function][contig].sequence)
                    genes_worksheet.write(row, col, gene_read_count, cell_numformat1)
                    col += 1
                    # Write RPKM
                    gene_rpkm = assembler.assembly.contigs[function][contig].get_rpkm(total_read_count) 
                    genes_worksheet.write(row, col, gene_rpkm, cell_numformat5)
                    col += 1
                    # Write coverage
                    genes_worksheet.write(row, col, assembler.assembly.contigs[function][contig].get_coverage(), cell_numformat1)
                    col += 1
                    # Write FAMA gene status
                    genes_worksheet.write(row, col, gene.get_status())
                    col += 1
                    if gene.get_status() == 'function':
                        # Write FAMA predicted functions
                        gene_functions = [y for x in gene.hit_list.get_hits() for y in x.functions]
                        genes_worksheet.write(row, col, ','.join(gene_functions))
                        col += 1
                        # Write FAMA identity
                        gene_identity = [x.get_identity() for x in gene.hit_list.get_hits()]
                        genes_worksheet.write(row, col, sum (gene_identity) / len (gene_identity), cell_numformat1)
                        col += 1
                        # Write CDS completeness
                        ref_lengths = [x.get_subject_length() for x in gene.hit_list.get_hits()]
                        genes_worksheet.write(row, col, len (gene.protein_sequence) * 100 * len(ref_lengths) / sum (ref_lengths), cell_numformat1)
                        col += 1
                        # Write FAMA best hits
                        fama_hits = [cleanup_protein_id(x.get_subject_id()) for x in gene.hit_list.get_hits()]
                        genes_worksheet.write(row, col, ','.join(fama_hits))
                        col += 1
                        # Write FAMA taxonomy ID
                        gene_taxonomy = [assembler.project.ref_data.lookup_protein_tax(cleanup_protein_id(x.get_subject_id())) for x in gene.hit_list.get_hits()]
                        genes_worksheet.write(row, col, ','.join(gene_taxonomy))
                        col += 1

                        # Write Fama best hit organism
                        gene_organism = [assembler.project.taxonomy_data.names[x]['name'] for x in gene_taxonomy]
                        genes_worksheet.write(row, col, ','.join(gene_organism))
                        col += 1
                        # Write Fama best hit taxonomy
                        best_hit_taxonomy = [assembler.project.taxonomy_data.get_lineage_string(x) for x in gene_taxonomy]
                        genes_worksheet.write(row, col, '|'.join(best_hit_taxonomy))
                        col += 1
                        
                        # Write Fama LCA taxonomy ID
                        lca_taxonomy_id = gene.taxonomy
                        genes_worksheet.write(row, col, lca_taxonomy_id)
                        col += 1
                        # Write Fama LCA organism
                        lca_organism = assembler.project.taxonomy_data.names[lca_taxonomy_id]['name']
                        genes_worksheet.write(row, col, lca_organism)
                        col += 1
                        # Write Fama LCA taxonomy
                        lca_taxonomy = assembler.project.taxonomy_data.get_lineage_string(lca_taxonomy_id)
                        genes_worksheet.write(row, col, lca_taxonomy)

                    else:
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')
                        col += 1
                        genes_worksheet.write(row, col, 'N/A')

                    #~ # Write UniRef taxon name
                    #~ taxonomy_id = gene.get_taxonomy_id()
                    #~ if taxonomy_id in taxonomy_data.names:
                        #~ genes_worksheet.write(row, col, taxonomy_data.names[taxonomy_id]['name'])
                    #~ else:
                        #~ genes_worksheet.write(row, col, 'N/A')
                    #~ col += 1

                    #~ # Write UniRef taxonomy ID
                    #~ if taxonomy_id:
                        #~ genes_worksheet.write(row, col, taxonomy_id)
                    #~ else:
                        #~ genes_worksheet.write(row, col, 'N/A')
                    #~ col += 1
                    
                    #~ if gene.uniref_hit:
                        #~ uniref_id = gene.uniref_hit.get_subject_id()
                        #~ if uniref_id:
                            #~ # Write Uniref best hit ID
                            #~ genes_worksheet.write(row, col, uniref_id)
                            #~ col += 1
                            #~ # Write Uniref best hit identity
                            #~ genes_worksheet.write(row, col, gene.uniref_hit.get_identity())
                            #~ col += 1
                            #~ # Write Uniref best hit description
                            #~ genes_worksheet.write(row, col, assembler.uniprot.get_uniprot_description(uniref_id))
                    #~ else:
                        #~ genes_worksheet.write(row, col, 'N/A')
                        #~ col += 1
                        #~ genes_worksheet.write(row, col, 'N/A')
                        #~ col += 1
                        #~ genes_worksheet.write(row, col, 'N/A')


                    for sample in samples_list:
                        col += 1
                        gene_read_count = assembler.assembly.contigs[function][contig].get_read_count(sample) * len(gene.protein_sequence) * 3 / len(assembler.assembly.contigs[function][contig].sequence)
                        genes_worksheet.write(row, col, gene_read_count, cell_numformat1)
                        col += 1
                        gene_rpkm = assembler.assembly.contigs[function][contig].get_rpkm(assembler.project.options.get_fastq1_readcount(sample),sample) #* len(gene.protein_sequence) * 3 / len(assembler.assembly.contigs[function][contig].sequence)
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

def have_similar_functions(read1, read2):
    ret_val = False
    read1_functions = read1.get_functions()
    for function in read1.get_functions():
        if function in read1_functions:
            ret_val = True
    return ret_val
