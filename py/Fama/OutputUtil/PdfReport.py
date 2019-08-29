from fpdf import FPDF,HTMLMixin
import os
from collections import defaultdict,Counter,OrderedDict
from Fama.DiamondParser.hit_utils import get_paired_end
from Fama.utils import cleanup_protein_id
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData


class MyFPDF(FPDF, HTMLMixin):
    pass
    
def generate_pdf_report(parser):
    pdf = MyFPDF('P', 'mm', 'Letter')
    pdf.add_page()

    pdf.set_font('Arial', 'B', 16)
    pdf.cell(40, 10, 'Run info')
    pdf.ln(h = '')

    pdf.set_font('Arial', 'B', 12)
    pdf.cell(60, 10, 'Sequence data')
    pdf.ln(h = '')
    pdf.set_font('Arial', '', 12)
    pdf.cell(75, 10, 'Sample ID: ' + parser.options.get_sample_name(parser.sample.sample_id))
    pdf.ln(h = 5)
    pdf.cell(90, 10, 'Paired end: ' + parser.end)
    pdf.ln(h = 5)
    pdf.cell(105, 10, 'FASTQ file: ' + parser.options.get_fastq_path(parser.sample.sample_id, parser.end))
    pdf.ln(h = 5)
    if parser.options.get_fastq_path(parser.sample.sample_id, get_paired_end(parser.end)):
        pdf.cell(120, 10, 'Paired-end FASTQ file: ' + parser.options.get_fastq_path(parser.sample.sample_id, get_paired_end(parser.end)))
        pdf.ln(h = 5)
    pdf.cell(135, 10, 'Total number of reads: ' + str(parser.options.get_fastq1_readcount(parser.sample.sample_id)))
    pdf.ln(h = 20)

    pdf.set_font('Arial', 'B', 12)
    pdf.cell(160, 10, 'Reference data ')
    pdf.ln(h = '')
    pdf.set_font('Arial', '', 12)
    pdf.cell(175, 10, 'Reference collection: ' + parser.collection)
    pdf.ln(h = 5)
    pdf.cell(190, 10, 'Number of functions in reference collection: ' + str(len(parser.ref_data.functions_dict)))
    pdf.ln(h = 5)
    pdf.cell(205, 10, 'Number of proteins in reference collection: ' + str(len(parser.ref_data.proteins_dict)))
    pdf.ln(h = 5)
    pdf.cell(220, 10, 'Reference DB size (aa): ' + str(parser.config.get_reference_db_size(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(235, 10, 'Background DB size (aa): ' + str(parser.config.get_background_db_size(parser.collection)))
    pdf.ln(h = 20)

    pdf.set_font('Arial', 'B', 12)
    pdf.cell(260, 10, 'Parameters ')
    pdf.ln(h = '')
    pdf.set_font('Arial', '', 12)
    pdf.cell(275, 10, 'Protein identity threshold (%): ' + str(parser.config.get_identity_cutoff(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(290, 10, 'Alignment length threshold (aa): ' + str(parser.config.get_length_cutoff(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(305, 10, 'Hits overlap threshold (aa): ' + str(parser.config.get_overlap_cutoff(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(320, 10, 'Bitscore range threshold (%): ' + str(100 * parser.config.get_biscore_range_cutoff(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(335, 10, 'e-value threshold for reference DB search: ' + '{:.2e}'.format(parser.config.get_evalue_cutoff(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(335, 10, 'e-value threshold for background DB search: ' + '{:.2e}'.format(parser.config.get_background_db_size(parser.options.get_collection(parser.sample.sample_id)) 
                        * parser.config.get_evalue_cutoff(parser.options.get_collection(parser.sample.sample_id))
                        / parser.config.get_reference_db_size(parser.options.get_collection(parser.sample.sample_id))))
    pdf.ln(h = 20)

    # Write search stats
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(40, 10, 'Search stats')
    pdf.ln(h = '')
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(65, 10, 'Number of reads found in reference DB search: ' + str(len(parser.reads)))
    pdf.ln(h = 10)

    read_stats = Counter()
    for read in sorted(parser.reads.keys()):
        read_stats[parser.reads[read].get_status()] += 1
    pdf.cell(65, 10, 'Number of reads found in background DB search: ' + str(len(parser.reads) - read_stats['unaccounted']))
    pdf.ln(h = 5)
    pdf.set_font('Arial', '', 12)
    
    table_rows = ['<table border="0" align="center" width="100%">',
                    '<thead><tr><th width="90%">Status</th><th width="10%">Read count</th></tr></thead>',
                    '<tbody>']
    for status in OrderedDict(read_stats.most_common()):
        if status == 'unaccounted':
            pass
        elif status == 'nofunction':
            table_rows.append('<tr><td>Reads not mapped to any function</td><td>' 
                                + str(read_stats[status]) + '</td></tr>')
        elif status == 'function':
            table_rows.append('<tr><td>Reads mapped to a function of interest</td><td>' 
                                + str(read_stats[status]) + '</td></tr>')
#        elif status == 'function,besthit':
#            table_rows.append('<tr><td>Reads mapped to a function of interest AND functional protein</td><td>' 
#                                + str(read_stats[status]) + '</td></tr>')
        else:
            table_rows.append('<tr><td>' + status + '</td><td>' 
                                + str(read_stats[status]) + '</td></tr>')
        
    table_rows.append('</tbody>\n</table>' )
    pdf.write_html('\n'.join(table_rows))
     
    # Write group scores
    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(40, 10, 'Function statistics by category')
    pdf.ln(h = 5)

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
                for function in hit.functions:
                    func_identity[parser.ref_data.lookup_function_group(function)] += hit.identity
                    func_hit_counts[parser.ref_data.lookup_function_group(function)] += 1
    for function in func_identity:
        func_identity[function] = func_identity[function]/func_hit_counts[function]

    table_rows = ['<table border="0" align="center" width="100%">',
                    '<thead><tr><<th width="55%">Definition</th><th width="15%">RPKM Score</th><th width="15%">Read count</th><th width="15%">Avg. identity</th></tr></thead>',
                    '<tbody>']
    for function in sorted(func_stats.keys()):
        table_rows.append('<tr><td>' + function
                        + '</td><td>' + '{0:.2f}'.format(func_stats[function])
                        + '</td><td>' + '{0:g}'.format(func_counts[function])
                        + '</td><td>' + '{0:.1f}'.format(func_identity[function])
                        + '</td></tr>')
    table_rows.append('</tbody></table>' )
    pdf.write_html('\n'.join(table_rows))

    # Write function scores
    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(40, 10, 'Function statistics')
    pdf.ln(h = 5)

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
                for function in hit.functions:
                    func_identity[function] += hit.identity
                    func_hit_counts[function] += 1
    for function in func_identity:
        func_identity[function] = func_identity[function]/func_hit_counts[function]

    table_rows = ['<font size="8"><table border="1" align="center" width="100%">',
                    '<thead><tr><th width="12%">ID</th><th width="61%">Definition</th><th width="9%">RPKM Score</th><th width="9%">Read count</th><th width="9%">Avg. identity</th></tr></thead>',
                    '<tbody>']
    for function in sorted(func_stats.keys()):
        function_definition = parser.ref_data.lookup_function_name(function)
        if len(function_definition) > 80:
            function_definition = function_definition[:81] + '...'
        table_rows.append('<tr><td>' + function 
                        + '</td><td>' + function_definition
                        + '</td><td>' + '{0:.2f}'.format(func_stats[function])
                        + '</td><td>' + '{0:g}'.format(func_counts[function])
                        + '</td><td>' + '{0:.1f}'.format(func_identity[function])
                        + '</td></tr>')
    table_rows.append('</tbody></table></font>' )
    pdf.write_html('\n'.join(table_rows))

    # Write taxonomy stats
    tax_stats = Counter()
    identity_stats = Counter()
    rpkm_stats = defaultdict(float)
    for read in parser.reads.keys():
        if parser.reads[read].get_status() == 'function':
            hits = parser.reads[read].get_hit_list().get_hits()
            for hit in hits:
                protein_taxid = parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.subject_id))
                tax_stats[protein_taxid] += 1
                identity_stats[protein_taxid] += hit.identity
            if len(hits) == 1:
                read_functions = parser.reads[read].get_functions()
                for function in read_functions:
                    rpkm_stats[parser.ref_data.lookup_protein_tax(cleanup_protein_id(hits[0].subject_id))] += read_functions[function]
            else:
                read_functions = parser.reads[read].get_functions()
                protein_taxids = {}
                for hit in hits:
                    hit_taxid = parser.ref_data.lookup_protein_tax(cleanup_protein_id(hit.subject_id))
                    for hit_function in hit.functions:
                        protein_taxids[hit_taxid] = hit_function
                for taxid in protein_taxids:
                    if protein_taxids[taxid] in read_functions:
                        rpkm_stats[taxid] += read_functions[protein_taxids[taxid]]
                    
#    tax_data = TaxonomyData(parser.config)
#    tax_data.load_taxdata(parser.config)
    tax_data = parser.taxonomy_data
    counts_per_rank, identity_per_rank, rpkm_per_rank = tax_data.get_taxonomy_profile(counts=tax_stats, identity=identity_stats, scores = rpkm_stats)
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for rank in ranks:

        pdf.add_page()
        pdf.set_font('Arial', 'B', 16)
        pdf.cell(40, 10, 'Taxonomy statistics for best hits (rank: ' + rank + ')')
        pdf.ln(h = 6)
        pdf.set_font('Arial', '', 10)
        pdf.cell(50, 10, '*top 100 entries are shown')

        table_rows = ['<table border="0" align="center" width="100%">',
                    '<thead><tr><<th width="60%">Taxon</th><th width="15%">RPKM score</th><th width="15%">Read count</th><th width="15%">Avg. identity</th></tr></thead>',
                    '<tbody>']
        for tax in OrderedDict(Counter(rpkm_per_rank[rank]).most_common(100)):
            table_rows.append('<tr><td>' + tax 
                            + '</td><td>' + '{0:.2f}'.format(rpkm_per_rank[rank][tax])
                            + '</td><td>' + '{0:g}'.format(counts_per_rank[rank][tax])
                            + '</td><td>' + '{0:.1f}'.format(identity_per_rank[rank][tax])
                            + '</td></tr>')
        table_rows.append('</tbody>\n</table>' )
        pdf.write_html('\n'.join(table_rows))


    outfile = os.path.join(parser.options.get_project_dir(parser.sample.sample_id), parser.options.get_output_subdir(parser.sample.sample_id),parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_report_name()+'.pdf')
    print('Writing PDF output to ', outfile)
    pdf.output(outfile, 'F')

def generate_protein_pdf_report(parser):
    pdf = MyFPDF('P', 'mm', 'Letter')
    pdf.add_page()

    pdf.set_font('Arial', 'B', 16)
    pdf.cell(40, 10, 'Run info')
    pdf.ln(h = '')

    pdf.set_font('Arial', 'B', 12)
    pdf.cell(60, 10, 'Sequence data')
    pdf.ln(h = '')
    pdf.set_font('Arial', '', 12)
    pdf.cell(75, 10, 'Sample ID: ' + parser.options.get_sample_id(parser.sample.sample_id))
    pdf.ln(h = 5)
    pdf.cell(105, 10, 'FASTA file: ' + parser.options.get_fastq_path(parser.sample.sample_id, parser.end))
    pdf.ln(h = 5)
    pdf.cell(135, 10, 'Number of reads: ' + str(parser.options.get_fastq1_readcount(parser.sample.sample_id)))
    pdf.ln(h = 20)

    pdf.set_font('Arial', 'B', 12)
    pdf.cell(160, 10, 'Reference data ')
    pdf.ln(h = '')
    pdf.set_font('Arial', '', 12)
    pdf.cell(175, 10, 'Reference collection: ' + parser.collection)
    pdf.ln(h = 5)
    pdf.cell(190, 10, 'Number of functions in reference collection: ' + str(len(parser.ref_data.functions_dict)))
    pdf.ln(h = 5)
    pdf.cell(205, 10, 'Number of proteins in reference collection: ' + str(len(parser.ref_data.proteins_dict)))
    pdf.ln(h = 5)
    pdf.cell(220, 10, 'Reference DB size (aa): ' + str(parser.config.get_reference_db_size(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(235, 10, 'Background DB size (aa): ' + str(parser.config.get_background_db_size(parser.collection)))
    pdf.ln(h = 20)

    pdf.set_font('Arial', 'B', 12)
    pdf.cell(260, 10, 'Parameters ')
    pdf.ln(h = '')
    pdf.set_font('Arial', '', 12)
    pdf.cell(275, 10, 'Protein identity threshold (%): ' + str(parser.config.get_identity_cutoff(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(290, 10, 'Alignment length threshold (aa): ' + str(parser.config.get_length_cutoff(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(305, 10, 'Hits overlap threshold (aa): ' + str(parser.config.get_overlap_cutoff(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(320, 10, 'Bitscore range threshold (%): ' + str(100 * parser.config.get_biscore_range_cutoff(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(335, 10, 'e-value threshold for reference DB search: ' + '{:.2e}'.format(parser.config.get_evalue_cutoff(parser.collection)))
    pdf.ln(h = 5)
    pdf.cell(335, 10, 'e-value threshold for background DB search: ' + '{:.2e}'.format(parser.config.get_background_db_size(parser.options.get_collection(parser.sample.sample_id)) 
                        * parser.config.get_evalue_cutoff(parser.options.get_collection(parser.sample.sample_id))
                        / parser.config.get_reference_db_size(parser.options.get_collection(parser.sample.sample_id))))
    pdf.ln(h = 20)

    # Write search stats
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(40, 10, 'Search stats')
    pdf.ln(h = '')
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(65, 10, 'Number of proteins found in reference DB search: ' + str(len(parser.reads)))
    pdf.ln(h = 10)

    read_stats = Counter()
    for read in sorted(parser.reads.keys()):
        read_stats[parser.reads[read].get_status()] += 1
    pdf.cell(65, 10, 'Number of proteins found in background DB search: ' + str(len(parser.reads) - read_stats['unaccounted']))
    pdf.ln(h = 5)
    pdf.set_font('Arial', '', 12)
    
    table_rows = ['<table border="0" align="center" width="100%">',
                    '<thead><tr><th width="90%">Status</th><th width="10%">Read count</th></tr></thead>',
                    '<tbody>']
    for status in OrderedDict(read_stats.most_common()):
        if status == 'unaccounted':
            pass
        elif status == 'nofunction':
            table_rows.append('<tr><td>Proteins not mapped to any function</td><td>' 
                                + str(read_stats[status]) + '</td></tr>')
        elif status == 'function':
            table_rows.append('<tr><td>Proteins mapped to a function of interest</td><td>' 
                                + str(read_stats[status]) + '</td></tr>')
        elif status == 'function,besthit':
            table_rows.append('<tr><td>Proteins mapped to a function of interest AND functional protein</td><td>' 
                                + str(read_stats[status]) + '</td></tr>')
        else:
            table_rows.append('<tr><td>' + status + '</td><td>' 
                                + str(read_stats[status]) + '</td></tr>')
        
    table_rows.append('</tbody>\n</table>' )
    pdf.write_html('\n'.join(table_rows))
     
    # Write group scores
    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(40, 10, 'Function statistics by category')
    pdf.ln(h = 5)

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

    table_rows = ['<table border="0" align="center" width="100%">',
                    '<thead><tr><<th width="55%">Definition</th><th width="15%">Norm. abundance</th><th width="15%">Prot. count</th><th width="15%">Avg. identity</th></tr></thead>',
                    '<tbody>']
    for function in sorted(func_stats.keys()):
        table_rows.append('<tr><td>' + function 
                        + '</td><td>' + '{0:.2f}'.format(func_stats[function])
                        + '</td><td>' + '{0:g}'.format(func_counts[function])
                        + '</td><td>' + '{0:.1f}'.format(func_identity[function])
                        + '</td></tr>')
    table_rows.append('</tbody></table>' )
    pdf.write_html('\n'.join(table_rows))

    # Write function scores
    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(40, 10, 'Function statistics')
    pdf.ln(h = 5)

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

    table_rows = ['<font size="8"><table border="1" align="center" width="100%">',
                    '<thead><tr><th width="12%">ID</th><th width="61%">Definition</th><th width="9%">Norm. abundance</th><th width="9%">Prot. count</th><th width="9%">Avg. identity</th></tr></thead>',
                    '<tbody>']
    for function in sorted(func_stats.keys()):
        function_definition = parser.ref_data.lookup_function_name(function)
        if len(function_definition) > 80:
            function_definition = function_definition[:81] + '...'
        table_rows.append('<tr><td>' + function 
                        + '</td><td>' + function_definition
                        + '</td><td>' + '{0:.2f}'.format(func_stats[function])
                        + '</td><td>' + '{0:g}'.format(func_counts[function])
                        + '</td><td>' + '{0:.1f}'.format(func_identity[function])
                        + '</td></tr>')
    table_rows.append('</tbody></table></font>' )
    pdf.write_html('\n'.join(table_rows))

    # Write taxonomy stats
    tax_stats = Counter()
    identity_stats = Counter()
    rpkm_stats = defaultdict(float)
    for read in parser.reads.keys():
        if parser.reads[read].get_status() == 'function,besthit' or parser.reads[read].get_status() == 'function':
            hits = parser.reads[read].get_hit_list().get_hits()
            for hit in hits:
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

        pdf.add_page()
        pdf.set_font('Arial', 'B', 16)
        pdf.cell(40, 10, 'Taxonomy statistics for best hits (rank: ' + rank + ')')
        pdf.ln(h = 6)
        pdf.set_font('Arial', '', 10)
        pdf.cell(50, 10, '*top 100 entries are shown')

        table_rows = ['<table border="0" align="center" width="100%">',
                    '<thead><tr><<th width="60%">Taxon</th><th width="15%">Norm. abundance</th><th width="15%">Prot. count</th><th width="15%">Avg. identity</th></tr></thead>',
                    '<tbody>']
        for tax in OrderedDict(Counter(rpkm_per_rank[rank]).most_common(100)):
            table_rows.append('<tr><td>' + tax 
                            + '</td><td>' + '{0:.2f}'.format(rpkm_per_rank[rank][tax])
                            + '</td><td>' + '{0:g}'.format(counts_per_rank[rank][tax])
                            + '</td><td>' + '{0:.1f}'.format(identity_per_rank[rank][tax])
                            + '</td></tr>')
        table_rows.append('</tbody>\n</table>' )
        pdf.write_html('\n'.join(table_rows))


    outfile = os.path.join(parser.options.get_project_dir(parser.sample.sample_id), parser.options.get_output_subdir(parser.sample.sample_id),parser.sample.sample_id + '_' + parser.end + '_'+ parser.options.get_report_name()+'.pdf')
    print('Writing PDF output to ', outfile)
    pdf.output(outfile, 'F')
