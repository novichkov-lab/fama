"""Various functions for Excel table generation"""
import os
import xlsxwriter
import pandas as pd

from lib.utils.const import STATUS_GOOD
from lib.utils.utils import autovivify, cleanup_protein_id, sanitize_file_name
from lib.taxonomy.taxonomy_profile import TaxonomyProfile


def make_function_sample_xlsx(project, scores, metric, sample_id=None):
    """Generates XLSX file for function scores for one or more samples.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        scores (dict[str, dict[str, dict[str, float]]]): outer key is function
        identifier, middle-level key is sample identifier,
        inner key is metric, value id float
        metric (str, optional): acceptable values are 'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'
        sample_id (str, optional): sample identifier
    """
    if sample_id is None:
        xlsxfile = sanitize_file_name(
            os.path.join(
                project.options.work_dir,
                project.options.project_name + '_' + metric + '_functions.xlsx'
            )
        )
    else:
        xlsxfile = sanitize_file_name(
            os.path.join(
                project.options.work_dir,
                sample_id + '_' + metric + '_functions.xlsx'
            )
        )

    print('Writing', xlsxfile)
    workbook = xlsxwriter.Workbook(xlsxfile)
    bold = workbook.add_format({'bold': True})

    functions_list = sorted(project.ref_data.functions_dict.keys())
    categories_list = sorted(
        list(
            set(
                [project.ref_data.functions_dict[x]['group']for x in
                 project.ref_data.functions_dict.keys()]
            )
        )
    )

    scores_cat = autovivify(2, float)

    # generate tables for functions
    scores_worksheet = workbook.add_worksheet('Functions ' + metric)

    row = 0
    col = 0
    scores_worksheet.write(row, col, 'Function', bold)
    for sample in project.list_samples():
        if sample_id is not None and sample != sample_id:
            continue
        col += 1
        scores_worksheet.write(row, col, sample, bold)

    col += 1
    scores_worksheet.write(row, col, 'Definition', bold)

    for function in functions_list:
        category = project.ref_data.lookup_function_group(function)
        row += 1
        col = 0
        scores_worksheet.write(row, col, function, bold)
        for sample in project.list_samples():
            if sample_id is not None and sample != sample_id:
                continue
            col += 1
            if function in scores and sample in scores[function]:
                scores_worksheet.write(row, col, scores[function][sample][metric])
                scores_cat[category][sample] += scores[function][sample][metric]
            else:
                scores_worksheet.write(row, col, 0.0)

        col += 1
        scores_worksheet.write(row, col, project.ref_data.lookup_function_name(function))

    # adjust column width
    scores_worksheet.set_column(0, 0, 10)
    scores_worksheet.set_column(col, col, 50)

    # Write worksheet for categories
    scores_cat_worksheet = workbook.add_worksheet('Categories ' + metric)
    row = 0
    col = 0
    scores_cat_worksheet.write(row, col, 'Categories', bold)

    for sample in project.list_samples():
        if sample_id is not None and sample != sample_id:
            continue
        col += 1
        scores_cat_worksheet.write(row, col, sample, bold)

    for category in categories_list:
        row += 1
        col = 0
        scores_cat_worksheet.write(row, col, category, bold)
        for sample in project.list_samples():
            if sample_id is not None and sample != sample_id:
                continue
            col += 1
            if category in scores_cat and sample in scores_cat[category]:
                scores_cat_worksheet.write(row, col, scores_cat[category][sample])
            else:
                scores_cat_worksheet.write(row, col, 0.0)
    # adjust column width
    scores_cat_worksheet.set_column(0, 0, 50)

    workbook.close()


def make_func_tax_sample_xlsx(project, scores, metric, sample_id=None, rank=None):
    """Generates XLSX file for function and taxon scores for one or more samples.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        scores (dict[str, dict[str, dict[str, float]]]): outer key is function
        identifier, middle-level key is sample identifier,
        inner key is metric, value id float
        metric (str, optional): acceptable values are 'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'
        sample_id (str, optional): sample identifier
        rank (str, optional): taxonomic rank. if rank parameter is not None, the
            resulting XLSX file will contain only entries for this rank.
    """
    if sample_id is None:
        if rank is None:
            xlsxfile = sanitize_file_name(
                os.path.join(
                    project.options.work_dir,
                    project.options.project_name + '_' + metric + '_functions_taxonomy.xlsx'
                )
            )
        else:
            xlsxfile = sanitize_file_name(
                os.path.join(
                    project.options.work_dir,
                    project.options.project_name + '_' + metric + '_functions_'
                    + rank + '_taxonomy.xlsx'
                )
            )
    else:
        if rank is None:
            xlsxfile = sanitize_file_name(
                os.path.join(
                    project.options.work_dir,
                    sample_id + '_' + metric + '_functions_taxonomy.xlsx'
                )
            )
        else:
            xlsxfile = sanitize_file_name(
                os.path.join(
                    project.options.work_dir,
                    sample_id + '_' + metric + '_functions_' + rank + '_taxonomy.xlsx'
                )
            )

    print('Writing', xlsxfile)
    writer = pd.ExcelWriter(xlsxfile, engine='xlsxwriter')

    for sample in project.list_samples():
        if sample_id is not None and sample != sample_id:
            continue

        # Subsetting scores
        sample_scores = autovivify(3, float)
        for taxonomy_id in scores.keys():
            for function_id in scores[taxonomy_id].keys():
                if sample in scores[taxonomy_id][function_id]:
                    for key, val in scores[taxonomy_id][function_id][sample].items():
                        sample_scores[taxonomy_id][function_id][key] = val

        tax_profile = TaxonomyProfile()
        tax_profile.make_function_taxonomy_profile(project.taxonomy_data, sample_scores)

        if sample_id is None:
            taxonomy_df = tax_profile.convert_profile_into_score_df(metric=metric)
        else:
            taxonomy_df = tax_profile.convert_profile_into_df(metric=metric)

        if rank is None:
            taxonomy_df.to_excel(writer, sheet_name=sample, merge_cells=False)
        else:
            filtered_df = taxonomy_df[taxonomy_df[('', 'Rank')] == rank]
            filtered_df.to_excel(writer, sheet_name=sample, merge_cells=False)

        format_taxonomy_worksheet(writer, sample)

    writer.save()


def format_taxonomy_worksheet(xlsx_writer, worksheet_label):
    """Applies formatting to a worksheet in Excel workbook

    Args:
        xlsx_writer (xlsxwriter): xlsxwriter instance with existing Workbook
        worksheet_label (str): labal of worksheet to format
    """
    workbook = xlsx_writer.book
    worksheet = xlsx_writer.sheets[worksheet_label]
    superkingdom_format = workbook.add_format({'bg_color': '#FF6666'})
    phylum_format = workbook.add_format({'bg_color': '#FF9900'})
    class_format = workbook.add_format({'bg_color': '#FFCC99'})
    order_format = workbook.add_format({'bg_color': '#FFFFCC'})
    family_format = workbook.add_format({'bg_color': '#99FFCC'})
#    genus_format = workbook.add_format({'bg_color': '#99FFFF'})
    worksheet.conditional_format('B4:B1048560', {'type': 'text',
                                                 'criteria': 'containing',
                                                 'value': 'superkingdom',
                                                 'format': superkingdom_format})
    worksheet.conditional_format('B4:B1048560', {'type': 'text',
                                                 'criteria': 'containing',
                                                 'value': 'phylum',
                                                 'format': phylum_format})
    worksheet.conditional_format('B4:B1048560', {'type': 'text',
                                                 'criteria': 'containing',
                                                 'value': 'class',
                                                 'format': class_format})
    worksheet.conditional_format('B4:B1048560', {'type': 'text',
                                                 'criteria': 'containing',
                                                 'value': 'order',
                                                 'format': order_format})
    worksheet.conditional_format('B4:B1048560', {'type': 'text',
                                                 'criteria': 'containing',
                                                 'value': 'family',
                                                 'format': family_format})
    worksheet.set_column(1, 1, 15)
    worksheet.set_column(2, 2, 35)


def make_sample_tax_func_xlsx(project, scores, metric, function_id=None, rank=None):
    """Generates XLSX file for taxa scores for one or all functions in all samples.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        scores (dict[str, dict[str, dict[str, float]]]): outer key is function
        identifier, middle-level key is sample identifier,
        inner key is metric, value id float
        metric (str, optional): acceptable values are 'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'
        function_id (str, optional): function identifier. If function_id is None, all
            functions will be included into workbook.
        rank (str, optional): taxonomic rank. if rank parameter is not None, the
            resulting XLSX file will contain only entries for this rank.
    """
    if function_id is None:
        if rank is None:
            xlsxfile = sanitize_file_name(
                os.path.join(
                    project.options.work_dir,
                    project.options.project_name + '_' + metric + '_samples_taxonomy.xlsx'
                )
            )
        else:
            xlsxfile = sanitize_file_name(
                os.path.join(
                    project.options.work_dir,
                    project.options.project_name + '_' + metric + '_samples_'
                    + rank + '_taxonomy.xlsx'
                )
            )

    else:
        if rank is None:
            xlsxfile = sanitize_file_name(
                os.path.join(
                    project.options.work_dir,
                    function_id + '_' + metric + '_samples_taxonomy.xlsx'
                )
            )
        else:
            xlsxfile = sanitize_file_name(
                os.path.join(
                    project.options.work_dir,
                    function_id + '_' + metric + '_samples_' + rank + '_taxonomy.xlsx'
                )
            )

    print('Writing', xlsxfile)
    writer = pd.ExcelWriter(xlsxfile, engine='xlsxwriter')

    for function in sorted(project.ref_data.functions_dict.keys()):
        if function_id is not None and function != function_id:
            continue

        # Subsetting scores
        sample_scores = autovivify(3, float)
        for taxonomy_id in scores.keys():
            if function in scores[taxonomy_id].keys():
                for sample in project.list_samples():
                    if sample in scores[taxonomy_id][function]:
                        for key, val in scores[taxonomy_id][function][sample].items():
                            sample_scores[taxonomy_id][sample][key] = val
                    else:
                        sample_scores[taxonomy_id][sample][metric] = 0.0

        tax_profile = TaxonomyProfile()
        tax_profile.make_function_taxonomy_profile(project.taxonomy_data, sample_scores)

        taxonomy_df = tax_profile.convert_profile_into_score_df(metric=metric)

        if rank is None:
            taxonomy_df.to_excel(writer, sheet_name=function, merge_cells=False)
        else:
            filtered_df = taxonomy_df[taxonomy_df[('', 'Rank')] == rank]
            filtered_df.to_excel(writer, sheet_name=function, merge_cells=False)
        format_taxonomy_worksheet(writer, function)

    # Make 'Average' sheet
    if function_id is None:
        sample_scores = autovivify(3, float)
        for taxonomy_id in scores:
            for function in sorted(project.ref_data.functions_dict.keys()):
                if function in scores[taxonomy_id]:
                    for sample in project.list_samples():
                        if sample in scores[taxonomy_id][function]:
                            for key, val in scores[taxonomy_id][function][sample].items():
                                sample_scores[taxonomy_id][sample][key] += val
                        else:
                            sample_scores[taxonomy_id][sample][metric] += 0.0
        for taxonomy_id in sample_scores:
            for sample in sample_scores[taxonomy_id]:
                sample_scores[taxonomy_id][sample][metric] = \
                    sample_scores[taxonomy_id][sample][metric] \
                    / len(project.ref_data.functions_dict.keys())

        tax_profile = TaxonomyProfile()
        tax_profile.make_function_taxonomy_profile(project.taxonomy_data, sample_scores)

        taxonomy_df = tax_profile.convert_profile_into_score_df(metric=metric)

        if rank is None:
            taxonomy_df.to_excel(writer, sheet_name='Average', merge_cells=False)
        else:
            filtered_df = taxonomy_df[taxonomy_df[('', 'Rank')] == rank]
            filtered_df.to_excel(writer, sheet_name='Average', merge_cells=False)

        format_taxonomy_worksheet(writer, 'Average')

    writer.save()


def make_assembly_xlsx(assembler):
    """Generates XLSX file for assembly.

    Args:
        assembler (:obj:'GeneAssembler'): gene assembler object
    """
    xlsxfile = sanitize_file_name(
        os.path.join(
            assembler.project.options.assembly_dir,
            'out',
            assembler.project.options.project_name + '_assembly.xlsx'
        )
    )
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
    function_read_counts = autovivify(2, float)  # function_read_counts[function][sample]
    gene_rpkm = autovivify(3, float)  # gene_rpkm[function][gene][sample],
    # parameters are RPKM, coverage, identity

    # count reads per function, per sample
    for function in assembler.assembly.reads:
        functions_list.add(function)
        for read in assembler.assembly.reads[function]:
            function_read_counts[function][assembler.assembly.reads[function][read]] += 1

    # collect RPKM scores for contigs per function, per sample (for contigs? for genes?)
    # calculate total read count
    total_read_count = 0

    for sample in samples_list:
        total_read_count += assembler.project.options.get_fastq1_readcount(sample)
        total_read_count += assembler.project.options.get_fastq2_readcount(sample)

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
                reads_worksheet.write(
                    row, col, function_read_counts[function][sample]*2, cell_numformat0
                )
            else:
                reads_worksheet.write(row, col, 0, cell_numformat0)
        col += 1
        all_reads = sum(function_read_counts[function].values())*2
        reads_worksheet.write(row, col, all_reads, cell_numformat0)
        col += 1
        assembled_reads = 0
        if function in assembler.assembly.contigs:
            assembled_reads = sum(
                [len(c.reads) for c in assembler.assembly.contigs[function].values()]
                )
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
                contigs_worksheet.write(
                    row, col, len(assembler.assembly.contigs[function][contig].sequence)
                    )
                col += 1
                contigs_worksheet.write(
                    row, col, assembler.assembly.contigs[function][contig].get_read_count()
                    )
                col += 1
                contigs_worksheet.write(
                    row, col, assembler.assembly.contigs[function][contig].get_rpkm(
                        total_read_count
                        ),
                    cell_numformat5
                    )
                col += 1
                contigs_worksheet.write(
                    row, col, assembler.assembly.contigs[function][contig].get_coverage(),
                    cell_numformat1
                    )
                col += 1
                contigs_worksheet.write(
                    row, col, len(assembler.assembly.contigs[function][contig].genes)
                    )
                col += 1

                for sample in samples_list:
                    contigs_worksheet.write(
                        row, col,
                        assembler.assembly.contigs[function][contig].get_read_count(sample)
                        )
                    col += 1
                    contigs_worksheet.write(
                        row, col,
                        assembler.assembly.contigs[function][contig].get_rpkm(
                            assembler.project.options.get_fastq1_readcount(sample),
                            sample
                            ),
                        cell_numformat5
                        )
                    col += 1
                    contigs_worksheet.write(
                        row, col,
                        assembler.assembly.contigs[function][contig].get_coverage(sample),
                        cell_numformat1
                        )
                    col += 1
                contigs_worksheet.write(
                    row, col, assembler.project.ref_data.lookup_function_name(function)
                    )

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
        if function not in assembler.assembly.contigs:
            continue
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
                # Write read count (calculated from read count of contig,
                # adjusted by gene length)
                gene_read_count = assembler.assembly.contigs[function][contig].get_read_count()\
                    * gene_length \
                    / len(assembler.assembly.contigs[function][contig].sequence)
                genes_worksheet.write(row, col, gene_read_count, cell_numformat1)
                col += 1
                # Write RPKM
                gene_rpkm = assembler.assembly.contigs[function][contig].get_rpkm(
                    total_read_count
                    )
                genes_worksheet.write(row, col, gene_rpkm, cell_numformat5)
                col += 1
                # Write coverage
                genes_worksheet.write(
                    row,
                    col,
                    assembler.assembly.contigs[function][contig].get_coverage(),
                    cell_numformat1
                    )
                col += 1
                # Write FAMA gene status
                genes_worksheet.write(row, col, gene.status)
                col += 1
                if gene.status == STATUS_GOOD:
                    # Write FAMA predicted functions
                    gene_functions = [y for x in gene.hit_list.hits for y in x.functions]
                    genes_worksheet.write(row, col, ','.join(gene_functions))
                    col += 1
                    # Write FAMA identity
                    gene_identity = [x.identity for x in gene.hit_list.hits]
                    genes_worksheet.write(
                        row, col,
                        sum(gene_identity) / len(gene_identity), cell_numformat1
                        )
                    col += 1
                    # Write CDS completeness
                    ref_lengths = [x.s_len for x in gene.hit_list.hits]
                    genes_worksheet.write(
                        row, col,
                        len(gene.protein_sequence) * 100 * len(ref_lengths) / sum(ref_lengths),
                        cell_numformat1
                        )
                    col += 1
                    # Write FAMA best hits
                    fama_hits = [cleanup_protein_id(x.subject_id) for x in gene.hit_list.hits]
                    genes_worksheet.write(row, col, ','.join(fama_hits))
                    col += 1
                    # Write FAMA taxonomy ID
                    gene_taxonomy = [assembler.project.ref_data.lookup_protein_tax(
                        cleanup_protein_id(x.subject_id)
                        ) for x in gene.hit_list.hits]
                    genes_worksheet.write(row, col, ','.join(gene_taxonomy))
                    col += 1

                    # Write Fama best hit organism
                    gene_organism = [assembler.project.taxonomy_data.names[x]['name'] for x
                                     in gene_taxonomy]
                    genes_worksheet.write(row, col, ','.join(gene_organism))
                    col += 1
                    # Write Fama best hit taxonomy
                    best_hit_taxonomy = [assembler.project.taxonomy_data.get_lineage_string(x)
                                         for x in gene_taxonomy]
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
                    lca_taxonomy = assembler.project.taxonomy_data.get_lineage_string(
                        lca_taxonomy_id
                        )
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

                for sample in samples_list:
                    col += 1
                    gene_read_count = assembler.assembly.contigs[function][contig].get_read_count(
                        sample
                        ) * len(
                            gene.protein_sequence
                            ) * 3 / len(
                                assembler.assembly.contigs[function][contig].sequence
                                )

                    genes_worksheet.write(row, col, gene_read_count, cell_numformat1)
                    col += 1
                    gene_rpkm = assembler.assembly.contigs[function][contig].get_rpkm(
                        assembler.project.options.get_fastq1_readcount(sample),
                        sample
                        )
                    genes_worksheet.write(row, col, gene_rpkm, cell_numformat5)
                    col += 1
                    genes_worksheet.write(
                        row, col,
                        assembler.assembly.contigs[function][contig].get_coverage(sample),
                        cell_numformat1
                        )
                col += 1
                genes_worksheet.write(
                    row, col,
                    assembler.project.ref_data.lookup_function_name(function)
                    )

    # adjust column width
    genes_worksheet.set_column(0, 0, 20)
    genes_worksheet.set_column(1, 1, 10)
    genes_worksheet.set_column(7, 9, 15)
    genes_worksheet.set_column(col, col, 50)
    workbook.close()
