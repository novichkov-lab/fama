"""Various functions for Report generation"""

import os
from collections import defaultdict, Counter, OrderedDict

from lib.utils.const import RANKS, STATUS_CAND, STATUS_GOOD, STATUS_BAD
from lib.utils.utils import autovivify, sanitize_file_name
from lib.diamond_parser.hit_utils import get_efpk_score, get_fpk_score
from lib.taxonomy.taxonomy_profile import TaxonomyProfile
from lib.output.xlsx_util import make_function_sample_xlsx, make_func_tax_sample_xlsx, \
    make_sample_tax_func_xlsx
from lib.output.krona_xml_writer import make_taxonomy_series_chart


def generate_fastq_report(parser):
    """Writes report for a single FASTQ file into a text file

    Args:
        parser (:obj:DiamondParser): DiamondParser object with annotated reads
    """
    outfile = os.path.join(
        parser.options.get_project_dir(parser.sample.sample_id),
        parser.options.get_output_subdir(parser.sample.sample_id),
        parser.sample.sample_id + '_' + parser.end + '_' + parser.options.report_name
        )
    with open(outfile, 'w') as out_f:
        # Write general info
        out_f.write('\nRun info\n\n')
        out_f.write('Sample ID:\t' + parser.options.get_sample_name(parser.sample.sample_id) + '\n')
        out_f.write('Paired end:\t' + parser.end + '\n')
        out_f.write(
            'FASTQ file:\t' + parser.options.get_fastq_path(
                parser.sample.sample_id, parser.end
                ) + '\n'
            )
        out_f.write(
            'Total number of reads:\t' + str(
                parser.options.get_fastq1_readcount(parser.sample.sample_id)
                )
            + '\n'
            )
        out_f.write('*****************************************\n\n')

        # Write read statistics
        out_f.write('\nRead statistics\n\n')
        read_stats = Counter()
        for read in sorted(parser.reads.keys()):
            read_stats[parser.reads[read].status] += 1
        for status in OrderedDict(read_stats.most_common()):
            if status == STATUS_CAND:
                out_f.write(
                    'Reads missing from background DB search result\t'
                    + str(read_stats[status]) + '\n'
                    )
            elif status == STATUS_BAD:
                out_f.write('Reads not mapped to any function\t' + str(read_stats[status]) + '\n')
            elif status == STATUS_GOOD:
                out_f.write(
                    'Reads mapped to a function of interest\t' + str(read_stats[status]) + '\n'
                    )
            else:
                out_f.write(status + '\t' + str(read_stats[status]) + '\n')

        out_f.write('*****************************************\n\n')

        # Write function scores
        out_f.write('\nFunction statistics\n')
        func_stats = defaultdict(float)
        func_counts = Counter()
        func_identity = defaultdict(float)
        func_hit_counts = Counter()
        for read in parser.reads.keys():
            if parser.reads[read].status == STATUS_GOOD:
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
        out_f.write('\nFunction\tDefinition\tRPK score\tRead count\tAvg. identity\n')
        for function in sorted(func_stats.keys()):
            out_f.write(
                function + '\t' + parser.ref_data.lookup_function_name(function) + '\t'
                + str(func_stats[function]) + '\t' + str(func_counts[function]) + '\t'
                + str(func_identity[function]) + '\n'
                )
        out_f.write('*****************************************\n\n')

        # Write group scores
        out_f.write('\nFunction statistics by category\n')
        func_stats = defaultdict(float)
        func_counts = Counter()
        func_identity = defaultdict(float)
        func_hit_counts = Counter()
        for read in parser.reads.keys():
            if parser.reads[read].status == STATUS_GOOD:
                functions = parser.reads[read].functions
                for function in functions:
                    func_stats[parser.ref_data.lookup_function_group(
                        function
                    )] += functions[function]
                    func_counts[parser.ref_data.lookup_function_group(
                        function
                    )] += 1/len(functions)
                for hit in parser.reads[read].hit_list.hits:
                    for function in hit.functions:
                        func_identity[parser.ref_data.lookup_function_group(
                            function
                            )] += hit.identity
                        func_hit_counts[parser.ref_data.lookup_function_group(function)] += 1
        for function in func_identity:
            func_identity[function] = func_identity[function]/func_hit_counts[function]
        out_f.write('\nCategory\tRPK score\tRead count\tAvg. identity\n')
        for function in sorted(func_stats.keys()):
            out_f.write(
                function + '\t' + str(func_stats[function]) + '\t'
                + str(func_counts[function]) + '\t' + str(func_identity[function]) + '\n'
                )

        out_f.write('*****************************************\n\n')
        # Write taxonomy stats
        out_f.write('\nTaxonomy statistics for best hits\n\n')
        tax_stats = Counter()
        identity_stats = defaultdict(float)
        rpkm_stats = defaultdict(float)
        for read in parser.reads.keys():
            if parser.reads[read].status == STATUS_GOOD:
                taxonomy = parser.reads[read].taxonomy
                if taxonomy is None:
                    print('No taxonomy ID assigned to ', read)
                    continue
                tax_stats[taxonomy] += 1
                read_functions = parser.reads[read].functions
                for funcion_id in read_functions:
                    rpkm_stats[taxonomy] += read_functions[funcion_id]
                identity_stats[taxonomy] += sum(
                    list(hit.identity for hit in parser. reads[read].hit_list.hits)
                    )
        counts_per_rank, identity_per_rank, rpkm_per_rank = get_scores_per_tax_rank(
            counts=tax_stats, identity=identity_stats, scores=rpkm_stats,
            taxonomy_data=parser.taxonomy_data
            )

        ranks = RANKS[1:]
        for rank in ranks:
            out_f.write('Taxonomy report for rank ' + rank + '\n\n')
            out_f.write('Taxon\tRead count\tRPKM score\tAverage identity\n')
            for tax in OrderedDict(Counter(counts_per_rank[rank]).most_common()):
                out_f.write(
                    rank + '\t' + tax + '\t' + str(counts_per_rank[rank][tax]) + '\t'
                    + str(rpkm_per_rank[rank][tax]) + '\t'
                    + str(identity_per_rank[rank][tax]) + '\n'
                    )
            out_f.write('*****************************************\n\n')

        out_f.write('\nList of reads\n')
        # Write list of reads
        for read in sorted(parser.reads.keys()):
            out_f.write(
                read + ': ' + parser.reads[read].status + ': ' + ','.join(
                    sorted(parser.reads[read].functions.keys())
                    )
                )
            if parser.reads[read].taxonomy is not None:
                out_f.write(' Taxonomy :' + parser.reads[read].taxonomy + '\n')
            else:
                out_f.write(' No taxonomy\n')
            hit_index = 0
            for hit in parser.reads[read].hit_list.hits:
                hit_index += 1
                out_f.write('hit' + str(hit_index) + '\t' + str(hit) + '\n')
        out_f.write('\n\n*** End of report ***\n')


def generate_fasta_report(parser):
    """Writes report for a single FASTA file into a text file

    Args:
        parser (:obj:DiamondParser): DiamondParser object with annotated reads
    """
    outfile = os.path.join(
        parser.options.get_project_dir(parser.sample.sample_id),
        parser.options.get_output_subdir(parser.sample.sample_id),
        parser.sample.sample_id + '_' + parser.end + '_' + parser.options.report_name
        )
    with open(outfile, 'w') as out_f:
        # Write general info
        out_f.write('\nRun info\n\n')
        out_f.write(
            'Sample ID:\t' + parser.options.get_sample_name(parser.sample.sample_id) + '\n'
            )
        out_f.write(
            'Sequence file:\t' + parser.options.get_fastq_path(
                parser.sample.sample_id, parser.end
                ) + '\n'
            )
        out_f.write(
            'Total number of reads:\t' + str(parser.sample.fastq_fwd_readcount) + '\n'
            )
        out_f.write('*****************************************\n\n')

        # Write read statistics
        out_f.write('\nRead statistics\n\n')
        read_stats = Counter()
        for read in sorted(parser.reads.keys()):
            read_stats[parser.reads[read].status] += 1
        for status in OrderedDict(read_stats.most_common()):
            if status == STATUS_CAND:
                out_f.write(
                    'Reads missing from background DB search result\t'
                    + str(read_stats[status]) + '\n'
                    )
            elif status == STATUS_BAD:
                out_f.write(
                    'Reads not mapped to any function\t' + str(read_stats[status]) + '\n'
                    )
            elif status == STATUS_GOOD:
                out_f.write(
                    'Reads mapped to a function of interest\t' + str(read_stats[status]) + '\n'
                    )
            else:
                out_f.write(status + '\t' + str(read_stats[status]) + '\n')
        out_f.write('*****************************************\n\n')

        # Write function scores
        out_f.write('\nFunction statistics\n')
        func_stats = defaultdict(float)
        func_counts = Counter()
        func_identity = defaultdict(float)
        func_hit_counts = Counter()
        for read in parser.reads.keys():
            if parser.reads[read].status == STATUS_GOOD:
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
        out_f.write('\nFunction\tDefinition\tRPK score\tRead count\tAvg. identity\n')
        for function in sorted(func_stats.keys()):
            out_f.write(
                function + '\t' + parser.ref_data.lookup_function_name(function) + '\t'
                + str(func_stats[function]) + '\t' + str(func_counts[function]) + '\t'
                + str(func_identity[function]) + '\n'
                )
        out_f.write('*****************************************\n\n')

        # Write group scores
        out_f.write('\nFunction statistics by category\n')
        func_stats = defaultdict(float)
        func_counts = Counter()
        func_identity = defaultdict(float)
        func_hit_counts = Counter()
        for read in parser.reads.keys():
            if parser.reads[read].status == STATUS_GOOD:
                functions = parser.reads[read].functions
                for function in functions:
                    func_stats[
                        parser.ref_data.lookup_function_group(function)
                        ] += functions[function]
                    func_counts[
                        parser.ref_data.lookup_function_group(function)
                        ] += 1/len(functions)
                for hit in parser.reads[read].hit_list.hits:
                    for function in hit.functions:
                        func_identity[
                            parser.ref_data.lookup_function_group(function)
                            ] += hit.identity
                        func_hit_counts[
                            parser.ref_data.lookup_function_group(function)
                            ] += 1
        for function in func_identity:
            func_identity[function] = func_identity[function]/func_hit_counts[function]
        out_f.write('\nCategory\tRPK score\tRead count\tAvg. identity\n')
        for function in sorted(func_stats.keys()):
            out_f.write(
                function + '\t' + str(func_stats[function]) + '\t'
                + str(func_counts[function]) + '\t' + str(func_identity[function]) + '\n'
                )
        out_f.write('*****************************************\n\n')
        # Write taxonomy stats
        out_f.write('\nTaxonomy statistics for best hits\n\n')
        tax_stats = Counter()
        identity_stats = defaultdict(float)
        rpkm_stats = defaultdict(float)
        for read in parser.reads.keys():
            if parser.reads[read].status == STATUS_GOOD:
                taxonomy = parser.reads[read].taxonomy
                if taxonomy is None:
                    print('No taxonomy ID assigned to ', read)
                    continue
                tax_stats[taxonomy] += 1
                read_functions = parser.reads[read].functions
                for function_id in read_functions:
                    rpkm_stats[taxonomy] += read_functions[function_id]
                identity_stats[taxonomy] += sum(
                    list(hit.identity for hit in parser.reads[read].hit_list.hits)
                )

        counts_per_rank, identity_per_rank, rpkm_per_rank = get_scores_per_tax_rank(
            counts=tax_stats, identity=identity_stats, scores=rpkm_stats,
            taxonomy_data=parser.taxonomy_data
            )

        ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus']
        for rank in ranks:
            out_f.write('Taxonomy report for rank ' + rank + '\n\n')
            out_f.write('Taxon\tRead count\tRPKM score\tAverage identity\n')

            for tax in OrderedDict(Counter(counts_per_rank[rank]).most_common()):
                out_f.write(
                    rank + '\t' + tax + '\t' + str(counts_per_rank[rank][tax]) + '\t'
                    + str(rpkm_per_rank[rank][tax]) + '\t'
                    + str(identity_per_rank[rank][tax]) + '\n'
                    )
            out_f.write('*****************************************\n\n')

        out_f.write('\nList of sequences and hits\n')
        # Write list of reads
        for read in sorted(parser.reads.keys()):
            out_f.write(
                read + ': ' + parser.reads[read].status + ': '
                + ','.join(sorted(parser.reads[read].functions.keys()))
                )
            if not parser.reads[read].taxonomy is None:
                out_f.write(' Taxonomy :' + parser.reads[read].taxonomy + '\n')
            else:
                out_f.write(' No taxonomy\n')
            for hit in parser.reads[read].hit_list.hits:
                out_f.write('\t' + str(hit) + '\n')
        out_f.write('\n\n*** End of report ***\n')


def get_function_scores(project, sample_id=None, metric=None):
    """Builds two-dimensional (function ID vs. sample ID) tables for read
    counts, hit counts, amino acid % identity (cumulative) and metric of choice.

    Supported metrics are:
        readcount: raw read count
        erpk: number of reads per kb of effective reference sequence
        fragmentcount: raw fragment count (for paired-end sequencing)
        fpk: number of fragments per kb of reference sequence
        efpk: number of fragments per kb of effective reference sequence
        fpkm: number of fragments per kb of reference sequence per million of fragments
        erpkm: number of reads per kb of effective reference sequence per million of reads
        efpkm: number of fragments per kb of effective reference sequence per million of fragments
        fpkg: number of fragments per kb of reference sequence per genome-equivalent
        erpkg: number of reads per kb of effective reference sequence per genome-equivalent
        efpkg: number of fragments per kb of effective reference sequence per genome-equivalent

    Note 1: Recommended metrics are erpkg for single-end sequencing and
    efpkg for paired-end sequencing.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        sample_id (str, optional): sample identifier
        metric (str, optional): acceptable values are 'readcount', 'erpk',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'erpkg', 'efpkg'

    Returns:
        three-level dictionary (defaultdict[str, defaultdict[str,defaultdict[str,float]]]):
            outer-level key is function identifier
            mid-level key is sample identifier
            inner-level keys are metric, 'count', 'hit_count', 'identity'
            values are: score for selected metric,
                        raw read count for 'count'
                        number of hits for 'hit_count'
                        cumulative amino acid % identity for 'identity'

    Note 2: 'identity' metric returned by this function is a sum of amino
    acid identity % of all hits for a function and a sample. To calculate
    average identity %, divide this value by number of hits stored under
    'hit_count' key. Cumulative identity % is additive, i.e. you can
    calculate average identity % for a group of function by summing up
    their cumulative identity % values and divide by total number of hits.

    """
    # This function actually returns read counts for readcount metric,
    # RPK for rpkm and rpkg
    # or FPK for fpkm and fpkg
    # it also always returns read count as 'count' metric, sum of identity % of
    # all best hits as 'identity' and number of best hits as 'hit_count'
    # for calculation of average identity %
    # resulting data structure is
    # ret_val[function_id][sample_id][metric|'count'|'identity'|'hit_count']

    ret_val = autovivify(3, float)
    for sample in project.list_samples():
        if sample_id is not None and sample != sample_id:
            continue
        length_cutoff = project.config.get_length_cutoff(project.options.get_collection(sample))
        average_read_length = project.samples[sample].get_avg_read_length('pe1')
        # Check if reads were processed or imported for this sample
        if project.samples[sample].reads is None or 'pe1' not in project.samples[sample].reads:
            raise KeyError('No reads data loaded for sample', sample, 'end pe1')
        if project.samples[sample].is_paired_end:
            if 'pe2' not in project.samples[sample].reads:
                raise KeyError('No reads data loaded for sample', sample, 'end pe2')

        norm_factor = 0.0
        if metric in ['readcount', 'erpk', 'fragmentcount', 'fpk', 'efpk']:
            norm_factor = 1.0
        elif metric in ['fpkm', 'erpkm', 'efpkm']:
            norm_factor = project.samples[sample].rpkm_scaling_factor
        elif metric in ['fpkg', 'erpkg', 'efpkg']:
            norm_factor = project.samples[sample].rpkg_scaling_factor
        else:
            raise ValueError('Unknown metric:' + metric)

        if norm_factor == 0.0:
            raise ValueError('Cannot get normalization factor')

        # Calculate scores
        if metric in ['readcount', 'erpk', 'erpkg', 'erpkm']:
            if project.samples[sample].is_paired_end:
                raise ValueError('No read count, RPKG and RPKM metric for paired-end input')

            for read_id, read in project.samples[sample].reads['pe1'].items():
                if read.status != STATUS_GOOD:  # Skip bad reads
                    continue
                for function in read.functions:
                    if metric == 'readcount':
                        ret_val[function][sample]['count'] += 1.0
                        ret_val[function][sample][metric] += 1.0
                    else:
                        ret_val[function][sample]['count'] += 1.0
                        ret_val[function][sample][metric] += \
                            norm_factor * read.functions[function]

                function_maxbitscores = defaultdict(dict)
                # Find max. bitscore for each function
                for hit in read.hit_list.hits:
                    for hit_function in [
                            function for function in hit.functions
                            if function in read.functions
                    ]:
                        if hit_function in function_maxbitscores:
                            if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                function_maxbitscores[hit_function]['identity'] = hit.identity
                        else:
                            function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                            function_maxbitscores[hit_function]['identity'] = hit.identity
                # Count hits and calculate cumulative %identity
                for function in read.functions:
                    if function in function_maxbitscores:
                        ret_val[function][sample]['hit_count'] += 1.0
                        ret_val[function][sample]['identity'] += function_maxbitscores[
                            function
                            ]['identity']
                    else:
                        print('Function', function, 'not found in hits of read', read_id)

        elif metric in ['fragmentcount', 'fpk', 'efpk', 'fpkg', 'fpkm', 'efpkg', 'efpkm']:
            if not project.samples[sample].is_paired_end:
                raise ValueError('Metrics based on fragment count require paired-end sequences')
            insert_size = project.get_insert_size(project.samples[sample])
            reads_processed = set()

            for read_id, read_pe1 in project.samples[sample].reads['pe1'].items():
                if read_pe1.status != STATUS_GOOD:
                    continue
                reads_processed.add(read_id)

                if 'pe2' in project.samples[sample].reads and (
                        read_id in project.samples[sample].reads['pe2']
                ):
                    read_pe2 = project.samples[sample].reads['pe2'][read_id]
                    if read_pe2.status == STATUS_GOOD:
                        # Both ends are mapped

                        fragment_functions = set()  # List of functions assigned to the current read
                        read1_functions = read_pe1.functions
                        read2_functions = read_pe2.functions

                        fragment_functions.update(read1_functions.keys())
                        fragment_functions.update(read2_functions.keys())

                        function_maxbitscores = defaultdict(dict)
                        hits1 = [hit for hit in read_pe1.hit_list.hits]
                        hits2 = [hit for hit in read_pe2.hit_list.hits]
                        # Find hit with max. bitscore for each function
                        for hit in hits1:
                            for hit_function in [
                                    function for function in hit.functions
                                    if function in fragment_functions
                            ]:
                                if hit_function in function_maxbitscores:
                                    if hit.bitscore > function_maxbitscores[hit_function][
                                            'bitscore'
                                    ]:
                                        function_maxbitscores[hit_function][
                                            'bitscore'
                                            ] = hit.bitscore
                                        function_maxbitscores[hit_function][
                                            'identity'
                                            ] = hit.identity
                                        function_maxbitscores[hit_function]['length'] = hit.s_len
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                        for hit in hits2:
                            for hit_function in [
                                    function for function in hit.functions
                                    if function in fragment_functions
                            ]:
                                if hit_function in function_maxbitscores:
                                    if hit.bitscore > function_maxbitscores[hit_function][
                                            'bitscore'
                                    ]:
                                        function_maxbitscores[hit_function][
                                            'bitscore'
                                            ] = hit.bitscore
                                        function_maxbitscores[hit_function][
                                            'identity'
                                            ] = hit.identity
                                        function_maxbitscores[hit_function]['length'] = hit.s_len
                                else:
                                    function_maxbitscores[hit_function][
                                        'bitscore'
                                        ] = hit.bitscore
                                    function_maxbitscores[hit_function][
                                        'identity'
                                        ] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                        # Count hits and calculate cumulative %identity
                        for function in fragment_functions:
                            ret_val[function][sample]['count'] += 1
                            if function in function_maxbitscores:
                                ret_val[function][sample]['hit_count'] += 1.0
                                ret_val[function][sample][
                                    'identity'
                                    ] += function_maxbitscores[function]['identity']
                                if metric == 'fragmentcount':
                                    ret_val[function][sample][metric] += 1
                                elif metric in ['fpk', 'fpkg', 'fpkm']:
                                    ret_val[function][sample][
                                        metric
                                        ] += norm_factor * get_fpk_score(
                                            function_maxbitscores[function]['length']
                                            )
                                else:
                                    ret_val[function][sample][
                                        metric
                                        ] += norm_factor * get_efpk_score(
                                            function_maxbitscores[function]['length'],
                                            average_read_length,
                                            length_cutoff,
                                            insert_size=insert_size
                                            )
                            else:
                                print('Function', function, 'not found in hits of read', read_id)

                else:  # Only end1 is mapped
                    read_functions = read_pe1.functions
                    fragment_functions = set(read_functions.keys())
                    # Count FPKM
                    # for function in fragment_functions:
                    #     ret_val[function][s]['count'] += 1
                    #     ret_val[function][s][metric] += norm_factor * read_functions[function]

                    function_maxbitscores = defaultdict(dict)
                    hits1 = [hit for hit in read_pe1.hit_list.hits]

                    # Find max. bitscore for each function
                    for hit in hits1:
                        for hit_function in [
                                function for function in hit.functions
                                if function in fragment_functions
                        ]:
                            if hit_function in function_maxbitscores:
                                if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                            else:
                                function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                function_maxbitscores[hit_function]['identity'] = hit.identity
                                function_maxbitscores[hit_function]['length'] = hit.s_len
                    # Count hits and calculate cumulative identity
                    for function in fragment_functions:
                        ret_val[function][sample]['count'] += 1
                        if function in function_maxbitscores:
                            ret_val[function][sample]['hit_count'] += 1.0
                            ret_val[function][sample][
                                'identity'
                                ] += function_maxbitscores[function]['identity']
                            if metric == 'fragmentcount':
                                ret_val[function][sample][metric] += 1
                            elif metric in ['fpk', 'fpkg', 'fpkm']:
                                ret_val[function][sample][metric] += norm_factor * get_fpk_score(
                                    function_maxbitscores[function]['length']
                                    )
                            else:
                                ret_val[function][sample][metric] += norm_factor * get_efpk_score(
                                    function_maxbitscores[function]['length'], average_read_length,
                                    length_cutoff, insert_size=insert_size
                                    )
                        else:
                            print('Function', function, 'not found in hits of read', read_id)

            for read_id, read_pe2 in project.samples[sample].reads['pe2'].items():
                if read_id in reads_processed:
                    continue  # Skip read if it was already counted
                if read_pe2.status != STATUS_GOOD:
                    continue

                fragment_functions = set()
                read_functions = read_pe2.functions
                fragment_functions.update(read_functions.keys())

                function_maxbitscores = defaultdict(dict)
                hits = [hit for hit in read_pe2.hit_list.hits]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [
                            function for function in hit.functions
                            if function in fragment_functions
                    ]:
                        if hit_function in function_maxbitscores:
                            if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                function_maxbitscores[hit_function]['identity'] = hit.identity
                                function_maxbitscores[hit_function]['length'] = hit.s_len
                        else:
                            function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                            function_maxbitscores[hit_function]['identity'] = hit.identity
                            function_maxbitscores[hit_function]['length'] = hit.s_len

                # Count hits and calculate cumulative identity
                for function in fragment_functions:
                    ret_val[function][sample]['count'] += 1
                    if function in function_maxbitscores:
                        ret_val[function][sample]['hit_count'] += 1.0
                        ret_val[function][sample]['identity'] += function_maxbitscores[function][
                            'identity'
                            ]
                        if metric == 'fragmentcount':
                            ret_val[function][sample][metric] += 1
                        elif metric in ['fpk', 'fpkg', 'fpkm']:
                            ret_val[function][sample][metric] += norm_factor * get_fpk_score(
                                function_maxbitscores[function]['length']
                                )
                        else:
                            ret_val[function][sample][metric] += norm_factor * get_efpk_score(
                                function_maxbitscores[function]['length'], average_read_length,
                                length_cutoff, insert_size=insert_size
                                )
                    else:
                        print('Function', function, 'not found in hits of read', read_id)
    return ret_val


def get_function_taxonomy_scores(project, sample_id=None, metric=None):
    """Builds three-dimensional (taxonomy ID/ function ID / sample ID) tables for read
    counts, hit counts, amino acid % identity (cumulative) and metric of choice.

    Supported metrics are:
        readcount: raw read count
        erpk: number of reads per kb of effective reference sequence
        rpkm: number of reads per kb of reference sequence per million of reads
        fragmentcount: raw fragment count (for paired-end sequencing)
        fpk: number of fragments per kb of reference sequence
        efpk: number of fragments per kb of effective reference sequence
        fpkm: number of fragments per kb of reference sequence per million of fragments
        rpkg: number of reads per kb of reference sequence per genome-equivalent
        erpkm: number of reads per kb of effective reference sequence per million of reads
        efpkm: number of fragments per kb of effective reference sequence per million of fragments
        fpkg: number of fragments per kb of reference sequence per genome-equivalent
        erpkg: number of reads per kb of effective reference sequence per genome-equivalent
        efpkg: number of fragments per kb of effective reference sequence per genome-equivalent

    Note 1: Recommended metrics are erpkg for single-end sequencing and
    efpkg for paired-end sequencing.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        sample_id (str, optional): sample identifier
        metric (str, optional): acceptable values are 'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'

    Returns:
        ret_val (defaultdict[str, defaultdict[str, defaultdict[str,defaultdict[str,float]]]]):
            four-level dictionary . where
            ret_val[taxonomy_id][function_id][sample_id][score type] = score value
            outermost level key is taxonomy identifier
            mid-outer level key is function identifier
            mid-inner level key is sample identifier
            innermost level keys are metric, 'count', 'hit_count', 'identity'
            values are: score for selected metric,
                        raw read count for 'count'
                        number of hits for 'hit_count'
                        cumulative amino acid % identity for 'identity'

    Note 2: 'identity' metric returned by this function is a sum of amino
    acid identity % of all hits for a function and a sample. To calculate
    average identity %, divide this value by number of hits stored under
    'hit_count' key. Cumulative identity % is additive, i.e. you can
    calculate average identity % for a group of function by summing up
    their cumulative identity % values and divide by total number of hits.
    """
    ret_val = autovivify(4, float)
    for sample in project.list_samples():
        if sample_id is not None and sample != sample_id:
            continue

        # Check if reads were processed or imported for this sample
        if project.samples[sample].reads is None or 'pe1' not in project.samples[sample].reads:
            raise ValueError('No reads data loaded for sample' + sample + 'end pe1')
        if project.samples[sample].is_paired_end:
            if 'pe2' not in project.samples[sample].reads:
                raise ValueError('No reads data loaded for sample' + sample + 'end pe2')

        norm_factor = 0.0
        if metric in ['readcount', 'erpk', 'fragmentcount', 'fpk', 'efpk']:
            norm_factor = 1.0
        elif metric in ['rpkm', 'fpkm', 'erpkm', 'efpkm']:
            norm_factor = project.samples[sample].rpkm_scaling_factor
        elif metric in ['rpkg', 'fpkg', 'erpkg', 'efpkg']:
            norm_factor = project.samples[sample].rpkg_scaling_factor

        if norm_factor == 0.0:
            raise ValueError('Cannot get normalization factor')

        # Calculate scores
        if metric in ['readcount', 'rpkg', 'rpkm', 'erpk', 'erpkg', 'erpkm']:
            if project.samples[sample].is_paired_end:
                raise ValueError('No Read count, RPKG and RPKM metrics for paired-end input')

            for read_id, read in project.samples[sample].reads['pe1'].items():
                if read.status != STATUS_GOOD:  # Filter unmapped reads
                    continue
                read_erpk_scores = read.functions
                for function in read_erpk_scores:
                    ret_val[read.taxonomy][function][sample]['count'] += 1.0
                    if metric == 'readcount':
                        ret_val[read.taxonomy][function][sample][metric] += 1
                    else:
                        ret_val[read.taxonomy][function][sample][
                            metric
                            ] += norm_factor * read_erpk_scores[function]

                function_maxbitscores = {}
                hits = [hit for hit in read.hit_list.hits]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [
                            function for function in hit.functions
                            if function in read.functions
                    ]:
                        if hit_function in function_maxbitscores:
                            if hit.bitscore > function_maxbitscores[hit_function]:
                                function_maxbitscores[hit_function] = hit.bitscore
                        else:
                            function_maxbitscores[hit_function] = hit.bitscore
                # Count hits and calculate cumulative identity
                for hit in hits:
                    for hit_function in [
                            function for function in hit.functions
                            if function in read.functions
                    ]:
                        if hit_function in function_maxbitscores and (
                                hit.bitscore == function_maxbitscores[hit_function]
                        ):
                            ret_val[read.taxonomy][function][sample]['hit_count'] += 1.0
                            ret_val[read.taxonomy][function][sample]['identity'] += hit.identity
                            del function_maxbitscores[hit_function]

        elif metric in ['fragmentcount', 'fpk', 'efpk', 'fpkg', 'fpkm', 'efpkg', 'efpkm']:
            if not project.samples[sample].is_paired_end:
                raise ValueError('FPKG and FPKM metric require paired-end sequences')
            insert_size = project.get_insert_size(project.samples[sample])
            reads_processed = set()
            length_cutoff = project.config.get_length_cutoff(project.options.get_collection(sample))
            average_read_length = project.samples[sample].get_avg_read_length('pe1')

            for read_id in project.samples[sample].reads['pe1'].keys():
                read_pe1 = project.samples[sample].reads['pe1'][read_id]
                if read_pe1.status != STATUS_GOOD:
                    continue
                reads_processed.add(read_id)

                if 'pe2' in project.samples[sample].reads and (
                        read_id in project.samples[sample].reads['pe2']
                ):
                    read_pe2 = project.samples[sample].reads['pe2'][read_id]
                    if read_pe2.status == STATUS_GOOD:
                        # Both ends are mapped
                        fragment_taxonomy = project.taxonomy_data.get_lca(
                            [read_pe1.taxonomy, read_pe2.taxonomy]
                            )

                        fragment_functions = set()
                        read1_functions = read_pe1.functions
                        read2_functions = read_pe2.functions

                        fragment_functions.update(read1_functions.keys())
                        fragment_functions.update(read2_functions.keys())

                        function_maxbitscores = defaultdict(dict)
                        hits1 = [hit for hit in read_pe1.hit_list.hits]
                        hits2 = [hit for hit in read_pe2.hit_list.hits]
                        # Find max. bitscore for each function
                        for hit in hits1:
                            for hit_function in [
                                    function for function in hit.functions
                                    if function in fragment_functions
                            ]:
                                if hit_function in function_maxbitscores:
                                    if hit.bitscore > function_maxbitscores[hit_function][
                                            'bitscore'
                                    ]:
                                        function_maxbitscores[hit_function][
                                            'bitscore'
                                            ] = hit.bitscore
                                        function_maxbitscores[hit_function][
                                            'identity'
                                            ] = hit.identity
                                        function_maxbitscores[hit_function]['length'] = hit.s_len
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                        for hit in hits2:
                            for hit_function in [
                                    function for function in hit.functions
                                    if function in fragment_functions
                            ]:
                                if hit_function in function_maxbitscores:
                                    if hit.bitscore > function_maxbitscores[hit_function][
                                            'bitscore'
                                    ]:
                                        function_maxbitscores[hit_function][
                                            'bitscore'
                                            ] = hit.bitscore
                                        function_maxbitscores[hit_function][
                                            'identity'
                                            ] = hit.identity
                                        function_maxbitscores[hit_function]['length'] = hit.s_len
                                else:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                        # Count hits and calculate cumulative identity of top hits for each function
                        for function in fragment_functions:
                            if function in function_maxbitscores:
                                ret_val[fragment_taxonomy][function][sample]['count'] += 1
                                ret_val[fragment_taxonomy][function][sample]['hit_count'] += 1.0
                                ret_val[fragment_taxonomy][function][sample][
                                    'identity'
                                    ] += function_maxbitscores[function]['identity']
                                if metric == 'fragmentcount':
                                    ret_val[fragment_taxonomy][function][sample][metric] += 1
                                elif metric in ['fpk', 'fpkg', 'fpkm']:
                                    ret_val[fragment_taxonomy][function][sample][
                                        metric
                                        ] += norm_factor * get_fpk_score(
                                            function_maxbitscores[function]['length']
                                            )
                                else:
                                    ret_val[fragment_taxonomy][function][sample][
                                        metric
                                        ] += norm_factor * get_efpk_score(
                                            function_maxbitscores[function]['length'],
                                            average_read_length, length_cutoff,
                                            insert_size=insert_size
                                            )
                            else:
                                print('Function', function, 'not found in hits of read', read_id)

                else:  # Only end1 is mapped
                    fragment_taxonomy = read_pe1.taxonomy
                    read_functions = read_pe1.functions
                    fragment_functions = set(read_functions.keys())
                    function_maxbitscores = defaultdict(dict)
                    hits1 = [hit for hit in read_pe1.hit_list.hits]
                    # Find max. bitscore for each function
                    for hit in hits1:
                        for hit_function in [
                                function for function in hit.functions
                                if function in fragment_functions
                        ]:
                            if hit_function in function_maxbitscores:
                                if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                    function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                    function_maxbitscores[hit_function]['identity'] = hit.identity
                                    function_maxbitscores[hit_function]['length'] = hit.s_len
                            else:
                                function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                function_maxbitscores[hit_function]['identity'] = hit.identity
                                function_maxbitscores[hit_function]['length'] = hit.s_len
                    # Count hits and calculate cumulative identity of top hit for each function
                    for function in fragment_functions:
                        if function in function_maxbitscores:
                            ret_val[fragment_taxonomy][function][sample]['count'] += 1
                            ret_val[fragment_taxonomy][function][sample]['hit_count'] += 1.0
                            ret_val[fragment_taxonomy][function][sample][
                                'identity'
                                ] += function_maxbitscores[function]['identity']
                            if metric == 'fragmentcount':
                                ret_val[fragment_taxonomy][function][sample][metric] += 1
                            elif metric in ['fpk', 'fpkg', 'fpkm']:
                                ret_val[fragment_taxonomy][function][sample][
                                    metric
                                ] += norm_factor * get_fpk_score(
                                    function_maxbitscores[function]['length']
                                )
                            else:
                                ret_val[fragment_taxonomy][function][sample][
                                    metric
                                ] += norm_factor * get_efpk_score(
                                    function_maxbitscores[function]['length'],
                                    average_read_length, length_cutoff,
                                    insert_size=insert_size
                                )
                        else:
                            print('Function', function, 'not found in hits of read', read_id)

            for read_id in project.samples[sample].reads['pe2'].keys():
                if read_id in reads_processed:
                    continue  # Skip read if it was already counted

                read_pe2 = project.samples[sample].reads['pe2'][read_id]
                if read_pe2.status != STATUS_GOOD:
                    continue

                fragment_functions = set()
                fragment_taxonomy = read_pe2.taxonomy
                read_functions = read_pe2.functions
                fragment_functions.update(read_functions.keys())

                # Build FPKM-based functional taxonomic profile
                function_maxbitscores = defaultdict(dict)
                hits = [hit for hit in read_pe2.hit_list.hits]
                # Find max. bitscore for each function
                for hit in hits:
                    for hit_function in [
                            function for function in hit.functions
                            if function in fragment_functions
                    ]:
                        if hit_function in function_maxbitscores:
                            if hit.bitscore > function_maxbitscores[hit_function]['bitscore']:
                                function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                                function_maxbitscores[hit_function]['identity'] = hit.identity
                                function_maxbitscores[hit_function]['length'] = hit.s_len
                        else:
                            function_maxbitscores[hit_function]['bitscore'] = hit.bitscore
                            function_maxbitscores[hit_function]['identity'] = hit.identity
                            function_maxbitscores[hit_function]['length'] = hit.s_len
                # Count hits and calculate cumulative %identity of top hits for each function
                for function in fragment_functions:
                    if function in function_maxbitscores:
                        ret_val[fragment_taxonomy][function][sample]['count'] += 1
                        ret_val[fragment_taxonomy][function][sample]['hit_count'] += 1.0
                        ret_val[fragment_taxonomy][function][sample][
                            'identity'
                        ] += function_maxbitscores[function]['identity']
                        if metric == 'fragmentcount':
                            ret_val[fragment_taxonomy][function][sample][metric] += 1
                        elif metric in ['fpk', 'fpkg', 'fpkm']:
                            ret_val[fragment_taxonomy][function][sample][
                                metric
                                ] += norm_factor * get_fpk_score(
                                    function_maxbitscores[function]['length']
                                    )
                        else:
                            ret_val[fragment_taxonomy][function][sample][
                                metric] += norm_factor * get_efpk_score(
                                    function_maxbitscores[function]['length'],
                                    average_read_length, length_cutoff,
                                    insert_size=insert_size
                                    )
                    else:
                        print('Function', function, 'not found in hits of read', read_id)
    return ret_val


def generate_sample_report(project, sample_id, metric=None):
    """ Creates report files for a FASTQ sample in project's directory.
    Output files include report in text format, tables in XLSX format, Krona charts.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        sample_id (str, optional): sample identifier
        metric (str, optional): acceptable values are 'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'
    """

    if metric is None:
        if project.samples[sample_id].is_paired_end:
            metric = 'efpkg'
        else:
            metric = 'erpkg'
    outfile = sanitize_file_name(os.path.join(project.options.work_dir, sample_id + '_report.txt'))
    with open(outfile, 'w') as out_f:
        out_f.write('Report for ' + sample_id + '\n\n')
        out_f.write('Sample ID:\t' + project.options.get_sample_name(sample_id) + '\n')
        out_f.write('Replicate:\t' + str(project.samples[sample_id].replicate) + '\n\n')
        out_f.write('Project name:\t' + project.options.project_name + '\n\n')
        out_f.write('Reference data set:\t' + project.options.get_collection(sample_id) + '\n\n')

        out_f.write('Source files\n')
        out_f.write('FASTQ file 1:\t' + project.options.get_fastq_path(sample_id, 'pe1') + '\n')
        out_f.write(
            '    Number of reads:\t' + str(project.samples[sample_id].fastq_fwd_readcount) + '\n'
            )
        out_f.write(
            '    Number of bases:\t' + str(project.samples[sample_id].fastq_fwd_basecount) + '\n'
            )

        if project.samples[sample_id].is_paired_end:
            out_f.write('FASTQ file 2:\t' + project.options.get_fastq_path(sample_id, 'pe2') + '\n')
            out_f.write('    Number of reads:\t'
                        + str(project.samples[sample_id].fastq_rev_readcount) + '\n')
            out_f.write('    Number of bases:\t'
                        + str(project.samples[sample_id].fastq_rev_basecount) + '\n')

        out_f.write('\nNormalization\n')
        if not project.samples[sample_id].rpkm_scaling_factor is None:
            out_f.write('RPKM normalization factor:\t'
                        + str(project.samples[sample_id].rpkm_scaling_factor) + '\n')
        if not project.samples[sample_id].rpkg_scaling_factor is None:
            out_f.write('Average genome size:\t' + format(
                project.samples[sample_id].rpkg_scaling_factor
                * project.samples[sample_id].fastq_fwd_basecount, "0.0f"
                ) + '\n')
            out_f.write('RPKG normalization factor:\t'
                        + str(project.samples[sample_id].rpkg_scaling_factor) + '\n')
        if project.samples[sample_id].is_paired_end:
            out_f.write('Average insert size:\t' + str(
                project.get_insert_size(project.samples[sample_id])
                ) + '\n')

        out_f.write('\nNumber of mapped reads\n')
        out_f.write(
            'FASTQ file 1:\t' + str(len(project.samples[sample_id].reads['pe1'])) + '\n'
            )
        if project.samples[sample_id].is_paired_end:
            out_f.write(
                'FASTQ file 2:\t' + str(len(project.samples[sample_id].reads['pe2'])) + '\n'
                )

    if project.samples[sample_id].rpkg_scaling_factor is None and (
            metric in ['efpkg', 'fpkg', 'erpkg', 'rpkg']
    ):
        raise ValueError('Not enough data to normalize by average genome size')
    if project.samples[sample_id].rpkm_scaling_factor is None and (
            metric in ['efpkm', 'fpkm', 'erpkm', 'rpkm']
    ):
        raise ValueError('Not enough data to normalize by sample size')

    scores_function = get_function_scores(project, sample_id, metric=metric)
    scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id, metric=metric)
    make_function_sample_xlsx(
        project, scores_function, metric=metric, sample_id=sample_id
        )
    make_func_tax_sample_xlsx(
        project, scores_function_taxonomy, metric=metric, sample_id=sample_id
        )
    sample_scores_taxonomy = autovivify(3, float)
    for tax in scores_function_taxonomy.keys():
        for function_id in scores_function_taxonomy[tax].keys():
            if sample_id in scores_function_taxonomy[tax][function_id]:
                sample_scores_taxonomy[tax][function_id] = \
                    scores_function_taxonomy[tax][function_id][sample_id]
                #~ for key, val in scores_function_taxonomy[tax][function_id][sample_id].items():
                    #~ sample_scores_taxonomy[tax][function_id][key] = val

    tax_profile = TaxonomyProfile()
    outfile = sanitize_file_name(
        os.path.join(
            project.options.work_dir,
            sample_id + '_' + metric + '_functional_taxonomy_profile.xml'
            )
        )
    tax_profile.make_function_taxonomy_profile(
        project.taxonomy_data, sample_scores_taxonomy)
    functions = sorted(scores_function.keys())  # project.ref_data.functions_dict.keys())
    make_taxonomy_series_chart(
        tax_profile, functions, outfile, project.config.krona_path, metric=metric
        )


def generate_protein_sample_report(project, sample_id, metric=None):
    """ Creates report files for a protein sample in project's directory.
    Output files include report in text format, tables in XLSX format, Krona charts.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        sample_id (str, optional): sample identifier
        metric (str, optional): scoring metric. Acceptable values are
            'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'
    """
    outfile = sanitize_file_name(os.path.join(project.options.work_dir,
                                              sample_id + '_report.txt'))
    with open(outfile, 'w') as out_f:
        out_f.write('Report for ' + sample_id + '\n\n')
        out_f.write('Sample ID:\t' + project.options.get_sample_name(sample_id) + '\n')
        out_f.write('Replicate:\t' + str(project.samples[sample_id].replicate) + '\n\n')
        out_f.write('Project name:\t' + project.options.project_name + '\n\n')
        out_f.write('Reference data set:\t' + project.options.get_collection(sample_id) + '\n\n')

        out_f.write('Source files\n')
        out_f.write('Input file:\t' + project.options.get_fastq_path(sample_id, 'pe1') + '\n')
        out_f.write(
            '    Number of proteins:\t' + str(project.samples[sample_id].fastq_fwd_readcount) + '\n'
            )
        out_f.write(
            '    Number of amino acids:\t' + str(
                project.samples[sample_id].fastq_fwd_basecount) + '\n'
            )

        out_f.write(
            '\nNumber of mapped proteins:' + str(
                len(project.samples[sample_id].reads['pe1'])) + '\n'
            )

        for read in sorted(project.samples[sample_id].reads['pe1'].keys()):
            protein = project.samples[sample_id].reads['pe1'][read]
            if protein.status == STATUS_BAD:
                continue
            out_f.write(read + ': ' + ','.join(sorted(protein.functions.keys())))
            tax_id = protein.taxonomy
            if tax_id is None:
                out_f.write(' No taxonomy\n')
            else:
                out_f.write(' Taxonomy :' + project.taxonomy_data.names[tax_id]['name']
                            + '(' + tax_id + ')\n')
            out_f.write(' Hits:\n')
            for hit in protein.hit_list.hits:
                out_f.write('\t' + str(hit) + '\n')

    if not project.samples[sample_id].reads['pe1']:
        return

    scores_function = get_function_scores(project, sample_id, metric=metric)
    scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id, metric=metric)
    make_function_sample_xlsx(
        project, scores_function, metric=metric, sample_id=sample_id
        )
    make_func_tax_sample_xlsx(
        project, scores_function_taxonomy, metric=metric, sample_id=sample_id
        )
    sample_scores_taxonomy = autovivify(3, float)
    for tax in scores_function_taxonomy.keys():
        for function_id in scores_function_taxonomy[tax].keys():
            if sample_id in scores_function_taxonomy[tax][function_id]:
                for key, val in scores_function_taxonomy[tax][function_id][sample_id].items():
                    sample_scores_taxonomy[tax][function_id][key] = val

    tax_profile = TaxonomyProfile()
    outfile = sanitize_file_name(
        os.path.join(
            project.options.work_dir,
            sample_id + '_' + metric + '_functional_taxonomy_profile.xml'
            )
        )
    tax_profile.make_function_taxonomy_profile(
        project.taxonomy_data, sample_scores_taxonomy
        )
    function_list = sorted(project.ref_data.functions_dict.keys())
    make_taxonomy_series_chart(
        tax_profile, function_list, outfile, project.config.krona_path, metric=metric
        )


def generate_project_report(project, metric=None):
    """ Creates project report files for nucleotide project.
    Output files include report in text format, tables in XLSX format, Krona charts.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        metric (str, optional): acceptable values are 'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'
    """
    is_paired_end = None
    for sample_id in project.list_samples():
        if project.samples[sample_id].is_paired_end:
            if is_paired_end is None:
                is_paired_end = True
            elif not is_paired_end:
                print('Project contains both single-end and paired-end sequences.',
                      'No comparative tables will be generated.')
                return
        else:
            if is_paired_end is None:
                is_paired_end = False
            elif is_paired_end:
                print('Project contains both single-end and paired-end sequences.',
                      'No comparative tables will be generated.')
                return

        # If there are no read data, try to load them from JSON
        if not project.samples[sample_id].reads:
            project.import_reads_json(sample_id, project.ENDS)
    print('Generating spreadsheets and interactive diagrams for the project...')

    if metric is None:
        if is_paired_end:
            metric = 'efpkg'
        else:
            metric = 'erpkg'

    scores = get_function_scores(project, sample_id=None, metric=metric)
    make_function_sample_xlsx(project, scores, metric=metric)
    scores_function_taxonomy = get_function_taxonomy_scores(
        project, sample_id=None, metric=metric
        )
    make_func_tax_sample_xlsx(
        project, scores_function_taxonomy, metric=metric, sample_id=None
        )
    make_sample_tax_func_xlsx(
        project, scores_function_taxonomy, metric=metric, function_id=None
        )
    generate_project_mkdocs(project, scores, metric=metric)


def generate_protein_project_report(project):
    """ Creates project report files for protein project.
    Output files include report in text format, tables in XLSX format, Krona charts.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
    """
    for sample_id in project.list_samples():
        # If there are no read data, try to load them from JSON
        if not project.samples[sample_id].reads:
            project.import_reads_json(sample_id, project.ENDS)
    print('Generating spreadsheets and interactive diagrams for the project...')
    outfile = os.path.join(project.options.work_dir, 'project_report.txt')
    with open(outfile, 'w') as out_f:
        out_f.write(project.options.project_name + '\n\n')
        for sample_id in project.list_samples():
            out_f.write(
                sample_id + ':\t' + project.samples[sample_id].sample_name
                + '\tproteins mapped: ' + str(len(project.samples[sample_id].reads['pe1']))
                )
            out_f.write('\n')
        out_f.write('\nList of mapped proteins\n')

        for sample_id in project.list_samples():
            for protein_id in sorted(project.samples[sample_id].reads['pe1'].keys()):
                protein = project.samples[sample_id].reads['pe1'][protein_id]
                if protein.status == STATUS_BAD:
                    continue
                for hit in protein.hit_list.hits:
                    out_f.write(
                        '\t'.join([
                            sample_id,
                            protein_id,
                            project.taxonomy_data.names[protein.taxonomy]['name'],
                            str(hit)
                            ]) + '\n'
                        )
    metric = 'readcount'
    scores = get_function_scores(project, sample_id=None, metric=metric)
    make_function_sample_xlsx(project, scores, metric=metric)
    scores_function_taxonomy = get_function_taxonomy_scores(
        project, sample_id=None, metric=metric
        )
    make_func_tax_sample_xlsx(project, scores_function_taxonomy, metric=metric, sample_id=None)
    make_sample_tax_func_xlsx(project, scores_function_taxonomy, metric=metric, function_id=None)


def generate_project_mkdocs(project, scores, sample_id=None, metric=None):
    """ Creates project markdown documents.
    Output files include report in text format, tables in XLSX format, Krona charts.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        scores (dict[str, dict[str, dict[str, float]]]): outer key is function
        identifier, middle-level key is sample identifier,
        inner key is metric, value id float
        sample_id (str, optional): sample identifier
        metric (str, optional): acceptable values are 'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'
    """
    functions_list = sorted(scores.keys())
    # write output
    if sample_id is None:
        samples_list = sorted(project.list_samples())
    else:
        samples_list = [sample_id]

    outfile = sanitize_file_name(os.path.join(project.options.work_dir, 'index.md'))

    with open(outfile, 'w') as out_f:
        out_f.write('# ')
        out_f.write(project.options.project_name)
        out_f.write('\n\n')

        out_f.write('\n## Functions\n\n' + metric + ' scores of individual functions\n\n')

        out_f.write('Function|')
        out_f.write('|'.join(samples_list))
        out_f.write('|Definition|\n')
        out_f.write('|---'*(len(samples_list) + 2))
        out_f.write('|\n')

        for function in sorted(functions_list):
            out_f.write(function)
            for sample in samples_list:
                out_f.write('|')
                if function in scores and sample in scores[function]:
                    out_f.write('{0:.3f}'.format(scores[function][sample][metric]))
                else:
                    out_f.write('0.000')
            out_f.write('|')
            out_f.write(project.ref_data.lookup_function_name(function))
            out_f.write('|\n')

        targetfile = sanitize_file_name(os.path.join(
            'data', project.options.project_name + '_functions.xlsx'
        ))

        out_f.write('\n<a href="')
        out_f.write(targetfile)
        out_f.write(
            '" target="_blank">Download table of RPKM scores and read counts in XLSX format</a>\n\n'
            )

        out_f.write('## Taxonomy profile for all reads mapped to nitrogen cycle genes\n\n')

        targetfile = sanitize_file_name(
            os.path.join('data', project.options.project_name + '_taxonomy_profile.xml.html')
            )
        out_f.write('<div class="krona-wrapper">\n<iframe src="')
        out_f.write(targetfile)
        out_f.write('" height="800" width="100%">')
        out_f.write(project.options.project_name)
        out_f.write(' taxonomy profile</iframe>\n<br>\n<a href="')
        out_f.write(targetfile)
        out_f.write('" target="_blank">Open chart in a new window</a>\n</div>\n\n')

        out_f.write('## Taxonomy profiles for individual functions\n\n')
        targetfile = sanitize_file_name(
            os.path.join('data', project.options.project_name + '_functions_taxonomy.xlsx')
            )
        out_f.write('<a href="')
        out_f.write(targetfile)
        out_f.write(
            '" target="_blank">Download detailed taxonomic profile for all functions '
            + 'in all samples (XLSX format)</a>\n\n<div>\n'
            )
        for sample in samples_list:
            targetfile = sanitize_file_name(
                os.path.join('data', sample + '_functional_taxonomy_profile.xml.html')
                )
            out_f.write('<a href="')
            out_f.write(targetfile)
            out_f.write('" target="_blank">Taxonomy profile for individual functions in sample ')
            out_f.write(sample)
            out_f.write(' (interactive chart)</a><br>\n')
        out_f.write('</div>\n')

    # write files for samples
    for sample in samples_list:
        outfile = sanitize_file_name(os.path.join(project.options.work_dir, sample + '.md'))
        with open(outfile, 'w') as out_f:
            out_f.write('# Sample ')
            out_f.write(sample)
            out_f.write('\n\n## Functional taxonomy profile for reads mapped to ')
            out_f.write(os.path.join(project.options.get_collection()))
            out_f.write(' dataset\n\n')

            targetfile = sanitize_file_name(os.path.join(
                '..', 'data', sample + '_functional_taxonomy_profile.xml.html'
            ))
            out_f.write('<div class="krona-wrapper">\n<iframe src="')
            out_f.write(targetfile)
            out_f.write('" height="800" width="100%">')
            out_f.write(project.options.project_name)
            out_f.write(' taxonomy profile</iframe>\n<br>\n<a href="')
            out_f.write(targetfile)
            out_f.write('" target="_blank">Open chart in a new window</a>\n</div>\n\n')

            out_f.write('## Reports for FASTQ files\n\n')
            if os.path.exists(
                    os.path.join(
                        project.options.get_project_dir(sample),
                        project.options.get_output_subdir(sample),
                        sample + '_pe1_' + project.options.report_name + '.pdf'
                    )
            ):
                targetfile = os.path.join(
                    '..', 'data', sample + '_pe1_' + project.options.report_name + '.pdf'
                    )
                out_f.write('<a href="')
                out_f.write(targetfile)
                out_f.write('">Download report for read end 1 (PDF format)</a>\n<br>\n')
            if os.path.exists(
                    os.path.join(
                        project.options.get_project_dir(sample),
                        project.options.get_output_subdir(sample),
                        sample + '_pe2_' + project.options.report_name + '.pdf'
                    )
            ):
                targetfile = os.path.join(
                    '..', 'data', sample + '_pe2_' + project.options.report_name + '.pdf'
                    )
                out_f.write('<a href="')
                out_f.write(targetfile)
                out_f.write('">Download report for read end 2 (PDF format)</a>\n<br>\n')


def generate_functions_stamp_input(project, scores, metric):
    """ Creates functional profile in text format for downstream analysis with STAMP.
    There are two output files with data and metadata.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        scores (dict[str, dict[str, dict[str, float]]]): outer key is function
        identifier, middle-level key is sample identifier,
        inner key is metric, value id float
        metric (str, optional): acceptable values are 'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'
    """
    scores_outfile = sanitize_file_name(
        os.path.join(
            project.options.work_dir,
            project.options.project_name + '_' + metric + '_functions.stamp.tsv'
            )
        )
    metadata_outfile = sanitize_file_name(
        os.path.join(
            project.options.work_dir,
            project.options.project_name + '_functions.metadata.stamp.tsv'
            )
        )

    with open(scores_outfile, 'w') as out_f:
        out_f.write('Category\tFunction\t' + '\t'.join(project.list_samples()) + '\n')
        for function in sorted(scores.keys()):
            out_f.write(project.ref_data.lookup_function_group(function) + '\t' + function)
            for sample_id in project.list_samples():
                if sample_id in scores[function]:
                    out_f.write('\t' + str(scores[function][sample_id][metric]))
                else:
                    out_f.write('\t0')
            out_f.write('\n')

    with open(metadata_outfile, 'w') as out_f:
        out_f.write('SampleID\tSample_name\tReplicate\n')
        for sample_id in project.list_samples():
            out_f.write(
                sample_id + '\t' + project.samples[sample_id].sample_name
                + '\t' + project.samples[sample_id].replicate + '\n'
                )


def print_stamp_dataseries_taxonomy(tax_profile, sample_list, metric,
                                    taxonomy_id='1', prefix='', taxonomy_level=0):
    """Makes a taxon entry for export as STAMP input files. Recursively calls
    itself for all children taxa.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        sample_list (list of str): sample identifiers
        metric (str, optional): acceptable values are 'readcount', 'erpk',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'erpkg', 'efpkg'
        taxonomy_id (str): taxonomy identifier
        prefix (str): taxonomy lineage, tab-separated
        taxonomy_level (int): number of taxonomy level

    Returns:
        ret_val (str): taxon data
        attribute_values (defaultdict[str, defaultdict[str, float]]): outer key
            is sample identifier, inner key is score type, value is score value.
    """
    attribute_values = autovivify(2, float)

    if taxonomy_id not in tax_profile.tree.data:
        raise ValueError(taxonomy_id + 'not found in the tree!!!')

    if taxonomy_id == '1':  # top level
        new_prefix = prefix
    elif prefix == '':  # superkingdom level: do not start with tab
        new_prefix = tax_profile.tree.data[taxonomy_id].name
    else:
        new_prefix = prefix + '\t' + tax_profile.tree.data[taxonomy_id].name
    new_taxonomy_level = taxonomy_level + 1
    ret_val = ''
    if tax_profile.tree.data[taxonomy_id].children:
        for child_taxid in tax_profile.tree.data[taxonomy_id].children:
            child_node, child_values = print_stamp_dataseries_taxonomy(
                tax_profile, sample_list, metric, taxonomy_id=child_taxid,
                prefix=new_prefix, taxonomy_level=new_taxonomy_level
                )
            ret_val += child_node
            for datapoint in child_values.keys():
                for key, val in child_values[datapoint].items():
                    attribute_values[datapoint][key] += val

        unclassified_flag = False
#        print(tax_profile.tree.data[taxid].attributes)
        for sample_id in sample_list:
            if sample_id in tax_profile.tree.data[taxonomy_id].attributes:
                if attribute_values[sample_id][metric] < tax_profile.tree.data[
                        taxonomy_id
                ].attributes[sample_id][metric]:
                    unclassified_flag = True
                    break

        if unclassified_flag:
            if taxonomy_id == '1':  # line sohuld not start with tab symbol
                ret_val += 'unclassified' + '\tunclassified' * (len(
                    tax_profile.RANKS
                    ) - new_taxonomy_level - 1)
            else:
                ret_val += new_prefix + '\tUnclassified ' + tax_profile.tree.data[
                    taxonomy_id
                    ].name + '\tunclassified' * (len(tax_profile.RANKS) - new_taxonomy_level - 1)
            for sample_id in sample_list:
                if sample_id in tax_profile.tree.data[taxonomy_id].attributes and (
                        attribute_values[sample_id][metric]
                        < tax_profile.tree.data[taxonomy_id].attributes[sample_id][metric]
                ):
                    ret_val += '\t' + str(
                        tax_profile.tree.data[taxonomy_id].attributes[sample_id][
                            metric
                            ] - attribute_values[sample_id][metric]
                        )
                else:
                    ret_val += '\t0'
            ret_val += '\n'
    else:
        ret_val += new_prefix + '\tunclassified' * (len(tax_profile.RANKS) - new_taxonomy_level)
        if tax_profile.tree.data[taxonomy_id].attributes:
            for sample_id in sample_list:
                if sample_id in tax_profile.tree.data[taxonomy_id].attributes and (
                        metric in tax_profile.tree.data[taxonomy_id].attributes[sample_id]
                ):
                    ret_val += '\t' + str(
                        tax_profile.tree.data[taxonomy_id].attributes[sample_id][metric]
                        )
                else:
                    ret_val += '\t0'
        else:
            ret_val += '\t0' * len(sample_list)
        ret_val += '\n'

    attribute_values = autovivify(1)
    for sample_id in sample_list:
        if sample_id in tax_profile.tree.data[taxonomy_id].attributes and (
                metric in tax_profile.tree.data[taxonomy_id].attributes[sample_id]
        ):
            attribute_values[sample_id][metric] = \
                tax_profile.tree.data[taxonomy_id].attributes[sample_id][metric]
    return ret_val, attribute_values


def generate_func_tax_stamp_input(project, scores, metric):
    """ Creates functional-taxonomic profile in text format for downstream analysis with STAMP.
    There are two output files with data and metadata.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        scores (dict[str, dict[str, dict[str, float]]]): outer key is function
        identifier, middle-level key is sample identifier,
        inner key is metric, value id float
        metric (str, optional): acceptable values are 'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'
    """
    metadata_outfile = sanitize_file_name(
        os.path.join(
            project.options.work_dir,
            project.options.project_name + '_functions_taxonomy.metadata.stamp.tsv'
        )
    )
    sample_list = project.list_samples()
    root_id = '1'

    for function in project.ref_data.list_functions():
        # Subsetting scores
        sample_scores = autovivify(3, float)
        for tax in scores.keys():
            if function in scores[tax].keys():
                for sample_id in project.list_samples():
                    if sample_id in scores[tax][function] and (
                            metric in scores[tax][function][sample_id]
                    ):
                        sample_scores[tax][sample_id][metric] = \
                            scores[tax][function][sample_id][metric]
                    else:
                        sample_scores[tax][sample_id][metric] = 0.0

        if not sample_scores:
            continue
        scores_outfile = sanitize_file_name(
            os.path.join(
                project.options.work_dir,
                project.options.project_name + '_' + metric + '_'
                + function + '_taxonomy.stamp.tsv'
                )
            )

        with open(scores_outfile, 'w') as out_f:
            out_f.write(
                '\t'.join(
                    project.taxonomy_data.RANKS[1:-1]
                    )
                + '\t' + '\t'.join(project.list_samples()) + '\n'
                )
            tax_profile = TaxonomyProfile()
            tax_profile.make_function_taxonomy_profile(project.taxonomy_data, sample_scores)
            child_nodes, _ = print_stamp_dataseries_taxonomy(
                tax_profile, sample_list, metric, taxonomy_id=root_id, prefix='', taxonomy_level=0
                )
            out_f.write(child_nodes)

    with open(metadata_outfile, 'w') as out_f:
        out_f.write('SampleID\tSample_name\tReplicate\n')
        for sample_id in project.list_samples():
            out_f.write(
                sample_id + '\t' + project.samples[sample_id].sample_name + '\t'
                + project.samples[sample_id].replicate + '\n'
                )


def generate_sample_text_report(project, sample_id, metric=None):
    """ Creates sample report in text format.

    Args:
        project (:obj:'Project'): Project object that stores all annotated reads
        sample_id (str, optional): sample identifier
        metric (str, optional): acceptable values are 'readcount', 'erpk', 'rpkm',
            'fragmentcount', 'fpk', 'efpk', 'fpkm', 'erpkm', 'efpkm',
            'fpkg', 'rpkg', 'erpkg', 'efpkg'
    """
    # This function creates output files only in project's working directory
    if metric is None:
        if project.samples[sample_id].is_paired_end:
            metric = 'efpkg'
        else:
            metric = 'erpkg'
    # outfile = os.path.join(project.samples[sample_id].work_directory, sample_id + '_report.txt')
    outfile = sanitize_file_name(os.path.join(project.options.work_dir, sample_id + '_report.txt'))
    with open(outfile, 'w') as out_f:
        out_f.write('Report for ' + sample_id + '\n\n')
        out_f.write('Sample ID:\t' + project.options.get_sample_name(sample_id) + '\n')
        out_f.write('Replicate:\t' + str(project.samples[sample_id].replicate) + '\n\n')
        out_f.write('Project name:\t' + project.options.project_name + '\n\n')
        out_f.write('Reference data set:\t' + project.options.get_collection(sample_id) + '\n\n')

        out_f.write('Source files\n')
        out_f.write('FASTQ file 1:\t' + project.options.get_fastq_path(sample_id, 'pe1') + '\n')
        out_f.write(
            '    Number of reads:\t' + str(project.samples[sample_id].fastq_fwd_readcount) + '\n'
            )
        out_f.write(
            '    Number of bases:\t' + str(project.samples[sample_id].fastq_fwd_basecount) + '\n'
            )

        if project.samples[sample_id].is_paired_end:
            out_f.write('FASTQ file 2:\t' + project.options.get_fastq_path(sample_id, 'pe2') + '\n')
            out_f.write('    Number of reads:\t'
                        + str(project.samples[sample_id].fastq_rev_readcount) + '\n')
            out_f.write('    Number of bases:\t'
                        + str(project.samples[sample_id].fastq_rev_basecount) + '\n')

        out_f.write('\nNormalization\n')
        if not project.samples[sample_id].rpkm_scaling_factor is None:
            out_f.write(
                'RPKM normalization factor:\t' + str(
                    project.samples[sample_id].rpkm_scaling_factor
                    )
                + '\n'
                )
        if not project.samples[sample_id].rpkg_scaling_factor is None:
            out_f.write(
                'Average genome size:\t' + format(
                    project.samples[sample_id].rpkg_scaling_factor
                    * project.samples[sample_id].fastq_fwd_basecount, "0.0f"
                    )
                + '\n'
                )
            out_f.write(
                'RPKG normalization factor:\t' + str(
                    project.samples[sample_id].rpkg_scaling_factor
                    )
                + '\n'
                )
        if project.samples[sample_id].is_paired_end:
            out_f.write(
                'Average insert size:\t' + str(
                    project.get_insert_size(project.samples[sample_id])
                    )
                + '\n'
                )

        out_f.write('\nNumber of mapped reads\n')
        out_f.write(
            'FASTQ file 1:\t' + str(
                len(project.samples[sample_id].reads['pe1'])
                )
            + '\n'
            )
        if project.samples[sample_id].is_paired_end:
            out_f.write(
                'FASTQ file 2:\t' + str(
                    len(project.samples[sample_id].reads['pe2'])
                    )
                + '\n'
                )

    if project.samples[sample_id].rpkg_scaling_factor is None and (
            metric in ['efpkg', 'fpkg', 'erpkg', 'rpkg']
    ):
        raise ValueError('Not enough data to normalize by average genome size')
    if project.samples[sample_id].rpkm_scaling_factor is None and (
            metric in ['efpkm', 'fpkm', 'erpkm', 'rpkm']
    ):
        raise ValueError('Not enough data to normalize by sample size')

    scores_function = get_function_scores(project, sample_id, metric=metric)
    scores_function_taxonomy = get_function_taxonomy_scores(project, sample_id, metric=metric)
    make_function_sample_xlsx(
        project, scores_function, metric=metric, sample_id=sample_id
        )
    make_func_tax_sample_xlsx(
        project, scores_function_taxonomy, metric=metric, sample_id=sample_id
        )
    sample_scores_taxonomy = autovivify(3, float)
    for tax in scores_function_taxonomy.keys():
        for function_id in scores_function_taxonomy[tax].keys():
            if sample_id in scores_function_taxonomy[tax][function_id]:
                for key, val in scores_function_taxonomy[tax][function_id][sample_id].items():
                    sample_scores_taxonomy[tax][function_id][key] = val

    tax_profile = TaxonomyProfile()
    outfile = sanitize_file_name(
        os.path.join(
            project.options.work_dir,
            sample_id + '_' + metric + '_functional_taxonomy_profile.xml'
            )
        )
    tax_profile.make_function_taxonomy_profile(project.taxonomy_data, sample_scores_taxonomy)
    function_list = sorted(project.ref_data.functions_dict.keys())
    make_taxonomy_series_chart(
        tax_profile, function_list, outfile, project.config.krona_path, metric=metric
        )


def generate_assembly_report(assembler):
    """ Creates assemly report in text format.

    Args:
        assembler (:obj:'GeneAssembler'): gene assembler object
    """
    outfile = os.path.join(assembler.assembly_dir, 'out', 'assembly_report.txt')
    with open(outfile, 'w') as out_f:
        out_f.write('\nFunction statistics\n\n')
        for function in assembler.assembly.contigs:
            out_f.write(function + ':\n')
            for contig in assembler.assembly.contigs[function]:
                for gene_id in assembler.assembly.contigs[function][contig].genes:
                    gene = assembler.assembly.contigs[function][contig].genes[gene_id]
                    if gene.status == STATUS_GOOD:
                        out_f.write(
                            function + '\t' + contig + '\t' + str(
                                len(assembler.assembly.contigs[function][contig].sequence)
                                )
                            + '\t' + gene_id + '\t' + ','.join(gene.functions.keys())
                            )
                        for sample in assembler.assembly.contigs[function][contig].read_count:
                            out_f.write(
                                '\t' + sample + ':read count ' + str(
                                    assembler.assembly.contigs[function][contig].read_count[sample]
                                    )
                                + ';coverage ' + str(
                                    assembler.assembly.contigs[function][contig].get_coverage(
                                        sample
                                        )
                                    )
                                + ';rpkm ' + str(
                                    assembler.assembly.contigs[function][contig].get_rpkm(
                                        sample,
                                        assembler.project.options.get_fastq1_readcount(sample))
                                    )
                                )
                        out_f.write('\n')
        out_f.write('\n\n')
        for function in assembler.assembly.reads:
            out_f.write(function + '\t' + str(len(assembler.assembly.reads[function])) + '\n')
        out_f.write('\n\n*** End of report ***\n')


def get_scores_per_tax_rank(counts, identity, scores, taxonomy_data):
    """Calculates taxonomy profile for a number of taxonomy identifiers.
    This function takes three dictionaries, assuming that all of
    them have equal size and identical keys. This function is used only for
    generation of text and PDF reports for individual FASTQ/FASTA files.

    Args:
        counts (dict[str,int]): key is taxonomy identifier, value is count
        identity (dict[str,float]): key is taxonomy identifier, value is amino acid % identity
        scores (dict[str,float]): key is taxonomy identifier, value is score

    Returns:
        counts_per_rank (defaultdict[str, defaultdict[str, float]]):
            external key is rank, internal key is name, value is count
        identity_per_rank (defaultdict[str, defaultdict[str, float]]):
            external key is rank, internal key is name, value is amino acid % identity
        scores_per_rank (defaultdict[str, defaultdict[str, float]]): ():
            external key is rank, internal key is name, value is score
    """
    unknown_rank = 'superkingdom'
    cellular_organisms_taxid = '131567'
    non_cellular_organisms_name = 'Non-cellular'
    non_cellular_organisms_rank = 'superkingdom'

    rpkm_per_rank = defaultdict(lambda: defaultdict(float))
    counts_per_rank = defaultdict(lambda: defaultdict(int))
    identity_per_rank = defaultdict(lambda: defaultdict(float))

    for taxid in counts:
        current_id = taxid
        if taxid == UNKNOWN_TAXONOMY_ID:
            label = taxonomy_data.get_name(UNKNOWN_TAXONOMY_ID)
            rpkm_per_rank[unknown_rank][label] += scores[taxid]
            counts_per_rank[unknown_rank][label] += counts[taxid]
            identity_per_rank[unknown_rank][label] += identity[taxid]
            continue
        is_cellular = False
        not_found = False
        while current_id != ROOT_TAXONOMY_ID:
            if current_id == cellular_organisms_taxid:
                is_cellular = True
                break
            if not taxonomy_data.is_exist(current_id):
                print('A) ncbi_code not found in ncbi_nodes: \'' + current_id + '\'')
                not_found = True
                break
            current_id = taxonomy_data.get_parent(current_id)

        if not_found:
            continue

        if not is_cellular:
            rpkm_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] \
                += scores[taxid]
            counts_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] \
                += counts[taxid]
            identity_per_rank[non_cellular_organisms_rank][non_cellular_organisms_name] \
                += identity[taxid]
            continue

        current_id = taxid
        while current_id != ROOT_TAXONOMY_ID:
            if not taxonomy_data.is_exist(current_id):
                print('B) Got nothing for ncbi_code in ncbi_nodes: ' + current_id)
                break
            rank = taxonomy_data.get_rank(current_id)
            if rank in RANKS:
                name = taxonomy_data.get_name(current_id)
                rpkm_per_rank[rank][name] += scores[taxid]
                counts_per_rank[rank][name] += counts[taxid]
                identity_per_rank[rank][name] += identity[taxid]
            current_id = taxonomy_data.get_parent(current_id)

    for rank in identity_per_rank:
        for taxon in identity_per_rank[rank]:
            identity_per_rank[rank][taxon] = \
                identity_per_rank[rank][taxon]/counts_per_rank[rank][taxon]

    return counts_per_rank, identity_per_rank, rpkm_per_rank
