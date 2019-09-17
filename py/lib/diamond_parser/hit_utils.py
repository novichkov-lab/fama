"""Various functions working with DIAMOND hits"""
from collections import defaultdict, Counter
from lib.utils.const import ENDS, STATUS_GOOD, STATUS_BAD, ROOT_TAXONOMY_ID


def get_rpkm_score(hit, function_fraction, total_readcount, length_cutoff):
    """Calculates RPKM score of a single hit (number of reads per million
        reads per kilobase of reference gene)

    Args:
        hit (:obj:DiamondHit): DIAMOND hit
        function_fraction (float): fraction of the score assigned to a function
        total_readcount (int): total number of reads in the sample
        length_cutoff(int): minimal length of acceptable alignment (in amino acids)

    Returns:
        ret_val (float): RPKM score of the hit
    """
    ret_val = None
    if (hit.s_len - length_cutoff) > 0:
        ret_val = function_fraction*1000000000.0/((hit.s_len - length_cutoff)*3*total_readcount)
    else:
        ret_val = function_fraction*1000000000.0/(3*total_readcount)
    return ret_val


def get_rpk_score(protein_length):
    """Calculates RPK score of a single hit (number of reads per kilobase of reference gene)

    Args:
        protein_length (int): length of a reference protein

    Returns:
        ret_val (float): RPK score of the hit

    """
    return 1000/3/protein_length


def get_erpk_score(protein_length, average_read_length, length_cutoff):
    """Calculates ERPK score of a single hit (effective RPK)

    Note: For very short proteins (<<90 aa) and very short reads, effective gene
    length may be negative. In such cases, this function assumes effective
    gene length is 1 bp.

    Args:
        protein_length (int): length of a reference protein
        average_read_length (float): average length of reads in the sample
        length_cutoff(int): minimal length of acceptable alignment (in amino acids)

    Returns:
        ret_val (float): ERPK score of the hit

    """
    effective_gene_length = 3*protein_length - 6*length_cutoff + average_read_length + 1
    result = 1000
    if effective_gene_length > 0:
        result = 1000/effective_gene_length
    return result


def get_fpk_score(protein_length):
    """Calculates FPK score of a single hit

    Note: FPK calculation is like RPK, but for fragment count instead of read count

    Args:
        protein_length (int): length of a reference protein

    Returns:
        ret_val (float): FPK score of the hit

    """
    return get_rpk_score(protein_length)


def get_efpk_score(protein_length, average_read_length, length_cutoff, insert_size=None):
    """Calculates EFPK score of a single hit (effective FPK)

    Note: For very short proteins (<<90 aa) and very short reads, effective gene
    length may be negative. In such cases, this function assumes effective
    gene length is 1 bp.
    Note: if insert size is not known, this function assumes insert size is
    twice more than average gene length.

    Args:
        protein_length (int): length of a reference protein
        average_read_length (float): average length of reads in the sample
        length_cutoff(int): minimal length of acceptable alignment (in amino acids)
        insert size(float): average insert size in the library (fragment
            size excluding technical sequences)

    Returns:
        ret_val (float): EFPK score of the hit

    """
    if insert_size is None:
        effective_gene_length = 3*protein_length - 6*length_cutoff + 2*average_read_length + 1
    elif insert_size > 2*average_read_length + protein_length:
        effective_gene_length = 2 * (3*protein_length - 6*length_cutoff + average_read_length + 1)
    else:
        effective_gene_length = 3*protein_length - 6*length_cutoff + insert_size + 1
    result = 1000
    if effective_gene_length > 0:
        result = 1000/effective_gene_length
    return result


def compare_hits_erpk_lca(read, hit_start, hit_end, new_hit_list, bitscore_range_cutoff,
                          length_cutoff, average_read_length, taxonomy_data,
                          ref_data, rank_cutoffs=None):
    """Compares DiamondHit object assigned to AnnotatedRead object with list of new
    DiamondHit objects, assigns scores to functions and taxonomy

    Args:
        read (:obj:AnnotatedRead): sequence read being under analysis
        hit_start (int): start position of known hit
        hit_end (int): end position of known hit
        new_hit_list (:obj:'DiamondHitList'): list of hit from new search
        bitscore_range_cutoff (float): lowest acceptaple bitscore
            (relative to top bit-score, default 0.97)
        length_cutoff (int): length of shortest acceptable alignment, in
            amino acid residues
        average_read_length (float): average length of sequence reads in the sample
        taxonomy_data (:obj:TaxonomyData): taxonomic data
        ref_data (:obj:ReferenceData): functional reference data
        rank_cutoffs (:obj:dict[str, float]): key is taxonomy rank, value
            is of amino acid identity % threshold for this rank

    This function compares one hit assigned to an annotated read with a list
    of new hits. It looks through the hit list, finds hits with bitscore
    above cutoff and assigns their functions to the read. If any hits to
    functional protein of interest are found, the read gets status 'function'.
    Otherwise, it gets status 'nofunction'.

    This function does not return anything. It sets status of read and
    assigns RPKM score to each function of the read.

    """
    # Find best hit

    for hit in read.hit_list.hits:
        if hit.q_start != hit_start or hit.q_end != hit_end:
            continue
        best_bitscore = 0.0
        best_hit = None
        # Find top hit and highest bitscore
        for new_hit in new_hit_list.hits:
            if new_hit.bitscore > best_bitscore:
                best_hit = new_hit
                best_bitscore = new_hit.bitscore

        # Set status of read
        if best_hit is None or '' in best_hit.functions:
            read.set_status(STATUS_BAD)
            return
        else:
            read.set_status(STATUS_GOOD)

        # Filter list of hits by bitscore
        bitscore_lower_cutoff = best_bitscore * (1.0 - bitscore_range_cutoff)
        selected_hits = [new_hit for new_hit in new_hit_list.hits if
                         new_hit.bitscore > bitscore_lower_cutoff]
        # Add existing hit if it has acceptable bitscore
        if hit.subject_id not in [selected_hit.subject_id for selected_hit in selected_hits] and (
                hit.bitscore >= best_bitscore
        ):
            selected_hits.append(hit)

        # Collect taxonomy IDs of all hits for LCA inference
        taxonomy_ids = set()
        # If rank-specific AAI cutoffs are not set
        if not rank_cutoffs:
            taxonomy_ids = set([ref_data.lookup_protein_tax(selected_hit.subject_id)
                                for selected_hit in selected_hits])
        # If rank-specific AAI cutoffs were calculated for the reference dataset:
        else:
            for selected_hit in selected_hits:
                subject_taxon_id = ref_data.lookup_protein_tax(selected_hit.subject_id)
                subject_rank = taxonomy_data.get_taxonomy_rank(subject_taxon_id)
                while subject_taxon_id != ROOT_TAXONOMY_ID:
                    if subject_rank not in rank_cutoffs:
                        subject_taxon_id, subject_rank = \
                            taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                    elif selected_hit.identity < rank_cutoffs[subject_rank]:
                        subject_taxon_id, subject_rank = \
                            taxonomy_data.get_upper_level_taxon(subject_taxon_id)
                    else:
                        taxonomy_ids.add(subject_taxon_id)
                        break

        # Make non-redundant list of functions from hits after filtering
        selected_functions = {}
        selected_functions_counts = Counter()
        selected_functions_data = defaultdict(dict)
        # Find best hit for each function: only one hit with
        # highest bitscore to be reported for each function
        for selected_hit in selected_hits:
            for selected_hit_function in selected_hit.functions:
                selected_functions_counts[selected_hit_function] += 1
                if selected_hit_function in selected_functions_data:
                    if (
                            selected_hit.bitscore >
                            selected_functions_data[selected_hit_function]['bit_score']
                    ):
                        selected_functions_data[selected_hit_function]['bit_score'] = \
                            selected_hit.bitscore
                        selected_functions_data[selected_hit_function]['hit'] = \
                            selected_hit
                else:
                    selected_functions_data[selected_hit_function]['bit_score'] = \
                        selected_hit.bitscore
                    selected_functions_data[selected_hit_function]['hit'] = \
                        selected_hit

        # If the most common function in new hits is unknown, set
        # status STATUS_BAD and return
        if selected_functions_counts.most_common(1)[0][0] == '':
            read.set_status(STATUS_BAD)
            return

        # Calculate RPK scores for functions
        for function in selected_functions_data:
            if function == '':
                continue
            selected_functions[function] = \
                get_erpk_score(selected_functions_data[function]['hit'].s_len,
                               average_read_length, length_cutoff)

        read.append_functions(selected_functions)

        # Set new list of hits
        for func in selected_functions_data:
            if func == '':
                continue
            good_hit = selected_functions_data[func]['hit']
            good_hit.query_id = read.read_id
            good_hit.q_start = hit_start
            good_hit.q_end = hit_end
            good_hit.annotate_hit(ref_data)
            read.hit_list.add_hit(good_hit)
        # Set read taxonomy ID
        read.taxonomy = taxonomy_data.get_lca(taxonomy_ids)
        break


def get_paired_end(end):
    """ Returns end identifier of the opposite end for paired-end reads

    Args:
        end (str): end identifier

    """
    result = ''
    if end == ENDS[0]:
        result = ENDS[1]
    elif end == ENDS[1]:
        result = ENDS[0]
    return result


def get_paired_read_id(read_id):
    """ Returns read identifier of mate read for paired-end reads

    Args:
        read_id (str): read identifier

    """
    result = read_id
    if ' ' in read_id:
        line_tokens = read_id.split(' ')
        if len(line_tokens) == 2:
            # For Casava 1.8+ format, read_id should not contain end
            # number, so result = read_id
            pass
    elif read_id.endswith('/1'):
        # old Illumina format
        result = read_id[:-1] + '2'
    elif read_id.endswith('/2'):
        # old Illumina format
        result = read_id[:-1] + '1'
    elif read_id.endswith('.1'):
        # SRA format
        result = read_id[:-1] + '2'
    elif read_id.endswith('.2'):
        # SRA format
        result = read_id[:-1] + '1'
    return result


def parse_fastq_seqid(line):
    """Extracts read identifier and end identifier from different formats of FASTQ sequence IDs

    Args:
        line (str): sequence identifier from FASTQ file

    Returns:
        read_id (str): read identifier
        end_id (str): end identifier (if any) or empty string
    """
    if line.startswith('@'):
        line = line[1:]
    result = (line, '')
    if ' ' in line:
        line_tokens = line.split(' ')
        if len(line_tokens) == 2:
            # Casava 1.8+ format
            end = line_tokens[1]
            end = end[0]
            result = (line_tokens[0], end)
        elif len(line_tokens) == 3:
            # SRA format?
            if line_tokens[0].endswith('.1') or line_tokens[0].endswith('.2'):
                # SRA format
                result = (line_tokens[0][:-2], line_tokens[0][-1])
            else:
                # unknown format
                result = (line_tokens[0], '')
        else:
            # unknown format
            result = (line_tokens[0], '')
        # return (line.split('\s')[0], line.split('\s')[1][0])
    elif line.endswith('/1') or line.endswith('/2'):
        # Old Ilumina format
        result = (line[:-2], line[-1])
    elif line.endswith('.1') or line.endswith('.2'):
        # Converted SRA
        result = (line[:-2], line[-1])
    return result


def hits_do_overlap(new_hit, hit, overlap_cutoff):
    """Returns True if query sequence of new_hit DiamondHit overlap
    by at least 'overlap_cutoff' base pairs with query sequence of
    hit DiamondHit. Otherwise, returns False.

    Args:
        new_hit(:obj:'DiamondHit'): a hit to be tested
        hit(:obj:'DiamondHit'): a hit to be tested on
        overlap_cutoff(int): minimal length of a common area between
            hits to be considered overlapping

    Returns:
        result(bool)

    """
    new_hit_start = new_hit.q_start
    new_hit_end = new_hit.q_end
    start = hit.q_start
    end = hit.q_end
    result = True
    if new_hit_start < new_hit_end:
        # existing hit on + strand
        if start > end:
            # hits on different strands: no overlap
            result = False
        # hits on + strand
        elif end < (new_hit_start + overlap_cutoff):
            result = False
        elif start > (new_hit_end - overlap_cutoff):
            result = False
        else:
            return True  # overlap
    if new_hit_start > new_hit_end:
        # existing hit on - strand
        if start < end:
            # hits on different strands: no overlap
            result = False
        # hits on - strand
        elif start < (new_hit_end + overlap_cutoff):
            result = False
        elif end > (new_hit_start - overlap_cutoff):
            result = False
    return result


def hit_overlaps_any_hits(new_hit, hit_list, overlap_cutoff):
    """Returns False if query sequence of new_hit overlap
    by at least 'overlap_cutoff' base pairs with query sequence of
    at least one DiamondHit from hit_list list of DiamondHit objects.
    Otherwise, returns True

    Args:
        new_hit(:obj:'DiamondHit'): a hit to be tested
        hit_list(:obj:'list' of :obj:'DiamondHit'): list of hits to
            be tested on
        overlap_cutoff(int): minimal length of a common area between
            hits to be considered overlapping

    Returns:
        bool

    """
    for hit in hit_list:
        if hits_do_overlap(new_hit, hit, overlap_cutoff):
            return True
    return False


def has_higher_score(new_hit, hit_list, overlap_cutoff):
    """Returns True if bitscore of new_hit is higher than bitscore of
    all overlapping DiamondHit objects from hit_list list.
    Otherwise, returns False

    Args:
        new_hit(:obj:'DiamondHit'): a hit to be tested
        hit_list(:obj:'list' of :obj:'DiamondHit'): list of hits to
            be tested
        overlap_cutoff(int): minimal length of a common area between
            hits to be considered overlapping

    Returns:
        bool

    """
    for hit in hit_list:
        if hits_do_overlap(new_hit, hit, overlap_cutoff):
            if new_hit.bitscore <= hit.bitscore:
                return False
    return True


def replace_hit(new_hit, hit_list, overlap_cutoff):
    """Compares all DiamondHit objects from a list with a
    new_hit DiamondHit object. Removes all hits, which overlap by more
    than 'overlap_cutoff' base pairs with new_hit and have lower bit-score.
    If any hits were removed, new_hit is added to the list.

    Args:
        new_hit(:obj:'DiamondHit'): a hit to be tested
        hit_list(:obj:'list' of :obj:'DiamondHit'): list of hits to
            be tested on
        overlap_cutoff(int): minimal length of a common area between
            hits to be considered overlapping

    Returns:
        :obj:'list' of :obj:'DiamondHit': updated list of DIAMOND hits

    """
    ret_val = hit_list
    append_flag = False
    for index, hit in reversed(list(enumerate(hit_list))):
        if hits_do_overlap(new_hit, hit, overlap_cutoff):
            if new_hit.bitscore > hit.bitscore:
                del ret_val[index]
                append_flag = True
    if append_flag:
        ret_val.append(new_hit)
    return ret_val
