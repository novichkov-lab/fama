"""Functions for JSON-based serialization and deserialization of Fama objects"""
import os
import json
from lib.project.sample import Sample
from lib.diamond_parser.diamond_hit import DiamondHit
from lib.diamond_parser.diamond_hit_list import DiamondHitList
from lib.sequences.annotated_read import AnnotatedRead
from lib.gene_assembler.contig import Contig
from lib.gene_assembler.gene import Gene
from lib.gene_assembler.gene_assembly import GeneAssembly


def decode_reads(obj):
    """Custom JSON decoder for AnnotatedRead object

    Args:
        obj (obj): AnnotatedRead object to decode
    """
    if '__AnnotatedRead__' in obj:
        annotated_read = AnnotatedRead()
        annotated_read.__dict__.update(obj['__AnnotatedRead__'])
        return annotated_read
    elif '__DiamondHitList__' in obj:
        diamond_hit_list = DiamondHitList()
        diamond_hit_list.__dict__.update(obj['__DiamondHitList__'])
        return diamond_hit_list
    elif '__DiamondHit__' in obj:
        diamond_hit = DiamondHit()
        diamond_hit.__dict__.update(obj['__DiamondHit__'])
        return diamond_hit
    return obj


def decode_sample(obj):
    """Custom JSON decoder for Sample object

    Args:
        obj (obj): Sample object to decode
    """
    if '__Sample__' in obj:
        sample = Sample()
        sample.__dict__.update(obj['__Sample__'])
        return sample
    elif '__AnnotatedRead__' in obj:
        annotated_read = AnnotatedRead()
        annotated_read.__dict__.update(obj['__AnnotatedRead__'])
        return annotated_read
    elif '__DiamondHitList__' in obj:
        diamond_hit_list = DiamondHitList()
        diamond_hit_list.__dict__.update(obj['__DiamondHitList__'])
        return diamond_hit_list
    elif '__DiamondHit__' in obj:
        diamond_hit = DiamondHit()
        diamond_hit.__dict__.update(obj['__DiamondHit__'])
        return diamond_hit
    return obj


def decode_assembly(obj):
    """Custom JSON decoder for Sample object

    Args:
        obj (obj): GeneAssembly object to decode
    """
    if '__DiamondHitList__' in obj:
        diamond_hit_list = DiamondHitList()
        diamond_hit_list.__dict__.update(obj['__DiamondHitList__'])
        return diamond_hit_list
    elif '__DiamondHit__' in obj:
        diamond_hit = DiamondHit()
        diamond_hit.__dict__.update(obj['__DiamondHit__'])
        return diamond_hit
    elif '__Contig__' in obj:
        contig = Contig()
        contig.__dict__.update(obj['__Contig__'])
        return contig
    elif '__Gene__' in obj:
        gene = Gene()
        gene.__dict__.update(obj['__Gene__'])
        return gene
    elif '__GeneAssembly__' in obj:
        gene_assembly = GeneAssembly()
        gene_assembly.__dict__.update(obj['__GeneAssembly__'])
        return gene_assembly
    return obj


def import_annotated_reads(infile):
    """Imports annotated reads from JSON file

    Args:
        infile (str): JSON file name

    Returns:
        deserialized (dict[str, :obj:AnnotatedRead]): key is read identifier,
            value is annotated read
    """
    deserialized = None
    try:
        with open(infile, 'r') as file_handle:
            deserialized = json.load(file_handle, object_hook=decode_reads)
    except FileNotFoundError:
        deserialized = {}
    return deserialized


def import_sample(infile):
    """Imports sample from JSON file

    Args:
        infile (str): JSON file name

    Returns:
        deserialized (:obj:Sample): sample object
    """
    deserialized = None
    with open(infile, 'r') as file_handle:
        deserialized = json.load(file_handle, object_hook=decode_sample)
    return deserialized


def import_gene_assembly(infile):
    """Imports gene assembly from JSON file

    Args:
        infile (str): JSON file name

    Returns:
        deserialized (:obj:GeneAssembly): gene assembly object
    """
    deserialized = None
    with open(infile, 'r') as file_handle:
        deserialized = json.load(file_handle, object_hook=decode_assembly)
    return deserialized


def export_annotated_reads(parser):
    """Exports annotated reads in JSON format

    Args:
        parser (:obj:DiamondParser): parser with annotated reads
    """
    outfile = os.path.join(
        parser.options.get_project_dir(parser.sample.sample_id),
        parser.sample.sample_id + '_' + parser.end + '_' + parser.options.reads_json_name
    )
    # print pretty JSON: print(json.dumps(parser.reads,indent=4, cls=CustomEncoder))
    with open(outfile, 'w') as out:
        json.dump(parser.reads, out, cls=CustomEncoder)


def export_sample(sample):
    """Exports sample in JSON format

    Args:
        sample (:obj:Sample): sample object
    """
    outfile = os.path.join(sample.work_directory, sample.sample_id + '_sample.json')
    # print pretty JSON: print(json.dumps(parser.reads,indent=4, cls=CustomEncoder))
    with open(outfile, 'w') as out:
        # of.write(json.dumps(sample,indent=4, cls=CustomEncoder))
        json.dump(sample, out, cls=CustomEncoder)


def export_gene_assembly(assembly, outfile):
    """Saves GeneAssembly object as JSON file

    Args:
        assembly (:obj:GeneAssembly): gene assembly
        outfile (str): output file name
    """
    # print pretty JSON: print(json.dumps(assembly,indent=4, cls=CustomEncoder))
    with open(outfile, 'w') as out:
        json.dump(assembly, out, indent=4, cls=CustomEncoder)


class CustomEncoder(json.JSONEncoder):
    """Custom JSON encoder class"""
    def default(self, obj):
        """Returns encoded object as JSON-formatted string

        Args:
            obj (object): instance to encode
        """
        return {'__{}__'.format(obj.__class__.__name__): obj.__dict__}
