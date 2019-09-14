""" Utility functions"""
from collections import defaultdict

def autovivify(levels=1, final=dict):
    """Creates multi-level dictionary based on defaultdict

    Args:
        levels (int): number of levels
        final (type): type of innermost level
    """
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))

def cleanup_protein_id(protein):
    """ This function was added for back compatibility with early versions
    of Fama reference datasets.
    For compatibility with old format of protein IDs uncomment next 4 lines

    Args:
        protein (str): protein identifier

    Todo:
        remove
    """
    #if len(protein.split('_')) > 1:
    #    return "_".join(protein.split('_')[1:])
    #else:
    #    return protein
    return protein

def sanitize_file_name(filename):
    """ Replaces unsafe symbols in filenames

    Args:
        filename (str): file name
    """
    filename = filename.replace(' ', '_')
    filename = filename.replace("'", "")
    filename = filename.replace('"', '')
    return filename
