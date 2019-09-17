"""Fama assembly pipeline"""
import os
from lib.project.project import Project
from lib.gene_assembler.gene_assembler import GeneAssembler
from lib.output.json_util import export_gene_assembly


def assembly_pipeline(args):
    """Runs steps of assembly pipeline

    Args:
        args: ArgumentParser namespace with defined args.config (path to
            program config ini file) and args.project (path to project
            options ini file)
    """
    project = Project(config_file=args.config, project_file=args.project)
    assembler = GeneAssembler(project, assembler=args.assembler)
    assembler.export_reads(do_coassembly=args.coassembly)
    assembler.assemble_contigs()
    assembler.predict_genes()
    assembler.annotate_genes()
    assembler.generate_output()
    export_gene_assembly(
        assembler.assembly,
        os.path.join(project.options.assembly_dir, 'all_contigs_assembly.json')
        )
