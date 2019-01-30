import os,sys
from Fama.Project import Project
from Fama.GeneAssembler.GeneAssembler import GeneAssembler
from Fama.OutputUtil.JSONUtil import export_gene_assembly


def assembly_pipeline(args):
    project = Project(config_file=args.config, project_file=args.project)
    assembler = GeneAssembler(project)
    if args.coassembly:
        assembler.coassemble_contigs()
    else:
        assembler.assemble_contigs()
    assembler.map_genes()
    assembler.map_genes2uniprot()
    assembler.generate_output()
    export_gene_assembly(assembler.assembly,os.path.join(project.options.get_assembly_dir(), 'all_contigs_assembly.json'))

