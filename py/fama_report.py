#!/usr/bin/python3
"""This script re-creates Fama output files for existing read annotations.

Note 1: Project must have ouput JSON files in samples' directories.

Note 2: Reports may be created for all samples in the project or only for one sample.
For one sample, project reports (comparative tables etc.) will not be
re-created.

Note 3: Use this script to create new project reports with non-default metrics.

"""
import sys
import argparse
from lib.utils.const import ENDS
from lib.project.project import Project
from lib.functional_profiling_pipeline import run_microbecensus
import lib.output.report as fama_report
from lib.protein_functional_pipeline import generate_output


def get_args():
    """Returns command-line arguments"""
    desc = '''This program re-runs report generation for functional profiling pipeline.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', dest='config', type=str, help='Path to config.ini')
    parser.add_argument('-p', dest='project', type=str, help='Path to project.ini')
    parser.add_argument('-s', dest='sample', type=str, default=None,
                        help='Sample ID (optional)')
    parser.add_argument('-m', dest='metrics', type=str, default='efpkg',
                        help='Metrics (optional). Acceptable values: readcount, \
                        fragmentcount, fpk, fpkm, fpkg, erpk, erpkm, erpkg, \
                        efpk, efpkm, efpkg. Default values: efpkg for \
                        paired-end projects, erkpg for single-end projects.')
    parser.add_argument('--prot', dest='prot', action='store_true',
                        help='Process protein sequences in FASTA format (default: False)')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args


def main():
    """Loads annotated reads and calls report generators"""
    args = get_args()
    if args.prot:
        project = Project(config_file=args.config, project_file=args.project)
        project.load_project()
        for sample_id in project.list_samples():
            if args.sample is not None:
                if args.sample != sample_id:
                    continue
            project.options.set_sample_data(project.samples[sample_id])
            project.import_reads_json(sample_id, ENDS)
            fama_report.generate_protein_sample_report(project, sample_id, metrics='readcount')

        if args.sample is None:
            # Skip project report generation if the pipeline is running for only one sample
            fama_report.generate_protein_project_report(project)
        generate_output(project)

    else:
        project = Project(config_file=args.config, project_file=args.project)
        project.load_project()
        for sample_id in project.list_samples():
            if args.sample is not None and args.sample != sample_id:
                continue
            if project.samples[sample_id].rpkg_scaling_factor is None or project.samples[
                    sample_id
            ].rpkg_scaling_factor == 0.0:
                project.samples[sample_id].import_rpkg_scaling_factor()
                if project.samples[
                        sample_id
                    ].rpkg_scaling_factor is None or project.samples[
                        sample_id
                        ].rpkg_scaling_factor == 0.0:
                    run_microbecensus(
                        sample=project.samples[sample_id],
                        config=project.config
                        )
                    project.samples[sample_id].import_rpkg_scaling_factor()
            project.options.set_sample_data(project.samples[sample_id])
            project.import_reads_json(sample_id, ENDS)
            print('Generating report for', sample_id)
            fama_report.generate_sample_report(project, sample_id, metrics=args.metrics)
        print('Generating report for project')
        project.generate_report(metrics=args.metrics)
        project.save_project_options()
        print('Done!')

if __name__ == '__main__':
    main()
