#!/usr/bin/python3
"""This script renames a sample in Fama project in four easy steps:
1. Renames existing intermediate files associated with this sample.
2. Deletes output subdirectory in sample's directory
3. Replaces existing sample data in project options with new name and
saves project ini file.
4. Re-creates output files for the sample in the sample's output subdirectory.

Note: do not try to manually edit section names in project ini file. This may lead
to error messages about missing files.

ATTENTION: This script DOES NOT re-creates output files in the project's
working directory to avoid unneccessary actions if multiple samples has
to be renamed. After renaming a sample, files in the project's directory
still have old sample identifier. You must run fama_report.py after all
renamings has been done to re-create reports for project and samples.

"""
import os
import sys
import argparse
import shutil

from lib.utils.const import ENDS
from lib.project.project import Project
from lib.diamond_parser.diamond_parser import DiamondParser

from lib.output.report import generate_fastq_report
from lib.output.pdf_report import generate_pdf_report
from lib.output.krona_xml_writer import make_functions_chart


def get_args():
    """Returns command-line arguments"""
    desc = '''This program renames one sample in a project.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c', dest='config', type=str, help='Path to config.ini')
    parser.add_argument('-p', dest='project', type=str, help='Path to project.ini')
    parser.add_argument('-s', dest='sample', type=str, help='Sample ID ')
    parser.add_argument('-n', dest='name', type=str, help='New sample ID')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args


def check_id(project, old_id, new_id):
    """Checks if old identifier exists in the project and if the new
    identifier does not interfere with any existing sample
    """
    if old_id not in project.list_samples():
        raise ValueError('Sample ID ' + old_id + ' not found')
    if new_id in project.list_samples():
        raise ValueError('Sample ID ' + new_id
                         + ' already exists in this project. Choose other ID')


def rename_files(project, old_id, new_id):
    """Renames all intermediate files having old sample identifier with new sample identifier"""
    for end in ENDS:
        if not project.samples[old_id].is_paired_end and end == 'pe2':
            continue
        if os.path.exists(os.path.join(project.options.get_project_dir(old_id),
                                       old_id + '_' + end + '_' + project.options.ref_output_name)):
            os.rename(
                os.path.join(
                    project.options.get_project_dir(old_id),
                    old_id + '_' + end + '_' + project.options.ref_output_name
                    ),
                os.path.join(
                    project.options.get_project_dir(old_id),
                    new_id + '_' + end + '_' + project.options.ref_output_name
                    )
                )
        if os.path.exists(os.path.join(project.options.get_project_dir(old_id),
                                       old_id + '_' + end + '_'
                                       + project.options.ref_hits_fastq_name)):
            os.rename(
                os.path.join(
                    project.options.get_project_dir(old_id),
                    old_id + '_' + end + '_' + project.options.ref_hits_fastq_name
                    ),
                os.path.join(
                    project.options.get_project_dir(old_id),
                    new_id + '_' + end + '_' + project.options.ref_hits_fastq_name
                    )
                )
        if os.path.exists(os.path.join(project.options.get_project_dir(old_id),
                                       old_id + '_' + end + '_'
                                       + project.options.background_output_name)):
            os.rename(
                os.path.join(
                    project.options.get_project_dir(old_id),
                    old_id + '_' + end + '_' + project.options.background_output_name
                    ),
                os.path.join(
                    project.options.get_project_dir(old_id),
                    new_id + '_' + end + '_' + project.options.background_output_name
                    )
                )
        if os.path.exists(os.path.join(project.options.get_project_dir(old_id),
                                       old_id + '_' + end + '_'
                                       + project.options.reads_fastq_name + '.gz')):
            os.rename(
                os.path.join(
                    project.options.get_project_dir(old_id),
                    old_id + '_' + end + '_' + project.options.reads_fastq_name + '.gz'
                    ),
                os.path.join(
                    project.options.get_project_dir(old_id),
                    new_id + '_' + end + '_' + project.options.reads_fastq_name + '.gz'
                    )
                )
        if os.path.exists(os.path.join(project.options.get_project_dir(old_id),
                                       old_id + '_' + end + '_'
                                       + project.options.pe_reads_fastq_name + '.gz')):
            os.rename(
                os.path.join(
                    project.options.get_project_dir(old_id),
                    old_id + '_' + end + '_' + project.options.pe_reads_fastq_name + '.gz'
                    ),
                os.path.join(
                    project.options.get_project_dir(old_id),
                    new_id + '_' + end + '_' + project.options.pe_reads_fastq_name + '.gz'
                    )
                )
        if os.path.exists(os.path.join(project.options.get_project_dir(old_id),
                                       old_id + '_' + end + '_'
                                       + project.options.ref_hits_list_name)):
            os.rename(
                os.path.join(
                    project.options.get_project_dir(old_id),
                    old_id + '_' + end + '_' + project.options.ref_hits_list_name
                    ),
                os.path.join(
                    project.options.get_project_dir(old_id),
                    new_id + '_' + end + '_' + project.options.ref_hits_list_name
                    )
                )
        if os.path.exists(os.path.join(project.options.get_project_dir(old_id),
                                       old_id + '_' + end + '_'
                                       + project.options.reads_json_name)):
            os.rename(
                os.path.join(
                    project.options.get_project_dir(old_id),
                    old_id + '_' + end + '_' + project.options.reads_json_name
                    ),
                os.path.join(
                    project.options.get_project_dir(old_id),
                    new_id + '_' + end + '_' + project.options.reads_json_name
                    )
                )
    if os.path.exists(os.path.join(project.options.get_project_dir(old_id),
                                   old_id + '_sample.json')):
        os.rename(
            os.path.join(
                project.options.get_project_dir(old_id),
                old_id + '_sample.json'
                ),
            os.path.join(
                project.options.get_project_dir(old_id),
                new_id + '_sample.json'
                )
            )


def rename_sample(config_file, project_file, old_id, new_id):
    """This function is doing all the work to rename a sample

    Args:
        config_file (str): path to program config ini file
        project_file (str): path to project options ini file
        old_id (str): existing sample identifier
        new_id (str): new sample identifier

    """
    project = Project(config_file=config_file, project_file=project_file)
    project.load_project()
    check_id(project, old_id, new_id)

    # Rename files
    rename_files(project, old_id, new_id)

    # Change samples

    new_sample = project.samples[old_id]
    new_sample.sample_id = new_id
    # Delete sample output directory
    shutil.rmtree(
        os.path.join(
            project.options.get_project_dir(old_id),
            project.options.get_output_subdir(old_id)
            )
        )

    items = []
    for option in project.options.parser.options(old_id):
        if option not in project.options.parser.defaults():
            items.append([option, project.options.parser.get(old_id, option)])
    project.options.parser.add_section(new_id)
    for item in items:
        project.options.parser.set(new_id, item[0], item[1])
    project.options.parser.remove_section(old_id)

    project.samples.pop(old_id, None)
    project.samples[new_id] = new_sample
    project.options.set_sample_data(project.samples[new_id])
    project.save_project_options()

    # Re-open project with new version of project file
    project = Project(config_file=config_file, project_file=project_file)
    project.load_project()
    os.mkdir(
        os.path.join(
            project.options.get_project_dir(new_id),
            project.options.get_output_subdir(new_id)
            )
        )
    project.import_reads_json(new_id, ENDS)

    for end in ENDS:
        if not project.samples[new_id].is_paired_end and end == 'pe2':
            continue

        parser = DiamondParser(config=project.config,
                               options=project.options,
                               taxonomy_data=project.taxonomy_data,
                               ref_data=project.ref_data,
                               sample=project.samples[new_id],
                               end=end)
        parser.reads = project.samples[new_id].reads[end]
        # Re-create output files
        # Generate output
        generate_fastq_report(parser)
        generate_pdf_report(parser)
        make_functions_chart(parser)


def main():
    """Main function: gets args and calls rename_sample"""
    args = get_args()
    rename_sample(args.config, args.project, args.sample, args.name)
    print('Run fama_report.py after you finish with renamings!')
    print('Done!')

if __name__ == '__main__':
    main()
