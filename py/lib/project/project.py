"""Describes Project class"""
import os
from collections import defaultdict

from lib.utils.const import ENDS
from lib.project.program_config import ProgramConfig
from lib.project.project_options import ProjectOptions
from lib.project.sample import Sample

from lib.reference_library.reference_data import ReferenceData
from lib.reference_library.taxonomy_data import TaxonomyData

from lib.output.report import generate_project_report
from lib.output.json_util import import_sample, import_annotated_reads


class Project(object):
    """Project is an umbrella object for all samples currently analyzed.

    Attributes:
        samples (:obj:'dict'[str,:obj:'Sample']): dictionary with sample
            identifiers as keys and Sample objects as values
        config (:obj:'ProgramConfig'): Fama configuration parameters
        options  (:obj:'ProjectOptions'): Fama project options
        collection (str): collection identifier
        collection (str): reference collection identifier
        ref_data (:obj:'ReferenceData'): reference dataset for the
            collection (list of functions, list of proteins etc.)
        taxonomy_data (:obj:'TaxonomyData'): NCBI taxonomy dataset for
            the collection
    """
    def __init__(self, config_file, project_file):
        """
        Args:
            config_file (str): full path to program configuration ini file.
            project_file (str): full path to project ini file.
        """
        self.samples = {}
        self.config = ProgramConfig(config_file)
        self.options = ProjectOptions(project_file)
        collection = self.options.get_collection()
        if collection not in self.config.collections:
            raise Exception(
                'Collection ' + collection + ' not found. Available collections are: '
                + (',').join(self.config.collections)
                )
        self.collection = collection
        self.ref_data = ReferenceData(self.config, self.collection)
        # self.ref_data.load_reference_data(self.collection)
        self.taxonomy_data = TaxonomyData(self.config)
        # self.taxonomy_data.load_taxdata(self.config)
        if not os.path.exists(self.options.work_dir) and not os.path.isdir(self.options.work_dir):
            os.makedirs(self.options.work_dir, exist_ok=True)

    def list_samples(self):
        """Returns list of sample identifiers"""
        return self.options.list_samples()

    def save_project_options(self):
        """Saves project options as new version of project ini file"""
        for sample_id in self.samples:
            self.options.set_sample_data(self.samples[sample_id])
        self.options.save_options()

    def load_project(self):
        """Populates reads attribute with samples data"""
        for sample_id in self.list_samples():
            sample = Sample(sample_id=sample_id)
            sample.load_sample(self.options)
            self.samples[sample_id] = sample

    def load_sample(self, sample):
        """Loads sample data from JSON file into memory

        Args:
            sample (:obj:'Sample'): Sample object
        """
        self.samples[sample.sample_id] = \
            import_sample(os.path.join(sample.work_directory,
                                       sample.sample_id + '_' + self.options.reads_json_name))

    def import_reads_json(self, sample_id, ends):
        """Loads annotated reads from one or two JSON files into memory
        Args:
            sample_id (str): sample identifier
            ends (:obj:'list' of str): either ['pe1','pe2'] or ['pe1'] or ['pe2']
        """
        for end_id in ends:
            if end_id == 'pe2' and not self.samples[sample_id].is_paired_end:
                continue
            self.samples[sample_id].reads[end_id] = \
                import_annotated_reads(os.path.join(self.options.get_project_dir(sample_id),
                                                    sample_id + '_' + end_id + '_'
                                                    + self.options.reads_json_name))

    def get_insert_size(self, sample):
        """Returns average insert size for paired-end sample. If calculation of
        insert size is not possible, returns None.
        """
        result = None
        if not sample.is_paired_end or sample.insert_size is None:
            pass
        elif sample.insert_size == 0:
            insert_size = sample.estimate_average_insert_size(
                self.config.get_length_cutoff(self.options.get_collection(sample.sample_id)))
            sample.insert_size = insert_size
            result = insert_size
        elif sample.insert_size > 0:
            result = sample.insert_size
        return result

    def generate_report(self, metrics=None):
        """Writes project report in text format. Also, calls XLSX report
        generation.

        Args:
            metrics (str, optional): metrics for report score calculation
        """
        outfile = os.path.join(self.options.work_dir, 'project_report.txt')
        with open(outfile, 'w') as outfile:
            outfile.write(self.options.project_name + '\n\n')
            for sample_id in self.list_samples():
                outfile.write('\t'.join([sample_id + ':',
                                         self.samples[sample_id].sample_name,
                                         'pe1 reads: '
                                         + str(len(self.samples[sample_id].reads['pe1']))]))
                if self.samples[sample_id].is_paired_end:
                    outfile.write('\tpe2 reads: '
                                  + str(len(self.samples[sample_id].reads['pe2'])))
                outfile.write('\n')
        generate_project_report(self, metrics)

    def check_project(self):
        """Checks if all files and directories of a project do exist.

        Todo:
            Ensure that it looks up the last version of file names
        """
        problems = defaultdict(list)
        print('Checking project', self.options.project_name)
        for sample in self.list_samples():
            print('Checking sample', sample)
            for end in ENDS:
                skip_output_check = False
                if not os.path.exists(self.options.get_fastq_path(sample, end)):
                    problems[sample].append('Input FASTQ file not found for sample',
                                            sample, ',end', end, ':',
                                            self.options.get_fastq_path(sample, end))
                outdir = self.options.get_project_dir(sample)
                if not os.path.isdir(outdir):
                    problems[sample].append('Directory not found for sample', sample, ':',
                                            outdir, 'CHECK INTERRUPTED')
                    continue
                if not os.path.isdir(os.path.join(outdir, self.options.get_output_subdir(sample))):
                    problems[sample].append('Output directory not found for sample', sample, ':',
                                            os.path.join(outdir,
                                                         self.options.get_output_subdir(sample)),
                                            'OUTPUT FILES NOT CHECKED')
                    skip_output_check = True
                if not os.path.exists(os.path.join(outdir,
                                                   sample + '_' + end + '_'
                                                   + self.options.ref_output_name)):
                    problems[sample].append('Reference DB search output not found for sample '
                                            + sample + ', end ' + end + ':' +
                                            os.path.join(outdir,
                                                         sample + '_' + end + '_'
                                                         + self.options.ref_output_name))
                if not os.path.exists(os.path.join(outdir,
                                                   sample + '_' + end + '_'
                                                   + self.options.background_output_name)):
                    problems[sample].append('Background DB search output not found \
                                            for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir,
                                                         sample + '_' + end + '_'
                                                         + self.options.background_output_name))
                if not os.path.exists(os.path.join(outdir,
                                                   sample + '_' + end + '_'
                                                   + self.options.reads_fastq_name + '.gz')):
                    problems[sample].append('Output FASTQ file with reads not found for sample '
                                            + sample + ', end ' + end + ':' +
                                            os.path.join(outdir,
                                                         sample + '_' + end + '_'
                                                         + self.options.reads_fastq_name + '.gz'))
                if not os.path.exists(os.path.join(outdir,
                                                   sample + '_' + end + '_'
                                                   + self.options.ref_hits_fastq_name)):
                    problems[sample].append('Reference hits FASTQ file not found for sample '
                                            + sample + ', end ' + end + ':' +
                                            os.path.join(outdir,
                                                         sample + '_' + end + '_'
                                                         + self.options.ref_hits_fastq_name))
                if not os.path.exists(os.path.join(outdir,
                                                   sample + '_' + end + '_'
                                                   + self.options.ref_hits_list_name)):
                    problems[sample].append('List of reference hits not found for sample '
                                            + sample + ', end ' + end + ':' +
                                            os.path.join(outdir,
                                                         sample + '_' + end + '_'
                                                         + self.options.ref_hits_list_name))
                if not os.path.exists(os.path.join(outdir,
                                                   sample + '_' + end + '_'
                                                   + self.options.pe_reads_fastq_name + '.gz')):
                    problems[sample].append('Output FASTQ file with paired-ends not found \
                                            for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir,
                                                         sample + '_' + end + '_'
                                                         + self.options.pe_reads_fastq_name
                                                         + '.gz'))
                if not os.path.exists(os.path.join(outdir,
                                                   sample + '_' + end + '_'
                                                   + self.options.reads_json_name)):
                    problems[sample].append('Output JSON file with annotated reads not found \
                        for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir,
                                                         sample + '_' + end + '_'
                                                         + self.options.reads_json_name))
                if skip_output_check:
                    continue
                if not os.path.exists(os.path.join(outdir,
                                                   self.options.get_output_subdir(sample),
                                                   sample + '_' + end + '_'
                                                   + self.options.report_name)):
                    problems[sample].append('Text report file not found for sample ' + sample
                                            + ', end ' + end + ':' +
                                            os.path.join(outdir,
                                                         self.options.get_output_subdir(sample),
                                                         sample + '_' + end + '_'
                                                         + self.options.report_name))
                if not os.path.exists(os.path.join(outdir,
                                                   self.options.get_output_subdir(sample),
                                                   sample + '_' + end + '_'
                                                   + self.options.xml_name)):
                    problems[sample].append('Krona XML file for functional profile \
                                            not found for sample ' + sample + ', end '
                                            + end + ':' +
                                            os.path.join(outdir,
                                                         self.options.get_output_subdir(sample),
                                                         sample + '_' + end + '_'
                                                         + self.options.xml_name))
                if not os.path.exists(os.path.join(outdir,
                                                   self.options.get_output_subdir(sample),
                                                   sample + '_' + end + '_'
                                                   + self.options.html_name)):
                    problems[sample].append('HTML file for functional profile not found for sample '
                                            + sample + ', end ' + end + ':' +
                                            os.path.join(outdir,
                                                         self.options.get_output_subdir(sample),
                                                         sample + '_' + end + '_'
                                                         + self.options.html_name))
        if not problems:
            print('No problems found in your project. Could be worse.')
        else:
            print('Problems found:')
            for sample in problems:
                print('********* ' + sample + ' *********')
                print('\n'.join(problems[sample]))
                print('*********************************\n\n')
