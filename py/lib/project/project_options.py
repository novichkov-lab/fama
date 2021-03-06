"""Descibes ProjectOptions class"""
import os
import configparser
from lib.utils.utils import singleton


@singleton
class ProjectOptions(object):
    """Project parameters for Fama run. Wraps a ConfigParser instance,
    which reads project ini file.

    Attributes:
        project_file (str): path to project ini file
        parser (:obj:ConfigParser): underlying ConfigParser instance
    """
    def __init__(self, project_file):
        """
        Args:
            project_file (str): path to project ini file
        """
        self.project_file = project_file
        self.parser = configparser.ConfigParser()
        self.load_project(project_file)

    # Properties
    @property
    def project_name(self):
        """Text label for project name"""
        return self.parser['DEFAULT']['project_name']

    @property
    def work_dir(self):
        """Path to project's working directory"""
        return self.parser['DEFAULT']['work_dir']

    @property
    def ref_output_name(self):
        """Ending of file name for reference DB search results"""
        return self.parser['DEFAULT']['ref_output_name']

    @property
    def background_output_name(self):
        """Ending of file name for background DB search results"""
        return self.parser['DEFAULT']['background_output_name']

    @property
    def ref_hits_list_name(self):
        """Ending of file name for list of hits"""
        return self.parser['DEFAULT']['ref_hits_list_name']

    @property
    def ref_hits_fastq_name(self):
        """Ending of file name for FASTQ file of hits"""
        return self.parser['DEFAULT']['ref_hits_fastq_name']

    @property
    def reads_fastq_name(self):
        """Ending of file name for FASTQ file of selected reads"""
        return self.parser['DEFAULT']['reads_fastq_name']

    @property
    def pe_reads_fastq_name(self):
        """Ending of file name for FASTQ file of paired-ends of selected reads"""
        return self.parser['DEFAULT']['pe_reads_fastq_name']

    @property
    def report_name(self):
        """Ending of file name for report"""
        return self.parser['DEFAULT']['report_name']

    @property
    def xml_name(self):
        """Ending of file name for Krona XML"""
        return self.parser['DEFAULT']['xml_name']

    @property
    def html_name(self):
        """Ending of file name for Krona HTML"""
        return self.parser['DEFAULT']['html_name']

    @property
    def reads_json_name(self):
        """Ending of JSON file name for results"""
        return self.parser['DEFAULT']['reads_json_name']

    @property
    def assembly_dir(self):
        """Name of subdirectory for assembly"""
        return os.path.join(self.work_dir, self.parser['DEFAULT']['assembly_subdir'])

    # Methods
    def load_project(self, project_file):
        """Loads(reloads) project ini file

        Args:
            project_file (str): path to project ini file
        """
        self.parser.read(project_file)

    def get_collection(self, sample=None):
        """Returns collection identifier for a sample or for the project

        Args:
            sample (str, optional): sample identifier
        """
        result = self.parser['DEFAULT']['collection']
        if sample and self.parser.has_option(sample, 'collection'):
            result = self.parser[sample]['collection']
        return result

    def get_sample_name(self, sample):
        """Returns text label for a sample

        Args:
            sample (str): sample identifier
        """
        return self.parser[sample]['sample_id']

    def list_samples(self):
        """Returns list of sample identifiers"""
        return self.parser.sections()

    def get_fastq1_readcount(self, sample):
        """Returns number of reads in FASTQ file 1 of a sample

        Args:
            sample (str): sample identifier
        """
        return self.parser.getint(sample, 'fastq_pe1_readcount', fallback=0)

    def get_fastq2_readcount(self, sample):
        """Returns number of reads in FASTQ file 2 of a sample

        Args:
            sample (str): sample identifier
        """
        return self.parser.getint(sample, 'fastq_pe2_readcount', fallback=0)

    def get_fastq1_basecount(self, sample):
        """Returns number of bases in FASTQ file 1 of a sample

        Args:
            sample (str): sample identifier
        """
        return self.parser.getint(sample, 'fastq_pe1_basecount', fallback=0)

    def get_fastq2_basecount(self, sample):
        """Returns number of bases in FASTQ file 2 of a sample

        Args:
            sample (str): sample identifier
        """
        return self.parser.getint(sample, 'fastq_pe2_basecount', fallback=0)

    def get_project_dir(self, sample):
        """Returns path to sample's working directory. If the name is
        not set (not recommended), returns default path to working subdirectory.

        Args:
            sample (str): sample identifier
        """
        result = self.parser['DEFAULT']['work_dir']
        if self.parser.has_option(sample, 'sample_dir'):
            result = self.parser[sample]['sample_dir']
        return result

    def get_output_subdir(self, sample):
        """Returns name of sample's output subdirectory. If the name is
        not set, returns default name for output subdirectory.

        Args:
            sample (str): sample identifier
        """
        result = self.parser['DEFAULT']['output_subdir']
        if self.parser.has_option(sample, 'output_subdir'):
            return self.parser[sample]['output_subdir']
        return result

    def get_fastq_path(self, sample, end):
        """Returns path to input FASTQ(FASTA) file for sample and end identifiers.
        If the path is not set, returns None.

        Args:
            sample (str): sample identifier
            end (str): end identifier
        """
        result = None
        if end == 'pe1':
            if self.parser.has_option(sample, 'fastq_pe1'):
                result = self.parser[sample]['fastq_pe1']
        elif end == 'pe2':
            if self.parser.has_option(sample, 'fastq_pe2'):
                result = self.parser[sample]['fastq_pe2']
        return result

    def get_coverage_path(self, sample):
        """Returns path to file with contig coverage data (for protein
        input files). If the path is not set, returns default path to
        contig coverage data file.
        If the default path is not set, returns None.

        Args:
            sample (str): sample identifier
        """
        result = None
        if self.parser.has_option(sample, 'coverage'):
            result = self.parser[sample]['coverage']
        elif self.parser.has_option('', 'coverage'):
            result = self.parser['DEFAULT']['coverage']
        return result

    def get_rpkg_scaling_factor(self, sample):
        """Returns normalization coefficient for RPKG normalization
        (normalization by sample size and average genome size). If the
        coefficient is not set, returns None.

        Args:
            sample (str): sample identifier
        """
        return self.parser.getfloat(sample, 'rpkg_scaling', fallback=0.0)

    def set_sample_data(self, sample):
        """Writes parameters of sample from Sample object to the corresponding
        section of underlying ConfigParser. Call it for each sample
        before calling save_options()

        Args:
            sample (:obj:Sample): sample object
        """
        self.parser[sample.sample_id]['sample_id'] = sample.sample_name
        self.parser[sample.sample_id]['fastq_pe1'] = sample.fastq_fwd_path
        self.parser[sample.sample_id]['fastq_pe1_readcount'] = str(sample.fastq_fwd_readcount)
        self.parser[sample.sample_id]['fastq_pe1_basecount'] = str(sample.fastq_fwd_basecount)
        if sample.is_paired_end:
            self.parser[sample.sample_id]['fastq_pe2'] = sample.fastq_rev_path
            self.parser[sample.sample_id]['fastq_pe2_readcount'] = str(sample.fastq_rev_readcount)
            self.parser[sample.sample_id]['fastq_pe2_basecount'] = str(sample.fastq_rev_basecount)
        self.parser[sample.sample_id]['sample_dir'] = sample.work_directory
        if sample.rpkg_scaling_factor is not None:
            self.parser[sample.sample_id]['rpkg_scaling'] = str(sample.rpkg_scaling_factor)
        self.parser[sample.sample_id]['replicate'] = sample.replicate
        if sample.insert_size is not None:
            self.parser[sample.sample_id]['insert_size'] = str(sample.insert_size)

    def save_options(self):
        """Writes previous version of project ini file under different name.
        Saves current state of ProjectOptions object as project ini file.
        """
        i = 0
        while True:
            if os.path.exists(self.project_file + '.old.' + str(i)):
                i += 1
            else:
                backup_file = self.project_file + '.old.' + str(i)
                break
        os.rename(self.project_file, backup_file)
        with open(self.project_file, 'w') as outfile:
            self.parser.write(outfile)
        print('Old version of the project file copied to', backup_file)
