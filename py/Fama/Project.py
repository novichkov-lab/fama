import os
from collections import defaultdict
from subprocess import Popen, PIPE, CalledProcessError

from Fama.ProjectUtil.ProgramConfig import ProgramConfig
from Fama.ProjectUtil.ProjectOptions import ProjectOptions
from Fama.Sample import Sample

from Fama.ReferenceLibrary.ReferenceData import ReferenceData
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.GeneAssembler.GeneAssembler import GeneAssembler

from Fama.OutputUtil.Report import generate_project_report
from Fama.OutputUtil.PdfReport import generate_pdf_report
from Fama.OutputUtil.KronaXMLWriter import generate_functions_chart
from Fama.OutputUtil.JSONUtil import export_annotated_reads
from Fama.OutputUtil.JSONUtil import import_sample
from Fama.OutputUtil.JSONUtil import import_annotated_reads


class Project(object):
    
    def __init__(self, config_file, project_file):
        self.ENDS = ['pe1','pe2']
        self.samples = {}
        self.config = ProgramConfig(config_file)
        self.options = ProjectOptions(project_file)
        collection = self.options.get_collection()
        if collection not in self.config.list_collections():
            raise Exception ('Collection ' + collection + ' not found. Available collections are: ' + (',').join(self.config.list_collections()))
        self.collection = collection
        self.ref_data = ReferenceData(self.config)
        self.ref_data.load_reference_data(self.collection)
        self.taxonomy_data = TaxonomyData(self.config)
        self.taxonomy_data.load_taxdata(self.config)
        if not os.path.exists(self.options.get_work_dir()) and not os.path.isdir(self.options.get_work_dir()):
            os.makedirs(self.options.get_work_dir(), exist_ok=True)

    def list_samples(self):
        return self.options.list_samples()

    def save_project_options(self):
        for sample_id in self.samples:
            self.options.set_sample_data(self.samples[sample_id])
        self.options.save_options()

    def load_project(self):
        for sample_id in self.list_samples():
            sample = Sample(sample_id=sample_id)
            sample.load_sample(self.options)
            self.samples[sample_id] = sample

    def load_sample(self, sample):
        self.samples[sample.sample_id] = import_sample(os.path.join(sample.work_directory, sample.sample_id + '_' + self.options.get_reads_json_name()))

    def import_reads_json(self, sample_id, ends):
        for end_id in ends:
            if end_id == 'pe2' and not self.samples[sample_id].is_paired_end:
                continue
            self.samples[sample_id].reads[end_id] = import_annotated_reads(os.path.join(self.options.get_project_dir(sample_id), sample_id + '_' + end_id + '_' + self.options.get_reads_json_name()))

    def get_insert_size(self, sample):
        if not sample.is_paired_end or sample.insert_size is None:
            return None
        elif sample.insert_size == 0:
            insert_size = sample.estimate_average_insert_size(self.config.get_length_cutoff(self.options.get_collection(sample.sample_id)))
            sample.insert_size = insert_size
            return insert_size
        elif sample.insert_size > 0:
            return sample.insert_size
        else: # that was not expected
            return None

    def generate_report(self, metrics = None):
        outfile = os.path.join(self.options.get_work_dir(), 'project_report.txt')
        with open (outfile, 'w') as of:
            of.write(self.options.get_name() + '\n\n')
            for sample_id in self.list_samples():
                of.write(sample_id + ':\t' + self.samples[sample_id].sample_name + '\tpe1 reads: ' + str(len(self.samples[sample_id].reads['pe1'])))
                if self.samples[sample_id].is_paired_end:
                    of.write('\tpe2 reads: ' + str(len(self.samples[sample_id].reads['pe2'])))
                of.write('\n')
            of.closed
        generate_project_report(self, metrics)

    def check_project(self):
        problems = defaultdict(list)
        print ('Checking project', self.options.get_name())
        for sample in self.list_samples():
            print ('Checking sample', sample)
            for end in self.ENDS:
                skip_output_check = False
                if not os.path.exists(self.options.get_fastq_path(sample,end)):
                    problems[sample].append('Input FASTQ file not found for sample', sample, ',end', end, ':', 
                                            self.options.get_fastq_path(sample,end))
                outdir = self.options.get_project_dir(sample)
                if not os.path.isdir(outdir):
                    problems[sample].append('Directory not found for sample', sample, ':', 
                                            outdir, 'CHECK INTERRUPTED')
                    continue
                if not os.path.isdir(os.path.join(outdir,self.options.get_output_subdir(sample))):
                    problems[sample].append('Output directory not found for sample', sample, ':', 
                                            os.path.join(outdir,self.options.get_output_subdir(sample)), 'OUTPUT FILES NOT CHECKED')
                    skip_output_check = True
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_'+ self.options.get_ref_output_name())):
                    problems[sample].append('Reference DB search output not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_'+ self.options.get_ref_output_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_'+ self.options.get_background_output_name())):
                    problems[sample].append('Background DB search output not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_'+ self.options.get_background_output_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_' + self.options.get_reads_fastq_name())):
                    problems[sample].append('Output FASTQ file with reads not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_' + self.options.get_reads_fastq_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_' + self.options.get_ref_hits_fastq_name())):
                    problems[sample].append('Reference hits FASTQ file not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_' + self.options.get_ref_hits_fastq_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_' + self.options.get_ref_hits_list_name())):
                    problems[sample].append('List of reference hits not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_' + self.options.get_ref_hits_list_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_' + self.options.get_pe_reads_fastq_name())):
                    problems[sample].append('Output FASTQ file with paired-ends not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_' + self.options.get_pe_reads_fastq_name()))
                if not os.path.exists(os.path.join(outdir, sample + '_' + end + '_' + self.options.get_reads_json_name())):
                    problems[sample].append('Output JSON file with annotated reads not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, sample + '_' + end + '_' + self.options.get_reads_json_name()))
                if skip_output_check:
                    continue
                if not os.path.exists(os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_report_name())):
                    problems[sample].append('Text report file not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_report_name()))
                if not os.path.exists(os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_report_name() + '.pdf')):
                    problems[sample].append('PDF report file not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_report_name() + '.pdf'))
                if not os.path.exists(os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_xml_name())):
                    problems[sample].append('Krona XML file for functional profile not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_xml_name()))
                if not os.path.exists(os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_html_name())):
                    problems[sample].append('HTML file for functional profile not found for sample ' + sample + ', end ' + end + ':' +
                                            os.path.join(outdir, self.options.get_output_subdir(sample),sample + '_' + end + '_'+ self.options.get_html_name()))
        if not problems:
            print('No problems found in your project. Could be worse.')
        else:
            print('Problems found:')
            for sample in problems:
                print ('********* ' + sample + ' *********')
                print ('\n'.join(problems[sample]))
                print ('*********************************\n\n')

        
    
