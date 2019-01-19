#!/usr/bin/python3
import os, csv, operator
import unittest
from collections import Counter
from fpdf import FPDF

from context import Fama
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.Project import Project
from Fama.OutputUtil.PdfReport import generate_pdf_report
from Fama.OutputUtil.KronaXMLWriter import generate_functions_chart
from Fama.OutputUtil.PdfReport import generate_pdf_report
from Fama.OutputUtil.XlsxUtil import create_functions_xlsx
from Fama.OutputUtil.Report import create_functions_markdown_document


data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_nitrogen9.ini')
sample_id = 'test_sample'
#sample_id = 'sample1'
end = 'pe1'

class FamaReportTest(unittest.TestCase):

    def setUp(self):
        self.project = Project(config_file=config_path, project_file=project_path)
        self.project.load_project()
        self.parser = DiamondParser(config = self.project.config, 
                            options=self.project.options, 
                            taxonomy_data=self.project.taxonomy_data,
                            ref_data=self.project.ref_data,
                            sample=self.project.samples[sample_id], 
                            end=end)


#    @unittest.skip("for faster testing")
    def test_1_collection_pdf_output(self):
        outfile = os.path.join(data_dir,'collection_list.pdf')
        if (os.path.exists(outfile)):
            os.remove(outfile)
        pdf = FPDF('P', 'mm', 'Letter')
        pdf.add_page()
        pdf.add_font('DejaVu', '', '/usr/share/fonts/truetype/dejavu/DejaVuSansCondensed.ttf', uni=True)
        pdf.add_font('DejaVu', 'B', '/usr/share/fonts/truetype/dejavu/DejaVuSansCondensed-Bold.ttf', uni=True)
        pdf.set_font('DejaVu', 'B', 16)
        pdf.cell(40, 10, 'List of functions in collection ' + self.parser.collection)
        pdf.ln(h = '')
        pdf.set_font('DejaVu', '', 10)
        
        v_limit = 40
        for function in self.parser.ref_data.functions_dict:
            v_limit += 10
            pdf.cell(v_limit, 10, function + '  ' + self.parser.ref_data.functions_dict[function]['name'] 
                    + ': group ' + self.parser.ref_data.functions_dict[function]['group'] )
            pdf.ln(h = 10)

        outfile = os.path.join(data_dir,'collection_list.pdf')
        pdf.output(outfile, 'F')
        self.assertTrue(os.path.exists(outfile))
        # if this test fails, function names may contain 'bad' symbols 

#    @unittest.skip("for faster testing")
    def test_2_get_functions_in_group(self):
        urease_list = self.parser.ref_data.get_functions_in_group('Urease')
        self.assertEqual(len(urease_list), 3)
        self.assertEqual(sorted(urease_list)[0], 'UreA')
        
#    @unittest.skip("for faster testing")
    def test_3_generate_pdf_report(self):
        self.parser.parse_background_output()
        generate_pdf_report(self.parser)

#    @unittest.skip("for faster testing")
    def test_4_generate_functions_chart(self):
        self.parser.parse_background_output()
        generate_functions_chart(self.parser)

#    @unittest.skip("for faster testing")
    def test_5_generate_functions_xlsx(self):
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        create_functions_xlsx(self.project, 'Read count')
        #create_functions_xlsx(self.project, 'RPKM')

#    @unittest.skip("for faster testing")
    def test_6_generate_functions_markdown(self):
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        create_functions_markdown_document(self.project)

#    @unittest.skip("for faster testing")
    def test_7_generate_function_table(self):
        function = 'AmoA'
        for sample_id in self.project.list_samples():
            self.project.import_reads_json(sample_id, self.project.ENDS)
        outfile = 'test_output_20181130.txt'
        with open(outfile, 'w') as of:
            for read_id,read in self.project.samples[sample_id].reads[end].items():
#                if function in read.get_functions():
                for hit in read.get_hit_list().get_hits():
                    if function in hit.get_functions():
                        print (read_id, function, read.status, hit.get_subject_id())
                        of.write(('\t').join([read_id, function, read.status, hit.get_subject_id()]) + '\n')
            of.closed
            

    def tearDown(self):
        self.parser = None

if __name__=='__main__':
    unittest.main()
