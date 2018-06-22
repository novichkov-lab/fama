#!/usr/bin/python
import os, csv, operator
import unittest
from collections import Counter
from fpdf import FPDF

from context import Fama
from Fama.DiamondParser.DiamondParser import DiamondParser
from Fama.OutputUtil.PdfReport import generate_pdf_report
from Fama.OutputUtil.KronaXMLWriter import generate_functions_chart



data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
sample = 'test_sample'
end = 'pe1'

class DiamondParserTest(unittest.TestCase):

    def setUp(self):
        self.parser = DiamondParser(config_file=config_path, project_file=project_path, sample=sample, end=end)

    def test_1_collection_pdf_output(self):
        outfile = os.path.join(data_dir,'collection_list.pdf')
        if (os.path.exists(outfile)):
            os.remove(outfile)
        pdf = FPDF('P', 'mm', 'Letter')
        pdf.add_page()
        pdf.set_font('Arial', 'B', 16)
        pdf.cell(40, 10, 'List of functions in collection ' + self.parser.collection)
        pdf.ln(h = '')
        pdf.set_font('Arial', '', 10)
        
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

    @unittest.skip("for faster testing")
    def test_3_generate_pdf_report(self):
        self.parser.parse_background_output()
        generate_pdf_report(self.parser)

    def test_2_get_functions_in_group(self):
        urease_list = self.parser.ref_data.get_functions_in_group('Urease')
        self.assertEqual(len(urease_list), 4)
        self.assertEqual(sorted(urease_list)[0], 'UreA')
        
    def test_4_generate_functions_chart(self):
        self.parser.parse_background_output()
        generate_functions_chart(self.parser)

    def tearDown(self):
        self.parser = None

if __name__=='__main__':
    unittest.main()
