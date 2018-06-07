#!/usr/bin/python
import os, csv, operator
import unittest
from context import lib
from collections import Counter

from lib.DiamondParser.DiamondParser import DiamondParser
from lib.OutputUtil.PdfReport import generate_pdf_report


data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project.ini')
sample = 'test_sample'
end = 'pe1'

class DiamondParserTest(unittest.TestCase):

    def setUp(self):
        self.parser = DiamondParser(config_path, project_path, sample, end)

    def test_1_generate_pdf_report(self):
        self.parser.parse_background_output()
        generate_pdf_report(self.parser)

    def tearDown(self):
        self.parser = None

if __name__=='__main__':
    unittest.main()
