#!/usr/bin/python3
import os, csv, operator
import unittest
import json
from collections import Counter

from context import lib

from lib.diamond_parser.diamond_hit import DiamondHit
from lib.diamond_parser.diamond_hit_list import DiamondHitList
from lib.diamond_parser.diamond_parser import DiamondParser

from lib.sequences.annotated_read import AnnotatedRead
from lib.output.json_util import import_annotated_reads,export_annotated_reads

from lib.protein_functional_pipeline import functional_profiling_pipeline


data_dir = 'data'
config_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'config_test.ini')
project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'projects', 'protein_project.ini')
sample = 'test_sample1'
end = 'pe1'

class DiamondParserTest(unittest.TestCase):

    def setUp(self):
        self.parser = DiamondParser(config_file=config_path, project_file=project_path, sample=sample, end=end)
        
    def test_1_test_pipeline(self):
        #parser = functional_profiling_pipeline(config_file=config_path, project_file=project_path, sample=sample, end=end)
        pass
        
    def tearDown(self):
        self.parser = None


if __name__=='__main__':
    unittest.main()
