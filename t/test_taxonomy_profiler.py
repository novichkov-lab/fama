#!/usr/bin/python3
import os
import unittest
from collections import Counter, defaultdict

from context import lib

from lib.utils.const import ENDS, STATUS_GOOD, ROOT_TAXONOMY_ID
from lib.utils.utils import autovivify, cleanup_protein_id, sanitize_file_name
from lib.diamond_parser.diamond_parser import DiamondParser
from lib.output.json_util import import_annotated_reads
from lib.output.report import generate_fastq_report, get_function_taxonomy_scores, get_function_scores, slice_function_taxonomy_scores
from lib.taxonomy.tree import Node, Tree
from lib.taxonomy.taxonomy_profile import TaxonomyProfile
from lib.output.krona_xml_writer import make_taxonomy_chart, make_taxonomy_series_chart, make_function_taxonomy_chart
from lib.project.project import Project
from lib.reference_library.taxonomy_data import TaxonomyData

data_dir = 'data'
config_path = os.path.join(
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    'py',
    'config.ini'
)
#~ config_path = os.path.join(
    #~ os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    #~ 'py',
    #~ 'config_rpL6_singleDB.ini'
#~ )
project_path = os.path.join(
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
    'py',
    'project.ini'
)
sample_id = 'test_sample'
end = 'pe1'


#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_Hans_sulfate_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_EB271_nitrogen_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW306_cazy_t.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'py', 'project_FW3062M_universal1.ini')
#project_path = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')), 'p', 'project_FW3062M_rpl6.ini')
#sample_id = 'FW306-4701'

class TaxonomyProfilingTest(unittest.TestCase):

    def setUp(self):
        self.project = Project(config_file=config_path, project_file=project_path)
        self.project.load_project()
        self.parser = DiamondParser(config = self.project.config, 
                            options=self.project.options, 
                            taxonomy_data=self.project.taxonomy_data,
                            ref_data=self.project.ref_data,
                            sample=self.project.samples[sample_id], 
                            end=end)
        self.taxonomy_data = TaxonomyData(self.project.config, '')
        self.node = Node(rank='genus', name='Escherichia', taxid='561', parent='543', children=None)
        self.tree = Tree()

# Tests of TaxonomyData class
    def test_0010_init_taxonomy_data(self):
        self.assertEqual(len(self.taxonomy_data.data), 6568)  # Nitrogen collection

    def test_0020_is_exist(self):
        self.assertTrue(self.taxonomy_data.is_exist('1'))
        self.assertTrue(self.taxonomy_data.is_exist('0'))
        self.assertTrue(self.taxonomy_data.is_exist('562'))
        self.assertFalse(self.taxonomy_data.is_exist('fake_identifier'))

    def test_0030_get_name(self):
        self.assertEqual(self.taxonomy_data.get_name('1'), 'root')
        self.assertEqual(self.taxonomy_data.get_name('0'), 'Unknown')
        self.assertEqual(self.taxonomy_data.get_name('562'), 'Escherichia coli')

    def test_0040_get_rank(self):
        self.assertEqual(self.taxonomy_data.get_rank('1'), 'norank')  # Test root
        self.assertEqual(self.taxonomy_data.get_rank('0'), 'norank')  # Test Unknown
        self.assertEqual(self.taxonomy_data.get_rank('2'), 'superkingdom')  # Test Bacteria
        self.assertEqual(self.taxonomy_data.get_rank('562'), 'species')  # Test E. coli

    def test_0050_get_parent(self):
        self.assertEqual(self.taxonomy_data.get_parent('1'), '1')  # Test root
        self.assertEqual(self.taxonomy_data.get_parent('0'), '1')  # Test Unknown
        self.assertEqual(self.taxonomy_data.get_parent('2'), '131567')  # Test Bacteria
        self.assertEqual(self.taxonomy_data.get_parent('562'), '561')  # Test E. coli

    def test_0060_get_lca(self):
        self.assertEqual(self.taxonomy_data.get_lca(['1']), '1')  # Test root
        self.assertEqual(self.taxonomy_data.get_lca(['']), '0')  # Test empty string
        self.assertEqual(self.taxonomy_data.get_lca([]), '0')  # Test empty list
        self.assertEqual(self.taxonomy_data.get_lca(['0']), '0')  # Test Unknown
        # Anything with root goes to Unknown
        self.assertEqual(self.taxonomy_data.get_lca(['571','1']), '0')  # K. oxytoca, root
        # Anything with Unknown is ignored
        self.assertEqual(self.taxonomy_data.get_lca(['571','0']), '571')  # K. oxytoca, Unknown
        self.assertEqual(self.taxonomy_data.get_lca(['2','2157']), '131567')  # Bacteria, Archaea
        self.assertEqual(self.taxonomy_data.get_lca(['571','573']), '570')  # K. oxytoca, K. pneumoniae

    def test_0070_get_upper_level_taxon(self):
        self.assertEqual(self.taxonomy_data.get_upper_level_taxon('1'), ('1', 'norank'))  # Test root
        self.assertEqual(self.taxonomy_data.get_upper_level_taxon('0'), ('1', 'norank'))  # Test Unknown
        # For Bacteria, returns 1 (root), not 131567 ('cellular organisms')
        self.assertEqual(self.taxonomy_data.get_upper_level_taxon('2'), ('1', 'norank'))  # Test Bacteria
        # For 651137 (Thaumarchaeota), returns 2157 (Archaea), not 1783275 ('TACK group')
        self.assertEqual(self.taxonomy_data.get_upper_level_taxon('651137'), ('2157', 'superkingdom'))
        # Test 44260 (Moorella): skip parent taxon 42857 and report 186814 (Thermoanaerobacteriaceae family)
        self.assertEqual(self.taxonomy_data.get_upper_level_taxon('44260'), ('186814', 'family'))

    def test_0080_get_upper_level_taxon(self):
        self.assertEqual(self.taxonomy_data.get_taxonomy_lineage('1'), '')  # Test root
        self.assertEqual(self.taxonomy_data.get_taxonomy_lineage('0'), 'Unknown')  # Test Unknown
        # For Bacteria, returns 1 (root), not 131567 ('cellular organisms')
        self.assertEqual(self.taxonomy_data.get_taxonomy_lineage('2'), 'Bacteria')  # Test Bacteria
        # Test 651137 (Thaumarchaeota), returns 2157 (Archaea), not 1783275 ('TACK group')
        self.assertEqual(self.taxonomy_data.get_taxonomy_lineage('651137'), 'Archaea_Thaumarchaeota')
        # Test 1660251 (Acidobacteria bacterium Mor1): skip two taxa and report 57723 (Acidobacteria)
        self.assertEqual(self.taxonomy_data.get_taxonomy_lineage('1525'),
            'Bacteria_Firmicutes_Clostridia_Thermoanaerobacterales_Thermoanaerobacteraceae_Moorella_Moorella_thermoacetica'
            )

# Tests of Node class
    def test_0090_init_node(self):
        node = Node(rank='genus', name='Escherichia', taxid='561', parent='543', children=None)
        self.assertEqual(node.rank, 'genus')
        self.assertEqual(node.name, 'Escherichia')
        self.assertEqual(node.taxid, '561')
        self.assertEqual(node.parent, '543')

    def test_0100_add_child(self):
        self.node.add_child('562')
        self.assertEqual(len(self.node.children),1)
        self.assertTrue(self.node.has_child('562'))
        self.assertFalse(self.node.has_child('561'))

    def test_0110_set_parent(self):
        node = Node(rank='species')
        node.set_parent('')
        self.assertIsNone(node.parent)
        # Set parent if None
        node.set_parent('561')
        self.assertEqual(node.parent,'561')
        # Change parent
        node.set_parent('2')
        self.assertEqual(node.parent,'2')

    def test_0120_set_taxid(self):
        node = Node(rank='species')
        # Set taxid if None
        node.set_taxid('561')
        self.assertEqual(node.taxid,'561')
        # Changing taxid is not possible
        node.set_taxid('2')
        self.assertEqual(node.taxid,'561')

    def test_0130_set_rank(self):
        self.assertEqual(self.node.rank,'genus')
        # Change rank
        self.assertTrue(self.node.set_rank('species'))
        self.assertEqual(self.node.rank,'species')
        # Rank must be defined in RANKS
        self.assertFalse(self.node.set_rank('#^$%@!'))
        self.assertEqual(self.node.rank,'species')

    def test_0140_set_attribute(self):
        self.node.set_attribute('score', 0.001)
        self.assertEqual(self.node.attributes['score'], 0.001)
        self.node.set_attribute('score', 1)
        self.assertEqual(self.node.attributes['score'], 1)

    def test_0150_add_attribute(self):
        self.node.set_attribute('score_float', 0.5)
        self.node.add_attribute('score_float', 0.2)
        self.assertEqual(self.node.attributes['score_float'], 0.7)
        self.node.set_attribute('score_int', 1)
        self.node.add_attribute('score_int', 999)
        self.assertEqual(self.node.attributes['score_int'], 1000)

    def test_0160_get_attribute(self):
        self.node.set_attribute('score_float', 0.5)
        self.assertEqual(self.node.get_attribute('score_float'), 0.5)
        self.assertIsNone(self.node.get_attribute('nonexisting_key'))

    def test_0170_is_in_children(self):
        self.node.add_child('562')
        self.assertEqual(len(self.node.children),1)
        self.assertTrue(self.node.has_child('562'))
        self.assertFalse(self.node.has_child('561'))

# Tests of Tree class
    def test_0180_init_tree(self):
        tree = Tree()
        self.assertEqual(tree.root.taxid, '1')
        self.assertEqual(len(tree.data), 1)
        self.assertEqual(tree.data[ROOT_TAXONOMY_ID].taxid, ROOT_TAXONOMY_ID)

    def test_0190_add_node(self):
        # Adding node with non-existing parent must fail
        self.assertFalse(self.tree.add_node(
            Node(rank='genus', name='Escherichia', taxid='561', parent='543', children=None)
            ))
        # Adding empty node must fail
        self.assertFalse(self.tree.add_node(
            Node(rank='species')
            ))
        # Adding second root must fail
        self.assertFalse(self.tree.add_node(
            Node(rank='norank', name='root', taxid='1', parent='1', children=None)
            ))
        # Adding node with existing parent must succeed
        self.assertTrue(self.tree.add_node(
            Node(rank='superkingdom', name='Bacteria', taxid='2', parent='1', children=None)
            ))
        self.assertEqual(len(self.tree.get_node('1').children), 1)
        # Adding node if its parent has children must succeed
        self.assertTrue(self.tree.add_node(
            Node(rank='superkingdom', name='Archaea', taxid='2157', parent='1', children=None)
            ))
        self.assertEqual(len(self.tree.get_node('1').children), 2)

    def test_0200_get_node(self):
        self.tree.add_node(
            Node(rank='superkingdom', name='Bacteria', taxid='2', parent='1', children=None)
            )
        # Getting existing node must succeed
        self.assertEqual(self.tree.get_node('2').parent, '1')
        # Getting non-existing node must fail
        self.assertIsNone(self.tree.get_node('561'))

    def test_0210_is_in_tree(self):
        self.tree.add_node(
            Node(rank='superkingdom', name='Bacteria', taxid='2', parent='1', children=None)
            )
        self.assertTrue(self.tree.is_in_tree('2'))
        self.assertFalse(self.tree.is_in_tree('562'))
        self.assertFalse(self.tree.is_in_tree(True))

    def test_0220_add_node_recursively(self):
        self.assertFalse(self.tree.is_in_tree('562'))
        self.assertTrue(self.tree.add_node_recursively(
            Node(rank='species', name='Escherichia coli', taxid='562', parent='561', children=None),
            self.taxonomy_data
            ))
        self.assertTrue(self.tree.is_in_tree('562'))
        self.assertTrue(self.tree.is_in_tree('561'))
        self.assertTrue(self.tree.is_in_tree('2'))
        # Adding second root must fail
        self.assertFalse(self.tree.add_node_recursively(
            Node(rank='norank', name='root', taxid='1', parent='1', children=None),
            self.taxonomy_data
            ))
        # Adding node with existing parent must succeed
        self.assertFalse(self.tree.is_in_tree('564'))
        self.assertTrue(self.tree.add_node_recursively(
            Node(rank='species', name='Escherichia fergusonii', taxid='564', parent='561', children=None),
            self.taxonomy_data
            ))
        self.assertTrue(self.tree.is_in_tree('564'))

    def test_0230_add_attribute(self):
        self.tree.add_node(
            Node(rank='superkingdom', name='Bacteria', taxid='2', parent='1', children=None)
            )
        # Add attribute to existing node
        self.assertFalse(self.tree.get_node('2').attributes)
        self.tree.add_attribute('2', 'score_float', 0.42, self.taxonomy_data)
        self.assertTrue(self.tree.get_node('2').attributes)
        self.assertEqual(self.tree.get_node('2').attributes['score_float'], 0.42)
        # Add attribute to non-existing node
        self.assertFalse(self.tree.get_node('1').attributes)
        self.tree.add_attribute('2157', 'score_float', 0.42, self.taxonomy_data)
        self.assertTrue(self.tree.get_node('1').attributes)
        self.assertEqual(self.tree.get_node('1').attributes['score_float'], 0.42)

    def test_0240_add_attribute_recursively(self):
        self.tree.add_node(
            Node(rank='superkingdom', name='Bacteria', taxid='2', parent='1', children=None)
            )
        # Add attribute to existing node
        self.assertFalse(self.tree.get_node('2').attributes)
        self.tree.add_attribute_recursively('2', 'score_float', 0.42, self.taxonomy_data)
        self.assertTrue(self.tree.get_node('2').attributes)
        self.assertEqual(self.tree.get_node('2').attributes['score_float'], 0.42)
        self.assertEqual(self.tree.get_node('1').attributes['score_float'], 0.42)
        # Add attribute to non-existing node
        self.tree.add_attribute_recursively('2157', 'score_float', 0.42, self.taxonomy_data)
        self.assertEqual(self.tree.get_node('1').attributes['score_float'], 0.84)

    def test_0250_get_parent(self):
        self.assertTrue(self.tree.add_node_recursively(
            Node(rank='species', name='Escherichia coli', taxid='562', parent='561', children=None),
            self.taxonomy_data
            ))
        # Getting parent of existing node must succeed
        node = self.tree.get_node('561')
        self.assertEqual(node.parent, '543')
        parent_node = self.tree.get_parent(node, self.taxonomy_data)
        self.assertEqual(parent_node.taxid, '543')
        # Getting parent of non-existing node must fail
        node = Node(rank='species')
        parent_node = self.tree.get_parent(node, self.taxonomy_data)
        self.assertIsNone(parent_node)

# Tests of TaxonomyProfile class
    def test_0260_init_taxonomy_profile(self):
        taxonomy_profile = TaxonomyProfile()
        self.assertIsNotNone(taxonomy_profile.tree)
        self.assertEqual(len(taxonomy_profile.tree.data), 1)
        self.assertEqual(taxonomy_profile.tree.root.taxid, ROOT_TAXONOMY_ID)
        self.assertEqual(taxonomy_profile.tree.data[ROOT_TAXONOMY_ID].taxid, ROOT_TAXONOMY_ID)

    def test_0270_make_function_taxonomy_profile(self):
        self.project.import_reads_json(sample_id, ENDS)
        scores = get_function_taxonomy_scores(self.project, sample_id=sample_id, metric='fpk')
        sample_scores_taxonomy = slice_function_taxonomy_scores(scores, sample_id)
        taxonomy_profile = TaxonomyProfile()
        taxonomy_profile.make_function_taxonomy_profile(self.project.taxonomy_data, sample_scores_taxonomy)
        print(taxonomy_profile)
        self.assertEqual(len(taxonomy_profile.tree.data), 22)
        self.assertEqual(taxonomy_profile.tree.data[ROOT_TAXONOMY_ID].attributes['NirK']['count'], 1.0)
        self.assertEqual(taxonomy_profile.tree.data[ROOT_TAXONOMY_ID].attributes['UreA']['count'], 3.0)
        self.assertEqual(taxonomy_profile.tree.data['118883'].attributes['UreC']['count'], 1.0)

    def test_0280_str(self):
        taxonomy_profile = TaxonomyProfile()
        self.assertEqual(str(taxonomy_profile),
            '1\tnorank\troot\tParent:None\tChildren:None\tScore:N/A\tIdentity:N/A\tRead count:N/A\n'
        )

    def test_0290_stringify_node(self):
        self.project.import_reads_json(sample_id, ENDS)
        scores = get_function_taxonomy_scores(self.project, sample_id=sample_id, metric='fpk')
        sample_scores_taxonomy = slice_function_taxonomy_scores(scores, sample_id)
        taxonomy_profile = TaxonomyProfile()
        taxonomy_profile.make_function_taxonomy_profile(self.project.taxonomy_data, sample_scores_taxonomy)
        self.assertEqual(taxonomy_profile.stringify_node('118883',0),
            "118883\tfamily\tSulfolobaceae\tParent:2281\tChildren:None\tUreC:{'count': 1.0, 'hit_count': 1.0"
            ", 'identity': 68.8, 'fpk': 0.5984440454817475}\n"
        )

    def test_0300_convert_profile_into_df(self):
        self.project.import_reads_json(sample_id, ENDS)
        scores = get_function_taxonomy_scores(self.project, sample_id=sample_id, metric='fpk')
        sample_scores_taxonomy = slice_function_taxonomy_scores(scores, sample_id)
        taxonomy_profile = TaxonomyProfile()
        taxonomy_profile.make_function_taxonomy_profile(self.project.taxonomy_data, sample_scores_taxonomy)
        result = taxonomy_profile.convert_profile_into_df(metric='fpk')
        self.assertEqual(result.iloc[0][1], 'root')
        self.assertEqual(result.iloc[1][0], 'superkingdom')

    def test_0310_convert_node_into_dict(self):
        self.project.import_reads_json(sample_id, ENDS)
        scores = get_function_taxonomy_scores(self.project, sample_id=sample_id, metric='fpk')
        sample_scores_taxonomy = slice_function_taxonomy_scores(scores, sample_id)
        print(sample_scores_taxonomy)
        taxonomy_profile = TaxonomyProfile()
        taxonomy_profile.make_function_taxonomy_profile(self.project.taxonomy_data, sample_scores_taxonomy)
        result, attributes = taxonomy_profile.convert_node_into_dict('118883', ['UreC', 'UreA'], 1, metric='fpk')
        self.assertEqual(result[1][('', 'Taxon name')], 'Sulfolobaceae')
        self.assertEqual(result[1][('UreC', '1.Score')], 0.5984440454817475)
        self.assertEqual(result[1][('UreA', '1.Score')], 0.0)
        self.assertEqual(attributes['UreC']['fpk'], 0.5984440454817475)
        self.assertEqual(attributes['UreA']['fpk'], 0.0)

    def test_0320_convert_profile_into_score_df(self):
        self.project.import_reads_json(sample_id, ENDS)
        scores = get_function_taxonomy_scores(self.project, sample_id=sample_id, metric='fpk')
        sample_scores_taxonomy = slice_function_taxonomy_scores(scores, sample_id)
        taxonomy_profile = TaxonomyProfile()
        taxonomy_profile.make_function_taxonomy_profile(self.project.taxonomy_data, sample_scores_taxonomy)
        result = taxonomy_profile.convert_profile_into_score_df(metric='fpk')
        self.assertEqual(result.iloc[0][1], 'root')
        self.assertEqual(result.iloc[1][0], 'superkingdom')

    def test_0330_convert_node_into_dict(self):
        self.project.import_reads_json(sample_id, ENDS)
        scores = get_function_taxonomy_scores(self.project, sample_id=sample_id, metric='fpk')
        sample_scores_taxonomy = slice_function_taxonomy_scores(scores, sample_id)
        taxonomy_profile = TaxonomyProfile()
        taxonomy_profile.make_function_taxonomy_profile(self.project.taxonomy_data, sample_scores_taxonomy)
        result, attributes = taxonomy_profile.convert_node_into_values_dict('118883', ['UreC', 'UreA'], 1, metric='fpk')
        self.assertEqual(result[1][('', 'Taxon name')], 'Sulfolobaceae')
        self.assertEqual(result[1][('UreC', 'fpk')], 0.5984440454817475)
        self.assertEqual(result[1][('UreA', 'fpk')], 0.0)
        self.assertEqual(attributes['UreC']['fpk'], 0.5984440454817475)
        self.assertNotIn('UreA', attributes.keys())

    def tearDown(self):
        self.parser = None


if __name__ == '__main__':
    unittest.main()
