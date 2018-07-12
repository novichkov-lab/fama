#!/usr/bin/python
import os
from collections import defaultdict,Counter,OrderedDict
from subprocess import Popen, PIPE, CalledProcessError
from Fama.ReferenceLibrary.TaxonomyData import TaxonomyData


def generate_functions_chart(parser):
    
    
    outfile = os.path.join(parser.project.get_project_dir(parser.sample), parser.project.get_output_subdir(parser.sample),parser.sample + '_' + parser.end + '_'+ parser.project.get_xml_name())
    
    with open(outfile, 'w') as of:
        # Write header
        of.write('<krona key="false">\n')
        of.write('\t<attributes magnitude="rpkm">\n')
        of.write('\t\t<attribute display="Read count">readcount</attribute>\n')
        of.write('\t\t<attribute display="RPKM">rpkm</attribute>\n')
        of.write('\t\t<attribute display="Best hit identity %" mono="true">identity</attribute>\n')
        of.write('\t</attributes>\n')
        of.write('\t<color attribute="identity" valueStart="50" valueEnd="100" hueStart="0" hueEnd="240" default="true"></color>\n')
        
        # Write dataset
        of.write('\t<datasets>\n')
        of.write('\t\t<dataset>' + parser.sample + '</dataset>\n')
        of.write('\t</datasets>\n')
        
        read_count = 0
        total_rpkm = 0.0
        groups_rpkm = defaultdict(float)
        groups_counts = defaultdict(set)
        groups_identity = defaultdict(list)
        functions_counts = defaultdict(set)
        functions_rpkm = defaultdict(float)
        functions_identity = defaultdict(list)
        for read in parser.reads.keys():
            if parser.reads[read].get_status() == 'function,besthit' or parser.reads[read].get_status() == 'function':
                read_count += 1
                read_functions = parser.reads[read].get_functions()
                for function in read_functions:
                    total_rpkm += read_functions[function]
                    groups_counts[parser.ref_data.lookup_function_group(function)].add(read)
                    functions_rpkm[function] += read_functions[function]
                    groups_rpkm[parser.ref_data.lookup_function_group(function)] += read_functions[function]
                    functions_counts[function].add(read)
                for hit in parser.reads[read].get_hit_list().get_hits():
                    for function in hit.get_functions():
                        functions_identity[function].append(hit.get_identity())
                        groups_identity[parser.ref_data.lookup_function_group(function)].append(hit.get_identity())
        
        # Write nodes
        # Write top-level node
        of.write('\t<node name="' + parser.sample + '_' + parser.end + '">\n')
        of.write('\t\t<readcount><val>' + str(read_count) + '</val></readcount>\n')
        of.write('\t\t<rpkm><val>' + str(total_rpkm) + '</val></rpkm>\n')
        
        for group in groups_rpkm:
            # Write group-level node
            of.write('\t\t<node name="' + group + '">\n')
            of.write('\t\t\t<readcount><val>' + str(len(groups_counts[group])) + '</val></readcount>\n')
            of.write('\t\t\t<rpkm><val>' + str(groups_rpkm[group]) + '</val></rpkm>\n')
            if group in groups_identity:
                of.write('\t\t\t<identity><val>' + str(sum(groups_identity[group])/len(groups_identity[group])) + '</val></identity>\n')
            else:
                of.write('\t\t\t<identity><val>0.0</val></identity>\n')
            for function in parser.ref_data.get_functions_in_group(group):
                if function in functions_rpkm:
                    # Write function-level node
                    of.write('\t\t\t<node name="' + function + '">\n')
                    of.write('\t\t\t\t<readcount><val>' + str(len(functions_counts[function])) + '</val></readcount>\n')
                    of.write('\t\t\t\t<rpkm><val>' + str(functions_rpkm[function]) + '</val></rpkm>\n')
                    if function in functions_identity:
                        of.write('\t\t\t\t<identity><val>' + str(sum(functions_identity[function])/len(functions_identity[function])) + '</val></identity>\n')
                    else:
                        of.write('\t\t\t\t<identity><val>0.0</val></identity>\n')
                    of.write('\t\t\t</node>\n')
            # Close group-level node
            of.write('\t\t</node>\n')
        # Close top-level node
        of.write('\t</node>\n')
        
        # Close XML
        of.write('</krona>')
        of.closed

    # Run Krona
    html_file = os.path.join(parser.project.get_project_dir(parser.sample), parser.project.get_output_subdir(parser.sample),parser.sample + '_' + parser.end + '_'+ parser.project.get_html_name())
    krona_cmd = ['ktImportXML', '-o', html_file, outfile]

    with Popen(krona_cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

def print_tax_xml(tax_profile, taxid, offset):
    #print(taxid)
    if taxid not in tax_profile.tree.data:
        raise Exception (taxid,'not found in the tree!!!')
    ret_val = '\t'*offset + '<node name="' + tax_profile.tree.data[taxid].name + '">\n'
    offset += 1
    if tax_profile.tree.data[taxid].attributes:
        ret_val += '\t'*offset + '<readcount><val>' + format(tax_profile.tree.data[taxid].attributes['count'], "0.0f") + '</val></readcount>\n'
        ret_val += '\t'*offset + '<rpkm><val>' + format((tax_profile.tree.data[taxid].attributes['rpkm']), "0.2f") + '</val></rpkm>\n'
        ret_val += '\t'*offset + '<identity><val>' + format((tax_profile.tree.data[taxid].attributes['identity']/tax_profile.tree.data[taxid].attributes['count']), "0.1f") + '</val></identity>\n'
    else:
        ret_val += '\t'*offset + '<readcount><val>0</val></readcount>\n'
        ret_val += '\t'*offset + '<rpkm><val>0.0</val></rpkm>\n'
        ret_val += '\t'*offset + '<identity><val>0.0</val></identity>\n'
        
    if tax_profile.tree.data[taxid].children:
        for child_taxid in tax_profile.tree.data[taxid].children:
            #print('Called print_tax_xml', child_taxid, offset)
            ret_val += print_tax_xml(tax_profile, child_taxid, offset)
    offset -= 1
    ret_val += '\t'*offset + '</node>\n'
    return ret_val

def generate_taxonomy_chart(tax_profile, sample, outfile):
    with open(outfile, 'w') as of:
        # Write header
        of.write('<krona key="false">\n')
        of.write('\t<attributes magnitude="rpkm">\n')
        of.write('\t\t<attribute display="Read count">readcount</attribute>\n')
        of.write('\t\t<attribute display="RPKM">rpkm</attribute>\n')
        of.write('\t\t<attribute display="Best hit identity %" mono="true">identity</attribute>\n')
        of.write('\t</attributes>\n')
        of.write('\t<color attribute="identity" valueStart="50" valueEnd="100" hueStart="0" hueEnd="240" default="true"></color>\n')
        
        # Write dataset
        of.write('\t<datasets>\n')
        of.write('\t\t<dataset>' + sample + '</dataset>\n')
        of.write('\t</datasets>\n')
        
        # Write nodes
        root_id = '1'
        offset = 1
        of.write(print_tax_xml(tax_profile, root_id, offset))
        
        # Close XML
        of.write('</krona>')
        of.closed

    # Run Krona
    html_file = outfile + '.html'
    krona_cmd = ['ktImportXML', '-o', html_file, outfile]

    with Popen(krona_cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

def generate_taxonomy_series_chart(tax_profile, sample_list, outfile):
    with open(outfile, 'w') as of:
        # Write header
        of.write('<krona key="false">\n')
        of.write('\t<attributes magnitude="rpkm">\n')
        of.write('\t\t<attribute display="Read count">readcount</attribute>\n')
        of.write('\t\t<attribute display="RPKM">rpkm</attribute>\n')
        of.write('\t\t<attribute display="Best hit identity %" mono="true">identity</attribute>\n')
        of.write('\t</attributes>\n')
        of.write('\t<color attribute="identity" valueStart="50" valueEnd="100" hueStart="0" hueEnd="240" default="true"></color>\n')
        
        # Write dataset
        of.write('\t<datasets>\n')
        for sample in sample_list:
            of.write('\t\t<dataset>' + sample + '</dataset>\n')
        of.write('\t</datasets>\n')

        
        # Write nodes
        root_id = '1'
        offset = 1
        of.write(print_dataseries_tax_xml(tax_profile, sample_list, root_id, offset))
        
        # Close XML
        of.write('</krona>')
        of.closed

    # Run Krona
    html_file = outfile + '.html'
    krona_cmd = ['ktImportXML', '-o', html_file, outfile]

    with Popen(krona_cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

def print_dataseries_tax_xml(tax_profile, dataseries, taxid, offset):
    #print(taxid)
    if taxid not in tax_profile.tree.data:
        raise Exception (taxid,'not found in the tree!!!')
    ret_val = '\t'*offset + '<node name="' + tax_profile.tree.data[taxid].name + '">\n'
    offset += 1
    if tax_profile.tree.data[taxid].attributes:
        ret_val += '\t'*offset + '<readcount>'
        for datapoint in dataseries:
            if datapoint in tax_profile.tree.data[taxid].attributes:
                ret_val += '<val>' + format(tax_profile.tree.data[taxid].attributes[datapoint]['count'], "0.0f") + '</val>'
            else:
                ret_val += '<val>0</val>'
        ret_val += '</readcount>\n' + '\t'*offset + '<rpkm>'
        for datapoint in dataseries:
            if datapoint in tax_profile.tree.data[taxid].attributes:
                ret_val += '<val>' + format((tax_profile.tree.data[taxid].attributes[datapoint]['rpkm']), "0.2f") + '</val>'
            else:
                ret_val += '<val>0.0</val>'
        ret_val += '</rpkm>\n' + '\t'*offset + '<identity>'
        for datapoint in dataseries:
            if datapoint in tax_profile.tree.data[taxid].attributes:
                ret_val += '<val>' + format((tax_profile.tree.data[taxid].attributes[datapoint]['identity']/tax_profile.tree.data[taxid].attributes[datapoint]['count']), "0.1f") + '</val>'
            else:
                ret_val += '<val>0.0</val>'
        ret_val += '</identity>\n'
    else:
        ret_val += '\t'*offset + '<readcount>'
        ret_val += '<val>0</val>'*len(dataseries)
        ret_val += '</readcount>\n' + '\t'*offset + '<rpkm>'
        ret_val += '<val>0.0</val>'*len(dataseries)
        ret_val += '</rpkm>\n' + '\t'*offset + '<identity>'
        ret_val += '<val>0.0</val>'*len(dataseries)
        ret_val += '</identity>\n'
        
    if tax_profile.tree.data[taxid].children:
        for child_taxid in tax_profile.tree.data[taxid].children:
            #print('Called print_tax_xml', child_taxid, offset)
            ret_val += print_dataseries_tax_xml(tax_profile, dataseries, child_taxid, offset)
    offset -= 1
    ret_val += '\t'*offset + '</node>\n'
    return ret_val

def generate_functional_taxonomy_chart(tax_profile, function_list, outfile):
    with open(outfile, 'w') as of:
        # Write header
        of.write('<krona key="false">\n')
        of.write('\t<attributes magnitude="rpkm">\n')
        of.write('\t\t<attribute display="Read count">readcount</attribute>\n')
        of.write('\t\t<attribute display="RPKM">rpkm</attribute>\n')
        of.write('\t\t<attribute display="Best hit identity %" mono="true">identity</attribute>\n')
        of.write('\t</attributes>\n')
        of.write('\t<color attribute="identity" valueStart="50" valueEnd="100" hueStart="0" hueEnd="240" default="true"></color>\n')
        
        # Write dataset
        of.write('\t<datasets>\n')
        for function in function_list:
            of.write('\t\t<dataset>' + function + '</dataset>\n')
        of.write('\t</datasets>\n')
        
        # Write nodes
        root_id = '1'
        offset = 1
        of.write(print_dataseries_tax_xml(tax_profile, function_list, root_id, offset))
        
        # Close XML
        of.write('</krona>')
        of.closed

    # Run Krona
    html_file = outfile + '.html'
    krona_cmd = ['ktImportXML', '-o', html_file, outfile]

    with Popen(krona_cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='')
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

