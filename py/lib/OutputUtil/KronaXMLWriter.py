#!/usr/bin/python
import os
from collections import defaultdict,Counter,OrderedDict
from lib.ReferenceLibrary.TaxonomyData import TaxonomyData
from subprocess import Popen, PIPE, CalledProcessError


def generate_xml(parser):
    
    
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

