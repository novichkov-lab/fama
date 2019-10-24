"""Various functions for Krona chart generation

Functions starting with "make_" are actual generators. Functions starting with "get_"
return XML internal nodes for taxa.

make_functions_chart: makes Krona chart of function abundances in one sample
make_taxonomy_chart: makes Krona chart of taxon abundances in one sample
make_taxonomy_series_chart: makes Krona chart of taxon abundances in multiple samples
make_function_taxonomy_chart: makes Krona chart of taxon abundances for multiple
    functions in one sample
make_assembly_taxonomy_chart: makes Krona chart of taxon/gene abundances in assemby
"""
import os
from collections import defaultdict
from lib.utils.const import STATUS_GOOD, ROOT_TAXONOMY_ID
from lib.utils.utils import autovivify, run_external_program


def make_functions_chart(parser, metric='efpkg'):
    """Writes XML file for functions chart and generates Krona plot from it

    Args:
        parser (:obj:DiamondParser): parser object with annotated reads
        metric (str): scoring metric (efpkg by default)
    """
    outfile = os.path.join(
        parser.options.get_project_dir(parser.sample.sample_id),
        parser.options.get_output_subdir(parser.sample.sample_id),
        parser.sample.sample_id + '_' + parser.end + '_' + parser.options.xml_name
    )
    with open(outfile, 'w') as out:
        # Write header
        out.write('<krona key="false">\n' +
                  '\t<attributes magnitude="' + metric + '">\n' +
                  '\t\t<attribute display="Read count">readcount</attribute>\n')
        if metric != 'readcount':
            out.write('\t\t<attribute display="Score:' + metric + '">' + metric + '</attribute>\n')
        out.write('\t\t<attribute display="AAI %" mono="true">identity</attribute>\n'
                  + '\t</attributes>\n'
                  + ' '.join(['\t<color attribute="identity"', 'valueStart="50"', 'valueEnd="100"',
                              'hueStart="0"', 'hueEnd="240"', 'default="true"></color>\n'])
                  # Write dataset
                  + '\t<datasets>\n\t\t<dataset>' + parser.sample.sample_id
                  + '</dataset>\n\t</datasets>\n')

        read_count = 0
        total_rpkm = 0.0
        groups_rpkm = defaultdict(float)
        groups_counts = defaultdict(set)
        groups_identity = defaultdict(list)
        functions_counts = defaultdict(set)
        functions_rpkm = defaultdict(float)
        functions_identity = defaultdict(list)
        for _, read in parser.reads.items():
            if read.status == STATUS_GOOD:
                read_count += 1
                for function in read.functions:
                    total_rpkm += read.functions[function]
                    groups_counts[parser.ref_data.lookup_function_group(function)].add(read.read_id)
                    functions_rpkm[function] += read.functions[function]
                    groups_rpkm[parser.ref_data.lookup_function_group(function)] += \
                        read.functions[function]
                    functions_counts[function].add(read.read_id)
                for hit in read.hit_list.hits:
                    for function in hit.functions:
                        functions_identity[function].append(hit.identity)
                        groups_identity[parser.ref_data.lookup_function_group(function)]\
                            .append(hit.identity)

        # Write nodes
        # Write top-level node
        out.write('\t<node name="' + parser.sample.sample_id + '_' + parser.end + '">\n')
        if metric != 'readcount':
            out.write('\t\t<readcount><val>' + str(read_count) + '</val></readcount>\n')
        out.write('\t\t<' + metric + '><val>' + str(total_rpkm) + '</val></' + metric + '>\n')

        for group in groups_rpkm:
            # Write group-level node
            out.write('\t\t<node name="' + group + '">\n')
            if metric != 'readcount':
                out.write('\t\t\t<readcount><val>'
                          + str(len(groups_counts[group]))
                          + '</val></readcount>\n')
            out.write('\t\t\t<' + metric + '><val>' + str(groups_rpkm[group])
                      + '</val></' + metric + '>\n')
            if group in groups_identity:
                out.write('\t\t\t<identity><val>'
                          + str(sum(groups_identity[group]) / len(groups_identity[group]))
                          + '</val></identity>\n')
            else:
                out.write('\t\t\t<identity><val>0.0</val></identity>\n')
            for function in parser.ref_data.get_functions_in_group(group):
                if function in functions_rpkm:
                    # Write function-level node
                    out.write('\t\t\t<node name="' + function + '">\n')
                    if metric != 'readcount':
                        out.write('\t\t\t\t<readcount><val>'
                                  + str(len(functions_counts[function]))
                                  + '</val></readcount>\n')
                    out.write('\t\t\t\t<' + metric + '><val>' + str(functions_rpkm[function])
                              + '</val></' + metric + '>\n')
                    if function in functions_identity:
                        out.write(
                            '\t\t\t\t<identity><val>' + str(sum(
                                functions_identity[function]
                            ) / len(functions_identity[function])) + '</val></identity>\n'
                        )
                    else:
                        out.write('\t\t\t\t<identity><val>0.0</val></identity>\n')
                    out.write('\t\t\t</node>\n')
            # Close group-level node
            out.write('\t\t</node>\n')
        # Close top-level node
        out.write('\t</node>\n</krona>')
    # Run Krona
    html_file = os.path.join(parser.options.get_project_dir(parser.sample.sample_id),
                             parser.options.get_output_subdir(parser.sample.sample_id),
                             parser.sample.sample_id + '_' + parser.end + '_'
                             + parser.options.html_name)
    run_external_program([parser.config.krona_path, '-o', html_file, outfile])


def get_taxon_xml(tax_profile, taxid, offset, metric='efpkg'):
    """Returns XML node for a phylogenetic tree node and all its children

    Args:
        tax_profile (:obj:TaxonomyProfile): taxonomy profile
        taxid (str): taxonomy identifier of a node of interest
        offset (int): number of starting tabs
        metric (str): scoring metric (default value 'efpkg')

    Returns:
        ret_val (str): XML node
    """
    if taxid not in tax_profile.tree.data:
        raise KeyError(taxid, 'not found in the tree!!!')
    ret_val = '\t'*offset + '<node name="' + tax_profile.tree.data[taxid].name + '">\n'
    offset += 1
    if tax_profile.tree.data[taxid].attributes:
        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount><val>' + format(
                tax_profile.tree.data[taxid].attributes['count'], "0.0f"
            ) + '</val></readcount>\n'
        ret_val += '\t'*offset + '<' + metric + '><val>' + format(
            tax_profile.tree.data[taxid].attributes[metric], "0.2f"
        ) + '</val></' + metric + '>\n'
        ret_val += '\t'*offset + '<identity><val>' + format((
            tax_profile.tree.data[taxid].attributes['identity']
            / tax_profile.tree.data[taxid].attributes['hit_count']
        ), "0.1f") + '</val></identity>\n'
    else:
        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount><val>0</val></readcount>\n'
        ret_val += '\t'*offset + '<' + metric + '><val>0.0</val></' + metric + '>\n'
        ret_val += '\t'*offset + '<identity><val>0.0</val></identity>\n'

    if tax_profile.tree.data[taxid].children:
        for child_taxid in tax_profile.tree.data[taxid].children:
            ret_val += get_taxon_xml(tax_profile, child_taxid, offset, metric)
    offset -= 1
    ret_val += '\t'*offset + '</node>\n'
    return ret_val


def get_lca_tax_xml(tax_profile, taxid, offset, metric='efpkg'):
    """Returns XML node for a phylogenetic tree node and all its children.
    Creates additional child node for a fictional "Unclassified..." taxon
    if not all reads of the current node are mapped to children nodes.

    Args:
        tax_profile (:obj:TaxonomyProfile): taxonomy profile
        taxid (str): taxonomy identifier of a node of interest
        offset (int): number of starting tabs
        metric (str): scoring metric (default value 'efpkg')

    Returns:
        ret_val (str): XML node
    """
    attribute_values = defaultdict(float)
    try:
        ret_val = '\t'*offset + '<node name="' + tax_profile.tree.data[taxid].name + '">\n'
    except KeyError:
        print(taxid, 'not found in the tree data!!!')
        raise KeyError
    offset += 1
    if tax_profile.tree.data[taxid].attributes:
        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount><val>' + format(
                tax_profile.tree.data[taxid].attributes['count'], "0.0f"
            ) + '</val></readcount>\n'
        ret_val += '\t'*offset + '<' + metric + '><val>' + format(
            tax_profile.tree.data[taxid].attributes[metric], "0.2f"
        ) + '</val></' + metric + '>\n'
        ret_val += '\t'*offset + '<identity><val>' + format((
            tax_profile.tree.data[taxid].attributes['identity']
            / tax_profile.tree.data[taxid].attributes['hit_count']
        ), "0.1f") + '</val></identity>\n'
    else:
        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount><val>0</val></readcount>\n'
        ret_val += '\t'*offset + '<' + metric + '><val>0.0</val></' + metric + '>\n'
        ret_val += '\t'*offset + '<identity><val>0.0</val></identity>\n'

    if tax_profile.tree.data[taxid].children:
        for child_taxid in tax_profile.tree.data[taxid].children:
            child_node, child_values = get_lca_tax_xml(tax_profile, child_taxid, offset, metric)
            ret_val += child_node
            for key, val in child_values.items():
                attribute_values[key] += val

        # Add a fictional child node for unclassified reads, if any
        if tax_profile.tree.data[taxid].attributes and (
                attribute_values['count'] < tax_profile.tree.data[taxid].attributes['count']
        ):
            unknown_node = 'Unidentified ' + tax_profile.tree.data[taxid].name
            if offset == 2:
                unknown_node = 'Unknown'
            ret_val += '\t'*offset + '<node name="' + unknown_node + '">\n'
            offset += 1
            if metric != 'readcount':
                ret_val += '\t'*offset + '<readcount><val>' + format((
                    tax_profile.tree.data[taxid].attributes['count']
                    - attribute_values['count']
                ), "0.0f") + '</val></readcount>\n'
            ret_val += '\t'*offset + '<' + metric + '><val>' + format((
                tax_profile.tree.data[taxid].attributes[metric]
                - attribute_values[metric]
            ), "0.2f") + '</val></' + metric + '>\n'
            if tax_profile.tree.data[taxid].attributes['hit_count'] > attribute_values['hit_count']:
                ret_val += '\t'*offset + '<identity><val>' + format(((
                    tax_profile.tree.data[taxid].attributes['identity']
                    - attribute_values['identity']
                ) / (
                    tax_profile.tree.data[taxid].attributes['hit_count']
                    - attribute_values['hit_count']
                )), "0.1f") + '</val></identity>\n'
            else:
                ret_val += '\t'*offset + '<identity><val>0.0</val></identity>\n'
            offset -= 1
            ret_val += '\t'*offset + '</node>\n'

    offset -= 1
    ret_val += '\t'*offset + '</node>\n'
    attribute_values = defaultdict(float)
    attribute_values[metric] = tax_profile.tree.data[taxid].attributes[metric]
    attribute_values['count'] = tax_profile.tree.data[taxid].attributes['count']
    attribute_values['identity'] = tax_profile.tree.data[taxid].attributes['identity']
    attribute_values['hit_count'] = tax_profile.tree.data[taxid].attributes['hit_count']
    return ret_val, attribute_values


def make_taxonomy_chart(tax_profile, sample, outfile, krona_path, metric='efpkg'):
    """Writes XML file for taxonomy chart of one sample and generates Krona plot from it

    Args:
        tax_profile (:obj:TaxonomyProfile): taxonomy profile object
        sample (str): sample identifier
        outfile (str): path for XML output
        krona_path (str): Krona Tools command
        metric (str): scoring metric (efpkg by default)
    """
    with open(outfile, 'w') as out:
        # Write header
        out.write('<krona key="false">\n')
        out.write('\t<attributes magnitude="' + metric + '">\n')
        if metric != 'readcount':
            out.write('\t\t<attribute display="Read count">readcount</attribute>\n')
        out.write('\t\t<attribute display="Score:' + metric + '">' + metric + '</attribute>\n')
        out.write('\t\t<attribute display="AAI %" mono="true">identity</attribute>\n')
        out.write('\t</attributes>\n')
        out.write('\t<color attribute="identity" valueStart="50" valueEnd="100" hueStart="0"'
                  + ' hueEnd="240" default="true"></color>\n')
        # Write dataset
        out.write('\t<datasets>\n')
        out.write('\t\t<dataset>' + sample + '</dataset>\n')
        out.write('\t</datasets>\n')

        # Write nodes
        offset = 1
        child_node, _ = get_lca_tax_xml(tax_profile, ROOT_TAXONOMY_ID, offset, metric=metric)
        out.write(child_node)
        # Close XML
        out.write('</krona>')

    # Run Krona
    html_file = outfile + '.html'
    krona_cmd = [krona_path, '-o', html_file, outfile]
    run_external_program(krona_cmd)


def make_taxonomy_series_chart(tax_profile, sample_list, outfile, krona_path, metric='efpkg'):
    """Writes XML file for taxonomy chart of multiple samples and generates Krona plot for it.
    Taxonomy profile must have two-level attributes, with function identifier as outer key and
    a metric as inner key.

    Args:
        tax_profile (:obj:TaxonomyProfile): taxonomy profile object
        sample_list (list of str): sample identifiers
        outfile (str): path for XML output
        krona_path (str): Krona Tools command
        metric (str): scoring metric (efpkg by default)
    """
    with open(outfile, 'w') as out:
        # Write header
        out.write('<krona key="false">\n')
        out.write('\t<attributes magnitude="' + metric + '">\n')
        if metric != 'readcount':
            out.write('\t\t<attribute display="Read count">readcount</attribute>\n')
        out.write('\t\t<attribute display="Score:' + metric + '">' + metric + '</attribute>\n')
        out.write('\t\t<attribute display="AAI %" mono="true">identity</attribute>\n')
        out.write('\t</attributes>\n')
        out.write('\t<color attribute="identity" valueStart="50" valueEnd="100" hueStart="0" ' +
                  'hueEnd="240" default="true"></color>\n')
        # Write dataset
        out.write('\t<datasets>\n')
        for sample in sample_list:
            out.write('\t\t<dataset>' + sample + '</dataset>\n')
        out.write('\t</datasets>\n')
        # Write nodes
        offset = 1
        child_nodes, _ = get_lca_dataseries_tax_xml(
            tax_profile, sample_list, ROOT_TAXONOMY_ID, offset, metric=metric
            )
        out.write(child_nodes)
        # Close XML
        out.write('</krona>')

    # Run Krona
    html_file = outfile + '.html'
    krona_cmd = [krona_path, '-o', html_file, outfile]
    run_external_program(krona_cmd)


def get_lca_dataseries_tax_xml(tax_profile, dataseries, taxid, offset, metric='efpkg'):
    """Returns XML node for a phylogenetic tree node and all its children.
    Creates additional child node for a fictional "Unclassified..." taxon
    if not all reads of the current node were mapped to children nodes.

    Args:
        tax_profile (:obj:TaxonomyProfile): taxonomy profile
        dataseries (list of str): either sample identifiers or function identifiers,
            depending on profile type (functional or taxonomic)
        taxid (str): taxonomy identifier of a node of interest
        offset (int): number of starting tabs
        metric (str): scoring metric (default value 'efpkg')

    Returns:
        ret_val (str): XML node
        attribute_values (defaultdict[str,dict[str,float]]): outer key is
            one of dataseries members, inner key is in [metric, 'count', 'identity'
            'hit_count'], value is float.
    """
    attribute_values = autovivify(2, float)

    if taxid not in tax_profile.tree.data:
        raise KeyError(taxid, 'not found in the tree!!!')
    ret_val = '\t'*offset + '<node name="' + tax_profile.tree.data[taxid].name + '">\n'
    offset += 1
    if tax_profile.tree.data[taxid].attributes:
        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount>'
            for datapoint in dataseries:
                if (datapoint in tax_profile.tree.data[taxid].attributes) and (
                        'count' in tax_profile.tree.data[taxid].attributes[datapoint]
                ):
                    ret_val += '<val>' + format(
                        tax_profile.tree.data[taxid].attributes[datapoint]['count'], "0.0f"
                        ) + '</val>'
                else:
                    ret_val += '<val>0</val>'
            ret_val += '</readcount>\n'
        ret_val += '\t'*offset + '<' + metric + '>'
        for datapoint in dataseries:
            if datapoint in tax_profile.tree.data[taxid].attributes and (
                    metric in tax_profile.tree.data[taxid].attributes[datapoint]
            ):
                ret_val += '<val>' + format(
                    tax_profile.tree.data[taxid].attributes[datapoint][metric], "0.6f"
                ) + '</val>'
            else:
                ret_val += '<val>0.0</val>'
        ret_val += '</' + metric + '>\n' + '\t'*offset + '<identity>'
        for datapoint in dataseries:
            if datapoint in tax_profile.tree.data[taxid].attributes and (
                    'identity' in tax_profile.tree.data[taxid].attributes[datapoint]
            ):
                ret_val += '<val>' + format((
                    tax_profile.tree.data[taxid].attributes[datapoint]['identity']
                    / tax_profile.tree.data[taxid].attributes[datapoint]['hit_count']
                ), "0.1f") + '</val>'
            else:
                ret_val += '<val>0.0</val>'
        ret_val += '</identity>\n'
    else:
        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount>'
            ret_val += '<val>0</val>'*len(dataseries)
            ret_val += '</readcount>\n'
        ret_val += '\t'*offset + '<' + metric + '>'
        ret_val += '<val>0.0</val>'*len(dataseries)
        ret_val += '<' + metric + '>\n' + '\t'*offset + '<identity>'
        ret_val += '<val>0.0</val>'*len(dataseries)
        ret_val += '</identity>\n'

    if tax_profile.tree.data[taxid].children:
        for child_taxid in tax_profile.tree.data[taxid].children:
            child_node, child_values = get_lca_dataseries_tax_xml(tax_profile,
                                                                  dataseries,
                                                                  child_taxid,
                                                                  offset,
                                                                  metric=metric)
            ret_val += child_node
            for datapoint in child_values.keys():
                for key, val in child_values[datapoint].items():
                    attribute_values[datapoint][key] += val
        # Add a child node for unidentified child taxon, if needed
        unidentified_flag = False
        for datapoint in dataseries:
            if datapoint in tax_profile.tree.data[taxid].attributes:
                if (
                        attribute_values[datapoint]['count']
                        < tax_profile.tree.data[taxid].attributes[datapoint]['count']
                ):
                    unidentified_flag = True
                    break

        if unidentified_flag:
            if offset == 2:
                ret_val += '\t'*offset + '<node name="Unclassified">\n'
            else:
                ret_val += '\t'*offset + '<node name="Unclassified '\
                           + tax_profile.tree.data[taxid].name + '">\n'
            offset += 1
            if metric != 'readcount':
                ret_val += '\t'*offset + '<readcount>'
                for datapoint in dataseries:
                    if datapoint in tax_profile.tree.data[taxid].attributes and (
                            attribute_values[datapoint]['count']
                            < tax_profile.tree.data[taxid].attributes[datapoint]['count']
                    ):
                        ret_val += '<val>' + format((
                            tax_profile.tree.data[taxid].attributes[datapoint]['count']
                            - attribute_values[datapoint]['count']
                        ), "0.0f") + '</val>'
                    else:
                        ret_val += '<val>0</val>'
                ret_val += '</readcount>\n'
            ret_val += '\t'*offset + '<' + metric + '>'
            for datapoint in dataseries:
                if datapoint in tax_profile.tree.data[taxid].attributes and (
                        attribute_values[datapoint]['count']
                        < tax_profile.tree.data[taxid].attributes[datapoint]['count']
                ):
                    ret_val += '<val>' + format((
                        tax_profile.tree.data[taxid].attributes[datapoint][metric]
                        - attribute_values[datapoint][metric]
                    ), "0.6f") + '</val>'
                else:
                    ret_val += '<val>0.0</val>'
            ret_val += '</' + metric + '>\n'
            ret_val += '\t'*offset + '<identity>'
            for datapoint in dataseries:
                if datapoint in tax_profile.tree.data[taxid].attributes and (
                        'hit_count' in tax_profile.tree.data[taxid].attributes[datapoint]
                ) and (
                    attribute_values[datapoint]['hit_count']
                    < tax_profile.tree.data[taxid].attributes[datapoint]['hit_count']
                ):
                    ret_val += '<val>' + format(((
                        tax_profile.tree.data[taxid].attributes[datapoint]['identity']
                        - attribute_values[datapoint]['identity']
                    ) / (
                        tax_profile.tree.data[taxid].attributes[datapoint]['hit_count']
                        - attribute_values[datapoint]['hit_count']
                    )), "0.1f") + '</val>'
                else:
                    ret_val += '<val>0.0</val>'
            ret_val += '</identity>\n'
            offset -= 1
            ret_val += '\t'*offset + '</node>\n'
    offset -= 1
    ret_val += '\t'*offset + '</node>\n'
    attribute_values = autovivify(1)
    for datapoint in dataseries:
        if datapoint in tax_profile.tree.data[taxid].attributes:
            if metric in tax_profile.tree.data[taxid].attributes[datapoint]:
                attribute_values[datapoint][metric] = tax_profile.tree.data[taxid]\
                    .attributes[datapoint][metric]
            if 'count' in tax_profile.tree.data[taxid].attributes[datapoint]:
                attribute_values[datapoint]['count'] = tax_profile.tree.data[taxid]\
                    .attributes[datapoint]['count']
            if 'identity' in tax_profile.tree.data[taxid].attributes[datapoint]:
                attribute_values[datapoint]['identity'] = tax_profile.tree.data[taxid]\
                    .attributes[datapoint]['identity']
            if 'hit_count' in tax_profile.tree.data[taxid].attributes[datapoint]:
                attribute_values[datapoint]['hit_count'] = tax_profile.tree.data[taxid]\
                    .attributes[datapoint]['hit_count']
    return ret_val, attribute_values


def get_dataseries_tax_xml(tax_profile, dataseries, taxid, offset, metric='efpkg'):
    """Returns XML node for a phylogenetic tree node and all its children.

    Args:
        tax_profile (:obj:TaxonomyProfile): taxonomy profile
        dataseries (list of str): either sample identifiers or function identifiers,
            depending on profile type (functional or taxonomic)
        taxid (str): taxonomy identifier of a node of interest
        offset (int): number of starting tabs
        metric (str): scoring metric (default value 'efpkg')

    Returns:
        ret_val (str): XML node
    """
    if taxid not in tax_profile.tree.data:
        raise KeyError(taxid, 'not found in the tree!!!')
    ret_val = '\t'*offset + '<node name="' + tax_profile.tree.data[taxid].name + '">\n'
    offset += 1
    if tax_profile.tree.data[taxid].attributes:
        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount>'
            for datapoint in dataseries:
                if datapoint in tax_profile.tree.data[taxid].attributes:
                    ret_val += '<val>' + format(
                        tax_profile.tree.data[taxid].attributes[datapoint]['count'], "0.0f"
                    ) + '</val>'
                else:
                    ret_val += '<val>0</val>'
            ret_val += '</readcount>\n'
        ret_val += '\t'*offset + '<' + metric + '>'
        for datapoint in dataseries:
            if datapoint in tax_profile.tree.data[taxid].attributes:
                ret_val += '<val>' + format((
                    tax_profile.tree.data[taxid].attributes[datapoint][metric]
                ), "0.5f") + '</val>'
            else:
                ret_val += '<val>0.0</val>'
        ret_val += '</' + metric + '>\n' + '\t'*offset + '<identity>'
        for datapoint in dataseries:
            if datapoint in tax_profile.tree.data[taxid].attributes:
                ret_val += '<val>' + format((
                    tax_profile.tree.data[taxid].attributes[datapoint]['identity']
                    / tax_profile.tree.data[taxid].attributes[datapoint]['hit_count']
                ), "0.1f") + '</val>'
            else:
                ret_val += '<val>0.0</val>'
        ret_val += '</identity>\n'
    else:
        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount>'
            ret_val += '<val>0</val>'*len(dataseries)
            ret_val += '</readcount>\n'
        ret_val += '\t'*offset + '<' + metric + '>'
        ret_val += '<val>0.0</val>'*len(dataseries)
        ret_val += '<' + metric + '>\n' + '\t'*offset + '<identity>'
        ret_val += '<val>0.0</val>'*len(dataseries)
        ret_val += '</identity>\n'

    if tax_profile.tree.data[taxid].children:
        for child_taxid in tax_profile.tree.data[taxid].children:
            ret_val += get_dataseries_tax_xml(
                tax_profile, dataseries, child_taxid, offset, metric=metric
                )
    offset -= 1
    ret_val += '\t'*offset + '</node>\n'
    return ret_val


def make_function_taxonomy_chart(tax_profile, function_list, outfile, krona_path,
                                 metric='efpkg'):
    """Writes XML file for taxonomy chart of multiple functions in one sample
       and generates Krona plot from it

    Args:
        tax_profile (:obj:TaxonomyProfile): taxonomy profile object
        function_list (list of str): function identifiers
        outfile (str): path for XML output
        krona_path (str): Krona Tools command
        metric (str): scoring metric (efpkg by default)
    """
    with open(outfile, 'w') as out:
        # Write header
        out.write('<krona key="false">\n')
        out.write('\t<attributes magnitude="' + metric + '">\n')
        if metric != 'readcount':
            out.write('\t\t<attribute display="Read count">readcount</attribute>\n')
        out.write('\t\t<attribute display="Score:' + metric + '">' + metric + '</attribute>\n')
        out.write('\t\t<attribute display="Best hit identity %" mono="true">identity</attribute>\n')
        out.write('\t</attributes>\n')
        out.write('\t<color attribute="identity" valueStart="50" valueEnd="100" hueStart="0" '
                  + 'hueEnd="240" default="true"></color>\n')
        # Write dataset
        out.write('\t<datasets>\n')
        for function in function_list:
            out.write('\t\t<dataset>' + function + '</dataset>\n')
        out.write('\t</datasets>\n')
        # Write nodes
        offset = 1
        out.write(
            get_dataseries_tax_xml(
                tax_profile, function_list, ROOT_TAXONOMY_ID, offset, metric=metric
                )
        )
        # Close XML
        out.write('</krona>')

    # Run Krona
    html_file = outfile + '.html'
    krona_cmd = [krona_path, '-o', html_file, outfile]
    run_external_program(krona_cmd)


def get_genes_xml(gene_data, gene_ids, dataseries, offset, metric):
    """Returns XML nodes for all predicted gene from one taxon.

    Args:
        gene_data (defaultdict[str,defaultdict[str,dict[str,float]]]): outer key is
            gene identifier, middle key is function identifier, inner key is in
            [metric, 'count', 'identity', 'coverage', 'Length', 'Completeness'],
            value is float.
        gene_ids (list of str): gene identifiers
        dataseries (list of str): either sample identifiers or function identifiers,
            depending on profile type (functional or taxonomic)
        offset (int): number of starting tabs
        metric (str): scoring metric

    Returns:
        ret_val (str): XML node
        attribute_values (defaultdict[str,dict[str,float]]): outer key is
            one of dataseries members, inner key is in [metric, 'count', 'identity'
            'hit_count'], value is float.
    """
    # gene data: gene_data[gene_id][function][parameter] = parameter_value
    ret_val = ''
    for gene_id in gene_ids:
        ret_val += '\t'*offset + '<node name="' + gene_id + '">\n'
        offset += 1

        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount>'
            for datapoint in dataseries:
                if datapoint in gene_data[gene_id]:
                    ret_val += '<val>' + gene_data[gene_id][datapoint]['count'] + '</val>'
                else:
                    ret_val += '<val>0</val>'
            ret_val += '</readcount>\n'

        ret_val += '\t'*offset + '<' + metric + '>'
        for datapoint in dataseries:
            if datapoint in gene_data[gene_id]:
                ret_val += '<val>' + gene_data[gene_id][datapoint][metric] + '</val>'
            else:
                ret_val += '<val>0</val>'
        ret_val += '</' + metric + '>\n'
        ret_val += '\t'*offset + '<coverage>'
        for datapoint in dataseries:
            if datapoint in gene_data[gene_id]:
                ret_val += '<val>' + gene_data[gene_id][datapoint]['coverage'] + '</val>'
            else:
                ret_val += '<val>0</val>'
        ret_val += '</coverage>\n'

        ret_val += '\t'*offset + '<identity>'
        for datapoint in dataseries:
            if datapoint in gene_data[gene_id]:
                ret_val += '<val>' + gene_data[gene_id][datapoint]['identity'] + '</val>'
            else:
                ret_val += '<val>0</val>'
        ret_val += '</identity>\n'

        ret_val += '\t'*offset + '<Length>'
        for datapoint in dataseries:
            if datapoint in gene_data[gene_id]:
                ret_val += '<val>' + gene_data[gene_id][datapoint]['Length'] + '</val>'
            else:
                ret_val += '<val>0</val>'
        ret_val += '</Length>\n'

        ret_val += '\t'*offset + '<Completeness>'
        for datapoint in dataseries:
            if datapoint in gene_data[gene_id]:
                ret_val += '<val>' + gene_data[gene_id][datapoint]['Completeness'] + '</val>'
            else:
                ret_val += '<val>0</val>'
        ret_val += '</Completeness>\n'

        ret_val += '\t'*offset + '<best_hit>'
        for datapoint in dataseries:
            if datapoint in gene_data[gene_id] and 'Best hit' in gene_data[gene_id][datapoint]:
                ret_val += '<val href="' + gene_data[gene_id][datapoint]['Best hit'] + '">'\
                           + gene_data[gene_id][datapoint]['Best hit'] + '</val>'
            else:
                ret_val += '<val></val>'
        ret_val += '</best_hit>\n'

        offset -= 1
        ret_val += '\t'*offset + '</node>\n'
    return ret_val


def get_assembly_tax_xml(tax_profile, genes, dataseries, taxid, offset, metric='efpkg'):
    """Returns XML node for assembly phylogenetic tree node and all its children.

    Args:
        tax_profile (:obj:TaxonomyProfile): taxonomy profile
        genes (defaultdict[str,defaultdict[str,dict[str,float]]]): outer key is
            gene identifier, middle key is function identifier, inner key is in
            [metric, 'count', 'identity', 'coverage', 'Length', 'Completeness'],
            value is float (genes[gene_id][function_id][parameter_name] = parameter_value).
        dataseries (list of str): function identifiers
        taxid (str): taxonomy identifier of a node of interest
        offset (int): number of starting tabs
        metric (str): scoring metric (default value 'efpkg')

    Returns:
        ret_val (str): XML node
    """
    if taxid not in tax_profile.tree.data:
        raise KeyError(taxid, 'not found in the tree!!!')
    ret_val = '\t'*offset + '<node name="' + taxid + ':' + tax_profile.tree.data[taxid].name\
              + '">\n'
    offset += 1
    if tax_profile.tree.data[taxid].attributes:
        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount>'
            for datapoint in dataseries:
                if datapoint in tax_profile.tree.data[taxid].attributes:
                    ret_val += '<val>' + format(
                        tax_profile.tree.data[taxid].attributes[datapoint]['count'], "0.0f"
                    ) + '</val>'
                else:
                    ret_val += '<val>0</val>'
            ret_val += '</readcount>\n'
        ret_val += '\t' * offset + '<' + metric + '>'
        for datapoint in dataseries:
            if datapoint in tax_profile.tree.data[taxid].attributes:
                ret_val += '<val>' + format(
                    tax_profile.tree.data[taxid].attributes[datapoint][metric], "0.7f"
                ) + '</val>'
            else:
                ret_val += '<val>0.0</val>'
        ret_val += '</' + metric + '>\n' + '\t'*offset + '<identity>'
        for datapoint in dataseries:
            if (datapoint in tax_profile.tree.data[taxid].attributes) and (
                tax_profile.tree.data[taxid].attributes[datapoint]['hit_count'] > 0
            ):
                ret_val += '<val>' + format((
                    tax_profile.tree.data[taxid].attributes[datapoint]['identity']
                    / tax_profile.tree.data[taxid].attributes[datapoint]['hit_count']
                ), "0.1f") + '</val>'
            else:
                ret_val += '<val>0.0</val>'
        ret_val += '</identity>\n'
        # if not tax_profile.tree.data[taxid].children:
        gene_ids = set()
        for datapoint in tax_profile.tree.data[taxid].attributes:
            if 'genes' in tax_profile.tree.data[taxid].attributes[datapoint]:
                for gene_id in tax_profile.tree.data[taxid].attributes[datapoint]['genes'].split(
                        ' '
                ):
                    gene_ids.add(gene_id)
        ret_val += get_genes_xml(genes, sorted(gene_ids), dataseries, offset, metric)

    else:
        if metric != 'readcount':
            ret_val += '\t'*offset + '<readcount>'
            ret_val += '<val>0</val>'*len(dataseries)
            ret_val += '</readcount>\n'
        ret_val += '\t' * offset + '<' + metric + '>'
        ret_val += '<val>0.0</val>'*len(dataseries)
        ret_val += '</' + metric + '>\n' + '\t'*offset + '<identity>'
        ret_val += '<val>0.0%</val>'*len(dataseries)
        ret_val += '</identity>\n'

    if tax_profile.tree.data[taxid].children:
        for child_taxid in tax_profile.tree.data[taxid].children:
            ret_val += get_assembly_tax_xml(tax_profile, genes, dataseries, child_taxid, offset,
                                            metric)
    offset -= 1
    ret_val += '\t'*offset + '</node>\n'
    return ret_val


def make_assembly_taxonomy_chart(tax_profile, genes, function_list, outfile,
                                 krona_path, metric='efpkg'):
    """Writes XML file for taxonomy chart of assembly, one chart for all reads and separate charts
       for each function and generates Krona plot from it

    Args:
        tax_profile (:obj:TaxonomyProfile): taxonomy profile object
        genes (defaultdict[str,defaultdict[str,dict[str,float]]]): outer key is
            gene identifier, middle key is function identifier, inner key is in
            [metric, 'count', 'identity', 'coverage', 'Length', 'Completeness'],
            value is float (genes[gene_id][function_id][parameter_name] = parameter_value).
        function_list (list of str): function identifiers
        outfile (str): path for XML output
        krona_path (str): Krona Tools command
        metric (str): scoring metric (efpkg by default)
    """
    # genes contains gene data:, genes[gene_id][function][parameter] = parameter_value
    with open(outfile, 'w') as out:
        # Write header
        out.write('<krona key="false">\n')
        out.write('\t<attributes magnitude="' + metric + '">\n')
        if metric != 'readcount':
            out.write('\t\t<attribute display="Read count">readcount</attribute>\n')
        out.write('\t\t<attribute display="Score:' + metric + '">' + metric + '</attribute>\n')
        out.write('\t\t<attribute display="Coverage" mono="true">coverage</attribute>\n')
        out.write('\t\t<attribute display="Length" mono="true">Length</attribute>\n')
        out.write('\t\t<attribute display="CDS completeness %" mono="true">Completeness'
                  + '</attribute>\n')
        out.write('\t\t<attribute display="Best hit identity %" mono="true">identity</attribute>\n')
    # Obsolete
        out.write('\t\t<attribute display="UniRef hit" hrefbase="https://www.uniprot.org/uniref/" '
                  + 'target="uniref" mono="true">best_hit</attribute>\n')
        out.write('\t</attributes>\n')
        out.write('\t<color attribute="identity" valueStart="50" valueEnd="100" hueStart="0" '
                  + 'hueEnd="240" default="true"></color>\n')
        # Write dataset
        out.write('\t<datasets>\n')
        for function in function_list:
            out.write('\t\t<dataset>' + function + '</dataset>\n')
        out.write('\t</datasets>\n')
        # Write nodes
        offset = 1
        out.write(
            get_assembly_tax_xml(
                tax_profile, genes, function_list, ROOT_TAXONOMY_ID, offset, metric
                )
            )
        # Close XML
        out.write('</krona>')

    # Run Krona
    html_file = outfile + '.html'
    krona_cmd = [krona_path, '-o', html_file, outfile]
    run_external_program(krona_cmd)
