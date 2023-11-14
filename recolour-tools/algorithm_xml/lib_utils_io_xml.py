"""
Library Features:

Name:          lib_utils_io_xml
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20230719'
Version:       '1.0.0'
"""
import os.path

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import re
import xml.etree.ElementTree as ET

from copy import deepcopy
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to read file xml
def read_file_xml(file_name):
    xml_tree = ET.parse(file_name)
    xml_root = xml_tree.getroot()
    return xml_root, xml_tree
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to write xml file
def write_file_xml(file_name, xml_data):

    # delete previous file
    if os.path.exists(file_name):
        os.remove(file_name)

    # set xml namespace(s)
    ET.register_namespace("gmd", "http://www.isotc211.org/2005/gmd")
    ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")
    ET.register_namespace("gco", "http://www.isotc211.org/2005/gco")
    ET.register_namespace("srv", "http://www.isotc211.org/2005/srv")
    ET.register_namespace("gmx", "http://www.isotc211.org/2005/gmx")
    ET.register_namespace("gts", "http://www.isotc211.org/2005/gts")
    ET.register_namespace("gsr", "http://www.isotc211.org/2005/gsr")
    ET.register_namespace("gmi", "http://www.isotc211.org/2005/gmi")
    ET.register_namespace("gml", "http://www.opengis.net/gml")
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    ET.register_namespace(
        "schemaLocation",
        "http://www.isotc211.org/2005/gmd http://schemas.opengis.net/iso/19139/20060504/gmd/gmd.xsd")

    # dump xml data
    if xml_data is not None:
        xml_tree = ET.ElementTree(xml_data)
        xml_tree.write(file_name, encoding="UTF-8", xml_declaration=True)
    else:
        logging.error(' ===> XML data is defined by NoneType')
        raise RuntimeError('XML date must be defined by element objects')
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to update xml file
def update_file_xml(xml_root, xml_keys=None, xml_values=None, xml_envs=None,
                    xml_field_error_update=True,
                    xml_field_not_declared_text=False, xml_field_not_declared_attributes=False):

    # default xml keys and value
    if xml_keys is None:
        xml_keys = {'field_1': {
            'description': 'test', 'active': True,
            'node': [
                "MD_Metadata", "identificationInfo", "MD_DataIdentification", "citation",
                "CI_Citation",
                "date", "CI_Date", "date", 'CharacterString'], 'element': 1}}
    if xml_values is None:
        xml_values = {'field_1': {
            "element_1": {
                "type": "ENV", "variable": "ENV_FILENAME_PRODUCT",
                "value": "file_name_UPD.tiff", "attrib": {}}}}

    # manage root tag
    xml_root_tag, xml_root_attr = xml_root.tag, xml_root.attrib
    xml_root_string = re.search('{(.*)}', xml_root_tag).group(1)
    xml_root_string = '{' + xml_root_string + '}'

    # iterate over field(s)
    for xml_id, xml_fields in xml_keys.items():

        if xml_id in list(xml_values.keys()):
            xml_elements = xml_values[xml_id]
        else:
            logging.error(' ===> XML fields and values must be defined by the same key(s)')
            raise RuntimeError('The XML key "' + xml_id + '" is not available in the XML values obj')

        # get field info
        field_description, field_active = xml_fields['description'], xml_fields['active']
        field_node_list, field_element = xml_fields['node'], xml_fields['element']

        # info field start
        logging.info(' -----> Update field -- ID: "' + xml_id + '" -- Description: "' + field_description +
                     '" ... ')

        # flag field activation
        if field_active:

            xml_value_list, xml_attrs_list = [], []
            for element_id, element_fields in xml_elements.items():

                # get element info
                element_type, element_variable = element_fields['type'], element_fields['variable']
                element_value, element_attrib = element_fields['value'], element_fields['attrib']

                if element_type == 'ENV':
                    if element_variable in xml_envs.keys():
                        xml_value_step = xml_envs[element_variable]
                    else:
                        logging.warning(' ===> Variable "' + element_variable + '" is expected from envs')
                        xml_value_step = deepcopy(element_value)
                else:
                    xml_value_step = deepcopy(element_value)

                xml_value_list.append(xml_value_step)
                xml_attrs_list.append(element_attrib)

            # iterate over xml nodes
            xml_root = iterate_nodes(deepcopy(xml_root),
                                     xml_child_list=field_node_list,
                                     xml_value_list=xml_value_list, xml_attrs_list=xml_attrs_list,
                                     xml_element=field_element,
                                     xml_error_update=xml_field_error_update)

            '''
            # debug control
            tree_tmp = ET.ElementTree(xml_root)
            tree_tmp.write("debug.xml", encoding="UTF-8", xml_declaration=True)
            '''

            # info field end
            logging.info(' -----> Update field -- ID: "' + xml_id + '" -- Description: "' + field_description +
                         '" ... DONE')

        else:
            # info field end
            logging.info(' -----> Update field "' + field_description + '" ... SKIPPED. Field is not activated')

    return xml_root
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to iterate over the nodes
def iterate_nodes(xml_obj, xml_child_list=None, xml_value_list=None, xml_attrs_list=None, xml_element=1,
                  xml_tail=2,
                  xml_error_update=True, xml_not_declared_text=False, xml_not_declared_attributes=False):

    if xml_child_list is None:
        xml_child_list = []
    if xml_value_list is None:
        xml_value_list = []
    if xml_attrs_list is None:
        xml_attrs_list = []

    if not isinstance(xml_value_list, list):
        xml_value_list = [xml_value_list]
    if not isinstance(xml_attrs_list, list):
        xml_attrs_list = [xml_attrs_list]

    xml_update_flag, xml_update_iter = False, 0
    xlm_list_tag, xlm_list_n = [], 0
    for xml_child in xml_obj.iter():

        xml_child_tag, xml_child_attr, xml_child_txt = xml_child.tag, xml_child.attrib, xml_child.text
        xlm_string_tag = re.sub('{.*?}', '', xml_child_tag)

        if xlm_string_tag in xml_child_list:

            if not xlm_list_tag:
                xml_child_idx = xml_child_list.index(xlm_string_tag)
                xlm_list_tag = xml_child_list[:xml_child_idx]

            xlm_list_tag.append(xlm_string_tag)

            # check elements (unique and not order)
            if set(xlm_list_tag) == set(xml_child_list):

                # update iteration(s)
                xml_update_iter = xml_update_iter + 1
                logging.info(' ------- Iteration: "' + str(xml_update_iter) + '"')

                # elements equal
                if xlm_list_tag == xml_child_list:

                    logging.info(' ------- Tags and child list are the same - Update the field')

                    xml_child = update_nodes(
                        xml_child, xml_value_list, xml_attrs_list, xlm_list_n,
                        xml_not_declared_text=xml_not_declared_text,
                        xml_not_declared_attributes=xml_not_declared_attributes)

                    xlm_list_tag = []
                    xlm_list_n = xlm_list_n + 1
                    if xlm_list_n == xml_element:
                        xml_update_flag = True
                        break

                elif xlm_list_tag[-xml_tail:] == xml_child_list[-xml_tail:]:

                    logging.info(' ------- Tags and child list have the same tail --- try to update the field')

                    xml_child = update_nodes(
                        xml_child, xml_value_list, xml_attrs_list, xlm_list_n,
                        xml_not_declared_text=xml_not_declared_text,
                        xml_not_declared_attributes=xml_not_declared_attributes)

                    xlm_list_tag = []
                    xlm_list_n = xlm_list_n + 1
                    if xlm_list_n == xml_element:
                        xml_update_flag = True
                        break

                else:

                    string_list_tag = ','.join(xlm_list_tag)
                    string_list_child = ','.join(xml_child_list)
                    logging.warning(' ===> Elements are the same the nodes are not correctly ordered')
                    logging.warning(' XML SELECTED: "' + string_list_tag + '"')
                    logging.warning(' XML EXPECTED: "' + string_list_child + '"')

    if not xml_update_flag:
        if xml_error_update:
            logging.error(' ===> Update field failed due to the wrong list of the node. Check the correct list nodes')
            raise RuntimeError('Field is not updated; check the configuration file to correctly set the updating')
        else:
            logging.warning(' ===> Update field failed due to the wrong list of the node. Check the correct list nodes')

    return xml_obj
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to update node(s)
def update_nodes(xml_child, xml_value_list, xml_attrs_list, xlm_list_n,
                 xml_not_declared_text=False, xml_not_declared_attributes=False):

    # info element id
    logging.info(' -------- Element: ' + str(xlm_list_n + 1))

    # info text update start
    xlm_element_text_in = xml_child.text

    if xml_not_declared_text:
        xlm_element_text_in = ''

    # check text format
    if xlm_element_text_in is not None:
        logging.info(' --------- (1a) Get Text: "' + xlm_element_text_in + '"')
        # get text selected
        xml_value_n = xml_value_list[xlm_list_n]
        # set text selected
        xml_child.text = xml_value_n
        # info text update end
        xlm_element_text_out = xml_child.text
        logging.info(' --------- (2a) Update Text: "' + xlm_element_text_out + '"')
    else:
        logging.info(' --------- (1a) Get Text: "NoneType"')
        logging.info(' --------- (2a) Update Text: SKIPPED')

    # info attributes update start
    xlm_element_attrib_in = xml_child.attrib
    if xml_not_declared_attributes:
        xlm_element_attrib_in = {}

    # check attributes format
    if xlm_element_attrib_in is not None:

        logging.info(' --------- (1b) Get Attributes: "' + str(xlm_element_attrib_in) + '"')
        # get attributes selected
        xml_attrs_n = xml_attrs_list[xlm_list_n]
        # set attributes selected
        if (not xml_attrs_n) and xlm_element_attrib_in:
            logging.warning(' ===> Attributes by the settings file are empty; '
                            'attributes in the xml default file are not empty. Keep the attributes of the default xml. ')
            xml_attrs_n = deepcopy(xlm_element_attrib_in)

        xml_child.attrib = xml_attrs_n
        # info attributes update end
        xlm_element_attrib_out = xml_child.attrib
        logging.info(' --------- (2b) Update Attributes: "' + str(xlm_element_attrib_out) + '"')

    else:
        logging.info(' --------- (1a) Get Attributes: "NoneType"')
        logging.info(' --------- (2a) Update Attributes: SKIPPED')

    return xml_child
# ----------------------------------------------------------------------------------------------------------------------
