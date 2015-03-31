import xml.etree.cElementTree as et
import copy
import text_outputs as tOut
import constants
import optimizer
#import strength
import lamination as lam
import materials as matls

# xml-Core
# This module contains all the necessary functions to drive LaminateTools from XML files

def InputCollector(input_file=None, xml_out = True, fea_out = False, print_out = True):
	# Collects the root tag of a file and determines what module and analysis to use from thers.
	tree = et.parse(input_file)
	root = tree.getroot()

	# Figure out what to do with the input file, collect results
	if root.tag is 'constants':
		# Run the engineering constants calculator
		analysis = constantsAnalysis(root)
	elif root.tag is 'optimize':
		# Run the optimization module
		analysis = optimizeAnalysis(root)
	elif root.tag is 'strength':
		# Run the failure testing module
		analysis = optimizeAnalysis(root)
	else:
		raise NameError('Analysis not defined')

	# Figure out how to output the results
	if xml_out:
		outputCollector(input_file, results, analysis)
	if fea_out:
		tOut.FeaOut(analysis)
	if print_out:
		print(analysis)



def OutputCollector(input_file=None, results=None):
	# Collect results from analysis and append to the original XML file.

def ConstantsAnalysis(root=None):
	# Parse a given XML tree and run an engineering constants analysis.
	# Figure out the type and set the correct function calls
	if root.attrib['type'] is 'continuum':
		material = matls.continuum
		analysis_function = constants.continuum
	elif root.attrib['type'] is 'plate':
		material = matls.plate
		analysis_function = constants.plate
	elif root.attrib['type'] is 'beam':
		material = matls.beam
		analysis_function = constants.beam
	else:
		raise NameError('Constants analysis type not defined')

	# Initialize material dictionary, laminate dictionary and results dictionary
	material_dict = dict()
	laminate_dict = dict()
	results = dict()

	# Gather materials
	for matl in root.findall('material'):
		matl_entry = {matl.attrib['name'],material(matl.attrib)}
		material_dict.update(matl_entry)

	#Gather laminates, perform analysis and collect analysis objects
	for plybook in root.findall('plybook'):
		n_count = int(plybook.attrib['n'])
		symmetry = bool(plybook.attrib['s'])
		ply_stack = list()

		for ply in plybook.findall('ply'):
			material = Material_Dict[ply.attrib['material']]
			orientation = float(ply.attrib['orientation'])
			thickness = float(ply.attrib['thickness'])
			ply_stack.append(lam.Ply(material,orientation,thickness))

		laminate = lam.Laminate(copy.deepcopy(ply_stack),n_count,symmetry)
		analysis_entry = analysis_function(laminate)
		results.update({plybook.attrib['name']:analysis_entry})

	return results

def OptimizeAnalysis(root=None):
	# Parse a given XML tree and run an optimization analysis.
	analysis = None

	return analysis

def StrengthAnalysis(root=None):
	# Parse a given XML tree and run a failure analysis.
	analysis = None

	return analysis
