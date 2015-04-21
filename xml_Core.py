import xml.etree.cElementTree as et
import copy
import constants
#import optimizer
#import strength
import lamination
import materials

# xml-Core
# This module contains all the necessary functions to drive LaminateTools from XML files

def inputFromXML(input_file=None, xml_out=True, fea_out=False, text_out=True):
	# Collects the root tag of a file and builds the required analysis
	# objects.
	tree = et.parse(input_file)
	root = tree.getroot()

	# Figure out what to do with the input file, collect results
	if root.tag=='constants':
		# Run the engineering constants calculator
		analysis = ConstantsAnalysis(root)
		# text_out =
		# fea_out = constants.feaOut
		# xml_out = constants.xmlOut
	elif root.tag=='optimize':
		# Run the optimization module
		analysis = OptimizeAnalysis(root)
	elif root.tag=='strength':
		# Run the failure testing module
		analysis = StrengthAnalysis(root)
	else:
		raise NameError('Analysis not defined')

	return analysis

def outputToXML(input_fule=None, results=None):
	# Collect results (sub)trees and append them to the given file.
	pass

def ConstantsAnalysis(root=None):
	# Parse a given XML tree and run an engineering constants analysis.
	# Figure out the type and set the correct function calls
	if root.attrib['type']=='continuum':
		material = materials.Continuum
		analysis_function = constants.Continuum
	elif root.attrib['type']=='thinplate':
		material = materials.Plate
		analysis_function = constants.ThinPlate
	elif root.attrib['type']=='thickplate':
		material = materials.Continuum
		analysis_function = constants.ThickPlate
	elif root.attrib['type']=='beam':
		material = materials.Beam
		analysis_function = constants.Beam
	else:
		raise NameError('Constants analysis type not defined')

	# Initialize material, laminate, analysis dictionaries
	material_dict = dict()
	laminate_dict = dict()
	analysis = dict()

	# Gather materials
	for matl in root.findall('material'):
		matl_entry = {matl.attrib['name']:material(matl.attrib)}
		material_dict.update(matl_entry)

	#Gather laminates, perform analysis and collect analysis objects
	for plybook in root.findall('plybook'):
		n_count = int(plybook.attrib['n'])
		symmetry = bool(plybook.attrib['s'])
		ply_stack = list()

		for ply in plybook.findall('ply'):
			material = material_dict[ply.attrib['material']]
			orientation = float(ply.attrib['orientation'])
			thickness = float(ply.attrib['thickness'])
			ply_stack.append(lamination.Ply(material,orientation,thickness))

		laminate = lamination.Laminate(copy.deepcopy(ply_stack),n_count,symmetry)
		analysis_entry = analysis_function(laminate)
		analysis.update({plybook.attrib['name']:analysis_entry})

	return analysis

def OptimizeAnalysis(root=None):
	# Parse a given XML tree and run an optimization analysis.
	analysis = None

	return analysis

def StrengthAnalysis(root=None):
	# Parse a given XML tree and run a failure analysis.
	analysis = None

	return analysis

## BUILT IN TESTING
if __name__=='__main__':

	testRuns = inputFromXML(input_file='constantsTestInput.xml',
	xml_out = False, fea_out = False, text_out = False)

	for key in testRuns:
		print('---- '+key+' ----')
		print(str(testRuns[key]))
