import xml.etree.cElementTree as et
import copy
import numpy as np

import constants
#import optimizer
#import strength
import lamination
import materials
import plate_failure
import failure_envelope
import text_outputs

# xml-Core
# This module contains all the necessary functions to drive LaminateTools from XML files

def inputFromXML(input_file=None, xml_out=False, fea_out=False, text_out=True):
	# Collects the root tag of a file and builds the required analysis
	# objects.
	tree = et.parse(input_file)
	root = tree.getroot()

	# Figure out what to do with the input file, collect results
	if root.tag=='constants':
		# Run the engineering constants calculator
		analysis, envelopes, index = ConstantsAnalysis(root)
		if text_out:
			outfile = input_file+'_results.txt'
			text_outputs.FullOut(outfile,analysis, envelopes, index)
		if fea_out: pass
		if xml_out: pass
	elif root.tag=='optimize':
		# Run the optimization module
		analysis = OptimizeAnalysis(root)
	elif root.tag=='strength':
		# Run the failure testing module
		analysis = StrengthAnalysis(root)
	else:
		raise NameError('Analysis not defined')

	return [analysis, envelopes, index]

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
	loads_dict = dict()
	envelopes_dict = dict()
	analysis_dict = dict()
	index_dict = dict()


	# Gather materials
	for matl in root.findall('material'):
		matl_entry = {matl.attrib['name']:material(matl.attrib)}
		material_dict.update(matl_entry)
	for core in root.findall('weakcore'):
		matl_entry = {core.attrib['name']:materials.WeakCore()}
		material_dict.update(matl_entry)

	# Gather laminates, perform analysis and collect analysis objects
	for plybook in root.findall('plybook'):
		n_count = int(plybook.attrib['n'])
		symmetry = bool(plybook.attrib['s'])
		ply_stack = list()

		for ply in plybook.findall('ply'):
			plyDef = {'matl':material_dict[ply.attrib['material']],
				'orient':float(ply.attrib['orientation']),
				'thk':float(ply.attrib['thickness'])}
			ply_stack.append(lamination.Ply(plyDef))

		laminate = lamination.Laminate(copy.deepcopy(ply_stack),n_count,symmetry)
		print('Generating constants for: '+plybook.attrib['name'])
		analysis_entry = {plybook.attrib['name']:analysis_function(laminate)}
		analysis_dict.update(analysis_entry)

	if root.attrib['type']=='thinplate':
		# Gather load cases
		for case in root.findall('loads'):
			N = np.matrix([float(case.attrib['Nx']), float(case.attrib['Ny']), float(case.attrib['Nxy'])])
			M = np.matrix([float(case.attrib['Mx']), float(case.attrib['My']), float(case.attrib['Mxy'])])
			loads_entry = {case.attrib['name']:[N, M]}
			loads_dict.update(loads_entry)

		# Gather envelope cases, run and collect results
		for case in root.findall('failenvelope'):
			if case.attrib['failtype'] == 'MaxStressAny':
				kind = 'any'
				failtype = plate_failure.MaxStress
			elif case.attrib['failtype'] == 'MaxStressLong':
				kind = 'long'
				failtype = plate_failure.MaxStress
			elif case.attrib['failtype'] == 'MaxStressTrans':
				kind = 'trans'
				failtype = plate_failure.MaxStress
			elif case.attrib['failtype'] == 'MaxStressShear':
				kind = 'shear'
				failtype = plate_failure.MaxStress
			elif case.attrib['failtype'] == 'TsaiHill':
				kind = 'any'
				failtype = plate_failure.TsaiHill
			elif case.attrib['failtype'] == 'Hoffman':
				kind = 'any'
				failtype = plate_failure.Hoffman
			else:
				raise NameError('Failure index type not defined')

			if case.attrib['vartype']=='NM-LinVar':
				vartype = failure_envelope.makeNMLinVarEnvelope
			else:
				raise NameError('Failure envelope generator not defined')

			name = case.attrib['laminate']+'--'+case.attrib['loads']+'  '+case.attrib['failtype']+'/'+case.attrib['vartype']
			print('Generating envelope for: '+ name)
			N, M = loads_dict[case.attrib['loads']]
			N = np.matrix(N).T
			M = np.matrix(M).T
			analysis = analysis_dict[case.attrib['laminate']]
			result = vartype(N,M,analysis,failtype,kind)
			envelope_entry = {name:result}
			envelopes_dict.update(envelope_entry)

		for case in root.findall('failIndex'):
			if case.attrib['failtype'] == 'MaxStressAny':
				kind = 'any'
				failtype = plate_failure.MaxStress
			elif case.attrib['failtype'] == 'MaxStressLong':
				kind = 'long'
				failtype = plate_failure.MaxStress
			elif case.attrib['failtype'] == 'MaxStressTrans':
				kind = 'trans'
				failtype = plate_failure.MaxStress
			elif case.attrib['failtype'] == 'MaxStressShear':
				kind = 'shear'
				failtype = plate_failure.MaxStress
			elif case.attrib['failtype'] == 'TsaiHill':
				kind = 'any'
				failtype = plate_failure.TsaiHill
			elif case.attrib['failtype'] == 'Hoffman':
				kind = 'any'
				failtype = plate_failure.Hoffman
			else:
				raise NameError('Failure index type not defined')
			name = case.attrib['laminate']+'--'+case.attrib['loads']+'  '+case.attrib['failtype']+'Index'
			N, M = loads_dict[case.attrib['loads']]
			NM = np.matrix(N+M).T
			analysis = analysis_dict[case.attrib['laminate']]
			analysis.calculatePlyStressStrain(NM)
			result = failtype(analysis.laminate, kind)
			index_entry = {name:result}
			index_dict.update(index_entry)

	return [analysis_dict, envelopes_dict, index_dict]

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
