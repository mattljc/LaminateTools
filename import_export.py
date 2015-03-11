import xml.etree.ElementTree as et
import copy
from composite_materials import *
from ply_stack import *
from laminate_properties import *

def fromXML(laminate_file=None):
	tree = et.parse(laminate_file)
	root = tree.getroot()

	# Initialize dictionaries
	Material_Dict = dict()
	Laminate_Dict = dict()
	# Determine analysis type and build material dictionaries as appropriate
	if (root.tag == 'continuum'):
		# A 3d analysis with appropriate materials
		for material in root.findall('materialProps'):
			matlName = material.attrib["name"]
			E11 = float(material.attrib["E11"])
			E22 = float(material.attrib["E22"])
			E33 = float(material.attrib["E33"])
			nu12 = float(material.attrib["nu12"])
			nu13 = float(material.attrib["nu13"])
			nu23 = float(material.attrib["nu23"])
			G12 = float(material.attrib["G12"])
			G13 = float(material.attrib["G13"])
			G23 = float(material.attrib["G23"])
			CPT = float(material.attrib["CPT"])
			arealDens = float(material.attrib["arealDensity"])
			Material_Dict.update({matlName : ContinuumMaterial(name=matlName, E11_in=E11, E22_in=E22, E33_in=E33, Nu12_in=nu12, Nu13_in=nu13, Nu23_in=nu23, G12_in=G12, G13_in=G13, G23_in=G23, ArealDensity_in=arealDens, CPT=0)})

	elif (root.tag == 'plate'):
		# A plate analysis with appropriate materials
		for material in root.findall('materialProps'):
			matlName = material.attrib["name"]
			E11 = float(material.attrib["E11"])
			E22 = float(material.attrib["E22"])
			nu12 = float(material.attrib["nu12"])
			G12 = float(material.attrib["G12"])
			CPT = float(material.attrib["CPT"])
			arealDens = float(material.attrib["arealDensity"])
			Material_Dict.update({matlName : PlateMaterial(name=matlName, E11_in=E11, E22_in=E22, Nu12_in=nu12, G12_in=G12, ArealDensity_in=arealDens, CPT_in=CPT)})

	elif (root.tag == 'beam'):
		# A beam analysis, which is scheduled for future implementation
		raise TypeError('Beam analysis not implemented')

	else:
		raise TypeError('Invalid Analysis Type')

	for plybook in root.findall('plybook'):
		lamName = plybook.attrib['name']
		n_count = int(plybook.attrib['n_count'])
		sym = bool(int(plybook.attrib['sym']))
		ply_stack = list()
		for ply in plybook.findall('ply'):
			material = Material_Dict[ply.attrib['material']]
			orientation = float(ply.attrib['orientation'])
			thickness = float(ply.attrib['thickness'])
			ply_stack.append(Ply(material,orientation,thickness))
		Laminate_Dict.update({lamName:Laminate(copy.deepcopy(ply_stack),n_count,sym)})

	return {'type':root.tag, 'mats':Material_Dict, 'lams':Laminate_Dict}

def toXML(lamProps, properties_file='laminate.xml'):
	# Outputs the implementable details of a laminate with its properties written in a comment

	return None

def toText(lamProps, properties_file='laminate_properties.txt'):
	file = open(properties_file,'w')
	file.write(lamProps.__str__())
	file.write(lamProps.__repr__())
	file.close()

def toNASTRAN(type='continuum', properties_file='laminate_properties.mat9'):
	# Determines the appropriate material property deck to use and outputs the appropriate NASTRAN data file
	return None

def toABAQUS(type='continuum',properties_file='laminate_properties.ipt'):
	# Determines the appropriate material property deck to use and outputs the appropriate ABAQUS/Simulia data file
	return None
