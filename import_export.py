import xml.etree.ElementTree as et
import copy
from composite_materials import *
from ply_stack import *
from laminate_properties import *

def fromXML(laminate_file=None):
	tree = et.parse(laminate_file)
	root = tree.getroot()

	Material_Dict = dict()
	#Pull in materials and build material dictionary
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
		Material_Dict.update({matlName:RealCompositeMaterial(matlName, E11, E22, E33, nu12, nu13, nu23, G12, G13, G23)})

	Laminate_Dict = dict()
	for plybook in root.findall('plybook'):
		lamName = plybook.attrib['name']
		n_count = int(plybook.attrib['n_count'])
		sym = bool(plybook.attrib['sym'])
		ply_stack = list()
		for ply in plybook.findall('ply'):
			material = Material_Dict[ply.attrib['material']]
			orientation = float(ply.attrib['orientation'])
			thickness = float(ply.attrib['thickness'])
			ply_stack.append(Ply(material,orientation,thickness))

		Laminate_Dict.update({lamName:Laminate(copy.deepcopy(ply_stack),n_count,sym)})

	return (Material_Dict, Laminate_Dict)


def fromKeyboard():
	return None

def toXML(properties_file=None):
	return None

def toNASTRAN(properties_file=None):
	return None

def toABAQUS(properties_file=None):
	return None
