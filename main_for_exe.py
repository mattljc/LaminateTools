import import_export as inex
from laminate_properties import *


print('LaminateTools v.0.1.0\n (c) Matthew Cowan 2015 \n Released under the Apache v2 License \n')

fileIn = raw_input('Enter laminate XML file name:\n')
dataIn = inex.fromXML(laminate_file=fileIn)

if (dataIn['type'] == 'continuum'):
	def analysis(lam_in):
		return Laminate3D(lam=lam_in)
elif (dataIn['type'] == 'plate'):
	def analysis(lam_in):
		x = Laminate2D(lam=lam_in)
		try:
			x.getEffectiveProperties()
		except AssertionError:
			print('Laminate not symmetric. Effective properties will not be included.')
		return x
elif (dataIn['type'] == 'beam'):
	def analysis(lam_in):
		return Laminate1D(lam=lam_in)

for lam in dataIn['lams'].keys():
	print('Processing laminate: '+lam)
	#print(str(dataIn['lams'][lam]))
	nameout = lam + '_properties.txt'
	props = analysis(lam_in=dataIn['lams'][lam])
	inex.toText(props, nameout)

print('Analysis Complete')
