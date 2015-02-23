import import_export as inex
from laminate_properties import PlateLaminate


print('LaminateTools v.0.1.0\n (c) Matthew Cowan 2015 \n Released under the Apache v2 License \n')

fileIn = raw_input('Enter laminate XML file name:\n')
dataIn = inex.fromXML(laminate_file=fileIn)

for lam in dataIn['lams'].keys():
	print('Processing laminate: '+lam)
	#print(str(dataIn['lams'][lam]))
	nameout = lam + '_properties.txt'
	props = PlateLaminate(dataIn['lams'][lam])
	inex.toText(props, nameout)

print('Analysis Complete')
