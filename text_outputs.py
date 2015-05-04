import numpy as np

def FullOut(outfile='results.txt',analysis={}, envelopes={}, index={}):

	# Set global output options
	np.set_printoptions(precision=5, linewidth=1000)

	contents = outfile +'\n\n'
	contents += '==LAMINATE CONSTANTS==\n'
	for key in analysis:
		contents += key+'\n'
		contents += str(analysis[key])
		contents += '\n'

	contents += '==LAMINATE FAILURE INDICES==\n'
	for key in index:
		contents += key+'\n'
		contents += str(index[key])
		contents += '\n'

	contents += '==LAMINATE FAILURE ENVELOPES==\n'
	for key in envelopes:
		contents += key + '\n'
		for subkey in envelopes[key]:
			contents += subkey + ' = '+str(envelopes[key][subkey])
			contents += '\n'
		contents += '\n'

	file = open(outfile,'w')
	file.write(contents)
	file.close()
