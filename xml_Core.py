import xml.etree.cElementTree as et
import copy
	tree = et.parse(input_file)
	root = tree.getroot()

	# Figure out what to do with the input file, collect results
	if root.tag=='constants':
		# Run the engineering constants calculator
		analysis = ConstantsAnalysis(root)


