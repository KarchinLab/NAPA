# Row reading

# Column Reading #############################################################

# Small or intermediate size file
def readColumnFromShortFile(fileName, sep=",", colNum=0):
	""" Returns a list of nonempty strings in column number "colNum"
	in the file fileName"""
	columnStrings = []
	for (i, line) in enumerate(open(fileName,'rb')):
		try:
			columnString = each_line.split(sep)[colNum].strip()
		except IndexError:
			print 'Not enough columns on line %i of file.' % (i+1)
			continue
		if len(columnString): columnStrings.append(columnString)
	return columnStrings
		