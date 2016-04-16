from Bio import Entrez, SeqIO
import basicFileIO as bFIO

def parse_arguments():
	"""
	Builds the argument parser object
	"""
	parser = argparse.ArgumentParser(prog="genBank2dnaSequence.py", 
		description="get Genbank DNA sequences from a set of Genbank IDs", 
		usage="%(prog)s [options] GenbankIDfile FastaOutFile")
	parser.add_argument("-f", dest="gbid_file"
		help="The file containing the Genbank IDs")
	parser.add_argument("-d", dest="delim", 
		default=",", help="File delimiter for GenbankIDfile")
	parser.add_argument("-c", dest="gbid_col", 
		default=",", help="The column number where the gbIds are stored")
	parser.add_argument("-o", dest="fasta_out", 
		default=",", help="The output fasta file for the sequences obtained")
	args = parser.parse_args()
    return args	


def gbRecEfetch(gbId,DB='gb',RetType):
	gbId = gbId.split('.')[0].split('.')[0].strip()
	
	record = None
	for trialIdx in range(5):
		try:
			handle = Entrez.efetch(db=DB,id=gbId,rettype=RetType)
		except urllib2.HTTPError:
			print "HTTP Error 502: Bad Gateway -- Empty handle"
			continue
		handle_len = len(str(handle.read()).strip().split())
		if handle_len > 6: # Check that the handle was full
			handle = Entrez.efetch(db=DB,id=gbId, rettype=RetType)
			if handle: 
				try:
					record = SeqIO.read(handle,RetType)
				except ValueError:
					print "No records found in handle from user DB; Trial#", trialIdx
					continue
		if record: break
	return record
	
def genBankIds2fastaDictionary(genBankIds): 
	for gbId in genBankIds:
		record = gbRecEfetch(gbId,DB='gb',RetType)
	if fRec:
		return fRec.seq
	else:
		defRec = SeqRecord(Seq('', IUPAC.protein))
		return defRec.seq

gb_file = "NC_005213.gbk"
for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
    # now do something with the record
    print "Name %s, %i features" % (gb_record.name, len(gb_record.features))
    print str(gb_record.seq)
    
######################################################################
def main():
		args = parse_arguments()
		genBankIds = readColumnFromShortFile(args.gbid_file, sep=args.delim, colNum=args.gbid_col)
 		
if __name__ = "__main__": main()
