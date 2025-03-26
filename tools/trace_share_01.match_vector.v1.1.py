## Surya Nagaraja
## Harvard University
## 2021-11-09

import os
import argparse
import itertools
import gzip
import Levenshtein

def reverseComplement(seq):
        complement = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
        bases = list(seq)
        comp = [complement[base] for base in seq]
        return ''.join(comp[::-1])

def fuzzySeqSearch(search_seq,read_seq,mismatch):
	for i, base in enumerate(read_seq):
		read_subset = read_seq[i:i+len(search_seq)]
		dist = Levenshtein.distance(read_subset,search_seq)
		if dist <= mismatch:
			return i, dist
			break

def barcodeSet(barcode):
	bases = 'ATCGN'
	bcSet = set()
	bcSet.add(barcode)
	for ind, b in enumerate(barcode):
		if b in bases:
			for base in bases:
				if barcode[ind] != base:
					bcSet.add(barcode[:ind]+base+barcode[ind+1:])
	return bcSet

def convertShareBC(seq,bcInds,bcLen,dictP5,dictR1,dictR2,dictR3):
	barcodeMatch = 0
	seqsplit = seq.split('+')
	
	sample = seqsplit[1]
	if sample in dictP5:
		sampBC = dictP5[sample]
		barcodeMatch += 1
	else:
		sampBC = sample
	cell = seqsplit[0]
	# check for exact match cell barcode
	if cell[bcInds[0]:bcInds[0]+bcLen] in dictR1:
	        barcode = dictR1[cell[bcInds[0]:bcInds[0]+bcLen]]
	        barcodeMatch += 1
	else:
	        barcode = cell[bcInds[0]:bcInds[0]+bcLen]
	if cell[bcInds[1]:bcInds[1]+bcLen] in dictR2:
	        barcode = barcode + ',' + dictR2[cell[bcInds[1]:bcInds[1]+bcLen]]
	        barcodeMatch += 1
	else:
	        barcode = barcode + ',' + cell[bcInds[1]:bcInds[1]+bcLen:]
	if cell[-bcLen:] in dictR3:
	        barcode = barcode + ',' + dictR3[cell[-bcLen:]]
	        barcodeMatch += 1
	else:
	        barcode = barcode + ',' + cell[-bcLen:]
	barcode = barcode + "," + sampBC
	return barcode,barcodeMatch

def extractVector4Reads(r1,r2,r3,r4,vect_seq,mm,out,sampleBC,pT,doPos,seqType,keepQname):
	print("Searching for sequence: "+vect_seq)
	print("Number of mismatches allowed: "+str(mm))

	if pT == 'check':
		print('PolyT checking mode detected')
	else:
		print ('Not checking UMI for polyT')
	
	if not keepQname:
		print('Replacing query names in reads')
	
	if sampleBC:
		print('Forcing all extracted sequences to have sample barcode: '+sampleBC)
	
	Read1 = gzip.open(r1,'rt')
	tot_reads = int(len(Read1.readlines())/4)
	Read1.close()
	print('Total reads: '+str(tot_reads))

	Read1 = gzip.open(r1,'rb')	
	Read2 = gzip.open(r2,'rb')
	Read3 = gzip.open(r3,'rb')
	Read4 = gzip.open(r4,'rb')
	
	out1 = gzip.open((out + ".matched.R1.fastq.gz"),'wb')
	out2 = gzip.open((out + ".matched.R2.fastq.gz"),'wb')
	if doPos:
		posFile = gzip.open((out + ".match_pos.txt.gz"),'wb')
	
	i = 0
	n = 0
	match = 0
	tot_match = 0
	polyT = 'TTTTTT'
	correctT = barcodeSet(polyT)
	
	# set up barcode conversion
	bcIdx = [15,53,91]
	bcLength = 8
	mmR1 = dict()
	mmR2 = dict()
	mmR3 = dict()
	mmP5 = dict()
	for barcode, name in R1.items():
        	barcodes = barcodeSet(barcode)
	        for bc in barcodes:
        	        mmR1[bc] = name
	for barcode, name in R2.items():
        	barcodes = barcodeSet(barcode)
	        for bc in barcodes:
        	        mmR2[bc] = name
	for barcode, name in R3.items():
        	barcodes = barcodeSet(barcode)
	        for bc in barcodes:
        	        mmR3[bc] = name
	for barcode, name in P5.items():
        	if seqType == "NextSeq":
                	barcodes = barcodeSet(barcode)
	        else:
        	        barcodes = barcodeSet(reverseComplement(barcode))
	        for bc in barcodes:
        	        mmP5[bc] = name
	
	# R1 = bio read 1, R2 = cell barcode, R3 = sample barcode, R4 = UMI + TTTTT
	for line in zip(Read1,Read2,Read3,Read4):
		if i == 0:
			# sequence header
			head1 = line[0].decode('utf-8').replace(" 1:",":1:").encode('utf-8')
			head2 = line[1].decode('utf-8').replace(" :",":2:").replace(" 2:",":4:").encode('utf-8')
			head3 = line[2].decode('utf-8').replace(" :",":3:").replace(" 3:",":4:").encode('utf-8')
			head4 = line[3].decode('utf-8').replace(" 2:",":4:").replace(" 4:",":4:").encode('utf-8')
		if i == 1:
			# sequence read
			seqread=line[0].decode('utf8').strip('\n')
			search = fuzzySeqSearch(vect_seq,seqread,mm)
			if search:
				idx,mis = search
				if (mis <= mm):
					if pT == 'check':
						if line[3][11:17] in correctT:
							match1 = 1
					else:
						match1 = 1
					if match1 == 1:
						vectsplit = seqread[idx+len(vect_seq):]
						if doPos:
        	                                        posFile.write((str(idx)+'\n').encode('utf-8'))
						cell = line[1].decode('utf8').strip('\n')
						if sampleBC:
							sample = sampleBC
						else:
							sample = line[2].decode('utf8').strip('\n')
						convertBC = convertShareBC((cell + "+" + sample),bcIdx,bcLength,mmP5,mmR1,mmR2,mmR3)
						if convertBC[1] == 4:
							match = 1
							tot_match += 1
							seqread2 = line[3].decode('utf8').strip('\n')
							umi = seqread2[0:10]
							restR2 = seqread2[10::]
							if keepQname:
								newHead1 = (head1.decode('utf-8').strip('\n') + "_" + convertBC[0] + "_" + umi + '\n')
								newHead2 = (head4.decode('utf-8').strip('\n') + "_" + convertBC[0] + "_" + umi + '\n')
							else:
								newHead1 = ( "@Vector:1:N:0_" + convertBC[0] + "_" + umi + '\n')
								newHead2 = ( "@Vector:2:N:0_" + convertBC[0] + "_" + umi + '\n')
							out1.write(newHead1.encode('utf-8'))
							out1.write((vectsplit+'\n').encode('utf-8'))
							out2.write(newHead2.encode('utf-8'))
							out2.write((restR2+'\n').encode('utf-8'))
						
		if i == 2:
			# quality header
			if match == 1:
				out1.write(line[0])
				out2.write(line[3])
		if i == 3:
			# quality read
			if match == 1:
				trim1 = line[0].decode('utf-8').strip('\n')
				trim2 = line[3].decode('utf-8').strip('\n')
				out1.write((trim1[idx+len(vect_seq):] + '\n').encode('utf-8'))
				out2.write((trim2[10::] + '\n').encode('utf-8'))
		i += 1
		if i == 4:
			i = 0 
			n += 1
			match = 0
			if (n % 500000 == 0):
				print(str(n)+" out of "+str(tot_reads)+" parsed") 
	Read1.close()
	Read2.close()
	Read3.close()
	Read4.close()
	out1.close()
	out2.close()
	print(str(tot_match)+" total reads with sequence accepted")

def extractVector2Reads(r1,r2,vect_seq,mm,out,doPos):
	print("Two reads detected - assuming processed with SHARE pipeline. Not converting cell barcode or pulling UMIs.")
	print("Searching for sequence: "+vect_seq)
	print("Number of mismatches allowed: "+str(mm))

	Read1 = gzip.open(r1,'rt')
	tot_reads = int(len(Read1.readlines())/4)
	Read1.close()
	print('Total reads: '+str(tot_reads))

	Read1 = gzip.open(r1,'rb')	
	Read2 = gzip.open(r2,'rb')
    	
	out1 = gzip.open((out + ".matched.R1.fastq.gz"),'wb')
	out2 = gzip.open((out + ".matched.R2.fastq.gz"),'wb')
	if doPos:
		posFile = gzip.open((out + ".match_pos.txt.gz"),'wb')
	
	i = 0
	n = 0
	match = 0
	tot_match = 0
	polyT = 'TTTTTT'
	correctT = barcodeSet(polyT)
	
	for line in zip(Read1,Read2):
		if i == 0:
			# sequence header
			head1 = line[0]
			head2 = line[1]
		if i == 1:
			# sequence read
			seqread=line[0].decode('utf8').strip('\n')
			search = fuzzySeqSearch(vect_seq,seqread,mm)
			if search:
				idx,mis = search
				if (mis <= mm):
					match = 1
					tot_match += 1
					vectsplit = seqread[idx+len(vect_seq):]
					if doPos:
						posFile.write((str(idx)+'\n').encode('utf-8'))
						out1.write(head1)
						out1.write((vectsplit+'\n').encode('utf-8'))
						out2.write(head2)
						out2.write(line[1])
						
		if i == 2:
			# quality header
			if match == 1:
				out1.write(line[0])
				out2.write(line[1])
		if i == 3:
			# quality read
			if match == 1:
				trim1 = line[0].decode('utf-8').strip('\n')
				out1.write((trim1[idx+len(vect_seq):] + '\n').encode('utf-8'))
				out2.write(line[1])
		i += 1
		if i == 4:
			i = 0 
			n += 1
			match = 0
			if (n % 500000 == 0):
				print(str(n)+" out of "+str(tot_reads)+" parsed") 
	Read1.close()
	Read2.close()
	out1.close()
	out2.close()
	print(str(tot_match)+" total reads with sequence accepted")

def main():
	parser = argparse.ArgumentParser(
		description="Extract SHARE-seq fastq reads with vector sequence",
		epilog="Requires R1-R4 with biological read in R1")
	parser.add_argument(
		"-R1",
		metavar="Read 1",
		required=True,
		help="Path to the Read 1 fastq file")
	parser.add_argument(
		"-R2",
		metavar="Read 2",
		required=True,
		help="Path to the Read 2 fastq file")
	parser.add_argument(
		"-R3",
		metavar="Read 3",
		default=None,
		help="Path to the Read 3 fastq file")
	parser.add_argument(
		"-R4",
		metavar="Read 4",
		default=None,
		help="Path to the Read 4 fastq file")
	parser.add_argument(
		"-v",
		metavar="Vector sequence",
		required=True,
		help="String to search fastq for")
	parser.add_argument(
		"-m",
		metavar="Number of mismatches",
		required=True,
		type=int,
		help="Mismatch allowed in vector seq search")
	parser.add_argument(
		"--out",
		metavar="Output prefix",
		required=True,
		help="Label for vector read file")
	parser.add_argument(
		"--forceSample",
		metavar="Sample barcode",
		help="Force sample barcode ID on output files")
	parser.add_argument(
		"--checkPolyT",
		nargs='?',
		default="noCheck",
		const='check',
		help="Check for TTTTTT after UMI to accept read")
	parser.add_argument(
		"--writePosition",
		metavar="Write position",
		default=True,
		type=bool,
		help="Whether to write position of sequence match to file")
	parser.add_argument(
                "--sequencer",
                metavar="Type of sequencer",
                default="NovaSeq",
                help="If 'NextSeq', will check exact match for Index2. Otherwise, will check reverse complement.")
	parser.add_argument(
                "--keepQname",
                metavar="Keep query name in output",
		dest = 'keepQname',
                help="Giving --keepQname will keep query name, default is to remove")
	parser.set_defaults(keepQname=False)
	args = parser.parse_args()
	if ((args.R3 is None) or (args.R4 is None)):
		extractVector2Reads(args.R1,args.R2,args.v,args.m,args.out,args.writePosition)
	else:
		extractVector4Reads(args.R1,args.R2,args.R3,args.R4,args.v,args.m,args.out,args.forceSample,args.checkPolyT,args.writePosition,args.sequencer,args.keepQname)


from splitBarcodes2_1 import R1
from splitBarcodes2_1 import R2
from splitBarcodes2_1 import R3
from splitBarcodes2_1 import P5

main()
