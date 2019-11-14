# Author: Peter Finn C16407722
import sys

# this is a list of the start and stop codons
startCod=["ATG"]
stopCod=["TGA"]

# The following is how the program finds the start of the sequence
# it looks for anything with ATCG, adds it to an array and if another letter shows up
# that is not A,T,C or G then it will clear the string and start again
lookFor="ATCGatcg"

# This is a very simple function to translate the DNA to codons.
# It does this one at a time. It cannot take the entire string.
def translateAA(sequence):
	CodonDictionary={
	'ATT':'I',   'ATC':'I',  'ATA':'I',  'CTT':'L',  'CTC':'L',  
	'CTA':'L',  'CTG':'L',  'TTA':'L',  'TTG':'L',  'GTT':'V',  'GTC':'V',  
	'GTA':'V',  'GTG':'V',  'TTT':'F',  'TTC':'F',  'ATG':'M',  'TGT':'C',  
	'TGC':'C',  'GCT':'A',  'GCC':'A',  'GCA':'A',  'GCG':'A',  'GGT':'G',  
	'GGC':'G',  'GGA':'G',  'GGG':'G',  'CCT':'P',  'CCC':'P',  'CCA':'P',  
	'CCG':'P',  'ACT':'T',  'ACC':'T',  'ACA':'T',  'ACG':'T',  'TCT':'S',  
	'TCC':'S',  'TCA':'S',  'TCG':'S',  'AGT':'S',  'AGC':'S',  'TAT':'Y',  
	'TAC':'Y',  'TGG':'W',  'CAA':'Q',  'CAG':'Q',  'AAT':'N',  'AAC':'N',  
	'CAT':'H',  'CAC':'H',  'GAA':'E',  'GAG':'E',  'GAT':'D',  'GAC':'D',  
	'AAA':'K',  'AAG':'K',  'CGT':'R',  'CGC':'R',  'CGA':'R',  'CGG':'R',  
	'AGA':'R',  'AGG':'R',  'TAA':'X',  'TAG':'X',  'TGA':'X'}
	return(CodonDictionary[sequence])

# This calculates the results and prints them to stdio.
# It recieves the frame number and start and stop codons
def printResults(start,stop,frame, length):

	framePrint=" "+str(frame+1)

# This checks if its a + or - frame. The maths and printing are differant in both.
	if frame > 2:
		frame-=3
		framePrint="-"+str(frame+1)
		start =start * 3 -2
		stop  = stop * 3 + frame
		return("Frame number: "+framePrint+" Start: "+str(-(start-length*3))+" Stop: "+str(-(stop-length*3)))
	else:
		start =start * 3 + frame +1
		stop  = stop * 3 + frame +3
		return("Frame number: "+framePrint+" Start: "+str(start)+" Stop: "+str(stop))


def framescanner(start,sequence,startCod,stopCod,negitiveNum):
	minus=""
	if negitiveNum==True:
		minus="-"

	print(textColors.red+"Frame Number: "+minus+str(start+1)+textColors.end)
	codonStr=""

	for x in range(start,len(sequence)-2,3):
		current=""
		current+=sequence[x]
		current+=sequence[x+1]
		current+=sequence[x+2]
		trigger=False

		for start in startCod:
			if start == current:
				#This is a start codon. Will add + to the codon string
				sys.stdout.write(textColors.green+current+textColors.end+" ")
				codonStr+="+"
				trigger=True
				break

		if trigger==False:
			for stop in stopCod:
				if stop == current:
					#This is a stop codon. Will add - to the codon string
					sys.stdout.write(textColors.red+current+textColors.end+" ")
					codonStr+="-"
					trigger=True
			if trigger==False:
				sys.stdout.write("")
				#This is not a start or stop codon. Will add the codon to the string
				codonStr+=translateAA(current)
				sys.stdout.write(textColors.end+current+" ")
	return codonStr

#This function is used to get the compliment.
#It returns the compliment of the string
def getCompliment(DNAstr):
	result=""
	for charecter in DNAstr:
		if charecter == "A":
			result+="T"
		elif charecter == "T":
			result+="A"
		elif charecter == "C":
			result+="G"
		elif charecter == "G":
			result+="C"
	result=result[::-1]
	return result

# this class allows the terminal to print out the anster in color
class textColors:
    green = '\033[92m'
    red   = '\033[91m'
    end = '\033[0m'

# error checking to ensure the correct number of arguements are recieved
# if it is incorrect, it will stop and print out an error message
if len(sys.argv) != 2:
	print("Invalid number of arguements")
	print("Usage: python AssignmentQ5a.py [fiesta file] [path to file]")
	exit()
fiestaDirectory=str(sys.argv[1])

# DNSstr is the string that holds the ATCG string straight from the file
DNAstr=""

# This will open the fine and filter out all the irrellivant metadata
with open(fiestaDirectory,"r") as file1:	# opens the file
	while True:
		charecter = file1.read(1)
		if not charecter:
			break
		if charecter in lookFor:
			DNAstr+=charecter
		else:
			if charecter != "\n":
				DNAstr=""

DNAstr.replace(r"\n","")

rev=False
y=0
AminoAcidArray=["" for x in range(6)]
AAstr=""
for i in range(6):
	if i == 3:
		rev=True
		DNAstr=getCompliment(DNAstr)
		y=0
	AminoAcidArray[i]+=framescanner(y,DNAstr,startCod,stopCod,rev)
	y+=1
	print("")
	print("")

#This calculates the start and stop frames
for x in range(6):
	leftOff=0
	restartFrom=-1
	trigger=False
	# The varibles leftOff and startFrom are to ensure it doesnt go over sections that are already in a frame

	for i, start in enumerate(AminoAcidArray[x]):
		if start == "+" and i >= restartFrom:
			trigger=False
			for y, stop in enumerate(AminoAcidArray[x]):
				if i < y and stop == "-" and leftOff <= i:
					#the results are in function "printResults". This calculates a result and returns a formatted string
					print(printResults(i,y,x,len(AminoAcidArray[x]))+" Amino acid string: "+AminoAcidArray[x-1][i:y])
					leftOff=y
					restartFrom=i
					trigger=True

