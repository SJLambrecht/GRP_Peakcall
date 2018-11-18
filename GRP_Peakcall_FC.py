#! usr/env/bin/python
import sys, re, os


# 1. Name your Dataset
item='Name'

# 2. index of dataset to be used for peak calling, first column is 0
dataSet = 0

# 3. index count data for fold change calculation, first column is 0

countdata = 0

# 4. index of reference e.g. control dataset, first column is 0
referenceData = 0

# 5. Enter file name or path for forward reads .grp file
fwReader = open('Forward_Reads.grp','r')
fwHandle = fwReader.readlines()
fwReader.close()

# 6. Enter file name or path for reverse reads .grp file
rvReader = open('Reverse_Reads.grp','r')
rvHandle = rvReader.readlines()
rvReader.close()

# 7. Range for baseplate
baseRange = 5
# 8. Threshold for difference between baseplates for peak start
threshHold = 1.5
# 9. Threshold for difference between baseplates for peak end
downHold = 1.5
# 10. Range for scanning exact peak start
ScanRange = 100
# 11. Range for scanning exact peak end
ScanRangeEnd = 100
# 12. Space between baseplates for start peak scanning
spacer = 10
# 13. Space between baseplates for end peak scanning
spacerEnd = 10
# 14. Percent of coverage that defines the peak
percCov = 0.2
# 16. Frame
frame = '.'

# 17. choose color for annotation, 0 = white, 1 = dark grey, 2 = red, 3 = green, 4 = blue, 5 = cyan, 6 = magenta, 7 = yellow, 8 = pale green, 9 = light sky blue, 10 = orange, 11 = brown, 12 = pale pink, 13 = light grey, 14 = black, 15 = mid red, 16 = light red, 17 = pink 
color = 0

# 18. Enter genome size
genomeSize = 1


counter = 0
counter2 = 0
counter3 = 0
counter4=0
endCounter=0
numCov={}
comCov={}
numbers=[]
numbers1=[]
compare=[]
compare1=[]
scan=[]
scanEnd=[]
count=[]
count1=[]

#Reading forward GRP
for line in fwHandle:
	line = line.split('\t')
	numbers.append(float(line[dataSet].strip()))
	compare.append(float(line[referenceData].strip()))
	count.append(float(line[countdata].strip()))

#Reading Reverse GRP
for line in reversed(rvHandle):
	line = line.split('\t')
	numbers1.append(float(line[dataSet].strip())*-1)
	compare1.append(float(line[referenceData].strip())*-1)
	count1.append(float(line[countdata].strip())*-1)

#Calculates the average per nucleotide in the data set
coverage = sum(numbers)+sum(numbers1)
print coverage/(genomeSize*2), 'Average Covererage per Nucleotide'

#Sets the threshold for coverage of a peak
covThresh = (coverage/(genomeSize*2)*percCov)

TUgff = open('peakcall_'+item+'.gff','w')
#Writes settings to the result file
TUgff.writelines('##data_set '+str(dataSet)+', base_range '+str(baseRange)+', Threshold up '+str(threshHold)+', Threshold down '+str(downHold)+', scan range start '+str(ScanRange)+', scan range end '+str(ScanRangeEnd)+', spacer start '+str(spacer)+', spacer end '+str(spacerEnd)+', coverage threshold '+str(covThresh)+', percent coverage '+str(percCov)+'%\n')

#Finding forward peaks
FC=0
FC2=0
FC3=0
FC4=0

while counter <=genomeSize:
	FC=0
	FC2=0
	FC3=0
	FC4=0
	fineEnd=0
	a = sum(numbers[counter:counter+baseRange])
	b = sum(numbers[counter+baseRange+spacer:counter+baseRange+baseRange+spacer])
	if a==0 and b >0:
		FC = b/baseRange
	if a >0 and b>0:
		FC = (b/baseRange)/(a/baseRange)
	if a>0 and b==0:
		FC = 0
	if FC>threshHold:
		while counter3<=ScanRange:
			c=sum(numbers[counter+counter3:counter+counter3+baseRange])
			d=sum(numbers[counter+counter3+baseRange+spacer:counter+counter3+baseRange+baseRange+spacer])
			if c==0 and d >0:
				FC2 = d/baseRange
			if c >0 and d>0:
				FC2 = (d/baseRange)/(c/baseRange)
			if c>0 and d==0:
				FC2 = 0
			scan.append(FC2)
			counter3=counter3+1
		counter3=0
		endCounter=counter+scan.index(max(scan))+baseRange+spacer
		while endCounter<=genomeSize:
			d=sum(numbers[endCounter:endCounter+baseRange])
			e=sum(numbers[endCounter+baseRange+spacerEnd:endCounter+baseRange+baseRange+spacerEnd])
			if e==0 and d>0:
				break
			if d>0 and e>0:
				FC3=(d/baseRange)/(e/baseRange)
			if e>0 and d==0:
				break
			if FC3>downHold:
				fineEnd=1
				while counter3<=ScanRangeEnd:
					f=sum(numbers[endCounter:endCounter+baseRange]) 
					g=sum(numbers[endCounter+counter3+baseRange+spacerEnd:endCounter+counter3+baseRange+baseRange+spacerEnd])
					if f== 0 and g>0:
						FC=0
					if f>0 and g>0:
						FC4 = (f/baseRange)/(g/baseRange)
					if f>0 and g==0:
						FC4 = f/baseRange
					scanEnd.append(FC4)
					counter3=counter3+1
				counter3 =0
				break
			else:
				endCounter=endCounter+1
		if fineEnd==0:
			if sum(numbers[counter+scan.index(max(scan))+baseRange:endCounter+baseRange+spacerEnd])/len(numbers[counter+scan.index(max(scan))+baseRange:endCounter+baseRange+spacerEnd])>=covThresh:
				sample=sum(count[counter+scan.index(max(scan))+baseRange:endCounter+baseRange+spacerEnd])
				control =sum(compare[counter+scan.index(max(scan))+baseRange:endCounter+baseRange+spacerEnd])
				try:
					division=sample/control
				except ZeroDivisionError:
					if control ==0:
						division= str(sample)+' s*'
					if sample ==0:
						division= str(control)+' c*'
				TUgff.writelines('all_bases\tGenbank\tpeak\t'+str(counter+scan.index(max(scan))+baseRange)+'\t'+str(endCounter+baseRange+spacerEnd)+'\t'+str(division)+'\t+\t'+str(frame)+'\tcolour='+str(color)+' ;index "TU'+str(counter2)+'"'+'\n')
				counter2= counter2+1
			counter= endCounter
		if fineEnd==1:
			if sum(numbers[counter+scan.index(max(scan))+baseRange:endCounter+scanEnd.index(max(scanEnd))+baseRange+spacerEnd])/len(numbers[counter+scan.index(max(scan))+baseRange:endCounter+scanEnd.index(max(scanEnd))+baseRange+spacerEnd]) >=covThresh:
				sample=sum(count[counter+scan.index(max(scan))+baseRange:endCounter+scanEnd.index(max(scanEnd))+baseRange+spacerEnd])
				control=sum(compare[counter+scan.index(max(scan))+baseRange:endCounter+scanEnd.index(max(scanEnd))+baseRange+spacerEnd])
				try:
					division=sample/control
				except ZeroDivisionError:
					if control ==0:
						division= str(sample)+' s*'
					if sample ==0:
						division= str(control)+' c*'
				TUgff.writelines('all_bases\tGenbank\tpeak\t'+str(counter+scan.index(max(scan))+baseRange)+'\t'+str(endCounter+scanEnd.index(max(scanEnd))+baseRange+spacerEnd)+'\t'+str(division)+'\t+\t'+str(frame)+'\tcolour='+str(color)+' ;index "TU'+str(counter2)+'"'+'\n')
				counter2= counter2+1
			counter= endCounter
		del scan[:]
		del scanEnd[:]
	else:
		counter = counter+1


#Finding reverse peaks
counter = 0
counter3 = 0
counter4=0
endCounter=0
numCov={}
comCov={}
compare=[]
scan=[]
scanEnd=[]
FC=0
FC2=0
FC3=0
FC4=0

while counter <=genomeSize:
	FC=0
	FC2=0
	FC3=0
	FC4=0
	fineEnd=0
	a = sum(numbers1[counter:counter+baseRange])
	b = sum(numbers1[counter+baseRange+spacer:counter+baseRange+baseRange+spacer])
	if a==0 and b >0:
		FC = b/baseRange
	if a >0 and b>0:
		FC = (b/baseRange)/(a/baseRange)
	if a>0 and b==0:
		FC = 0
	if FC>threshHold:
		while counter3<=ScanRange:
			c=sum(numbers1[counter+counter3:counter+counter3+baseRange])
			d=sum(numbers1[counter+counter3+baseRange+spacer:counter+counter3+baseRange+baseRange+spacer])
			if c==0 and d >0:
				FC2 = d/baseRange
			if c >0 and d>0:
				FC2 = (d/baseRange)/(c/baseRange)
			if c>0 and d==0:
				FC2 = 0
			scan.append(FC2)
			counter3=counter3+1
		counter3=0
		endCounter=counter+scan.index(max(scan))+baseRange+spacer
		while endCounter<=genomeSize:
			d=sum(numbers1[endCounter:endCounter+baseRange])
			e=sum(numbers1[endCounter+baseRange+spacerEnd:endCounter+baseRange+baseRange+spacerEnd])
			if e==0 and d>0:
				break
			if d>0 and e>0:
				FC3=(d/baseRange)/(e/baseRange)
			if e>0 and d==0:
				break
			if FC3>downHold:
				fineEnd=1
				while counter3<=ScanRangeEnd:
					f=sum(numbers1[endCounter:endCounter+baseRange])
					g=sum(numbers1[endCounter+counter3+baseRange+spacerEnd:endCounter+counter3+baseRange+baseRange+spacerEnd])
					if f== 0 and g>0:
						FC4=0
					if f>0 and g>0:
						FC4 = (f/baseRange)/(g/baseRange)
					if f>0 and g==0:
						FC4 = f/baseRange
					scanEnd.append(FC4)
					counter3=counter3+1
				counter3 =0
				break
			else:
				endCounter=endCounter+1
		if fineEnd==0:
			if sum(numbers1[counter+scan.index(max(scan))+baseRange:endCounter+baseRange+spacerEnd])/len(numbers1[counter+scan.index(max(scan))+baseRange:endCounter+baseRange+spacerEnd])>=covThresh:
				sample=sum(count1[counter+scan.index(max(scan))+baseRange:endCounter+baseRange+spacerEnd])
				control=sum(compare1[counter+scan.index(max(scan))+baseRange:endCounter+baseRange+spacerEnd])
				try:
					division=sample/control
				except ZeroDivisionError:
					if control ==0:
						division= str(sample)+' s*'
					if sample ==0:
						division= str(control)+' c*'
				TUgff.writelines('all_bases\tGenbank\tpeak\t'+str(genomeSize-(endCounter+baseRange+spacerEnd))+'\t'+str(genomeSize-(counter+scan.index(max(scan))+baseRange))+'\t'+str(division)+'\t-\t'+str(frame)+'\tcolour='+str(color)+' ;index "TU'+str(counter2)+'"'+'\n')
				counter2= counter2+1
			counter= endCounter

		if fineEnd==1:
			if sum(numbers1[counter+scan.index(max(scan))+baseRange:endCounter+scanEnd.index(max(scanEnd))+baseRange+spacerEnd])/len(numbers1[counter+scan.index(max(scan))+baseRange:endCounter+scanEnd.index(max(scanEnd))+baseRange+spacerEnd]) >=covThresh:
				sample=sum(count1[counter+scan.index(max(scan))+baseRange:endCounter+scanEnd.index(max(scanEnd))+baseRange+spacerEnd])
				control=sum(compare1[counter+scan.index(max(scan))+baseRange:endCounter+scanEnd.index(max(scanEnd))+baseRange+spacerEnd])
				try:
					division=sample/control
				except ZeroDivisionError:
					if control ==0:
						division= str(sample)+' s*'
					if sample ==0:
						division= str(control)+' c*'
				TUgff.writelines('all_bases\tGenbank\tpeak\t'+str(genomeSize-(endCounter+scanEnd.index(max(scanEnd))+baseRange+spacerEnd))+'\t'+str(genomeSize-(counter+scan.index(max(scan))+baseRange))+'\t'+str(division)+'\t-\t'+str(frame)+'\tcolour='+str(color)+';index "TU'+str(counter2)+'"'+'\n')
				counter2= counter2+1
			counter= endCounter
		del scan[:]
		del scanEnd[:]
	else:
		counter = counter+1

TUgff.close()
print counter2, 'Transcriptional Units Defined'
print 'Finished'
