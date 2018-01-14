#!/usr/bin/env python
import argparse
import pysam 

parser = argparse.ArgumentParser(description="")
parser.add_argument("samfile", help="input sam or bam file" )
parser.add_argument("out", help="tab seperated output file of percent identities, can be passed to slideAlign.R" )
parser.add_argument('--step', type=int, default=100, help="step size for windowed identity, window must be divisable by step")
parser.add_argument('--window', type=int, default=1000, help = "window size for windowed identity, window must be divisable by step")
parser.add_argument("--header", action="store_true", default=False)
args = parser.parse_args()



M=0 #M	BAM_CMATCH	0
I=1 #I	BAM_CINS	1
D=2 #D	BAM_CDEL	2
N=3 #N	BAM_CREF_SKIP	3
S=4 #S	BAM_CSOFT_CLIP	4
H=5 #H	BAM_CHARD_CLIP	5
P=6 #P	BAM_CPAD	6
E=7 #=	BAM_CEQUAL	7
X=8 #X	BAM_CDIFF	8
B=9 #B	BAM_CBACK	9
NM=10 #NM	NM tag	10
conRef	=	[M, D, N, E, E] # these ones "consume" the reference
conQuery=	[M, I, S, E, X] # these ones "consume" the query
conAln	=	[M, I, D, N, S, E, X] # these ones "consume" the alignments
# format for the output table 
form = ( "{}\t"*13 + "{}\n" )

def makeHeader():
	rtn = form.format("Alignment_Coordinate", "Query_Coordinate", "Reference_Coordinate", 
			"perID_by_matches", "perID_by_events", "perID_by_all", "qname", "rname",
			"Matches", "Mismatches", "Insertions", "Deletions", "InsertionEvents", "DeletionEvents")
	return(rtn)

def perId(matches, missmatch, ins, dele, insEvent, delEvent):
	# avoid dividing by zero if matches and mismatches are 0
	if(matches + missmatch== 0):
		missmatch = 0.0000000000000001
	bymatches = (100.0 * matches)/(matches + missmatch)
	byevents = (100.0 * matches)/(matches + missmatch + insEvent + delEvent)
	byall = (100.0 * matches)/(matches + missmatch + ins + dele)
	return([bymatches, byevents, byall])

# step then the amound it should sum to 
def checkStep(step, test):
	counter = 0
	for aln_type, num in step:
		if( aln_type in conRef ): # these ones "consume" the reference
			counter += num
	assert counter == test, (step)


def stepStats(step):
	match = 0
	mismatch = 0
	ins = 0
	dele = 0
	insEvent = 0
	delEvent = 0
	for aln_type, num in step:
		if(aln_type == E): 
			match += num
		elif(aln_type == X):
			mismatch += num
		elif(aln_type == I):
			insEvent += 1
			ins += num
		elif(aln_type == D):
			delEvent += 1
			dele += num
	stats = {"match":match, "mismatch":mismatch, "ins":ins, "del":dele, "insEvent":insEvent, "delEvent":delEvent}
	return(stats)

	

def cigarSteps(cigar):
	allSteps = []
	counter = 0
	stepCigar = []
	for aln_type, num in cigar:
		# add tuple of cigar to the step
		stepCigar.append( ( aln_type, num) )
		# check is this tupple moves us along the reference
		if(aln_type in conRef): # these ones "consume" the reference
			counter += num
		# end this step if we have gone over the step size
		if(counter >= args.step):
			diff = counter - args.step
			lastType, lastNum = stepCigar[-1]
			#print(diff, counter, lastNum, lastType, stepCigar[-1])
			stepCigar[-1] = (lastType, lastNum - diff)
			# add step to all steps
			allSteps.append(stepCigar)
			# make sure this cigar has the correct window size
			checkStep(stepCigar, args.step)
			# reset the reference consumption counter
			counter = 0
			# reset the step cigar and if nessisary add in leftover from the last tuple
			stepCigar = []
			while(diff > 0):
				overflowNum = min(diff, args.step)
				overflow = (lastType, overflowNum)
				# if the overflow by itself is enough to fill a step do so
				if(overflowNum == args.step):
					allSteps.append( [overflow] )
				else: 
					stepCigar.append(overflow)
					counter += overflowNum
				diff -= args.step

	# add in the last step
	allSteps.append(stepCigar)
	return(allSteps)	

def alnPositions(window):
	qpos = 0
	rpos = 0
	alnpos = 0
	for aln_type, num in window:
		if(aln_type in conRef): # these ones "consume" the reference
			rpos += num
		if(aln_type in conQuery): # these ones "consume" the query
			qpos += num
		if(aln_type in conAln): # these ones "consume" the alignments
			alnpos += num
	return([alnpos, qpos, rpos])




def stepsToWindow(curSteps, steps):
	window = []
	for idx, toUse in enumerate(curSteps):
		if(idx == 0):
			window += steps[toUse]
		else:
			old = window[-1]
			cur = steps[toUse][0]
			# same type on the edge of a steps
			if(old[0] == cur[0]):
				window = window[:-1] + [(cur[0], cur[1] + old[1])] + steps[toUse][1:]
			else:
				window += steps[toUse]
	
	if(curSteps[-1] < len(steps) - args.window/args.step):
		checkStep(window, args.window)
	return(window)	


def windowID(curSteps, steps):
	window = stepsToWindow(curSteps, steps)
	# counts of types in window
	counts = stepStats(window)
	# calculate three types of %ID
	perid = perId( counts["match"], counts["mismatch"], counts["ins"], counts["del"], counts["insEvent"], counts["delEvent"] )
	#only the last step is adding new alignment positions becuase the window is asliging by step size
	poses = alnPositions(steps[curSteps[-1]] )
	# if it is the first window add the whole windows alignment postions 
	if(curSteps[0] == 0 ):
		poses = alnPositions( window )
	return( {"perIDs":perid, "window":window, "poses":poses, "counts":counts} )



def generateTable(windows, read):
	alnpos = 0 - args.window/2
	qpos = 0 - args.window/2
	rpos = read.reference_start - args.window/2 
	# if the beggining is clipped add to the qposition 
	cigarStart = read.cigartuples[0]
	if( cigarStart[0] in [S, H]):
		qpos += cigarStart[1]
	
	table = ""
	for idx, window in enumerate(windows):
		poses = window["poses"]
		perIDs = window["perIDs"]
		counts = window["counts"]
		# update cordiantes
		alnpos += poses[0]
		qpos += poses[1]
		rpos += poses[2]
		table += form.format(alnpos, qpos, rpos, perIDs[0], perIDs[1], perIDs[2], read.query_name, read.reference_name,
				counts["match"], counts["mismatch"], counts["ins"], counts["del"], counts["insEvent"], counts["delEvent"] )
	return(table)

def perIDwindow(read):
	cigar = read.cigartuples
	steps = cigarSteps(cigar)
	stepsInWindow = args.window / args.step
	# a list of dicts with perIDs and windows 
	windows = []
	# calcualtes the perIDs for all windows, also combines steps into windows 
	curSteps = []
	for idx, step in enumerate(steps):
		curSteps.append( idx )
		if(len(curSteps) == stepsInWindow ):
			windows.append( windowID(curSteps, steps) )
			curSteps.pop(0)
	return(generateTable(windows, read))


# make sure window size is divisable by step size
assert args.window % args.step == 0, "window size not divisable by step size"
assert args.window % 2 == 0, "window size not divisable by two"

out = ""
samfile = pysam.AlignmentFile(args.samfile)
for read in samfile.fetch():
	# only get the windowed per ID of the best alignments 
	if( read.flag in [0,16] ):
		out += perIDwindow(read)
samfile.close()

# write out 
if(args.header):
	out = makeHeader() + out
open(args.out, "w+").write(out)




