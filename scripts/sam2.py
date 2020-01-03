import sys
import pysam
import collections
import numpy as np


samFile = pysam.AlignmentFile(sys.argv[1], "rb")
in_fasta_path = sys.argv[2]
bb_outFile = sys.argv[3]
ins_outFile = sys.argv[4]
stats_outFile = sys.argv[5]
# Create a dict with lists of SimpleRead using referenceName as key
myDict = collections.defaultdict(list)
for read in samFile:
    # Remove matches that suck
    try:
        alnErr = sum(x[1] for x in read.cigartuples if x[0] in [1,2])
        if alnErr/float(read.reference_end) < 0.2:
        #read.reference_start + read.infer_query_length(always=True)
            myDict[read.reference_name].append((read.reference_start,read.reference_end,read.flag&(1<<4)>0,read.query_name))
    except TypeError as e:
        # if there is no matching reads in sam file
        touch(ins_outFile)
        touch(bb_outFile)
        touch(stats_outFile)
        exit() 
# Combine matches per read
sortedKeys =list( myDict.keys())
sortedKeys.sort()
for key in sortedKeys:
    myDict[key].sort()
    # print(key,myDict[key])

# exit()

min_gap = 50
max_gap = 500
# Work through the fastq file
with open(in_fasta_path,'r') as fasta_file, open(ins_outFile,'w') as dumpFile, \
        open(bb_outFile,'w') as dump_bb_File, open(stats_outFile,'w') as statFile:
    try:
        while True:
            # python3 = fqFile.readline().rstrip() # need testing. 
            readName=next(fasta_file).rstrip()
            readSeq=next(fasta_file).rstrip()
#            readPlus=fqFile.next().rstrip()
#            readPhred=fqFile.next().rstrip()
            readId=readName[1:].split()[0]


            cutId = 1
            if readId not in myDict:
                statFile.write(', '.join(
                    [str(y) for y in
                        [readId,
                        0,
                        0,
                        0,
                        0,
                        0]])+'\n')
                continue

            cur_bbs = myDict[readId]
            # Split sequence, dump backbone
            for i,x in enumerate(cur_bbs):
                # print(x[1]-x[0])
                # Extend the identifier
                dump_bb_File.write(">" + readId +  "_"+  str(i) + '\n')#+':'+str(x[0])+'-'+str(x[1])

                # Dump the actual sub sequence with primers
                dump_bb_File.write(readSeq[x[0]:x[1]] + '\n')
#                dump_bb_File.write(readPlus+'\n')

                # Add perfect phred scores for forced primers
#                dump_bb_File.write(readPhred[x[0]:x[1]] + '\n')

            # Mind the gaps, should be inserts
            gaps = []
            if min_gap < cur_bbs[0][0] < max_gap:
                gaps.append((0,cur_bbs[0][0]))

            for i,x in enumerate(cur_bbs[1:]):
                if min_gap < x[0] - cur_bbs[i][1] < max_gap:
                    gaps.append((cur_bbs[i][1],x[0]))

            if min_gap < len(readSeq) - cur_bbs[-1][1] < max_gap:
                gaps.append((cur_bbs[-1][1],len(readSeq)))

            for i,x in enumerate(gaps):
                # Extend the identifier
                dumpFile.write(">" + readId + "_" + str(i) + '\n')#+':'+str(x[0])+'-'+str(x[1])

                # Dump the actual sub sequence with primers
                dumpFile.write(readSeq[x[0]:x[1]] + '\n')
#                dumpFile.write(readPlus+'\n')

                # Add perfect phred scores for forced primers
#                dumpFile.write(readPhred[x[0]:x[1]] + '\n')

            statFile.write(', '.join(
                [str(y) for y in
                    [readId,
                    len(cur_bbs),
                    sum([x[2] for x in cur_bbs]),
                    np.median([x[1]-x[0] for x in cur_bbs]),
                    len(gaps),
                    (np.median([x[1]-x[0] for x in gaps]) if len(gaps)>0 else 0)]])+'\n')

    except StopIteration:
        pass

#print('Splitting data finished.')
