import sys
sys.path.append("/hpc/local/CentOS7/cog/lib/python2.7/site-packages/")
import pysam
import collections
import numpy as np


samFile = pysam.AlignmentFile(sys.argv[1], "rb")
inFile = sys.argv[2]
bb_outFile = sys.argv[3]
ins_outFile = sys.argv[4]
# Create a dict with lists of SimpleRead using referenceName as key
myDict = collections.defaultdict(list)
for read in samFile:
    # Remove matches that suck
    alnErr = sum(x[1] for x in read.cigartuples if x[0] in [1,2])
    if alnErr/float(read.reference_end) < 0.2:
        #read.reference_start + read.infer_query_length(always=True)
        myDict[read.reference_name].append((read.reference_start,read.reference_end,read.flag&(1<<4)>0,read.query_name))

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
with open(inFile,'r') as fqFile, open(ins_outFile,'w') as dumpFile, \
        open(bb_outFile,'w') as dump_bb_File, open(bb_outFile+'_stat.csv','w') as statFile:
    try:
        while True:
            readName=fqFile.next().rstrip()
            readSeq=fqFile.next().rstrip()
            readPlus=fqFile.next().rstrip()
            readPhred=fqFile.next().rstrip()
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
                dump_bb_File.write(readName + ' c_start=' + str(x[0]) + ' c_ref=' + x[3] + ' c_str=' + ('-' if x[2] else '+') + ' c_idx=' + str(i) + '\n')#+':'+str(x[0])+'-'+str(x[1])

                # Dump the actual sub sequence with primers
                dump_bb_File.write(readSeq[x[0]:x[1]] + '\n')
                dump_bb_File.write(readPlus+'\n')

                # Add perfect phred scores for forced primers
                dump_bb_File.write(readPhred[x[0]:x[1]] + '\n')

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
                dumpFile.write(readName + ' c_start=' + str(x[0]) + ' c_ref=gap' + ' c_idx=' + str(i) + '\n')#+':'+str(x[0])+'-'+str(x[1])

                # Dump the actual sub sequence with primers
                dumpFile.write(readSeq[x[0]:x[1]] + '\n')
                dumpFile.write(readPlus+'\n')

                # Add perfect phred scores for forced primers
                dumpFile.write(readPhred[x[0]:x[1]] + '\n')

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

print('Splitting data finished.')
