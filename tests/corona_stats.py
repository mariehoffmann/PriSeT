#python3
import argparse
from collections import Counter
import collections
import functools
import operator
import os
import re
import sys

# The minimal transcript length.
TRANSCRIPT_MIN_LEN = 60

# The minimal transcript length.
TRANSCRIPT_MAX_LEN = 150

# primers computed by PriSeT
primer_fname = 'primers_priset.csv'
# https://www.biosyn.com/tew/Primer-and-Probe-Collection-for-2019-Coronavirus-Detection.aspx#!
primer_corona_fname = 'primers_pub.csv'

lib_covid19 = '/Volumes/plastic_data/priset/library/corona/corona19/corona19.fasta'
libdir_corona_other = '/Volumes/plastic_data/priset/library/corona/corona_other/corona_other_all.fasta'
libdir_corona_other_txt = '/Volumes/plastic_data/priset/library/corona/corona_other/corona_other_all.txt'

work = '/Volumes/plastic_data/priset/work/covid19'
outfile = 'primer_priset_addfilter.csv'
corona19_rx = re.compile('M.+')
corona_oth_rx = re.compile('NC.+')
corona_oth_transcripts = set()
limit_primers = 1000

text = ''  # non covid19 genomes
with open(libdir_corona_other, 'r') as f:
    for line in f.readlines():
        text += line.strip()
# with open(libdir_corona_other_txt, 'w') as fw:
#     fw.write(text)

# build dictionary of primer sequences
def read_primers(filename):
    primers = {}
    with open(filename, 'r') as f:
        for line in f.readlines(): #[:limit_primers]:
            if (line.startswith('#') or len(line.strip()) is 0):
                continue
            ID, seq1, seq2 = line.strip().split(',')[:3]

            # ID, seq1, seq2, freq = line.strip().split(',')[:4]
            # if int(freq) > 15:
            primers[ID] = (seq1, seq2)
    return primers

def has_no_cooccurrence(transcript):
    if text.find(transcript) == -1:
        return True
    else:
        # print(text.find(transcript))
        return False

# given reference sequence, determine all primer positions via text search
def get_matches(acc, reference, primers):
    pID2seq = {}
    # print(primers)
    for primerID, seqs in primers.items():
        pos1 = reference.find(seqs[0])  # search forward primer
        if (pos1 > -1):
            pos2 = reference.find(seqs[1])   # search reverse primer
            if (pos2 > -1 and (pos2 - pos1 >= TRANSCRIPT_MIN_LEN) and (pos2 - pos1 <= TRANSCRIPT_MAX_LEN)):
                transcript = reference[pos1 + len(seqs[0]) : pos2]
                # check if transcript does not occur in other covid genomes
                if has_no_cooccurrence(transcript) is True:
                    pID2seq[primerID] = transcript
                    # print(acc, ', primer ', primerID, ' has no co-occurences')
                    # sys.exit()
                # else:
                #     print(acc, ', primer ', primerID, ' has co-occurences')
                #     print(transcript)
                #     sys.exit()
    return pID2seq

def merge_matches(src, dest):
    for primerID, seqHash in src.items():
        key = (primerID, seqHash)
        if key in dest:
            dest[key].append(fileID)
        else:
            dest[key] = [fileID]
    return dest

def intersect(p1, p2):
    ctrs = [0, 0, 0]
    for pID1, seqs1 in p1.items():
        for pID2, seqs2 in p2.items():
            if seqs1[0] == seqs2[0] and seqs1[1] == seqs2[1]:
                ctrs[0] += 1
            elif seqs1[0] == seqs2[0]:
                ctrs[1] += 1
            elif seqs1[1] == seqs2[1]:
                ctrs[2] += 1
    return ctrs

def has_four_run(sequence):
    return sequence.find('CCCC') > -1 or sequence.find('GGGG') > -1

'''
    1. no 4 runs of C or G
    2. primer length range 18-24
'''
def additional_constraints(sequence1, sequence2):
    if has_four_run(primers[key[0]][0]) is True or has_four_run(primers[key[0]][1]) is True:
        return False
    if len(sequence1) < 18 or len(sequence1) > 24 or len(sequence2) < 18 or len(sequence2) > 24:
        return False
    return True

if __name__ == "__main__":
    primers_priset = read_primers(os.path.join(work, primer_fname))  # PriSeT's computed primers
    primers_pub = read_primers(os.path.join(work, primer_corona_fname)) # published primers
    ctrs = intersect(primers_priset, primers_pub)
    primers = {**primers_priset, **primers_pub}
    pIDSeqHash2accs = {}
    acc = ''
    covid19_genomes = {}
    with open(lib_covid19, 'r') as r:
        sequence = ''
        for line in r.readlines():
            if (line.startswith('>')):
                if len(sequence) > 0:
                    pID2seq = get_matches(acc, sequence, primers)
                    print('handle acc = ', acc)
                    for pID, seqHash in pID2seq.items():
                        key = (pID, seqHash)
                        if key in pIDSeqHash2accs:
                            pIDSeqHash2accs[key].append(acc)
                        else:
                            pIDSeqHash2accs[key] = [acc]
                    sequence = ''
                acc = line[1:line.find(' ')]

            else:
                sequence += line.strip()
        if len(sequence) > 0:
            pID2seq = get_matches(acc, sequence, primers)
            print('handle acc = ', acc)
            for pID, seqHash in pID2seq.items():
                key = (pID, seqHash)
                if key in pIDSeqHash2accs:
                    pIDSeqHash2accs[key].append(acc)
                else:
                    pIDSeqHash2accs[key] = [acc]
            sequence = ''

    print("primerID\tprimer_fwd\tprimer_rev\ttranscript\tgenomes")
    print("-"*80)
    del_keys = set()
    # filter for unmixed primer-transcript combinations
    for key, val in pIDSeqHash2accs.items():
        if len(val) < 2:
            del_keys.add(key)
            continue
        ctr_corona19 = len([1 for acc in val if corona19_rx.match(acc) is not None])
        ctr_corona_oth = len([1 for acc in val if corona_oth_rx.match(acc) is not None])
        if ctr_corona19 < 18 and ctr_corona_oth > 5:
            del_keys.add(key)

    m_srt = [(key, val) for key, val in pIDSeqHash2accs.items() if key not in del_keys]
    m_srt = sorted(m_srt, key = lambda kv : len(kv[1]), reverse=True)
    ctr1 = 0  # differentiating at least one primer pair against others
    ctr2 = 0 # differentiating all 19 genomes against all others
    with open(os.path.join(work, outfile), 'w') as f:
        for key, val in m_srt:
            primerID = key[0]
            if len(val) > 18 and primerID in primers_priset:
                if additional_constraints(primers[primerID][0], primers[primerID][1]) is True:
                    ctr2 += 1
                    line = ','.join([key[0], primers[primerID][0], primers[primerID][1], pID2seq[primerID]])
                    print(line)
                    f.write(line + '\n')
            if key[0] in primers_priset:
                ctr1 += 1


    print('there are ', ctr1, ' primer pairs differentiating at least one genome against other corona virus genomes')
    print('there are ', ctr2, ' primer pairs differentiating all corona 2 genomes against other corona virus genomes')

    print('intersectioning pairs: ', ctrs[0], ', fwd: ', ctrs[1], ', rev: ', ctrs[2])
