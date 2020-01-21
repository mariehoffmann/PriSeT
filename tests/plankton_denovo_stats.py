#python3
import argparse
from collections import Counter
import collections
import functools
import sys

# The minimal transcript length.
TRANSCRIPT_MIN_LEN = 30

# The minimal transcript length.
TRANSCRIPT_MAX_LEN = 800

parser = argparse.ArgumentParser(description = 'Compute cluster statistics for primer pairs.')
parser.add_argument('primers', type = str, nargs=1,
                    help = 'file containing known primer sequences in csv format #ID,seq1,seq2')
parser.add_argument('primers_denovo', type = str, nargs=1,
                    help = 'file containing de novo computer primer sequences in csv format #ID,seq1,seq2')
parser.add_argument('library', type = str, nargs=1,
                    help = 'library file in fasta format with all accessions representing a clade')
parser.add_argument('taxfile', type = str, nargs=1,
                    help = 'give taxonomy file of clade #p_taxid,taxid,is_species')
parser.add_argument('accfile', type = str, nargs=1,
                    help = 'give accession file of clade #taxid,acc1,acc2,...')

class stats(object):
    def __init__(self):
        self.TP = 0
        self.TN = 0
        self.FP = 0
        self.FN = 0

    def to_string(self):
        return " ".join(["TP =", str(self.TP), ", TN =", str(self.TN), ", FP =", str(self.FP), ", FN =", str(self.FN)])

# build dictionary of primer sequences
def read_primers(filename):
    primers = {}
    with open(filename, 'r') as f:
        for line in f.readlines():
            if (line.startswith('#')):
                continue
            ID, seq1, seq2 = line.split(',')
            primers[ID] = (seq1, seq2)
    return primers

# store taxonomic IDs (sp., genus) for each accession, each key taxid is a species
def read_taxa(tax_file, acc_file):
    acc2tax = {}
    tax2ptax = {}
    with open(tax_file, 'r') as f:
        for line in f:
            p_taxid, taxid, is_species = line.strip().split(',')
            if is_species == '1':
                tax2ptax[int(taxid)] = int(p_taxid)
    with open(acc_file, 'r') as f:  #taxid,acc1,acc2,...
        for line in f:
            line = line.strip().split(',')
            taxid = int(line[0])
            for acc in line[1:]:
                if acc not in acc2tax:
                    acc2tax[acc] = taxid

    return tax2ptax, acc2tax

# given reference sequence, determine all primer positions via text search
def get_matches(reference, primers):
    found = set()
    for primer in primers:
        pos1 = reference.find(primer[1][0])
        if (pos1 > -1):
            pos2 = reference.find(primer[1][1])
            if (pos2 > -1 and (pos2 - pos1 >= TRANSCRIPT_MIN_LEN) and (pos2 - pos1 <= TRANSCRIPT_MAX_LEN)):
                found[primer[0]] = hash(reference[pos1 + len(primer[1][0]) : pos2])
    return found

# collect transcripts/amplicons accession-wise and store unique transcripts with
# associated accession (header line) in fasta format
def dereplication(libraryFile, primers):
    # outputFile = fastaFile.split('.')[0] + ".derep.fasta"
    seqhash2accs = {}  # hash(seq) : [acc1, acc2, ...]
    ref = (None,  "")
    with open(libraryFile, 'r') as f:
        if line.startswith(">"):
            if len(ref[1]) > 0: # nor more lines to add, check primer matches
                found = get_matches(reference, primers)
                for primerID, seqhash in found:
                    if seqhash in seqhash2accs:
                        seqhash2accs[seqhash].append(ref[0])
                    else:
                        seqhash2accs[seqhash] = [ref[0]]
            ref = (line.split(' ')[0][1:], "")
        else:
            ref[1] += line.strip()
    return seqhash2accs

# return counts from binary classification, consider multi-class problem as one taxID = class
# against all others. Level = 0 indicates species level and accessions referring to higher
# taxonomic levels are ignored, and level = 1 to genus level, where all species are replaced by
# their genus taxon, accessions assigned to higher than genus levels are ignored.
def bin_class_counts(acc2tax, seqhash2accs, tax2ptax, level = 0):
    # bin_class = collections.namedtuple('bin_class', 'TP TN FP FN')
    tax2cnts = {}
    labels = {}
    taxa = set()
    # iterate over all accessions with same sequence assigned, determine label and count TP, FP
    for seqhash, accs in seqhash2accs.items():
        tax_ctrs = Counter([acc2tax[acc] for acc in accs if acc2tax[acc] in tax2ptax]) # tax: ctr on species level!
        # majority vote for cluster, tax_major is cluster label
        cnt_major, tax_major = 0, None
        for tax, cnt in tax_ctrs.items():
            taxa.add(tax)
            if cnt > cnt_major:
                cnt_major, tax_major = cnt, tax
        if tax_major not in tax2cnts:
            tax2cnts[tax_major] = stats()
            labels[tax_major] = [seqhash]
        else:
            labels[tax_major].append(seqhash)
        for tax, cnt in tax_ctrs.items():
            if tax == tax_major:
                tax2cnts[tax_major].TP += cnt
            else:
                tax2cnts[tax_major].FP += cnt
    # now count FN and TN in other classes after labeling is done
    for tax in taxa:
        if tax not in tax2cnts:
            tax2cnts[tax] = stats()
        for label, seqhash_list in labels.items():
            for seqhash in seqhash_list:
                tax_ctrs = Counter([acc2tax[acc] for acc in seqhash2accs[seqhash] if acc2tax[acc] in tax2ptax])
                tax_ctr = tax_ctrs.get(tax, 0) # all accessions assigned to tax
                tax_all = functools.reduce(lambda acc, ctr: acc + ctr, tax_ctrs.values(), 0) # count all
                if tax != label:
                    tax2cnts[tax].FN += tax_ctr     # tax in cluster, but it is labeled differently
                    tax2cnts[tax].TN += tax_all - tax_ctr

    for tax in  labels.keys():
        print('tax = ', tax, ', num clusters = ', len(labels[tax]))
    return tax2cnts

# compute accuracy, i.e. sum of correctly labeled samples divided by total number of samples
def accuracy(stats):
    return 0 if (stats.TP + stats.TN == 0) else float((stats.TP + stats.TN)/(stats.TP + stats.TN + stats.FP + stats.FN))

# compute precision, i.e. true positive divided by all positively labeled items
# TP: assigned taxID of accession corresponds to cluster label
# FP: assigned taxID does not correspond to cluster label
def precision(stats):
    return 0 if (stats.TP == 0) else float(stats.TP/(stats.TP + stats.FP))

# compute recall, i.e. proportion of true positives out of all positive samples
# TP: assigned taxID of accession corresponds to cluster label
# FN: accessions in clusters with different class label
def recall(stats):
    return 0 if (stats.TP == 0) else float(stats.TP/(stats.TP + stats.FN))

# compute cluster distinguishability (count of label pure clusters)
# cluster labels assigned based on majority vote
def purity(acc2tax, seqhash2accs, tax2ptax):
    purity_species = [0, 0]    # pure, total on species level
    purity_genera = [0, 0]     # pure, total on genus level
    for seqhash, accs in seqhash2accs.items():
        species = set([acc2tax[acc] for acc in accs if acc2tax[acc] in tax2ptax])
        if len(species) == 1:
            purity_species[0] += 1
        purity_species[1] += 1
        # select those which are species and assign their parent taxa
        genera = [tax2ptax[acc2tax[acc]] for acc in accs if acc2tax[acc] in tax2ptax]
        # add those already assigned to genus level
        genera += [acc2tax[acc] for acc in accs if acc2tax[acc] in tax2ptax.values()]
        counts_genera = Counter(genera)
        if len(counts_genera) == 1:
            purity_genera[0] += 1
        purity_genera[1] += 1
    return purity_species, purity_genera

'''
expect for cluster1 = (A1, A2, C1) and cluster2 = (B1, B2)
    tax = 1: TP = 2 , TN = 2 , FP = 1 , FN = 0
    tax = 2: TP = 2 , TN = 3 , FP = 0 , FN = 0
    tax = 3: TP = 0 , TN = 4 , FP = 0 , FN = 1
                tax 1   tax 2   tax 3
    accuracy    0.8     1.0     0.8
    precision   0.67    1.0     0.0
    recall      1.0     1.0     0.0

    level       species     genus
    purity      0.5         1.0
'''
def test_bin_class_counts():
    acc2tax = {'A1': 1, 'A2': 1, 'B1': 2, 'B2': 2, 'C1': 3}
    seqhash2accs = {'X34_seq': ['A1', 'A2', 'C1'], 'Y56_seq': ['B1', 'B2']}
    tax2ptax = {1: 11, 2: 22, 3: 11}  # all species
    tax2cnts = bin_class_counts(acc2tax, seqhash2accs, tax2ptax)
    purity_species, purity_genera = purity(acc2tax, seqhash2accs, tax2ptax)
    print('tax = 1:\t', tax2cnts[1].to_string())
    print('tax = 2:\t', tax2cnts[2].to_string())
    print('tax = 3:\t', tax2cnts[3].to_string())
    print('accuracy =\t', round(accuracy(tax2cnts[1]), 2), "\t", round(accuracy(tax2cnts[2]), 2), "\t", round(accuracy(tax2cnts[3]), 2))
    print('precision =\t', round(precision(tax2cnts[1]), 2), "\t", round(precision(tax2cnts[2]), 2), "\t", round(precision(tax2cnts[3]), 2))
    print('recall =\t', round(recall(tax2cnts[1]), 2), "\t", round(recall(tax2cnts[2]), 2), "\t", round(recall(tax2cnts[3]), 2))
    print('purity species =\t', purity_species[0]/purity_species[1], "\tgenera =\t", purity_genera[0]/purity_genera[1])

if __name__ == "__main__":
    test_bin_class_counts()
    # args = parser.parse_args()
    # print(args)
    # primers = read_primers(args.primers)
    # seqhash2accs = dereplication(args.libraryFile, primers)
    # tax2ptax, acc2tax = read_taxa(args.primers)
