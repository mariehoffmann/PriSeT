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

'''
taxid=2825
lib=/Volumes/plastic_data/tactac/subset
python ../PriSeT/tests/plankton_denovo_stats.py --library $lib/$taxid/root_$taxid.fasta --accfile $lib/$taxid/root_$taxid.acc --taxfile $lib/$taxid/root_$taxid.tax --primers /Volumes/plastic_data/priset/primers.csv
'''

parser = argparse.ArgumentParser(description = 'Compute cluster statistics for primer pairs.')
parser.add_argument('--primers', type = str, nargs=1,
                    help = 'file containing known primer sequences in csv format #ID,seq1,seq2')
parser.add_argument('--primers_denovo', type = str, nargs='?',
                    help = 'file containing de novo computer primer sequences in csv format #ID,seq1,seq2')
parser.add_argument('--library', type = str, nargs=1,
                    help = 'library file in fasta format with all accessions representing a clade')
parser.add_argument('--taxfile', type = str, nargs=1,
                    help = 'give taxonomy file of clade #p_taxid,taxid,is_species')
parser.add_argument('--accfile', type = str, nargs=1,
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
            ID, seq1, seq2 = line.strip().split(',')
            primers[ID] = (seq1, seq2)
    return primers

# store taxonomic IDs (sp., genus) for each accession, each key taxid is a species
def read_taxa(taxfile, accfile):
    acc2tax = {}
    tax2ptax = {}
    with open(taxfile, 'r') as f:
        for line in f:
            p_taxid, taxid, is_species = line.strip().split(',')
            if is_species == '1':
                tax2ptax[int(taxid)] = int(p_taxid)
    with open(accfile, 'r') as f:  #taxid,acc1,acc2,...
        line = f.readline()
        line = f.readline()
        while line:
            if line.startswith("#"):
                continue
            print(line)
            line = line.strip().split(',')
            taxid = int(line[0])
            for acc in line[1:]:
                if acc not in acc2tax:
                    acc2tax[acc] = taxid
            line = f.readline()
    return tax2ptax, acc2tax

# given reference sequence, determine all primer positions via text search
def get_matches(reference, primers):
    found = {}
    print(primers)
    for primerID, seqs in primers.items():
        pos1 = reference.find(seqs[0])  # search forward primer
        if (pos1 > -1):
            pos2 = reference.find(seqs[1])   # search reverse primer
            if (pos2 > -1 and (pos2 - pos1 >= TRANSCRIPT_MIN_LEN) and (pos2 - pos1 <= TRANSCRIPT_MAX_LEN)):
                found[primerID] = hash(reference[pos1 + len(seqs[0]) : pos2])
    return found

# collect transcripts/amplicons accession-wise and store unique transcripts with
# associated accession (header line) in fasta format
def dereplication(libraryFile, primers):
    # outputFile = fastaFile.split('.')[0] + ".derep.fasta"
    seqhash2accs = {}  # (hash(seq), pID) : [acc1, acc2, ...]
    ref = [None,  ""]
    with open(libraryFile, 'r') as f:
        line = f.readline()
        while line:
            if line.startswith(">"):
                if len(ref[1]) > 0: # no more lines to add, check primer matches
                    found = get_matches(ref[1], primers)
                    for primerID, seqhash in found.items():
                        key = (primerID, seqhash)
                        if key in seqhash2accs:
                            seqhash2accs[key].append(ref[0])
                        else:
                            seqhash2accs[key] = [ref[0]]
                ref = [line.split(' ')[0][1:], ""]
                print(ref[0])
            else:
                ref[1] += line.strip()
            line = f.readline()
    return seqhash2accs

# return counts from binary classification, consider multi-class problem as one taxID = class
# against all others. Level = 0 indicates species level and accessions referring to higher
# taxonomic levels are ignored, and level = 1 to genus level, where all species are replaced by
# their genus taxon, accessions assigned to higher than genus levels are ignored.
def bin_class_counts(acc2tax, seqhash2accs, tax2ptax, level = 0):
    # bin_class = collections.namedtuple('bin_class', 'TP TN FP FN')
    tax2stats = {}
    labels = {}
    taxa = set()
    primerIDs = set()
    # iterate over all accessions with same sequence assigned, determine label and count TP, FP
    for key, accs in seqhash2accs.items():
        primerID, seqhash = key[0], key[1]
        primerIDs.add(primerID)
        tax_ctrs = Counter([acc2tax[acc] for acc in accs if acc2tax[acc] in tax2ptax]) # tax: ctr on species level!
        # majority vote for cluster, tax_major is cluster label
        cnt_major, tax_major = 0, None
        for tax, cnt in tax_ctrs.items():
            taxa.add(tax)
            if cnt > cnt_major:
                cnt_major, tax_major = cnt, tax
        if tax_major not in tax2stats:
            tax2stats[(primerID, tax_major)] = stats()
            labels[(primerID, tax_major)] = [seqhash]
        else:
            labels[(primerID, tax_major)].append(seqhash)
        for tax, cnt in tax_ctrs.items():
            if tax == tax_major:
                tax2stats[(primerID, tax_major)].TP += cnt
            else:
                tax2stats[(primerID, tax_major)].FP += cnt
    # now count FN and TN in other classes after labeling is done
    for tax in taxa:
        for primerID in primerIDs:
            key = (primerID, tax)
            if key not in tax2stats:
                tax2stats[key] = stats()
            for label, seqhash_list in labels.items():  # label = (primerID, tax)
                for seqhash in seqhash_list:
                    tax_ctrs = Counter([acc2tax[acc] for acc in seqhash2accs[(primerID, seqhash)] if acc2tax[acc] in tax2ptax])
                    tax_ctr = tax_ctrs.get(tax, 0) # all accessions assigned to tax
                    tax_all = functools.reduce(lambda acc, ctr: acc + ctr, tax_ctrs.values(), 0) # count all
                    if tax != label[1]:
                        tax2stats[key].FN += tax_ctr     # tax in cluster, but it is labeled differently
                        tax2stats[key].TN += tax_all - tax_ctr

    for label, seqhash_list in  labels.items():
        print('tax = ', label[1], ', num clusters = ', len(seqhash_list))
    return tax2stats

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
    purity_species = {}    # pure, total on species level
    purity_genera = {}     # pure, total on genus level
    for key, accs in seqhash2accs.items():
        PID, seqhash = key[0], key[1]
        if PID not in purity_species:
            purity_species[PID] = [0, 0]
        if PID not in purity_genera:
            purity_genera[PID] = [0, 0]
        species = set([acc2tax[acc] for acc in accs if acc2tax[acc] in tax2ptax])
        if len(species) == 1:
            purity_species[PID][0] += 1
        purity_species[PID][1] += 1
        # select those which are species and assign their parent taxa
        genera = [tax2ptax[acc2tax[acc]] for acc in accs if acc2tax[acc] in tax2ptax]
        # add those already assigned to genus level
        genera += [acc2tax[acc] for acc in accs if acc2tax[acc] in tax2ptax.values()]
        counts_genera = Counter(genera)
        if len(counts_genera) == 1:
            purity_genera[PID][0] += 1
        purity_genera[PID][1] += 1
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
    seqhash2accs = {('Primer1','X34_seq'): ['A1', 'A2', 'C1'], ('Primer1','Y56_seq'): ['B1', 'B2']}
    tax2ptax = {1: 11, 2: 22, 3: 11}  # all species
    tax2stats = bin_class_counts(acc2tax, seqhash2accs, tax2ptax)
    purity_species, purity_genera = purity(acc2tax, seqhash2accs, tax2ptax)
    key1, key2, key3 = ('Primer1', 1), ('Primer1', 2), ('Primer1', 3)
    print('tax = 1:\t', tax2stats[key1].to_string())
    print('tax = 2:\t', tax2stats[key2].to_string())
    print('tax = 3:\t', tax2stats[key3].to_string())
    print('accuracy =\t',  round(accuracy(tax2stats[key1]), 2), "\t",
                            round(accuracy(tax2stats[key2]), 2), "\t",
                            round(accuracy(tax2stats[key3]), 2))
    print('precision =\t',  round(precision(tax2stats[key1]), 2), "\t",
                            round(precision(tax2stats[key2]), 2), "\t",
                            round(precision(tax2stats[key3]), 2))
    print('recall =\t',     round(recall(tax2stats[key1]), 2), "\t",
                            round(recall(tax2stats[key2]), 2), "\t",
                            round(recall(tax2stats[key3]), 2))
    print('purity species =\t', purity_species['Primer1'][0]/purity_species['Primer1'][1], "\tgenera =\t",
                                purity_genera['Primer1'][0]/purity_genera['Primer1'][1])

if __name__ == "__main__":
    test_bin_class_counts()
    # args = parser.parse_args()
    # print(args)
    # primers = read_primers(args.primers[0])
    # # primers_denovo = read_primers(args.primers_denovo[0])
    # seqhash2accs = dereplication(args.library[0], primers)
    # tax2ptax, acc2tax = read_taxa(args.taxfile[0], args.accfile[0])
    # tax2stats = bin_class_counts(acc2tax, seqhash2accs, tax2ptax, 0)
    # for tax, stat in tax2stats.items():
    #     print(tax, ": ", stat.to_string())
