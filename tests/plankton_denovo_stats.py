#python3
import argparse
from collections import Counter

# The minimal transcript length.
TRANSCRIPT_MIN_LEN = 30

# The minimal transcript length.
TRANSCRIPT_MAX_LEN = 800

parser = argparse.ArgumentParser(description = 'Compute cluster statistics for primer pairs.')
parser.add_argument('primers', type = string, nargs=1,
                    help = 'file containing known primer sequences in csv format #ID,seq1,seq2')
parser.add_argument('primers_denovo', type = string, nargs=1,
                    help = 'file containing de novo computer primer sequences in csv format #ID,seq1,seq2')
parser.add_argument('library', type = string, nargs=1,
                    help = 'library file in fasta format with all accessions representing a clade')
parser.add_argument('taxfile', type = string, nargs=1,
                    help = 'give taxonomy file of clade #p_taxid,taxid,is_species')
parser.add_argument('accfile', type = string, nargs=1,
                    help = 'give accession file of clade #taxid,acc1,acc2,...')


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
    for primer : primers:
        pos1 = reference.find(primer[1][0])
        if (pos1 > -1):
            pos2 = reference.find(primer[1][1])
            if (pos2 > -1 and (pos2 - pos1 >= TRANSCRIPT_MIN_LEN) and (pos2 - pos1 <= TRANSCRIPT_MAX_LEN))
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

# compute cluster distinguishability (count of label pure clusters)
# cluster labels assigned based on majority vote
def purity(acc2tax):
    purity_species = (0,0)    # pure, total on species level
    purity_genera = (0,0)     # pure, total on genus level
    for item in seqhash2accs:
        counts = Counter([acc2tax[acc] for acc in item[1] if acc2tax[acc] in tax2ptax])
        if len(counts) == 1:
            purity_species[0] += 1
        purity_species[1] += 1

    return purity_species, purity_genera

# compute precision = TP/(TP + FP) based on clusters labeled with majority taxID
# TP: assigned taxID of accession corresponds to cluster label
# FP: assigned taxID does not correspond to cluster label
def precision():
    TP = 0
    FP = 0
    for seqhash, accs in seqhash2accs:
        counts = Counter([acc2tax[acc] for acc in accs if acc2tax[acc] in tax2ptax]) # tax: ctr
        majority_cnt = 0
        majority_tax = None
        for tax, cnt in counts:
            if cnt > majority_cnt:
                majority_cnt = cnt
                majority_tax = tax
        labels[seqhash] = majority_tax
        for tax, cnt in counts:
            if tax == labels[seqhash]:
                TP += cnt
            else:
                FP += cnt
    return float(TP)/float(TP + FP)

# compute recall = TP/(TP + FN) for each taxon
# FN: accessions in clusters with different class label
def recall():
    class2stat = {}  #tax : (TP, FN)
    for seqhash, accs in seqhash2accs:
        counts = Counter([acc2tax[acc] for acc in accs if acc2tax[acc] in tax2ptax]) # tax: ctr
        majority_cnt = 0
        majority_tax = None
        for tax, cnt in counts:
            if cnt > majority_cnt:
                majority_cnt = cnt
                majority_tax = tax
        labels[seqhash] = majority_tax
        for tax, cnt in counts:
            if tax == labels[seqhash]:
                if tax in class2stat:
                    class2stat[tax][0] += cnt
                else:
                    class2stat[tax] = (cnt, 0)
            else:  # false negative for alien taxa
                if tax in class2stat:
                    class2stat[tax][1] += cnt
                else:
                    class2stat[tax] = (0, cnt)
    recall = {}
    for tax, stat in class2stat:
        recall[tax] = float(stat[0])/float(stat[0] + stat[1])
    return recall

if __name__ == "__main__":
    args = parser.parse_args()
    print(args)
    primers = read_primers(args.primers)
    seqhash2accs = dereplication(args.libraryFile, primers)
    tax2ptax, acc2tax = read_taxa(args.primers)
