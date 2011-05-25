#!/usr/bin/env python
from optparse import OptionParser
import sys, string, pdb, sys, os

################################################################################
# extract_aa.py
#
# Make a fasta file of amino acid sequences from the Glimmer3 gene predictions.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    parser.add_option('-s', dest='seqs_file', help='Sequence file')
    parser.add_option('-p', dest='predict_file', help='Gene predictions file')
    parser.add_option('-o', dest='output_file', help='Output amino acid fasta file')
    (options,args) = parser.parse_args()

    if not options.seqs_file:
        parser.error('Must provide sequence file with -s')
    elif not options.predict_file:
        parser.error('Must provide prediction file with -p')

    if options.output_file:
        out_aa = open('%s.faa' % options.output_file, 'w')
        out_dna = open('%s.ffn' % options.output_file, 'w')
    else:
        (base,ext) = os.path.splitext(options.seqs_file)
        out_aa = open('%s.faa' % base, 'w')
        out_dna = open('%s.ffn' % base, 'w')

    frag_preds = get_preds(options.seqs_file, options.predict_file)

    header = ''
    for line in open(options.seqs_file):
        if line[0] == '>':
            if header:
                print_frag_genes(out_aa, out_dna, header, seq, frag_preds[header])

            header = line[1:].rstrip()
            seq = ''
        else:
            seq += line.rstrip()
    if header:
        print_frag_genes(out_aa, out_dna, header, seq, frag_preds[header])

    out_aa.close()
    out_dna.close()
        
   
################################################################################
# gene_compare
#
# For sorting genes by start
################################################################################
def gene_compare(x, y):
    return x.start - y.start


################################################################################
# get_preds
################################################################################
def get_preds(seqs_file, genepred_file):
    # get fragment lengths
    frag_lengths = {}
    for line in open(seqs_file):
        if line[0] == '>':
            header = line[1:].rstrip()
            frag_lengths[header] = 0
        else:
            frag_lengths[header] += len(line.rstrip())

    # process genes
    frag_preds = {}
    for line in open(genepred_file):
        if line[0] =='>':
            header = line[1:].rstrip()
            frag_preds[header] = []
            indel_plusminus = 0
        else:
            a = line.split()

            # indels
            insertions = []
            if len(a[5]) > 2:
                insertions = [int(x)-1 for x in a[5][2:].split(',')]
            deletions = []
            if len(a[6]) > 2:
                deletions = [int(x)-1 for x in a[6][2:].split(',')]
            substitutions = []
            if len(a[7]) > 2:
                substitutions = [int(x)-1 for x in a[7][2:].split(',')]

            if int(a[3]) > 0:
                # forward
                strand = 1
                start = int(a[1])-1+indel_plusminus
                indel_plusminus += len(deletions) - len(insertions)
                end = int(a[2])+indel_plusminus

                # partial on left
                start_codon = True
                if start < 0:
                    start_codon = False
                # partial on right
                stop_codon = True
                if end > frag_lengths[header]+indel_plusminus:
                    stop_codon = False
                
            else:
                # reverse
                strand = -1
                start = int(a[2])-1+indel_plusminus
                indel_plusminus += len(deletions) - len(insertions)
                end = int(a[1])+indel_plusminus

                # partial on left
                stop_codon = True
                if start < 0:
                    stop_codon = False

                # partial on right
                start_codon = True
                if end > frag_lengths[header]+indel_plusminus:
                    start_codon = False

            frag_preds[header].append(Pred(start, end, strand, start_codon, stop_codon, insertions, deletions, substitutions))

    for header in frag_preds:
        frag_preds[header].sort(gene_compare)

    return frag_preds


################################################################################
# predict_msa
#
# If there were predicted insertions or deletions, add them to the MSA
################################################################################
def predict_msa(preds, seq):
    frag_msa = [' ',' ',' '] + list(seq) + [' ',' ',' ']

    # combine indels from genes
    insertions = []
    deletions = []
    substitutions = []
    for p in preds:
        insertions += p.insertions
        deletions += p.deletions
        substitutions += p.substitutions
    insertions.sort()
    deletions.sort()
    substitutions.sort()

    del_len = len(deletions)
    ins_len = len(insertions)
    sub_len = len(substitutions)

    if del_len == ins_len == sub_len == 0:
        pred_msa = frag_msa

    else:
       i = 0 # index into insertions
       d = 0 # index into deletions
       s = 0 # index into substitutions
       p = 3 # index into pred msa
       # m = 0 # index into frag msa (unadjusted)
       f = 0 # index into frag seq 

       pred_msa = [' ']*(len(frag_msa)+del_len)
       old_msa_len = len(frag_msa)

       for m in range(3,old_msa_len-3):
        if i < ins_len and insertions[i] == f:
            # insertion
            pred_msa[p] = '-'

            if frag_msa[p] != '-':
                f += 1
            p += 1
            i += 1

        elif d < del_len and deletions[d] == f:
            # deletion
            frag_msa.insert(p,'-')
            pred_msa[p] = pred_msa[p-1] # assuming deletions were homopolymer

            p += 1
            d += 1

            pred_msa[p] = frag_msa[p]

            if frag_msa[p] != '-':
                f += 1
            p += 1

        elif s < sub_len and substitutions[s] == f:
            # substitution (change stop codon)
            if frag_msa[p] == '-':
                print >> sys.stderr, 'Hit a gap where a substitution should be:'
                print >> sys.stderr, seq
                exit(1)
            elif frag_msa[p] == 'C':                
                pred_msa[p] = 'G'
            else:
                pred_msa[p] = 'C'
                
            f += 1
            p += 1
            s += 1
                
        else:
            # normal
            pred_msa[p] = frag_msa[p]

            if frag_msa[p] != '-':
                f += 1
            p += 1

    return pred_msa


################################################################################
# print_frag_genes
#
# Print all genes from this fragment
################################################################################
def print_frag_genes(out_aa, out_dna, header, seq, preds):
    pred_msa = predict_msa(preds, seq)

    for g in preds:
        gene_seq = ''

        s = -3
        gene_frame = 0
        for m in range(len(pred_msa)):
            if pred_msa[m] != '-':
                # forward
                if g.strand == 1:
                    # start
                    if g.start <= s < g.start+3:
                        gene_frame = 1
                        # continue or add only if it's the first base of the start and within the sequence
                        if len(gene_seq) > 0 or (s == g.start and s >= 0):
                            gene_seq += pred_msa[m]

                    # end
                    elif g.end-3 <= s < g.end:
                        gene_frame = 0

                    # middle
                    elif gene_frame > 0:
                        gene_frame = 1 + (gene_frame % 3)
                        # continue or add only if it's the first base of a codon and within the sequence
                        if len(gene_seq) > 0 or (gene_frame == 2 and s >= 0):
                            gene_seq += pred_msa[m]

                # reverse
                else:
                    # end
                    if g.start <= s < g.start+3:
                        gene_frame = 9
                        # continue or add only if it's the first base of the start and within the sequence
                        #if len(gene_seq) > 0 or (s == g.start and s >= 0):
                        #    gene_seq += pred_msa[m]

                    # start
                    elif g.end-3 <= s < g.end:
                        gene_frame = 0
                        #if s < len(seq):
                        if pred_msa[m] != ' ':
                            gene_seq += pred_msa[m]
                    
                    # middle
                    elif gene_frame > 0:
                        gene_frame -= 1
                        if gene_frame == 6:
                            gene_frame = 9
                        # continue or add only if it's the first base of a codon and within the sequence
                        if len(gene_seq) > 0 or (gene_frame == 8 and s >= 0):
                            gene_seq += pred_msa[m]

                s += 1

        # trim end
        gene_seq = gene_seq[:3*(len(gene_seq)/3)]
        
        # orient
        if g.strand == 1:
            dna_seq = gene_seq
            strand = '+'
        else:
            dna_seq = rc(gene_seq)
            strand = '-'

        print >> out_aa, '>%s_%d,%d_%s\n%s' % (header, g.start, g.end, strand, translate(dna_seq))
        print >> out_dna, '>%s_%d,%d_%s\n%s' % (header, g.start, g.end, strand, dna_seq)



############################################################
# rc
#
# Reverse complement sequence
############################################################
def rc(seq):
    return seq.translate(string.maketrans("ATCGatcg","TAGCtagc"))[::-1]


############################################################
# translate
#
# Translate a dna sequence into an amino acid.  Attempts
# to maintain lowercase or uppercase.  If a codon contains
# both lowercase and uppercase, returns a lowercase codon.
############################################################
code = {     'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', \
             'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', \
             'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*', \
             'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', \
             'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', \
             'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', \
             'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', \
             'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', \
             'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', \
             'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', \
             'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', \
             'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', \
             'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', \
             'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', \
             'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', \
             'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G', \

             'ttt': 'f', 'tct': 's', 'tat': 'y', 'tgt': 'c', \
             'ttc': 'f', 'tcc': 's', 'tac': 'y', 'tgc': 'c', \
             'tta': 'l', 'tca': 's', 'taa': '*', 'tga': '*', \
             'ttg': 'l', 'tcg': 's', 'tag': '*', 'tgg': 'w', \
             'ctt': 'l', 'cct': 'p', 'cat': 'h', 'cgt': 'r', \
             'ctc': 'l', 'ccc': 'p', 'cac': 'h', 'cgc': 'r', \
             'cta': 'l', 'cca': 'p', 'caa': 'q', 'cga': 'r', \
             'ctg': 'l', 'ccg': 'p', 'cag': 'q', 'cgg': 'r', \
             'att': 'i', 'act': 't', 'aat': 'n', 'agt': 's', \
             'atc': 'i', 'acc': 't', 'aac': 'n', 'agc': 's', \
             'ata': 'i', 'aca': 't', 'aaa': 'k', 'aga': 'r', \
             'atg': 'm', 'acg': 't', 'aag': 'k', 'agg': 'r', \
             'gtt': 'v', 'gct': 'a', 'gat': 'd', 'ggt': 'g', \
             'gtc': 'v', 'gcc': 'a', 'gac': 'd', 'ggc': 'g', \
             'gta': 'v', 'gca': 'a', 'gaa': 'e', 'gga': 'g', \
             'gtg': 'v', 'gcg': 'a', 'gag': 'e', 'ggg': 'g' \
             }

def translate(dna):
    if len(dna) % 3 != 0:
        print 'DNA sequence is not have length divisible by 3.'
        return ''
    else:
        i = 0
        peptide = ''
        while i < len(dna):
            if code.has_key(dna[i:i+3]):
                peptide += code[dna[i:i+3]]
            else:
                peptide += 'X'
            i += 3
        return peptide


class Pred:
    def __init__(self, start, end, strand, start_codon, stop_codon, insertions, deletions, substitutions):
        self.start = start
        self.end = end
        self.strand = strand
        self.start_codon = start_codon
        self.stop_codon = stop_codon
        self.insertions = insertions
        self.deletions = deletions
        self.substitutions = substitutions

    def __str__(self):
        return '\t'.join([str(x) for x in [self.start,self.end,self.strand,int(self.start_codon),int(self.stop_codon)]])


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
