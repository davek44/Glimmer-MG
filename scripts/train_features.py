#!/usr/bin/env python
from optparse import OptionParser
from copy import deepcopy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import pdb, sys, string, os, subprocess

################################################################################
# train_features.py
#
# ...
#
# start:X, end:Y in the biopython gbk file refers to seq[X:Y].
#
# lengths.genes.txt/lengths.non.txt
# Gene length models in format '<length>\t<count>'.
#
# starts.genes.txt/starts.non.txt
# Start codon models in format '<start_codon>\t<count>'.
#
# adj_orients.genes.txt/adj_orients.non.txt
# Adjacent gene orientation model in format '+-1 +-1\t<count>'
#
# adj_dist.+-1.+-1.genes.txt/adj_dist.+-1.+-1.non.txt
# Adjacent gene distances given orientation model in format '<length>\t<count>'
################################################################################

scripts_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
bin_dir = os.path.abspath('%s/../bin' % scripts_dir)
elph_bin = os.path.abspath('%s/../ELPH/sources/elph' % scripts_dir)

forward_start_codons = ['ATG','GTG','TTG']
forward_stop_codons = ['TAG','TAA','TGA']

TESTING = False

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)
    parser.add_option('-l','--min_length', dest='min_length', type='int', default=75, help='Minimum length of gene (and ORF) in nucleotides to consider [default: %default]')
    parser.add_option('-o','--max_overlap', dest='max_overlap', type='int', default=50, help='Maximum overlap of two genes (or gene and ORF) to consider [default: %default]')
    parser.add_option('--gbk', dest='gbk_file', help='Gbk file to train on')
    parser.add_option('--predict', dest='predict_file', help='Glimmer predictions to train on')
    parser.add_option('--seq','--seqs', dest='seq_file', help='Sequence(s) file predictions were made on')
    parser.add_option('-f', dest='print_featurefile', action='store_true', help='Print a feature file for direct input to Glimmer [Default: %default]')
    parser.add_option('-z', dest='mycoplas', action='store_true', default=False, help='Ignore TGA stop codon a la Mycoplasma. Will check .gbk file for this automatically. [Default: %default]')
    parser.add_option('--rbs', dest='rbs', action='store_true', default=False, help='Just re-compute the RBS model [Default: %default]')
    parser.add_option('--icm', dest='icm', action='store_true', default=False, help='Just re-compute the gene ICM model [Default: %default]')
    parser.add_option('--indels', dest='indels', action='store_true', default=False, help='Gene predictions contain indels')
    parser.add_option('--min_icm', dest='min_icm', type='int', default=0, help='Minimum bp required to train an ICM [Default: %default]')
    (options,args) = parser.parse_args()

    if options.mycoplas:
        forward_stop_codons[2] = 'XXX'

    if options.gbk_file:
        (genes, seqs, hypothetical) = parse_gbk(options.gbk_file)
        out_prefix = os.path.splitext(options.gbk_file)[0]

    elif options.predict_file and options.seq_file:
        (genes, seqs) = parse_predict(options.predict_file, options.seq_file)
        hypothetical = {}
        out_prefix = os.path.splitext(options.predict_file)[0]

    else:
        parser.error('Must provide either a gbk file (--gbk) or a Glimmer ' +
                     'predict file (--predict) and sequence fasta file (--seq)')

    print >> sys.stderr, out_prefix    

    if options.icm:
        if options.indels:
            build_icm_indels(options.seq_file, options.predict_file, out_prefix, options.min_icm)
        else:
            build_icm(genes, seqs, hypothetical, out_prefix, options.min_icm)
        if options.rbs:
            rbs_model(genes, seqs, hypothetical, out_prefix)
        return
    elif options.rbs:
        rbs_model(genes, seqs, hypothetical, out_prefix)
        return

    # initialize
    gene_stats = init_stats()
    nongene_stats = init_stats()

    # count
    parse_genes(gene_stats, genes, seqs, hypothetical, options.min_length, options.max_overlap)
    parse_nongenes(nongene_stats, genes, seqs, options.min_length, options.max_overlap)

    # fix
    destrand_orientations(gene_stats)
    destrand_orientations(nongene_stats)

    # output
    if options.print_featurefile:
        feature_file = open('%s.features.txt' % out_prefix, 'w')
        output_featurefile(feature_file, gene_stats, 'GENE', options.min_length, options.max_overlap)
        output_featurefile(feature_file, nongene_stats, 'NON', options.min_length, options.max_overlap)
        feature_file.close()
        rbs_model(genes, seqs, hypothetical, out_prefix)
        if options.indels:
            build_icm_indels(options.seq_file, options.predict_file, out_prefix, options.min_icm)
        else:
            build_icm(genes, seqs, hypothetical, out_prefix, options.min_icm)
    else:
        output_stats(out_prefix, gene_stats, 'gene', options.min_length, options.max_overlap)
        output_stats(out_prefix, nongene_stats, 'nongene', options.min_length, options.max_overlap)
        rbs_model(genes, seqs, hypothetical, out_prefix)
        build_icm(genes, seqs, hypothetical, out_prefix, options.min_icm)
        open('%s.gc.txt' % out_prefix, 'w').write('%f\n' % compute_gc(seqs))


################################################################################
# parse_gbk
#
# Extract a list of genes and the genome sequence from a gbk file.
#
# Skipping pseudogenes and genes circling the origin.
#
# Noting hypothetical proteins.
################################################################################
def parse_gbk(gbk_file):
    gbk = SeqIO.read(open(gbk_file), 'genbank')

    genes = {gbk.id:[]}
    hypothetical = {}
    for seqf in gbk.features:
        if seqf.type == 'CDS' and seqf.location.nofuzzy_start < seqf.location.nofuzzy_end:
            # pseudogene
            if seqf.qualifiers.has_key('note') and seqf.qualifiers['note'][0].find('pseudo') != -1:
                continue

            # global TGA?
            if seqf.qualifiers.has_key('transl_table') and seqf.qualifiers['transl_table'] == ['4']:
                forward_stop_codons[2] = 'XXX'

            # make gene
            g = Gene(seqf.location.nofuzzy_start, seqf.location.nofuzzy_end, seqf.location.nofuzzy_start, seqf.location.nofuzzy_end, seqf.strand, True, True)
            genes[gbk.id].append(g)

            # hypothetical
            if seqf.qualifiers.has_key('product') and seqf.qualifiers['product'][0].find('hypothetical') != -1:
                hypothetical[g.start] = True

    seqs = {gbk.id:str(gbk.seq)}

    return genes, seqs, hypothetical


################################################################################
# parse_predict
#
# Extract a list of genes and the genome sequence from a Glimmer .predict
# file and a fasta file of the sequence(s).
#
# Skipping partial genes.
################################################################################
def parse_predict(predict_file, seq_file):
    seqs = {}
    for line in open(seq_file):
        if line[0] == '>':
            header = line[1:].rstrip()
            seqs[header] = ''
        else:
            seqs[header] += line.rstrip()

    genes = {}
    for line in open(predict_file):
        if line[0] == '>':
            header = line[1:].rstrip()
        else:
            a = line.split()
            
            if int(a[3]) > 0:
                strand = 1
                start = int(a[1])-1
                end = int(a[2])
                start_codon = (start >= 0)
                stop_codon = (end <= len(seqs[header]))
                frame_start = start + 3*(1-int(start_codon))
                frame_end = end - 3*(1-int(stop_codon))
            else:
                strand = -1
                start = int(a[2])-1
                end = int(a[1])                
                stop_codon = (start >= 0)
                start_codon = (end <= len(seqs[header]))
                frame_start = start + 3*(1-int(stop_codon))
                frame_end = end - 3*(1-int(start_codon))            

            g = Gene(max(0, start), min(end, len(seqs[header])), frame_start, frame_end, strand, start_codon, stop_codon)
            genes.setdefault(header, []).append(g)

    return genes, seqs


################################################################################
# init_stats
#
# Initialize gene model counts
################################################################################
def init_stats():
    lengths = {}
    start_codons = dict.fromkeys(forward_start_codons, 0)
    adj_orients = {(1,1):0, (1,-1):0, (-1,1):0, (-1,-1):0}
    adj_dist = {(1,1):{}, (1,-1):{}, (-1,1):{}, (-1,-1):{}}
    return {'start_codons':start_codons, 'lengths':lengths, 'adj_orients':adj_orients, 'adj_dist':adj_dist}


################################################################################
# parse_genes
#
# Collect information from the GenBank file about the annotated genes.
#
# -Not going to even consider the gene if it's < min_length.
# -Will consider if it's > max_overlap, but won't count 'prev_distance'.
################################################################################
def parse_genes(stats, genes, seqs, hypothetical, min_length, max_overlap):
    ##################################################
    # preliminaries
    ##################################################
    if TESTING:
        test_out = open('test.genes.txt', 'w')

    for header in genes:
        hgenes = genes[header]
        hseq = seqs[header]

        last_strand = ''
        last_end = ''

        for gene in hgenes:
            # is it too short?
            gene_len = (gene.end-3 - gene.start)/3
            #if 3*gene_len < min_length:
            #    continue

            # count length
            if not hypothetical.has_key(gene.start):
                stats['lengths'][gene_len] = stats['lengths'].get(gene_len, 0) + 1

            # print gene sequence to file
            if gene.strand == 1:
                gene_seq = hseq[gene.start:gene.end]
            elif gene.strand == -1:
                gene_seq = rc(hseq[gene.start:gene.end])
            else:
                continue

            # count start
            if gene.start_codon and str(gene_seq[:3]) in forward_start_codons:
                stats['start_codons'][str(gene_seq[:3])] += 1

            if last_strand != '':                
                # adjacent orientations
                orientation = (last_strand, gene.strand)
                stats['adj_orients'][orientation] += 1

                # does it overlap too much?
                prev_distance = gene.start - last_end   # end is one past stop
                if -prev_distance <= max_overlap:
                    # adjacent distance
                    stats['adj_dist'][orientation][prev_distance] = stats['adj_dist'][orientation].get(prev_distance,0) + 1

            if TESTING:
                print >> test_out, '%d - %d (%d)' % (gene.start,gene.end,3*gene_len)
                print >> test_out, str(gene_seq[:3])
                if last_strand:
                    print >> test_out, orientation
                    print >> test_out, prev_distance
                else:
                    print >> test_out, '-\n-'

            last_strand = gene.strand
            last_end = gene.end

    ##################################################
    # finish
    ##################################################
    if TESTING:
        test_out.close()


################################################################################
# parse_nongenes
#
# Collect information from the GenBank file about the non-genes.
################################################################################
def parse_nongenes(stats, genes, seqs, min_length, max_overlap):
    # parse forward
    forward_parse_nongenes(1, genes, seqs, min_length, max_overlap, stats['start_codons'], stats['lengths'], stats['adj_orients'], stats['adj_dist'])

    # reverse complement genes, seq
    (rgenes, rseqs) = reverse_complement_genes(genes, seqs)

    # parse (rc'd) forward
    forward_parse_nongenes(-1, rgenes, rseqs, min_length, max_overlap, stats['start_codons'], stats['lengths'], stats['adj_orients'], stats['adj_dist'])


################################################################################
# forward_parse_nongenes
#
# Parse GenBank file, collecting information about non-genes, only considering
# ORFs on the forward strand.
#
# I made a design choice here to say that I'm not going to differentiate
# between the true gene and the non-gene ORF being in the first or second
# position in an adjacency orientation.  If I did, I'd still have to hack
# things back together to make it usable for the dynamic programming gene
# parse algorithm because it's not clear in that case whether the first or
# second is the true gene.  Nevertheless, I still prefer this way of doing
# things to using random ORFs which would get quite complicated with respect
# to overlaps.
################################################################################
def forward_parse_nongenes(genome_strand, genes, seqs, min_length, max_overlap, start_codons, lengths, adj_orients, adj_dist):
    if TESTING:
        if genome_strand == 1:
            test_out = open('test.non.txt', 'w')
        else:
            seq_len = len(gbk.seq)
            test_out = open('test.non.txt', 'a')

    for header in genes:
        hseq = seqs[header]
        hgenes = genes[header]

        # reset adjacent genes
        preceeding_i = 0
        succeeding_i = 0

        # find stop codons
        stop_codons = [i for i in range(len(hseq)) if str(hseq[i:i+3]) in forward_stop_codons]
        stop_codons += [len(hseq), len(hseq)+1, len(hseq)+2]

        for stop_i in stop_codons:
            ###################################
            # update preceeding, succeeding
            ###################################
            # move up preceeding_i too far
            preceeding_i = max(preceeding_i, 0)
            while preceeding_i < len(hgenes) and hgenes[preceeding_i].end-3 < stop_i:   # end is one past stop
                preceeding_i += 1

            # if it exists, set it to succeeding_i
            if preceeding_i < len(hgenes):
                succeeding_i = preceeding_i
            else:
                succeeding_i = -1

            # move preceeding_i back by one (possibly to -1 if it doesn't exist)
            preceeding_i -= 1

            ###################################
            # succeeding overlap check
            ###################################
            if succeeding_i != -1:
                if hgenes[succeeding_i].end-3 == stop_i:
                    # ORF is a gene (will never be equal to preceeding)
                    continue
                succeeding_overlap = stop_i - hgenes[succeeding_i].start + 3
                # if overlap is reasonable and cds doesn't span origin which will always be no good
                if succeeding_overlap > max_overlap: # or hgenes[succeeding_i].start > hgenes[succeeding_i].end: (I dropped these above)
                    # succeeding overlap too great
                    continue

            ###################################
            # count # start codons
            ###################################
            num_starts = 0
            codon_i = stop_i
            while codon_i >= 0:
                # get next codon
                codon_i -= 3
                if codon_i >= 0:
                    codon = hseq[codon_i:codon_i+3]
                else:
                    codon = ''

                if str(codon) in forward_stop_codons:
                    # ORF is done
                    break

                elif codon == '' or str(codon) in forward_start_codons:
                    # check preceeding overlap
                    if preceeding_i != -1:
                        preceeding_overlap = hgenes[preceeding_i].end - codon_i # end is one past stop
                        if preceeding_overlap > max_overlap:
                            # preceeding overlap (and for all downstream starts) is too great
                            break

                    # is it long enough
                    nongene_len = (stop_i - codon_i)/3
                    if 3*nongene_len >= min_length:
                        num_starts += 1

            ###################################
            # search for nongene ORFs
            ###################################
            codon_i = stop_i
            while codon_i >= 0:
                # get next codon
                codon_i -= 3
                if codon_i >= 0:
                    codon = hseq[codon_i:codon_i+3]
                else:
                    codon = ''

                if str(codon) in forward_stop_codons:
                    # ORF is done
                    break

                elif codon == '' or str(codon) in forward_start_codons:
                    # check preceeding overlap
                    if preceeding_i != -1:
                        preceeding_overlap = hgenes[preceeding_i].end - codon_i # end is one past stop
                        if preceeding_overlap > max_overlap:
                            # preceeding overlap (and for all downstream starts) is too great
                            break

                    # is it too short
                    nongene_len = (stop_i - codon_i)/3
                    if 3*nongene_len < min_length:
                        # count length, but nothing else
                        lengths[nongene_len] = lengths.get(nongene_len, 0) + 1
                        continue

                    ###################################
                    # nongene ORF!
                    ###################################
                    # count length
                    lengths[nongene_len] = lengths.get(nongene_len, 0) + 1

                    # count start
                    if codon:
                        start_codons[str(codon)] += 1

                    if preceeding_i != -1:
                        # adjacent orientation
                        if genome_strand == 1:
                            pre_orientation = (hgenes[preceeding_i].strand, 1)
                        else:
                            pre_orientation = (-1, -1*hgenes[preceeding_i].strand)
                        adj_orients[pre_orientation] += 1.0/num_starts

                        # adjacent distance
                        pre_distance = codon_i - hgenes[preceeding_i].end # end is one past stop
                        adj_dist[pre_orientation][pre_distance] = adj_dist[pre_orientation].get(pre_distance,0) + 1.0/num_starts

                    if succeeding_i != -1:
                        if codon == '':
                            x = 7

                        # adjacent orientation
                        if genome_strand == 1:
                            suc_orientation = (1, hgenes[succeeding_i].strand)
                        else:
                            suc_orientation = (-1*hgenes[succeeding_i].strand, -1)
                        adj_orients[suc_orientation] += 1.0/num_starts

                        # adjacent distance
                        suc_distance = hgenes[succeeding_i].start - (stop_i+3)
                        adj_dist[suc_orientation][suc_distance] = adj_dist[suc_orientation].get(suc_distance,0) + 1.0/num_starts


                    if TESTING:
                        if genome_strand == 1:
                            print >> test_out, '[%d:%d]' % (codon_i,(stop_i+3))
                        else:
                            print >> test_out, '[%d:%d]' % (seq_len-(stop_i+3),seq_len-codon_i)
                        print >> test_out, nongene_len
                        print >> test_out, str(codon)
                        if preceeding_i != -1:
                            print >> test_out, pre_orientation
                            print >> test_out, pre_distance
                        else:
                            print >> test_out, '-\n-'
                        if succeeding_i != -1:                            
                            print >> test_out, suc_orientation
                            print >> test_out, suc_distance
                        else:
                            print >> test_out, '-\n-'

    if TESTING:
        test_out.close()


################################################################################
# reverse_complement_genes
#
# Reverse complement the sequence and fix the gene coordinates accordingly
################################################################################
def reverse_complement_genes(genes, seqs):
    rgenes = {}
    rseqs = {}
    for header in genes:
        # sequence
        rseqs[header] = rc(seqs[header])
        seq_len = len(rseqs[header])

        rgenes[header] = []
        for gene in genes[header][::-1]:
            g = Gene(seq_len - gene.end, seq_len - gene.start, seq_len - gene.frame_end, seq_len - gene.frame_start, -1*gene.strand, gene.start_codon, gene.stop_codon)
            rgenes[header].append(g)
    return rgenes, rseqs


############################################################
# rc
#
# Reverse complement sequence
############################################################
def rc(seq):
    return seq.translate(string.maketrans("ATCGatcg","TAGCtagc"))[::-1]


################################################################################
# compute_gc
################################################################################
def compute_gc(seqs):
    gc = 0
    at = 0
    for s in seqs.values():
        for i in range(len(s)):
            if s[i] == 'A' or s[i] == 'T':
                at += 1
            elif s[i] == 'C' or s[i] == 'G':
                gc += 1
    return (float(gc)/(float(at)+float(gc)))


################################################################################
# destrand_orientations
#
# Adjust adjacent orientation counts to account for both strands.  For example,
# anything that is (1,1) on the 5' strand is (-1,-1) on the 3' strand and vice
# versa.  Alternatively, (1,-1) and (-1,1) are the same on the 5' and 3'
# strands.  So I'm basically just averaging (1,1) and (-1,-1)
################################################################################
def destrand_orientations(stats):
    # adjacent orientations
    stats['adj_orients'][(1,1)] += stats['adj_orients'][(-1,-1)]
    stats['adj_orients'][(1,1)] /= 2.0
    stats['adj_orients'][(-1,-1)] = stats['adj_orients'][(1,1)]

    # adjacent distance
    for l in (stats['adj_dist'][(1,1)].keys() + stats['adj_dist'][(-1,-1)].keys()):
        stats['adj_dist'][(1,1)][l] = stats['adj_dist'][(1,1)].get(l,0) + stats['adj_dist'][(-1,-1)].get(l,0)
        stats['adj_dist'][(1,1)][l] /= 2.0
        stats['adj_dist'][(-1,-1)][l] = stats['adj_dist'][(1,1)][l]


################################################################################
# output_stats
#
# Print stats to files
################################################################################
def output_stats(outf, stats, orf_type, min_length, max_overlap):
    ##################################################
    # lengths
    ##################################################
    if orf_type == 'gene':
        out = open('%s.lengths.genes.txt' % outf, 'w')
    else:
        out = open('%s.lengths.non.txt' % outf, 'w')

    if stats['lengths']:
        #for l in range(int(.5+(min_length/3)), 1+max(stats['lengths'].keys())):
        for l in range(1+max(stats['lengths'].keys())):
            l_p = stats['lengths'].get(l, 0)
            print >> out, '%d\t%d' % (l, l_p)
    out.close()

    ##################################################
    # start_codons
    ##################################################
    if orf_type == 'gene':
        out = open('%s.starts.genes.txt' % outf, 'w')
    else:
        out = open('%s.starts.non.txt' % outf, 'w')

    for sc in forward_start_codons:
        sc_p = stats['start_codons'][sc]
        print >> out, '%s\t%d' % (sc, sc_p)
    out.close()

    ##################################################
    # adjacent orientation
    ##################################################
    if orf_type == 'gene':
        out = open('%s.adj_orients.genes.txt' % outf, 'w')
    else:
        out = open('%s.adj_orients.non.txt' % outf, 'w')

    for strand1 in [1,-1]:
        for strand2 in [1,-1]:
            ao_p = stats['adj_orients'][(strand1,strand2)]
            print >> out, '%d,%d\t%d' % (strand1,strand2,ao_p)
            #print >> out, '%d %d\t%f (%d)' % (strand1,strand2,ao_p,stats['adj_orients'][(strand1,strand2)])
    out.close()

    ##################################################
    # adjacent distances
    ##################################################
    for strand1 in [1,-1]:
        for strand2 in [1,-1]:
            if strand1 == -1 and strand2 == -1:
                continue

            if orf_type == 'gene':
                out = open('%s.adj_dist.%d.%d.genes.txt' % (outf,strand1,strand2), 'w')
            else:
                out = open('%s.adj_dist.%d.%d.non.txt' % (outf,strand1,strand2), 'w')

            if stats['adj_dist'][(strand1,strand2)]:
                for l in range(-max_overlap, 1+max(stats['adj_dist'][(strand1,strand2)].keys())):
                    l_p = stats['adj_dist'][(strand1,strand2)].get(l, 0)
                    print >> out, '%d\t%.1f' % (l,l_p)
            out.close()

################################################################################
# output_featurefile
#
# Print stats to file for direct input to Glimmer
################################################################################
def output_featurefile(out, stats, orf_type, min_length, max_overlap):
    ##################################################
    # lengths
    ##################################################
    print >> out, 'DIST LENGTH %s' % orf_type
    #for l in range(int(.5+(min_length/3)), 1+max(stats['lengths'].keys())):
    for l in range(1+max(stats['lengths'].keys())):
        l_p = stats['lengths'].get(l, 0)
        print >> out, '%d\t%d' % (l, l_p)
    print >> out, ''

    ##################################################
    # start_codons
    ##################################################
    print >> out, 'DIST START %s' % orf_type
    for sc in forward_start_codons:
        sc_p = stats['start_codons'][sc]
        print >> out, '%s\t%d' % (sc, sc_p)
    print >> out, ''

    ##################################################
    # adjacent orientation
    ##################################################
    print >> out, 'DIST ADJACENT_ORIENTATION %s' % orf_type
    for strand1 in [1,-1]:
        for strand2 in [1,-1]:
            ao_p = stats['adj_orients'][(strand1,strand2)]
            print >> out, '%d,%d\t%d' % (strand1,strand2,ao_p)
            #print >> out, '%d %d\t%f (%d)' % (strand1,strand2,ao_p,stats['adj_orients'][(strand1,strand2)])
    print >> out, ''

    ##################################################
    # adjacent distances
    ##################################################
    for strand1 in [1,-1]:
        for strand2 in [1,-1]:
            if strand1 == -1 and strand2 == -1:
                continue
            print >> out, 'DIST ADJACENT_DISTANCE_%d_%d %s' % (strand1,strand2,orf_type)
            if stats['adj_dist'][(strand1,strand2)]:
                for l in range(-max_overlap, 1+max(stats['adj_dist'][(strand1,strand2)].keys())):
                    l_p = stats['adj_dist'][(strand1,strand2)].get(l, 0)
                    print >> out, '%d\t%.1f' % (l,l_p)
            print >> out, ''


################################################################################
# rbs_model
#
# Make position weight matrix using ELPH for ribosomal binding sites
# upstream of genes.
################################################################################
def rbs_model(genes, seqs, hypothetical, out_prefix):
    rbs_len = 25

    # print RBS region to file
    rbs_out = open('%s.rbs.upstream' % out_prefix, 'w')
    for header in genes:
        hgenes = genes[header]
        hseq = seqs[header]

        for gene in hgenes:
            # skip hypothetical
            if hypothetical.has_key(gene.start):
                continue

            if gene.strand == 1:
                if gene.start >= rbs_len:
                    rbs_seq = hseq[gene.start-25:gene.start]
                    print >> rbs_out, '>%s\t%d %d\n%s' % (header,gene.start,gene.end,rbs_seq)
            elif gene.strand == -1:
                if gene.end <= len(hseq)-rbs_len:
                    rbs_seq = rc(hseq[gene.end:gene.end+25])
                    print >> rbs_out, '>%s\t%d %d\n%s' % (header,gene.start,gene.end,rbs_seq)
            else:
                print >> sys.stderr, 'Bad strand %s %s' % (header,strand)
    rbs_out.close()

    # elph
    if os.path.getsize('%s.rbs.upstream' % out_prefix) > 0:
        p = subprocess.Popen('%s %s.rbs.upstream LEN=6 2> /dev/null | %s/get-motif-counts.awk > %s.motif' % (elph_bin, out_prefix, scripts_dir, out_prefix), shell=True)
        os.waitpid(p.pid, 0)
    else:
        motif_out = open('%s.motif' % out_prefix, 'w')
        cols = (1, 1, 1, 1, 1, 1)
        print >> motif_out, '6'
        print >> motif_out, 'a %7d %7d %7d %7d %7d %7d' % cols
        print >> motif_out, 'c %7d %7d %7d %7d %7d %7d' % cols
        print >> motif_out, 'g %7d %7d %7d %7d %7d %7d' % cols
        print >> motif_out, 't %7d %7d %7d %7d %7d %7d' % cols
        motif_out.close()
    os.remove('%s.rbs.upstream' % out_prefix)


################################################################################
# build_icm
#
# Print gene sequence to file and train a 3-periodic ICM for gene prediction.
#
# Skip genes with frameshifts and hypothetical proteins
################################################################################
def build_icm(genes, seqs, hypothetical, out_prefix, icm_train_bp_t):
    gene_out = open('%s.gene.fasta' % out_prefix, 'w')
    bp_printed = 0
    for header in genes:
        hgenes = genes[header]
        hseq = seqs[header]

        for gene in hgenes:
            # skip hypothetical
            if hypothetical.has_key(gene.start):
                continue

            # check for frameshift?
            if gene.strand not in [-1,1]:
                print >> sys.stderr, 'Bad strand %s %s' % (header,strand)
            else:
                if gene.strand == 1:
                    gene_seq = hseq[gene.frame_start:gene.frame_end-3*int(gene.stop_codon)]
                else:
                    gene_seq = rc(hseq[gene.frame_start+3*int(gene.stop_codon):gene.frame_end])

                print >> gene_out, '>%s_%d-%d_%d%d\n%s' % (header,gene.start,gene.end,int(gene.start_codon),int(gene.stop_codon),gene_seq)
                bp_printed += len(gene_seq)
                
    gene_out.close()

    # filter by entropy density profile distance ratio
    #edp_t = 1.15
    #os.system('cat %s.gene.fasta | entropy-fasta > %s.gene.edp.fasta' % (out_prefix,out_prefix))
    #gene_out = open('%s.gene.fasta' % out_prefix, 'w')
    #for line in open('%s.gene.edp.fasta' % out_prefix):
    #    if line[0] == '>':
    #        a = line.split()
    #        if float(a[-1]) < edp_t:
    #            print_gene = True
    #        else:
    #            print_gene = False            

    #    if print_gene:
    #        print >> gene_out, line,
    #gene_out.close()

    # build-icm
    if bp_printed >= icm_train_bp_t:
        if os.path.isfile('%s.gicm' % out_prefix):
            os.remove('%s.gicm' % out_prefix)
        p = subprocess.Popen('%s/build-icm -r %s.gicm < %s.gene.fasta' % (bin_dir, out_prefix,out_prefix), shell=True)
        os.waitpid(p.pid,0)
    #os.remove('%s.gene.fasta' % out_prefix)


################################################################################
# build_icm_indels
#
# Print gene sequence to file and train a 3-periodic ICM for gene prediction.
################################################################################
def build_icm_indels(seq_file, predict_file, out_prefix, icm_train_bp_t):
    # print gene sequences to file
    p = subprocess.Popen('%s/extract_aa.py -s %s -p %s -o %s' % (scripts_dir,seq_file,predict_file,out_prefix), shell=True)
    os.waitpid(p.pid, 0)

    os.remove('%s.faa' % out_prefix)
    os.rename('%s.ffn' % out_prefix, '%s.gene.fasta' % out_prefix)

    bp_printed = 0
    for line in open('%s.gene.fasta' % out_prefix):
        if line[0] != '>':
            bp_printed += len(line.rstrip())

    if bp_printed >= icm_train_bp_t:
        p = subprocess.Popen('%s/build-icm -r %s.gicm < %s.gene.fasta' % (bin_dir,out_prefix,out_prefix), shell=True)
        os.waitpid(p.pid, 0) 


################################################################################
# Gene
################################################################################
class Gene:
    def __init__(self, start, end, frame_start, frame_end, strand, start_codon, stop_codon):
        self.start = start
        self.end = end
        self.frame_start = frame_start
        self.frame_end = frame_end
        self.strand = strand
        self.start_codon = start_codon
        self.stop_codon = stop_codon


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
