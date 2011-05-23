#!/usr/bin/env python
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
import os, sys, pdb, subprocess, glob, gzip

################################################################################
# glimmer-mg.py
#
# Run the full Glimmer-MG pipeline.
################################################################################

scripts_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
bin_dir = os.path.abspath('%s/../bin' % scripts_dir)
phymm_dir = os.path.abspath('%s/../phymm' % scripts_dir)

physcimm_bin = os.path.abspath('%s/../scimm/bin/physcimm.py' % scripts_dir)

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fasta file>'
    parser = OptionParser(usage)
    add_options(parser)
    (options,args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error('Must provide fasta file of sequences')
    else:
        sequence_file = args[0]

    if options.output_file:
        output_file = options.output_file
    else:
        output_file = os.path.splitext(sequence_file)[0]
    
    # classify
    if not options.class_done and not options.raw_done:
        p = subprocess.Popen('%s/phymm_par.py -b -p %d %s' % (scripts_dir, options.proc, sequence_file), shell=True)
        os.waitpid(p.pid, 0)

    # parse classifications
    class_file = '%s.class.txt' % output_file
    if not options.class_done:
        (sequence_classes,sequence_scores) = parse_phymm(sequence_file, options.top_hits, options.ignore)
        # output classifications
        class_out = open(class_file, 'w')
        for seq in sequence_classes:
            print >> class_out, '%s\t%s' % (seq,' '.join(sequence_classes[seq]))
        class_out.close()

    elif options.iterate != 0 and not options.single_cluster and not options.raw_done:
        parser.error('Cannot use --class for multiple iterations. We need the scores')

    # long-orfs
    if options.long_orfs:
        tt = get_transl_table(sequence_file, output_file)
        p = subprocess.Popen('%s/long-orfs -n -t 1.15 -z %s %s %s.longorfs' % (bin_dir, tt, sequence_file, output_file), shell=True)
        os.waitpid(p.pid, 0)
        p = subprocess.Popen('%s/extract -t %s %s.longorfs > %s.train' % (bin_dir, sequence_file, output_file, output_file), shell=True)
        os.waitpid(p.pid, 0)
        if options.iterate == 0:
            p = subprocess.Popen('%s/build-icm -r %s.icm < %s.train' % (bin_dir, output_file, output_file), shell=True)
        else:
            p = subprocess.Popen('%s/build-icm -r %s.run1.icm < %s.train' % (bin_dir, output_file, output_file), shell=True)
        os.waitpid(p.pid, 0)

    # determine glimmer options
    g3_cmd = glimmer_options(options)
    if options.quality_file:
        qual_str = '-q %s' % options.quality_file
    else:
        qual_str = ''

    # single iterations
    if options.iterate == 0:
        if options.long_orfs:
            p = subprocess.Popen('%s -m %s.icm -c %s %s %s %s' % (g3_cmd, output_file, class_file, qual_str, sequence_file, output_file), shell=True)
        else:
            p = subprocess.Popen('%s -c %s %s %s %s' % (g3_cmd, class_file, qual_str, sequence_file, output_file), shell=True)
        os.waitpid(p.pid, 0)    

    # multiple iterations
    else:
        # make initial predictions using classifications only
        if options.long_orfs:
            p = subprocess.Popen('%s -m %s.run1.icm -c %s %s %s %s.run1' % (g3_cmd, output_file, class_file, qual_str, sequence_file, output_file), shell=True)
        else:
            p = subprocess.Popen('%s -c %s %s %s %s.run1' % (g3_cmd, class_file, qual_str, sequence_file, output_file), shell=True)
        os.waitpid(p.pid, 0)

        # no clustering
        if options.single_cluster:
            repredict(g3_cmd, sequence_file, output_file, class_file, options.iterate, options.filter_t, options.all_features, options.indel, qual_str)

        # clustering
        else:
            if not options.clust_done:
                # cluster reads
                phymm_results_file = 'results.01.phymm_%s.txt' % sequence_file.replace('.','_')
                p = subprocess.Popen('%s -s %s -p %d -r %s --taxlevel %s --minbp_pct %f' % (physcimm_bin, sequence_file, options.proc, phymm_results_file, options.taxlevel, options.minbp_pct), shell=True)
                os.waitpid(p.pid, 0)

            elif len(glob.glob('cluster*fa')) == 0:
                print >> sys.stderr, 'Cluster fasta files not found. Exclude option --clust.'
                exit(1)                

            # repredict within clusters
            predict_out = open(output_file+'.predict', 'w')
            for clust_sequence_file in glob.glob('cluster*fa'):                
                cluster_repredict(g3_cmd, clust_sequence_file, class_file, output_file, options.iterate, options.filter_t, options.all_features, options.indel, options.quality_file)
                combine_predictions(predict_out, sequence_scores, clust_sequence_file, output_file)
            predict_out.close()


################################################################################
# add_options
################################################################################
def add_options(parser):
    ############################################
    # pipeline options
    ############################################
    parser.add_option('--iter', dest='iterate', type='int', default=1, help='Iterate Glimmer by re-training on the initial predictions, thus assuming the sequences came from the same organism. [Default=%default]')
    parser.add_option('--long_orfs', dest='long_orfs', default=False, action='store_true', help='Generate the ICM to be used for the initial prediction iteration using the glimmer long-orfs program.  [Default: %default]')
    parser.add_option('-o', dest='out', help='Prefix for output files. Default uses fasta file name.')
    parser.add_option('-p', dest='proc', type='int', default=1, help='Number of processes to run. [Default=%default]')
    parser.add_option('--single_cluster', dest='single_cluster', default=False, action='store_true', help='Rather than cluster the sequences using PhyScimm, treat the sequences as having come from a single genome/cluster. [Default=%default]')
    parser.add_option('-t', dest='top_hits', type='int', default=3, help='Number of top Phymm classifications to use for training.  [Default: %default]')

    # in iterative re-training, Glimmer score threshold
    parser.add_option('--filter', dest='filter_t', type='float', default=1.0, help=SUPPRESS_HELP)
    # add a suffix to the glimmer binary in order to run a different version
    parser.add_option('--glim_suffix', dest='glim_suffix', default='', help=SUPPRESS_HELP)
    # find the map.txt file and ignore each sequence's source ICM score
    parser.add_option('--ignore', dest='ignore', default=False, action='store_true', help=SUPPRESS_HELP)
    # retrain all features include length and adjacency models
    parser.add_option('--all_features', dest='all_features', default=False, action='store_true', help=SUPPRESS_HELP)

    ############################################
    # glimmer-mg options
    ############################################
    glim_group = OptionGroup(parser, 'Glimmer')
    #glim_group.add_option('-g', '--gene_len', dest='gene_len', type='int', default=75, help='Minimum gene length. Must match the training database. [Default: %default]')
    glim_group.add_option('-i','--indel', dest='indel', action='store_true', default=False, help='Predict genes in "indel-mode" where gene predictions may shift the coding frame, implicitly predicting an insertion or deletion in the sequence. [Default: %default]')
    glim_group.add_option('-q', dest='quality_file', help='Fasta file of Phred quality values matching up with the sequences fasta file to be used in "indel-mode" and "substitution-mode".')
    glim_group.add_option('-r', '--circular', dest='circular', action='store_true', default=False, help='Assume circular rather than linear genome, i.e., allow wraparound. [Default: %default]')
    glim_group.add_option('-s', '--sub', dest='sub', action='store_true', default=False, help='Predict genes in "substitution-mode" where gene predictions may predict a sequencing error in a stop codon and pass through it. [Default: %default]')
    glim_group.add_option('-u', '--fudge', dest='fudge', type='float', default=1.0, help='Value to be added to the log-likelihood ratio score of every ORF. [Default: %default]')
    parser.add_option_group(glim_group)

    ############################################
    # phymm options
    ############################################
    phymm_group = OptionGroup(parser, 'Phymm')
    phymm_group.add_option('--raw', dest='raw_done', default=False, action='store_true', help='Do not classify sequences with Phymm because output file already exists. [Default=%default]')
    phymm_group.add_option('--class', dest='class_done', default=False, action='store_true', help='Do not classify sequences with Phymm or parse raw Phymm output because class.txt file already exists. [Default=%default]')
    parser.add_option_group(phymm_group)

    ############################################
    # physcimm options
    ############################################
    scimm_group = OptionGroup(parser, 'PhyScimm')
    scimm_group.add_option('--clust', dest='clust_done', default=False, action='store_true', help='Do not cluster sequences with PhyScimm because output files already exists. [Default=%default]')
    scimm_group.add_option('--taxlevel', dest='taxlevel', default='family', help='Taxonomic level at which to cluster reads with Phymm. [Default=%default]')
    scimm_group.add_option('--minbp_pct', dest='minbp_pct', type='float', default=.01, help='Minimum proportion of bp assigned to a class to become a cluster. [Default=%default]')
    parser.add_option_group(scimm_group)


################################################################################
# classify
#
# Return a dict of the top hits for the sequence's scores
################################################################################
def classify(scores, genomes, top_hits):
    sequence_top_hits = []

    (max_score,max_i,min_score) = maximin(scores)

    # save top hit
    sequence_top_hits.append(genomes[max_i])
    # move score to worst
    scores[max_i] = min_score-1        

    for th in range(1,top_hits):
        # save top hit
        (max_score,max_i,min_tmp) = maximin(scores)

        if max_score >= min_score:
            # save top hit
            sequence_top_hits.append(genomes[max_i])
            # move score to worst
            scores[max_i] = min_score-1

    return sequence_top_hits


################################################################################
# cluster_repredict
#
# Retrain on initial gene predictions and re-predict genes within a cluster
################################################################################
def cluster_repredict(g3_cmd, sequence_file, all_class_file, all_output_file, iterations, filter_t, all_features, indels, quality_file):
    # hash cluster headers
    cluster_reads = set()
    for line in open(sequence_file):
        if line[0] == '>':
            cluster_reads.add(line[1:].split()[0])

    # make cluster class file
    output_file = '%s.%s' % (all_output_file,sequence_file[:-3])
    all_predict_file = '%s.run1.predict' % all_output_file

    class_file = '%s.class.txt' % output_file
    class_out = open(class_file, 'w')
    for line in open(all_class_file):
        a = line.split()
        if a[0] in cluster_reads:
            print >> class_out, line,
    class_out.close()

    # make cluster predict file
    predict_file = '%s.run1.predict' % output_file
    predict_out = open(predict_file, 'w')
    print_seq = False
    for line in open(all_predict_file):
        if line[0] == '>':
            header = line[1:].split()[0]
            if header in cluster_reads:
                print_seq = True
            else:
                print_seq = False
        if print_seq:
            print >> predict_out, line,
    predict_out.close()

    # make cluster quality file
    if quality_file:
        make_cluster_quality(cluster_reads, sequence_file, quality_file, output_file)
        qual_str = '-q %s.qual' % output_file
    else:
        qual_str = ''

    # repredict
    repredict(g3_cmd, sequence_file, output_file, class_file, iterations, filter_t, all_features, indels, qual_str)


################################################################################
# combine_predictions
#
# Combine predictions from the initial iteration and the final iteration based
# on the training data and fit of the cluster to the sequences.
################################################################################
def combine_predictions(predict_out, sequence_scores, clust_sequence_file, all_output_file):
    min_gene_bp = 80000
    min_clust_phymm_ratio = -.013

    output_file = '%s.%s' % (all_output_file, clust_sequence_file[:-3])

    # check for sufficient training data
    all_init = False
    gene_bp = 0
    for line in open('%s.run1.filt.gene.fasta' % output_file):
        if line[0] != '>':
            gene_bp += len(line.rstrip())
    if gene_bp < min_gene_bp:
        # print initial only
        for line in open('%s.run1.predict' % output_file):
            print >> predict_out, line,

    else:
        # get sequence lengths
        seq_lengths = {}
        for line in open(clust_sequence_file):
            if line[0] == '>':
                header = line[1:].rstrip()
                seq_lengths[header] = 0
            else:
                seq_lengths[header] += len(line.rstrip())

        # get cluster ICM to phymm ICM log likelihood ratios
        cluster = int(clust_sequence_file[clust_sequence_file.find('-')+1:clust_sequence_file.find('.')])
        sequence_ratios = {}
        for line in open('icm-%d.scores.tmp' % cluster):
            (header,score) = line.split('\t')
            header = header.rstrip()
            header_prefix = header.split()[0]
            if seq_lengths.has_key(header):
                sequence_ratios[header] = (float(score) - sequence_scores[header_prefix]) / seq_lengths[header]

        # get initial predictions
        init_preds = {}
        for line in open('%s.run1.predict' % output_file):
            if line[0] == '>':
                header = line[1:].rstrip()
                init_preds[header] = []
            else:
                init_preds[header].append(line)

        # get cluster predictions
        clust_preds = {}
        for line in open('%s.predict' % output_file):
            if line[0] == '>':
                header = line[1:].rstrip()
                clust_preds[header] = []
            else:
                clust_preds[header].append(line)

        # print final
        seq_headers = set(clust_preds.keys() + init_preds.keys())
        for header in seq_headers:
            print >> predict_out, '>%s' % header
            if sequence_ratios[header] < min_clust_phymm_ratio:
                for line in init_preds[header]:
                    print >> predict_out, line,
            else:
                for line in clust_preds[header]:
                    print >> predict_out, line,


############################################################
# data_integrity
#
# Check for uniqueness of headers.
############################################################
def data_integrity(fasta_file):
    sequences = {}
    for line in open(fasta_file):
        if line[0] == '>':
            r = line[1:].split()[0]
            if sequences.has_key(r):
                print 'Sorry, Phymm only considers fasta headers up to the first whitespace.  Please make these unique in your file.  E.g. %s is not unique' % r
                exit()
            sequences[r] = True


################################################################################
# filter_predictions
#
# Filter Glimmer predictions before a training set is constructed for the
# next iteration.
################################################################################
def filter_predictions(predict_file, filter_t):
    filt_out = open('%s.filt.predict' % os.path.splitext(predict_file)[0], 'w')
    for line in open(predict_file):
        if line[0] == '>':
            print >> filt_out, line,
        else:
            a = line.split()
            if float(a[4]) > filter_t:
                print >> filt_out, line,
    filt_out.close()


################################################################################
# get_transl_table
#
# Decide on a single GenBank translation table code for the sequences given
# based on their Phymm classifications.
################################################################################
def get_transl_table(sequence_file, output_file):
    # get sequence length
    seq_lens = {}
    for line in open(sequence_file):
        if line[0] == '>':
            header = line[1:].rstrip()
            seq_lens[header] = 0
        else:
            seq_lens[header] += len(line.rstrip())

    # get sequence classes
    class_file = '%s.class.txt' % output_file
    sequence_classes = {}
    for line in open(class_file):
        a = line.split()
        sequence_classes[a[0]] = a[1]

    # get class tables
    class_tables = {}
    for clas in [sequence_classes[seq] for seq in sequence_classes]:
        #print clas
        (strain,nc_num) = clas.split('|')
        tt_code = '11'
        gbk_file = '%s/.genomeData/%s/%s.gbk' % (phymm_dir,strain,nc_num)
        for line in open(gbk_file):
            tt_i = line.find('transl_table=')
            if tt_i != -1:
                tt_code = line[tt_i+13:].rstrip()
                break
        class_tables[clas] = tt_code

    # gather evidence for tables
    transl_tables = {}
    for seq in sequence_classes:
        transl_tables[class_tables[sequence_classes[seq]]] = transl_tables.get(class_tables[sequence_classes[seq]],0) + seq_lens[seq]

    # find best
    tt_items_flip = [(bp,tt) for (tt,bp) in transl_tables.items()]
    tt_items_flip.sort()

    return tt_items_flip[-1][1]


################################################################################
# glimmer_options
#
# Determine Glimmer options
################################################################################
def glimmer_options(options):
    cmd = '%s/glimmer-mg%s' % (bin_dir, options.glim_suffix)

    if options.indel or options.quality_file:
        cmd += ' -i'

    return cmd


################################################################################
# make_cluster_quality
#
# Make a quality value file for the given cluster fasta file, where the
# sequences match up in order.
################################################################################
def make_cluster_quality(cluster_reads, sequence_file, quality_file, output_file):
    # hash quality values
    quality_hash = {}
    for line in open(quality_file):
        if line[0] == '>':
            header = line[1:].split()[0]
            if header in cluster_reads:
                quality_hash[header] = ''
            else:
                header = ''
        elif header:
            quality_hash[header] += line

    # arrange in order to match sequence file
    quality_out = open('%s.qual' % output_file, 'w')
    for line in open(sequence_file):
        if line[0] == '>':
            header = line[1:].split()[0]
            if header in quality_hash:
                print >> quality_out, line + quality_hash[header],
            else:
                print >> sys.stderr, 'Missing quality values for %s in %s' % (header, sequence_file)
                exit(1)
    quality_out.close()


################################################################################
# maximin
#
# Return max, index of max, and min
################################################################################
def maximin(data):
    max_i = 0
    max_data = data[0]
    min_data = data[0]
    for i in range(1,len(data)):
        if data[i] > max_data:
            max_i = i
            max_data = data[i]
        if data[i] < min_data:
            min_data = data[i]
    return (max_data,max_i,min_data)


################################################################################
# parse_phymm
#
# Parse the Phymm raw scores for the top hits for each sequence
################################################################################
def parse_phymm(sequence_file, top_hits, ignore):
    sequence_organisms = {}
    if ignore:
        map_file = 'map.txt'
        if not os.path.isfile(map_file):
            map_file = 'map.err.txt'
        print >> sys.stderr, 'Using %s to ignore' % map_file
        for line in open(map_file):
            a = line.split()
            sequence_organisms[a[0]] = a[1].split('|')[0]
        
    informative_genomes = set()
    for line in open('%s/../data/informative_genomes.txt' % scripts_dir):
        informative_genomes.add(line.rstrip())

    raw_file = 'rawPhymmOutput_%s.txt' % sequence_file.replace('.','_')
    if os.path.isfile(raw_file):
        raw_open = open(raw_file)
    elif os.path.isfile(raw_file+'.gz'):
        raw_open = gzip.open(raw_file+'.gz', 'rb')
    else:
        print >> sys.stderr, 'Cannot find raw Phymm output file %s' % raw_file
        exit(1)

    # read all genome names
    line = raw_open.readline()
    line = raw_open.readline()
    genomes = []
    while not line.startswith('END_ICM_LIST'):
        a = line.split('/')
        genome = '%s|%s' % (a[-2],a[-1].split('.')[0])
        genomes.append(genome)
        line = raw_open.readline()

    # read all sequence names
    line = raw_open.readline()
    line = raw_open.readline()
    sequences = []
    sequence_scores = []
    while not line.startswith('END_READID_LIST'):
        sequences.append(line.rstrip())
        sequence_scores.append(['']*top_hits)
        line = raw_open.readline()

    # read informative scores and find top hits
    line = raw_open.readline()
    line = raw_open.readline()
    g = 0
    while not line.startswith('END_DATA_MATRIX'):
        if genomes[g] in informative_genomes:
            org = genomes[g].split('|')[0]
            genome_scores = line.split()
            for s in range(len(sequences)):
                if not ignore or not sequence_organisms.has_key(sequences[s]) or org != sequence_organisms[sequences[s]]:
                    score_insert(sequence_scores[s], float(genome_scores[s]), g)
        g += 1
        line = raw_open.readline()

    # get top genomes
    sequence_classes = {}
    sequence_top_scores = {}
    for s in range(len(sequences)):
        seq = sequences[s]
        sequence_top_scores[seq] = sequence_scores[s][0][0]
        sequence_classes[seq] = []
        for t in range(top_hits):
            g = sequence_scores[s][t][1]
            sequence_classes[seq].append(genomes[g])

    return sequence_classes, sequence_top_scores


################################################################################
# repredict
#
# Retrain on initial gene predictions and re-predict
################################################################################
def repredict(g3_cmd, sequence_file, output_file, class_file, iterations, filter_t, all_features, indels, qual_str):
    # train and re-predict
    for i in range(2,iterations+2):
        prev_iter = '%s.run%d' % (output_file,(i-1))
        if i < iterations:
            next_iter = '%s.run%d' % (output_file,i)
        else:
            next_iter = output_file

        retrain(sequence_file, prev_iter, filter_t, all_features, indels)
        p = subprocess.Popen('%s -b %s.filt.motif -m %s.filt.gicm -f %s.filt.features.txt -c %s %s %s %s' %
                  (g3_cmd, prev_iter, prev_iter, prev_iter, class_file, qual_str, sequence_file, next_iter), shell=True)
        os.waitpid(p.pid, 0)


################################################################################
# retrain
#
# Re-train gene prediction models on predictions from the previous iteration.
# Train ICM, RBS, and start codon models.
################################################################################
def retrain(sequence_file, prev_iter, filter_t, all_features, indels):
    indel_str = ''
    if indels:
        indel_str = '--indel'

    # filter
    filter_predictions('%s.predict' % prev_iter, filter_t)

    # train
    p = subprocess.Popen('%s/train_features.py -f %s --seq %s --predict %s.filt.predict' % (scripts_dir,indel_str,sequence_file,prev_iter), shell=True)
    os.waitpid(p.pid, 0)

    # only keep start codons
    if not all_features:
        feat_out = open('%s.filt.features.tmp' % prev_iter, 'w')
        be_printing = False
        for line in open('%s.filt.features.txt' % prev_iter):
            if line.startswith('DIST START'):
                be_printing = True
            elif line.startswith('DIST'):
                be_printing = False

            if be_printing:
                print >> feat_out, line,
        feat_out.close()
        os.rename('%s.filt.features.tmp' % prev_iter, '%s.filt.features.txt' % prev_iter)
             
 
################################################################################
# score_insert
#
# Consider inserting a (score,genome) tuple into the sorted list of prior 
# (score,genome) tuples. If it's lower than all of the scores, don't insert
################################################################################
def score_insert(score_list, score, g):
    # add if empty slot
    for i in range(len(score_list)):
        if score_list[i] == '':
            score_list[i] = (score,g)
            return

    # else find insertion point
    insert_point = 0
    while insert_point < len(score_list):
        if score > score_list[insert_point][0]:
            break
        else:
            insert_point += 1

    # if past end, no insertion
    if insert_point == len(score_list):
        return

    # else insert and slide others
    else:
        for i in range(len(score_list)-1,insert_point,-1):
            score_list[i] = score_list[i-1]
        score_list[insert_point] = (score,g)
      

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
