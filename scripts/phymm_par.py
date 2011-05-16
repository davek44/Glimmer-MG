#!/usr/bin/env python
from optparse import OptionParser, SUPPRESS_HELP
import os, random, subprocess, time, sys

################################################################################
# phymm_par.py
#
# Run phymm in parallel by parallelizing over ICMs rather than splitting up the
# input file.  This is better than splitting up the sequence file because
# building the ICMs requires some time, so building every ICM for every split
# is costly.
#
# Requires my version of Phymm that takes the '-i' option.
################################################################################

scripts_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
phymm_dir = os.path.abspath('%s/../phymm' % scripts_dir)

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <seqs_file>'
    parser = OptionParser(usage)
    parser.add_option('-p', dest='proc', type='int', default=2, help='Number of processes to run [default: %default]')
    parser.add_option('-b', dest='no_blast', action='store_true', default=False, help='Do not run BLAST')
    parser.add_option('-c', dest='chr_only', action='store_true', default=False, help='Score with only chromosomes, not plasmids')
    parser.add_option('-i', dest='ignore_file', help='File of IMMs to ignore')

    # help='Run on CONDOR grid'
    parser.add_option('--condor', dest='condor', action='store_true', default=False, help=SUPPRESS_HELP)

    (options,args) = parser.parse_args()

    if len(args) != 1 or not os.path.isfile(args[0]):
        parser.error('Please provide sequence fasta file')
    else:
        seqsf = os.path.abspath(args[0])

    # move to phymm directory
    origdir = os.getcwd()
    os.chdir(phymm_dir)

    # get ICMs
    icms = sorted(os.listdir('.genomeData'))

    if options.no_blast:
        # if no blast, imm-based parallelization

        # build parallel commands
        if options.condor:
            (tmp_cmds,pids) = build_cmds_imm(seqsf, options.ignore_file, icms, options)
            cmds = build_condor_cmds(tmp_cmds)
        else:
            (cmds,pids) = build_cmds_imm(seqsf, options.ignore_file, icms, options)

        # execute commands in parallel
        exec_par(cmds, options.proc)

        # combine results
        combine_imm(seqsf, pids, origdir)

        # clean up intermediate files
        clean(seqsf, pids)

    else:
        # if blast, sequence-based parallelization
        
        # build parallel commands
        if options.condor:
            (tmp_cmds,pids) = build_cmds_seq(seqsf, options.ignore_file, icms, options)
            cmds = build_condor_cmds(tmp_cmds)
        else:
            (cmds,pids) = build_cmds_seq(seqsf, options.ignore_file, icms, options)

        # execute commands in parallel
        exec_par(cmds, options.proc)
        
        # combine results
        combine_seq(seqsf, pids, origdir)

        # clean up intermediate files
        clean(seqsf, pids)


############################################################
# exec_par
#
# Execute the commands in the list 'cmds' in parallel, but
# only running 'max_proc' at a time.
############################################################
def exec_par(cmds, max_proc):
    total = len(cmds)
    finished = 0
    running = 0
    p = []

    while finished + running < total:
        # launch jobs up to max
        while running < max_proc and finished+running < total:
            p.append(subprocess.Popen(cmds[finished+running], shell=True))
            #print 'Running %d' % p[running].pid
            running += 1

        # are any jobs finished
        new_p = []
        for i in range(len(p)):
            if p[i].poll() != None:
                running -= 1
                finished += 1
            else:
                new_p.append(p[i])

        # if none finished, sleep
        if len(new_p) == len(p):
            time.sleep(1)
        p = new_p

    # wait for all to finish
    for i in range(len(p)):
        p[i].wait()


############################################################
# max_i
#
# Find max and return index and value
############################################################
def max_i(lis):
    max_index = 0
    max_val = lis[0]
    for i in range(1,len(lis)):
        if lis[i] > max_val:
            max_index = i
            max_val = lis[i]

    return (max_val,max_index)


################################################################################
# build_cmds_imm
#
# Build unix phymm commands to run, particularly with respect to telling
# the processes which ICMs to ignore.
################################################################################
def build_cmds_imm(seqsf, ignoref, icms, options):
    # options
    if options.no_blast:
        dash_b = '-b'
    else:
        dash_b = ''
    if options.chr_only:
        dash_c = '-c'
    else:
        dash_c = ''

    # ignored ICMs
    ignored_icms = set()
    if ignoref:
        for line in open(ignoref):
            ignored_icms.insert(line.rstrip())
 
    # reverse complement sequence
    (prefix,suffix) = os.path.splitext(seqsf)
    rc_seqsf = prefix + '.revComp' + suffix
    if not os.path.isfile(rc_seqsf):
        p = subprocess.Popen('.scripts/revCompFASTA.pl %s' % seqsf, shell=True)
        os.waitpid(p.pid, 0)

    # work out ICM ranges
    icms_per = len(icms) / options.proc
    icms_done = 0
    icm_starts = []
    icm_ends = []
    for p in range(options.proc):
        icm_starts.append(icms_done)
        icm_ends.append(icm_starts[-1] + icms_per)
        icms_done += icms_per
    icm_ends[-1] = len(icms)

    # make commands
    cmds = []
    pids = [0]*options.proc
    for p in range(options.proc):
        # choose a random id for the sequence sym link
        pids[p] = random.randint(0,1000000)
        while os.path.isfile('seqs_%d.fa'%pids[p]):
            pids[p] = random.randint(0,1000000)
        os.symlink(seqsf, 'seqs_%d.fa'%pids[p])
        os.symlink(rc_seqsf, 'seqs_%d.revComp.fa'%pids[p])

        # ICMs to ignore
        ignore_icms = [icms[i] for i in range(len(icms)) if i < icm_starts[p] or i >= icm_ends[p] or icms[i] in ignored_icms]
        ignore_out = open('ignore_%d.txt'%pids[p], 'w')
        print >> ignore_out, '\n'.join(ignore_icms)
        ignore_out.close()

        cmds.append('%s/scoreReadsGlim.pl seqs_%d.fa %s %s -i ignore_%d.txt' % (scripts_dir,pids[p],dash_b,dash_c,pids[p]))

    return (cmds,pids)



################################################################################
# build_cmds_seq
#
# Build unix phymm commands to run, particularly with respect to splitting up
# the sequence file
################################################################################
def build_cmds_seq(seqsf, ignoref, icms, options):
    # options
    if options.no_blast:
        dash_b = '-b'
    else:
        dash_b = ''
    if options.chr_only:
        dash_c = '-c'
    else:
        dash_c = ''

    # count sequences
    num_seqs = 0
    for line in open(seqsf):
        if line[0] == '>':
            num_seqs += 1

    # choose partition end points
    seqs_per = num_seqs / options.proc
    part_counts = []
    for i in range(options.proc-1):
        part_counts.append(seqs_per)
    part_counts.append(num_seqs - (options.proc-1)*seqs_per)

    # choose a unique random id for each partition
    pids = [0]*options.proc
    part_outs = []
    for p in range(options.proc):
        pids[p] = random.randint(0,1000000)
        while os.path.isfile('seqs_%d.fa'%pids[p]):
            pids[p] = random.randint(0,1000000)
        part_outs.append(open('seqs_%d.fa' % pids[p], 'w'))

    # partition sequences to new files
    p = 0
    seqs_current = 0
    for line in open(seqsf):
        if line[0] == '>':
            seqs_current += 1
            if seqs_current > part_counts[p]:
                p += 1
                seqs_current = 0
        print >> part_outs[p], line,

    for p in range(options.proc):
        part_outs[p].close()

    cmds = []
    for p in range(options.proc):
        cmds.append('%s/scoreReadsGlim.pl seqs_%d.fa %s %s' % (scripts_dir,pids[p],dash_b,dash_c))

    return cmds, pids

################################################################################
# build_condor_cmds
#
# Take unix phymm commands to run, and adjust to run on CONDOR grid.
################################################################################
def build_condor_cmds(tmp_cmds):
    cmds = []

    for tcmd in tmp_cmds:
        cmds.append('runCmd -c "%s"' % tcmd)

    return cmds


################################################################################
# combine_seq
#
# Combine the phymm output from the multiple processes running on different
# sets of sequences into single output files in the original directory.
################################################################################
def combine_seq(seqsf, pids, origdir):
    # replicate phymm's naming scheme
    seqsf_suf = os.path.split(seqsf)[1]
    seqsf_flat = seqsf_suf.replace('.','_')

    phymm_raw_file = '%s/rawPhymmOutput_%s.txt' % (origdir,seqsf_flat)
    combine_phymm_raw_seq(phymm_raw_file, pids)

    phymm_results_file = '%s/results.01.phymm_%s.txt' % (origdir, seqsf_flat)
    combine_results_seq('01.phymm', phymm_results_file, pids)

    blast_raw_file = '%s/rawBlastOutput_%s.txt' % (origdir,seqsf_flat)
    combine_blast_raw_seq(blast_raw_file, pids)

    blast_results_file = '%s/results.02.blast_%s.txt' % (origdir, seqsf_flat)
    combine_results_seq('02.blast', blast_results_file, pids)

    results_file = '%s/results.03.phymmBL_%s.txt' % (origdir, seqsf_flat)
    combine_results_seq('03.phymmBL', results_file, pids)


################################################################################
# combine_results_seq
#
# Combine the results output files defined by step_str for sequenced-based
# parallelization
################################################################################
def combine_results_seq(step_str, results_file, pids):
    results_out = open(results_file, 'w')

    results_opens = []
    for pid in pids:
        results_opens.append(open('results.%s_seqs_%d_fa.txt' % (step_str,pid)))
        line = results_opens[-1].readline()

    print >> results_out, line,

    for p in range(len(pids)):
        line = results_opens[p].readline()
        while line:
            print >> results_out, line,
            line = results_opens[p].readline()
        results_opens[p].close()

    results_out.close()


################################################################################
# combine_blast_raw_seq
#
# Combine the raw Blast output files for sequenced-based parallelization
################################################################################
def combine_blast_raw_seq(raw_file, pids):
    raw_out = open(raw_file, 'w')

    for pid in pids:
        for line in open('rawBlastOutput_seqs_%d_fa.txt' % pid):
            print >> raw_out, line,

    raw_out.close()


################################################################################
# combine_phymm_raw_seq
#
# Combine the raw Phymm output files for sequenced-based parallelization
################################################################################
def combine_phymm_raw_seq(raw_file, pids):
    # raw
    raw_files = []
    for pid in pids:
        raw_files.append(open('rawPhymmOutput_seqs_%d_fa.txt' % pid))
    
    raw_out = open(raw_file, 'w')

    # print icms
    print >> raw_out, 'BEGIN_ICM_LIST'
    lines = ['']*len(pids)
    imm_count = 0
    for p in range(len(pids)):
        raw_files[p].readline()
        lines[p] = raw_files[p].readline()
    while not lines[0].startswith('END_ICM_LIST'):
        imm_count += 1
        print >> raw_out, lines[0],
        for p in range(len(pids)):
            lines[p] = raw_files[p].readline()            
    print >> raw_out, 'END_ICM_LIST\nBEGIN_READID_LIST'

    # print reads
    for p in range(len(pids)):
        line = raw_files[p].readline()
        line = raw_files[p].readline()
        while not line.startswith('END_READID_LIST'):
            print >> raw_out, line,
            line = raw_files[p].readline()            
    print >> raw_out, 'END_READID_LIST\nBEGIN_DATA_MATRIX'

    # print scores
    imm_scores = []
    for i in range(imm_count):
        imm_scores.append([])

    for p in range(len(pids)):
        imm = 0
        line = raw_files[p].readline()
        line = raw_files[p].readline()
        while not line.startswith('END_DATA_MATRIX'):
            imm_scores[imm] += line.split()
            imm += 1
            line = raw_files[p].readline()

    for imm in range(len(imm_scores)):
        print >> raw_out, '\t'.join(imm_scores[imm])
    print >> raw_out, 'END_DATA_MATRIX'

    raw_out.close()

    


################################################################################
# combine_imm
#
# Combine the phymm output from the multiple processes running on different
# sets of ICMs into single output files in the original directory.
################################################################################
def combine_imm(seqsf, pids, origdir):
    # replicate phymm's naming scheme
    seqsf_suf = os.path.split(seqsf)[1]
    seqsf_flat = seqsf_suf.replace('.','_')

    # progress, who cares about these
    #progress_out = open('%s/%s_progress.txt' % (origdir,seqsf_flat), 'w')
    #for pid in pids:
    #    for line in open('seqs_%d_fa_progress.txt' % pid):
    #        print >> progress_out, line,
    #progress_out.close()

    # raw
    raw_files = []
    for pid in pids:
        raw_files.append(open('rawPhymmOutput_seqs_%d_fa.txt' % pid))
    
    raw_out = open('%s/rawPhymmOutput_%s.txt' % (origdir,seqsf_flat), 'w')

    print >> raw_out, 'BEGIN_ICM_LIST'
    for p in range(len(pids)):
        line = raw_files[p].readline()
        line = raw_files[p].readline()
        while not line.startswith('END_ICM_LIST'):
            print >> raw_out, line,
            line = raw_files[p].readline()

    print >> raw_out, 'END_ICM_LIST\nBEGIN_READID_LIST'

    for p in range(len(pids)):
        line = raw_files[p].readline()
        line = raw_files[p].readline()
        while not line.startswith('END_READID_LIST'):
            if p == 0:
                print >> raw_out, line,
            line = raw_files[p].readline()

    print >> raw_out, 'END_READID_LIST\nBEGIN_DATA_MATRIX'

    for p in range(len(pids)):
        line = raw_files[p].readline()
        line = raw_files[p].readline()
        while not line.startswith('END_DATA_MATRIX'):
            print >> raw_out, line,
            line = raw_files[p].readline()

    print >> raw_out, 'END_DATA_MATRIX'

    raw_out.close()

    # results
    results_out = open('%s/results.01.phymm_%s.txt' % (origdir,seqsf_flat), 'w')
    results_in = [open('results.01.phymm_seqs_%d_fa.txt' % pid) for pid in pids]

    lines = [ri.readline() for ri in results_in]
    while lines[0]:        
        if lines[0].startswith('QUERY_ID'):
            print >> results_out, lines[0],
        else:
            scores = [float(l.split('\t')[2]) for l in lines]
            (max_score,max_line) = max_i(scores)
            print >> results_out, lines[max_line],
        lines = [ri.readline() for ri in results_in]

    results_out.close()
                        

################################################################################
# clean
#
# Clean up phymm temporary files.
################################################################################
def clean(seqsf, pids):
    # reverse complement sequence
    (prefix,suffix) = os.path.splitext(seqsf)
    rc_seqsf = prefix + '.revComp' + suffix
    if os.path.isfile(rc_seqsf):
        os.remove(rc_seqsf)

    for pid in pids:
        os.remove('seqs_%d.fa'%pid)
        os.remove('seqs_%d.revComp.fa'%pid)
        os.remove('seqs_%d_fa_progress.txt'%pid)
        os.remove('rawPhymmOutput_seqs_%d_fa.txt'%pid)
        os.remove('results.01.phymm_seqs_%d_fa.txt'%pid)
        if os.path.isfile('rawBlastOutput_seqs_%d_fa.txt'%pid):
            os.remove('rawBlastOutput_seqs_%d_fa.txt'%pid)
            os.remove('results.02.blast_seqs_%d_fa.txt'%pid)
            os.remove('results.03.phymmBL_seqs_%d_fa.txt'%pid)
        if os.path.isfile('ignore_%d.txt'%pid):
            os.remove('ignore_%d.txt'%pid)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
