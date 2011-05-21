#!/usr/bin/env python
import os, subprocess, glob, sys, stat

############################################################
# install_glimmer.py
#
# Compile source and set paths for Glimmer-MG.
#
# Author: David Kelley
############################################################

prior_phymm_dir = ''

install_phymm = True
install_scimm = True
install_elph = True
install_glimmer = True

############################################################
# main
############################################################
def main():
    installdir = os.getcwd()

    ############################################
    # Phymm
    ############################################
    if install_phymm:
        if prior_phymm_dir:
            if not os.path.islink('phymm'):
                os.symlink(prior_phymm_dir, 'phymm')
        else:
            if not os.path.isdir('phymm'):
                os.mkdir('phymm')
            os.chdir('phymm')

            # download
            p = subprocess.Popen('curl -o phymmInstaller.tar.gz http://www.cbcb.umd.edu/software/phymm/phymmInstaller.tar.gz', shell=True)
            os.waitpid(p.pid, 0)

            # unzip
            p = subprocess.Popen('tar -xzvf phymmInstaller.tar.gz', shell=True)
            os.waitpid(p.pid, 0)

            # install
            p = subprocess.Popen('./phymmSetup.pl', shell=True)
            os.waitpid(p.pid, 0)

            os.chdir('..')

    ############################################
    # Scimm
    ############################################
    if install_scimm:
        # unzip
        if not os.path.isdir('scimm'):
            scimm_ver = glob.glob('scimm*tar.gz')[0][:-7]
            p = subprocess.Popen('tar -xzvf %s.tar.gz' % scimm_ver)
            os.waitpid(p.pid, 0)
            os.rename(scimm_ver, 'scimm')

        os.chdir('scimm')

        # compile IMM code
        os.chdir('glimmer3.02/src')
        p = subprocess.Popen('make clean; make', shell=True)
        os.waitpid(p.pid,0)

        os.chdir('../..')

        # set IMM links
        if not os.path.isfile('bin/simple-score'):
            os.symlink('../glimmer3.02/bin/simple-score','bin/simple-score')
        if not os.path.isfile('bin/build-icm'):
            os.symlink('../glimmer3.02/bin/build-icm', 'bin/build-icm')

        # set scimm bin variable
        p = subprocess.Popen('sed \'s,scimm_bin = "[a-zA-Z/]*",scimm_bin = "%s/bin",\' bin/scimm.py > sc.tmp' % installdir, shell=True)
        os.waitpid(p.pid, 0)
        os.rename('sc.tmp', 'bin/scimm.py')
        p = subprocess.Popen('chmod ug+x bin/scimm.py', shell=True)
        os.waitpid(p.pid,0)

        # set physcimm bin variable
        p = subprocess.Popen('sed \'s,phymmdir = "[a-zA-Z/]*",phymmdir = "%s/phymm",\' bin/physcimm.py ph.tmp' % installdir, shell=True)
        os.waitpid(p.pid, 0)
        os.rename('ph.tmp', 'bin/physcimm.py')
        p = subprocess.Popen('chmod ug+x bin/physcimm.py', shell=True)
        os.waitpid(p.pid,0)

        os.chdir('..')

    
    ############################################
    # ELPH
    ############################################
    if install_elph:
        # download
        p = subprocess.Popen('curl -o ELPH-1.0.1.tar.gz ftp://ftp.cbcb.umd.edu/pub/software/elph/ELPH-1.0.1.tar.gz', shell=True)
        os.waitpid(p.pid, 0)

        # unzip
        p = subprocess.Popen('tar -xzvf ELPH-1.0.1.tar.gz', shell=True)
        os.waitpid(p.pid, 0)

        # compile
        os.chdir('ELPH/sources')
        p = subprocess.Popen('make', shell=True)
        os.waitpid(p.pid, 0)
        os.chdir('../..')


    ############################################
    # Glimmer
    ############################################
    if install_glimmer:
        # compile
        os.chdir('src')
        p = subprocess.Popen('make clean; make', shell=True)
        os.waitpid(p.pid,0)
        os.chdir('..')

        set_awk_path()

        # build gene IMMs
        p = subprocess.Popen('scripts/train_all.py', shell=True)
        os.waitpid(p.pid, 0)

        # find informative genomes
        p = subprocess.Popen('scripts/informative_genomes.py', shell=True)
        os.waitpid(p.pid, 0)

        # make double icms
        p = subprocess.Popen('scripts/double_icms.py', shell=True)
        os.waitpid(p.pid, 0)


############################################################
# set_awk_path
############################################################
def set_awk_path():
    # find awk
    awk_path = ''
    for path in os.environ['PATH'].split(':'):
        if os.path.isfile(os.path.join(path,'awk')):
            awk_path = os.path.join(path,'awk')
            break

    if awk_path == '':
        print >> sys.stderr, 'Cannot find Awk'
        exit(1)
    else:
        # change all scripts
        for awkf in glob.glob('scripts/*.awk'):
            os.rename(awkf,awkf+'.tmp')

            newf = open(awkf, 'w')
            print >> newf, '#!%s -f' % awk_path

            oldf = open(awkf+'.tmp')
            line = oldf.readline()
            line = oldf.readline()
            while line:
                print >> newf, line,
                line = oldf.readline()
            oldf.close()

            newf.close()


############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
