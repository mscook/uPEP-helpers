"""
fabfile for doing uPEP things

Mitchell Stanton-Cook
m.stantoncook@gmail.com
"""

import os
import sys
import ftputil

from   fabric.api import task


def setup(outpath, home, db):
    """
    Build outpath and chdir to it
    """
    if not outpath:
        outpath = os.path.join(home, db)
    else:
        outpath = os.path.join(os.path.expanduser(outpath), db)
    try:
        os.mkdir(outpath)
    except OSError:
        print "Output directory exists"
        check = raw_input("Overwrite [y/n]?")
        if check != 'y':
            sys.exit(0)
    os.chdir(outpath)
    return outpath


def download(db):
    """
    Download a db file
    """
    HOST = "ftp.ncbi.nlm.nih.gov"
    BASE = "/refseq/release/"
    name = db.split('-')[-1]
    url = BASE+name
    host = ftputil.FTPHost(HOST, 'anonymous', 'beaton.lab@gmail.com')
    host.chdir(url)
    remote = host.listdir(host.curdir)
    for f in remote:
        if f.endswith("rna.gbff.gz"):
            host.download(f, f, 'b')




@task
def get_dbs(outpath=None, key=None):
    """
    Get the databases
    """
    home = os.getcwd()
    dbs = ['RefSeq-complete',
           'RefSeq-fungi'
           'RefSeq-invertebrate',
           'RefSeq-microbial',
           'RefSeq-mitochondrion',
           'RefSeq-plant',
           'RefSeq-plasmid',
           'RefSeq-plastid',
           'RefSeq-protozoa',
           'RefSeq-vertebrate_mammalian',
           'RefSeq-vertebrate_other',
           'RefSeq-viral']
    print key
    if key:
        if key in dbs:
            setup(outpath, home, key)
            download(key)
        else:
            print "Not a defined db"
            sys.exit(1)
    else:
        for db in dbs:
            setup(outpath, home, db)
            download(db)
    os.chdir(home)

