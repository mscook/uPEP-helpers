"""
fabfile for doing uPEP things

Mitchell Stanton-Cook
m.stantoncook@gmail.com
"""

import os
import sys
import gzip
import bisect

import ftputil
from   fabric.api import task


HOST = "ftp.ncbi.nlm.nih.gov"
BASE = "/refseq/release/"

def setup(outpath, home, db):
    """
    Build outpath and chdir to it

    :param outpath: the outpath
    :param home: the directory this fabfile is executed in
    :param db: the databse identifier

    :rtype: the fullpath to the outpath as a string
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


def download_db(db):
    """
    Download a db file

    :param db: the database identifier
    """
    name = db.split('-')[-1]
    url = BASE+name
    host = ftputil.FTPHost(HOST, 'anonymous', 'beaton.lab@gmail.com')
    host.chdir(url)
    remote = host.listdir(host.curdir)
    for f in remote:
        if f.endswith("rna.gbff.gz"):
            host.download(f, f, 'b')

@task
def get_NCBI_RefSeq_release():
    """
    Prints & returns the NCBI RefSeq release number
    
    :rtype: the release number as a string
    """
    url = BASE+"release-notes/"
    host = ftputil.FTPHost(HOST, 'anonymous', 'beaton.lab@gmail.com')
    host.chdir(url)
    remote = host.listdir(host.curdir)
    number = -1
    for f in remote:
        if f.startswith("RefSeq-"):
            number = f.split('release')[-1][:-4]
            break
    print 'Release %s' % (number)
    return number

@task
def get_uPEP_RefSeq_release(location='/var/RefSeq/.refseq_version'):
    """
    Prints @ returns the uPEP RefSeq release number
    
    :param location: [default=/var/RefSeq/.refseq_version] the file that the 
                     RefSeq version is stored in

    :rtype: the release number as a string
    """
    number = -1
    with open(location) as fin:
        number = fin.readlines()[0].strip()
    print 'Release %s' % (number)
    return number

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
    if key:
        print "Working with database "+key
        if key in dbs:
            wd = setup(outpath, home, key)
            download_db(key)
            compacted = compact_RefSeq(wd)
            compile_RefSeq(compacted)
            os.chdir(home)
        else:
            print "Not a defined db"
            sys.exit(1)
    else:
        for db in dbs:
            print "Working with database"+db
            setup(outpath, home, db)
            download_db(db)
            compacted = compact_RefSeq(wd)
            compile_RefSeq(compacted)
            os.chdir(home)

# The Follwoing are from Adam Skarshewski
#   compaction
#   compilation
#   uPEPfinder


@task
def compact_RefSeq(RefSeq_directory):
    """
    Removes 'junk' from the RefSeq files

    Originally:

    mkdir {directory}
    cd {directory}
    compaction.py ../{RefSeq_directory}/
    cd ..
    """
    stored = RefSeq_directory.split('-')[-1]
    try:
        os.mkdir("../"+stored)
    except OSError:
        print "Output directory exists"
        check = raw_input("Overwrite [y/n]?")
        if check != 'y':
            sys.exit(0)
    os.chdir("../"+stored)
    n=0;
    for root, dirs, files in os.walk(RefSeq_directory):
        for filename in files:
            filepath = root + '/' + filename
            print "Compacting "+filename
            if filepath[-3:] == '.gz':
                readfile = gzip.open(filepath,'rb', 9)
                writefile = gzip.open(filename,'wb', 9)
                #writefile = open(filename + '.txt','wb')
                skipflag = False
                line = readfile.readline()
                while True:
                    if not(line):
                        break
                    if line[0:5] == 'LOCUS':
                        n += 1
                        j = 0
                        templine = (line[5:]).lstrip()
                        for i in templine:
                            if i == ' ':
                                identifier = templine[0:j]
                                break
                            j += 1
                        if not(identifier[0:2] == 'XM' or                     \
                                    identifier[0:2] == 'NM'):
                            print identifier
                            while True:
                                line = readfile.readline()
                                if not line:
                                    break
                                if line.rstrip() == '//':
                                    break
                            continue
                        while True:
                            writefile.write(line)
                            line = readfile.readline()
                            if line[0:9] == 'REFERENCE' or                    \
                                    line[0:8] == 'FEATURES':
                                break
                        continue
                    if line[0:8] == 'FEATURES':
                        writefile.write(line)
                        line = readfile.readline()
                        while True:
                            if line[5:8]=='CDS' or line[5:11]=='source':
                                while True:
                                    writefile.write(line)
                                    line = readfile.readline()
                                    if not(line[0:6] == '      '):
                                        break
                                continue
                            if line[0:6] == 'ORIGIN':
                                while True:
                                    writefile.write(line)
                                    line = readfile.readline()
                                    if line == '//\n':
                                        writefile.write(line)
                                        break
                                break
                            line = readfile.readline()
                    line = readfile.readline()
                readfile.close()
                writefile.close()
    os.chdir("../")
    return stored


@task
def compile_RefSeq(compacted_dir):
    """
    Originally:
   
    compilation.py {compacted_dir}/
    """
    accessionlist = []
    GIlist = []
    n=0;
    for root, dirs, files in os.walk(compacted_dir):
        for filename in files:
            filepath = root + '/' + filename
            print "Compiling (compacted) "+filename
            if filepath[-3:] == '.gz':
                readfile = gzip.open(filepath,'rb', 9)
                skipflag = False
                pos = readfile.tell()
                while True:
                    line = readfile.readline()
                    if skipflag and not(line == '//\n'):
                        continue
                    #print line.rstrip()
                    if line[0:5] == 'LOCUS':
                        n+=1
                        temp = line[12:]
                        if (temp[:2] == 'NM') or (temp[:2] == 'XM'):
                            accession = temp[:temp.find(' ')]
                        else:
                            skipflag = True
                    elif line[0:7] == 'VERSION':
                        GI = line[line.find(':')+1:].rstrip()
                    #elif line[5:8] == 'CDS':
                        #temp = line[21:]
                        #print temp[:temp.find(' ')]
                    elif line == '//\n':
                        if not(skipflag):
                            bisect.insort(accessionlist, [accession, pos, 
                                                                    filepath])
                            bisect.insort(GIlist, [GI, accession])
                        skipflag = False
                        pos = readfile.tell()
                    if not(line):
                        break
                readfile.close()
    accessiondictionaryfilename = compacted_dir+                              \
                                    '-ACCcompletecompact.dict.rna.gbff'
    openfile = open(accessiondictionaryfilename, 'wb')
    for i in accessionlist:
        openfile.write('%s\t%s\t%s\n' %(i[0], i[1], i[2]))
    openfile.close()
    GIdictionaryfilename = compacted_dir+'-GIcompletecompact.dict.rna.gbff'
    openfile = open(GIdictionaryfilename,'wb')
    for i in GIlist:
        openfile.write('%s\t%s\n' %(i[0], i[1]))
    openfile.close()
    print "Working with %i elements" % (n)


def seqtocaps(seq):
    temp = ''
    for i in seq:
        if ord(i) > 96:
            temp = temp + chr(ord(i)-32)
        else: temp
    return temp

def returnjoins(string):
    temp = []
    for i in range(0,len(string)):
        if string[i] == '(':
            start = i+1
        if string[i] == ',':
            temp.append(string[start:i])
            start = i+1
        if string[i] == ')':
            temp.append(string[start:i])
    temp2 = []
    for i in temp:
        for j in range(0,len(i)):
            if i[j] == '.':
                temp2.append([i[0:j],i[j+2:]])
                break
    return temp2

def converttoint(string):
    if string[0] == '<':
        return 0
    elif string[0]=='>':
        return int(string[1:])
    else:
        return int(string)

def uPEP_finder():
    """
    To finish
    """
    specieslist = []
    for top, dirs, files in os.walk('.'):
        for directory in dirs:
            print directory
            if top == '.' and directory:
                savefile = open(directory + '.db', 'wb')
            for root, dirs, files in os.walk(directory):
                for filename in files:
                    print filename
                    try:
                        takejoins = True
                        scanfile = gzip.open(root + '/' +filename, 'rb', '9')
                        while True:
                            line = scanfile.readline()
                            if not line:
                                break
                            if line [0:5] == 'LOCUS':
                                n = 0
                                templine = (line[5:]).lstrip()
                                for i in templine:
                                    if i == ' ':
                                        identifier = templine[0:n]
                                        break
                                    n = n + 1
                            elif line[2:12] == 'ORGANISM  ':
                                if not(line[12:].rstrip() in specieslist):
                                    specieslist.append(line[12:].rstrip())
                                if not True: ##if not(line[12:].rstrip() in specieslist):
                                    while True:
                                        line = scanfile.readline()
                                        if not line:
                                            break
                                        if line.rstrip() == '//':
                                            break
                            elif line[5:8] == 'CDS':
                                if line[21:25] == 'join':
                                    if takejoins:
                                        more = returnjoins(line[25:].rstrip())
                                        CDSstart = converttoint(more[0][0])
                                        CDSend = converttoint(more[len(more)-1][1])
                                else:
                                    i = line[21:]
                                    for j in range(0,len(i)):
                                        if i[j] == '.':
                                            CDSstart = converttoint(i[0:j])
                                            CDSend = converttoint(i[j+2:].rstrip())
                                            break
                            elif line[0:6] == 'ORIGIN':
                                temp = ''
                                while True:
                                    line = scanfile.readline()
                                    if not line or line[0:2] == '//':
                                        break
                                    for i in line[10:].rstrip():
                                        if not i == ' ':
                                            temp = temp + i
                                record = [identifier, temp, CDSstart, CDSend]
                                uORFs = []
                                fail = False 
                                if record:
                                    CDSstart = record[2]
                                    for i in range(0, CDSstart):
                                        if not fail:
                                            if record[1][i] == 'a' and record[1][i+1] == 't' and record[1][i+2] == 'g':
                                                j = 0
                                                while (i + j) < CDSstart + min(100, (CDSend - CDSstart-5)):  #could be CDSstart + gracelength
                                                    j = j + 3
                                                    try:
                                                        if record[1][i+j] == 't' and (record[1][i+j+1] == 'a' and (record[1][i+j+2] == 'g' \
                                                            or record[1][i+j+2] == 'a') or (record[1][i+j+1] == 'g' and (record[1][i+j+2] == 'a'))):
                                                            uORFs.append([record[1][i:i+j],[i+1,i+j]])
                                                            break
                                                    except:
                                                        print identifier+" failed"
                                                        fail = True
                                                        break
                                if uORFs:
                                    j = 0
                                    while True:
                                        if len(uORFs[j][0]) < 15 or len(uORFs[j][0]) > 600:
                                            del uORFs[j]
                                            j = j - 1
                                        j = j + 1
                                        if j == len(uORFs):
                                            break
                                    if uORFs:
                                        for i in range(0, len(uORFs)):
                                            uPEPloc = uORFs[i][1]
                                            savefile.write('>gi|' + identifier + '_' + str(i+1)+ '|'+ str(uPEPloc)+'\n')
                                            uPEPtemp = seqtocaps(uORFs[i][0])
                                            while len(uPEPtemp) > 80:
                                                savefile.writelines(uPEPtemp[0:80] + '\n')
                                                uPEPtemp = uPEPtemp[80:]
                                            savefile.writelines(uPEPtemp + '\n')
                            if not line:
                                break
                    finally:
                        scanfile.close()
            savefile.close()
    return specieslist



@task 
def deploy(outpath=None, key=None):
    """
    """
    get_dbs(outpath, key)

