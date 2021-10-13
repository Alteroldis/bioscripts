#!/usr/bin/python3.6

from sys import argv, stdout
from os.path import getsize, abspath
from os import makedirs, chdir, listdir, remove
from Bio import SeqIO
import subprocess
import threading
import time

if (argv[1] == "-h" or argv[1] == "--help"):
    print("USAGE:   PKG-name.py cluster-file.clstr CDS-file.cds num-threads")
    exit(0)

file_homology = abspath(argv[1])
file_cds = abspath(argv[2])
THREADS = int(argv[3])

print("Processing %s file" % file_homology)
cmd_cap = "/home/alteroldis/Tools/CAP3/cap3 %s -r 0 -p 80 -o %s -y %s -f %s -t 500 -s 300 -i 32 -j 42 -x o%s > %s.o%s.log"

lib_prot2cdsid = {}
file = open(file_homology, "r")

count = 0
set_cds_id = set()
set_tmp = set()
for line in file:
        if ">Cluster" in line:
            if len( set_tmp ) > 1:
                lib_prot2cdsid[ "Cluster%s" % count ] = set_tmp
                lib_prot2cdsid[ "Cluster%s" % count ] = list( lib_prot2cdsid[ "Cluster%s" % count ] )
                set_cds_id.update( set_tmp )
            set_tmp = set()
            count += 1
        if "... " in line:
            set_tmp.add( line.split(">")[1].split("... ")[0] )
file.close()
print("Complete processing %s file" % file_homology)
print("Amount clusters: %s" % len(lib_prot2cdsid) )

print("Processing %s file" % file_cds)
lib_id2seq = {}
seq_singlets = []
seq_contigs = []
for seq in SeqIO.parse(file_cds, "fasta"):
    if seq.id in set_cds_id:
        lib_id2seq[ seq.id ] = seq
    else:
        seq_singlets.append( seq )

print("Complete processing %s file" % file_cds)
print("Amount sequences: %s" % len(lib_id2seq) )
print("Amount singlets: %s" % len(seq_singlets) )
        
makedirs("workdir_thread")
chdir("workdir_thread")

set_assembled_id = set()
completed = 0
set_processing_prot = set()
LOCK = threading.RLock()


def ASSEMBLE(lib_prot2cdsid, lib_id2seq):
    global seq_singlets, seq_contigs, set_assembled_id, set_processing_prot, completed
    for prot in lib_prot2cdsid:
        LOCK.acquire()
        if prot in set_processing_prot:
            LOCK.release()
            continue
        set_processing_prot.add( prot )
        print("\r    COMPLETED/ALL: %s/%s; Assembled: %s" % (completed, len(lib_prot2cdsid), len(seq_contigs) ), end='')
        stdout.flush()
        LOCK.release()
        seq_homology_cds = []
        list_lens = []
        file_seq = "%s.cds" % prot
        for cds_id in lib_prot2cdsid[ prot ]:
            seq_homology_cds.append( lib_id2seq[ cds_id ] )
            list_lens.append( len(lib_id2seq[ cds_id ]) )
        SeqIO.write( seq_homology_cds, file_seq, "fasta")
        list_lens.sort()
        overlap = sum( list_lens )/ len( list_lens )
        max_gap = overlap*0.03
        max_clip = list_lens[0]*0.2
        if max_clip > 120: max_clip = 120
        while overlap >= 40:
            PIPE=subprocess.PIPE
            #cmd_cap = "/home/alteroldis/HDD/Tools/CAP3/cap3 %s -r 0 -p 80 -o %s -y %s -f %s -t 500 -s 300 -i 32 -j 42 -x o%s > %s.o%s.log"
            process = subprocess.Popen(cmd_cap % (file_seq, overlap, max_clip, max_gap, overlap, file_seq, overlap), shell=True, stdout=PIPE)
            process.wait()
            if getsize("%s.o%s.contigs" % (file_seq, overlap)) != 0:
                count_i = 0
                seq_homology_cds = []
                for seq in SeqIO.parse("%s.o%s.contigs" % (file_seq, overlap), "fasta"):
                    seq.id = "Homolog_%s_i%s_length_%s" % (prot, count_i, len(seq) )
                    seq.name = ""
                    seq.description = ""
                    seq_homology_cds.append( seq )
                    count_i += 1
                if getsize("%s.o%s.singlets" % (file_seq, overlap)) == 0:
                    LOCK.acquire()
                    seq_contigs.extend( seq_homology_cds )
                    set_assembled_id.update( set(lib_prot2cdsid[ prot ]) )
                    completed += 1
                    LOCK.release()
                    for filename in listdir():
                        if file_seq in filename:
                            remove(filename)
                    break
                else:
                    for seq in SeqIO.parse("%s.o%s.singlets" % (file_seq, overlap), "fasta"):
                        if "Homolog" in seq.id:
                            seq.id = "Homolog_%s_i%s_length_%s" % (prot, count_i, len(seq) )
                            seq.name = ""
                            seq.description = ""
                            seq_homology_cds.append( seq )
                            count_i += 1
                        else:
                            seq_homology_cds.append( seq )
                SeqIO.write( seq_homology_cds, file_seq, "fasta")
            for filename in listdir():
                if "%s." % file_seq in filename:
                    remove(filename)
            overlap -= 100
        else:
            LOCK.acquire()
            for seq in seq_homology_cds:
                if "Homolog" in seq.id:
                    seq_contigs.append( seq )
                else:
                    if seq.id in set_assembled_id:
                        continue
                    seq_singlets.append( seq )
                    lib_prot2cdsid[ prot ].remove( seq.id )
            set_assembled_id.update( set(lib_prot2cdsid[ prot ]) )
            completed += 1
            LOCK.release()
            for filename in listdir():
                if file_seq in filename:
                    remove(filename)

for _ in range(THREADS):
    thread_ = threading.Thread(target=ASSEMBLE, args=(lib_prot2cdsid, lib_id2seq))
    thread_.start()
while threading.active_count() >1:
    time.sleep(1)

print("\nAssembling homology CDS finished")
seq_singlets_finish = []
set_singlets_writed = set()
for seq in seq_singlets:
    if seq.id in set_singlets_writed:
        continue
    if seq.id not in set_assembled_id:
        seq_singlets_finish.append( seq )
        set_singlets_writed.add( seq.id )
print("Count contigs = %s; count singlets = %s" % ( len(seq_contigs), len(seq_singlets_finish) ))
print("Writing singlets and contigs into files")
SeqIO.write(seq_singlets_finish, "%s.ThreadHomoloCAP3.singlets" % file_cds, "fasta")
SeqIO.write(seq_contigs, "%s.ThreadHomoloCAP3.contigs" % file_cds, "fasta")

print("Finish")
