#!/usr/bin/python3
import primer3
from sys import argv

if len(argv) < 2 or "-h" in argv or "--help" in argv:
    print( "USAGE default: OligoAnalyse.py ddPCR/qPCRmix-HS-SYBR file-with-primers 1Xdna_concentration" )
    print( "USAGE with multiplex reaction: ddPCR/qPCRmix-HS-SYBR OligoAnalyse.py file-with-primers \
        1Xdna_concentration multiplex list-primers-comma-separated" )
    exit()
if "file-example" in argv:
    print("any head line")
    print("SetName\tPrimerType\tSeq")
    print("PrimerType is line with type (forward, reverse, probe) of primer")
    exit()

def hairpin(seq, conc, dv_conc, dntp_conc):
    res = primer3.calcHairpin( seq, mv_conc=50, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=conc, \
        temp_c=30, max_loop=20 )
    gibbs = res.dg/1000
    return( gibbs )

def homodimer(seq, conc, dv_conc, dntp_conc):
    res = primer3.calcHomodimer( seq, mv_conc=50, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=conc, \
        temp_c=30, max_loop=20, output_structure=True )
    gibbs = res.dg/1000
    if res.ascii_structure != "":
        secstruct = res.ascii_structure.split("\n")[1].split("\t")[1].split()[0]
        if len(secstruct) > 2:
            if seq.rfind( secstruct ) == -1:
                gibbs += 1.5*len(secstruct)
    return( gibbs )

def heterodimer(seq1, seq2, conc, dv_conc, dntp_conc):
    res = primer3.calcHeterodimer( seq1, seq2, mv_conc=50, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=conc, \
        temp_c=30, max_loop=20, output_structure=True )
    gibbs = res.dg/1000
    if res.ascii_structure != "":
        seq1stroke = False; seq2stroke = False
        secstruct = res.ascii_structure.split("\n")[1].split("\t")[1].split()[0]
        if seq1.rfind( secstruct ) == -1:
            seq1stroke = True
        secstruct = res.ascii_structure.split("\n")[2].split("\t")[1].split()[0]
        if seq2.rfind( secstruct ) == -1:
            seq2stroke = True
        if len(secstruct) > 2:
            if seq1stroke and seq2stroke:
                gibbs += 1.5*len(secstruct)
            elif seq1stroke or seq2stroke:
                gibbs += 2*len(secstruct)
    return( gibbs )

conc = argv[3]
mp = False
if argv[1] == "qPCRmix-HS-SYBR":
    dv_conc=3
    dntp_conc=0.44
else:
    dv_conc=3.5
    dntp_conc=0.8
if "multiplex" in argv:
    refp = argv[5].split(",")
    mp = True
set2primers = {}
set2data = {}
set2gibbs = {}
file = open(argv[2], "r")
for line in file:
    header = line.split("\n")[0]
    break
for line in file:
    setp = line.split("\t")[0]
    seq = line.split("\t")[2]
    if setp not in set2primers:
        set2primers[ setp ] = []
        set2data[ setp ] = []
        set2gibbs[ setp ] = [[0], [0], [0]]
    if seq != "":
        set2primers[ setp ].append( seq )
    set2data[ setp ].append( line.split("\n")[0] )
file.close()

delete_set = set()
for setp in set2primers:
    for primer in set2primers[setp]:
        if hairpin(primer, conc, dv_conc, dntp_conc) < -2:
            delete_set.add( setp )
for setp in list(delete_set):
    set2primers.pop(setp)
    set2gibbs.pop(setp)
if len(set2primers) == 0:
    print("Hairpin stage: primers without dimers and hairpins does not exist")
    exit()
else:
    print( len(set2primers) )

delete_set = set()
for setp in set2primers:
    for primer in set2primers[setp]:
        gibbs = homodimer(primer, conc, dv_conc, dntp_conc)
        set2gibbs[ setp ][0].append( gibbs )
        if gibbs < -8:
            delete_set.add( setp )
for setp in list(delete_set):
    set2primers.pop(setp)
    set2gibbs.pop(setp)
if len(set2primers) == 0:
    print("Homodimer stage: primers without dimers and hairpins does not exist")
    exit()
else:
    print( len(set2primers) )

delete_set = set()
for setp in set2primers:
    length = len(set2primers[setp]) - 1
    for i in range(length):
        primer1 = set2primers[setp][i]
        for primer2 in set2primers[setp][i+1:]:
            gibbs = heterodimer(primer1, primer2, conc, dv_conc, dntp_conc)
            set2gibbs[ setp ][1].append( gibbs )
            if gibbs < -8:
                delete_set.add( setp )
for setp in list(delete_set):
    set2primers.pop(setp)
    set2gibbs.pop(setp)
if len(set2primers) == 0:
    print("Heterodimer stage: primers without dimers and hairpins does not exist")
    exit()
else:
    print( len(set2primers) )

delete_set = set()
if mp:
    for primer1 in refp:
        for setp in set2primers:
            for primer2 in set2primers[setp]:
                gibbs = heterodimer(primer1, primer2, conc, dv_conc, dntp_conc)
                set2gibbs[ setp ][2].append( gibbs )
                if gibbs < -10:
                    delete_set.add( setp )
for setp in list(delete_set):
    set2primers.pop(setp)
    set2gibbs.pop(setp)
if len(set2primers) == 0:
    print("MultiHeterodimer stage: primers without dimers and hairpins does not exist")
    exit()
else:
    print( len(set2primers) )

delete_set = set()
for setp in set2primers:
    if len(set2primers[ setp ] ) == 2:
        break
    Gnum = set2primers[ setp ][1].count("G")
    Cnum = set2primers[ setp ][1].count("C")
    if set2primers[ setp ][1][0] == "G":
        delete_set.add( setp )
    if Cnum <= Gnum:
        delete_set.add( setp )
for setp in list(delete_set):
    set2primers.pop(setp)
    set2gibbs.pop(setp)

sorted_sets = []
for setp in set2gibbs:
    tmp = []
    for gibbs in set2gibbs[ setp ]:
        tmp.append( min(gibbs) )
    sorted_sets.append([ min(tmp), str(round(tmp[0], 2)), \
        str(round(tmp[1], 2)), str(round(tmp[2], 2)), setp ])
sorted_sets = sorted( sorted_sets, reverse=True, key=lambda gibbs: gibbs[0] )

file = open(argv[2] + ".out", "w")
file.write( header + "\t" + \
    "\t".join([ "Homodimer", "Heterodimer", "MultiHeterodimer" ]) + "\n" )
for setp in sorted_sets:
    for data in set2data[setp[4]]:
        file.write( data + "\t" + "\t".join(setp[1:4]) + "\n" )
file.close()
