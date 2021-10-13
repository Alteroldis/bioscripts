#!/usr/bin/python3.6
from sys import argv
from Bio import SeqIO

if argv[1] == ("-h" or "-help" or "--help" or ""):
    print( "USAGE: %s full-blast-outfmt-results column-of-(qseqid, sseqid, bitscore) swiss(optional) gene-level(optional)" )
    print( "       Example: %s human.outfmt 0 1 2" )
    exit()

index_seqid = int(argv[2])
index_hit = int(argv[3])
index_score = int(argv[4])

file = open(argv[1], "r")
lib_full = {}
lib_seqid2hit = {}
lib_hit2seqid = {}
for line in file:
	seqid = line.split("\t")[index_seqid]
	if "swiss" in argv:
		if "GN=" not in line.split("\t")[index_hit]:
			hit = line.split("\t")[index_hit].split("|")[2].split("_")[0].upper()
		else:
			hit = line.split("\t")[index_hit].split("GN=")[1].split(" ")[0].upper()
	else:
		hit = line.split("\t")[index_hit]
	score = float( line.split("\t")[index_score] )
	if "gene-level" in argv:
		seqid = seqid.split("_i")[0]
		hit = hit.split("_i")[0]
	if seqid not in lib_seqid2hit:
		lib_seqid2hit[ seqid ] = [[ score, hit ]]
	else:
		lib_seqid2hit[ seqid ].append([ score, hit ])
	if hit not in lib_hit2seqid:
		lib_hit2seqid[ hit ] = [[ score, seqid ]]
	else:
		lib_hit2seqid[ hit ].append([ score, seqid ])
	if hit not in lib_full:
		lib_full[ hit ] = [[ score, seqid, line ]]
	else:
		lib_full[ hit ].append([ score, seqid, line ])
#фильтрация хитов по score
for seqid in lib_seqid2hit:
	lib_seqid2hit[ seqid ].sort( reverse=True )
	bestscore = lib_seqid2hit[ seqid ][0][0]
	tmp = []
	for hit in lib_seqid2hit[ seqid ]:
		if hit[1] in tmp: continue
		if hit[0] >= bestscore*0.8:
			tmp.append(hit[1])
	lib_seqid2hit[ seqid ] = tmp
#то же самое, но в обратную сторону
for hit in lib_hit2seqid:
	lib_hit2seqid[ hit ].sort( reverse=True )
	bestscore = lib_hit2seqid[ hit ][0][0]
	tmp = []
	for seqid in lib_hit2seqid[ hit ]:
		if seqid[1] in tmp: continue
		if seqid[0] >= bestscore*0.8:
			tmp.append(seqid[1])
	lib_hit2seqid[ hit ] = tmp
#фильтрация совсем уж плохих выравниваний в полном списке
for hit in lib_full:
	lib_full[ hit ].sort( reverse=True )
	bestscore = lib_full[ hit ][0][0]
	tmp = []
	for seqid in lib_full[ hit ]:
		if seqid[1] in tmp: continue
		if seqid[0] >= bestscore*0.7:
			tmp.append([seqid[1], seqid[2]])
	lib_full[ hit ] = tmp

def main(lib_hit2seqid, lib_seqid2hit, remove_hits, remove_seqids, lock_pair):
	for hit in lib_hit2seqid:
        #если нет seqid для хита - удаления хита
		#if len( lib_hit2seqid[hit] ) == 0:
		#	remove_hits.add(hit)
		#	continue
        #если для лучшего seqid, соответствующему hit, нет хитов - удаление seqid
		seqid = lib_hit2seqid[hit][0]
		if seqid not in lib_seqid2hit:
			continue
		#if len( lib_seqid2hit[seqid] ) == 0:
            #если такому seqid соответствует только этот hit, то удаление и hit
		#	if len( lib_hit2seqid[hit] ) == 1:
		#		lib_hit2seqid[hit] = []
		#		remove_hits.add(hit)
		#	remove_seqids.add( seqid )
		#	continue
        #сохранение пары seqid-hit, если она безусловно лучшая
		if hit == lib_seqid2hit[ seqid ][0]:
			#lib_hit2seqid[hit] = [ lib_hit2seqid[hit][0] ]
			lock_pair.add("%s:::%s" % (hit, seqid))
			remove_seqids.add( seqid )
			remove_hits.add(hit)

def clean(lib_hit2seqid, lib_seqid2hit, remove_hits, remove_seqids, lock_pair):
    while remove_seqids != set() or remove_hits != set():
        remove_keys = set()
        for seqid in lib_seqid2hit:
            for hit in list( remove_hits & set( lib_seqid2hit[seqid] )):
                if "%s:::%s" % (hit, seqid) in lock_pair: continue
                lib_seqid2hit[seqid].remove(hit)
                #try: lib_hit2seqid[hit].remove(seqid)
                #except ValueError: continue
            if lib_seqid2hit[seqid] == []:
                remove_keys.add( seqid )
                remove_seqids.add( seqid )
        for seqid in remove_keys:
            lib_seqid2hit.pop(seqid)
        remove_hits = set()
        remove_keys = set()
        for hit in lib_hit2seqid:
            for seqid in list( remove_seqids & set( lib_hit2seqid[hit] )):
                if "%s:::%s" % (hit, seqid) in lock_pair: continue
                lib_hit2seqid[hit].remove(seqid)
                #try: lib_seqid2hit[seqid].remove(hit)
                #except ValueError: continue
            if lib_hit2seqid[hit] == []:
                remove_keys.add( hit )
                remove_hits.add( hit )
        for hit in remove_keys:
            lib_hit2seqid.pop(hit)
        remove_seqids = set()

def clean1(lib_hit2seqid, lib_seqid2hit, remove_hits, remove_seqids, lock_pair):
    remove_keys = set()
    for seqid in lib_seqid2hit:
        for hit in list( remove_hits & set( lib_seqid2hit[seqid] )):
            if "%s:::%s" % (hit, seqid) in lock_pair: continue
            lib_seqid2hit[seqid].remove(hit)
                    #try: lib_hit2seqid[hit].remove(seqid)
                    #except ValueError: continue
        if lib_seqid2hit[seqid] == []:
            remove_keys.add( seqid )
            remove_seqids.add( seqid )
    for seqid in remove_keys:
        lib_seqid2hit.pop(seqid)
    remove_hits = set()
    remove_keys = set()
    for hit in lib_hit2seqid:
        for seqid in list( remove_seqids & set( lib_hit2seqid[hit] )):
            if "%s:::%s" % (hit, seqid) in lock_pair: continue
            lib_hit2seqid[hit].remove(seqid)
                    #try: lib_seqid2hit[seqid].remove(hit)
                    #except ValueError: continue
        if lib_hit2seqid[hit] == []:
            remove_keys.add( hit )
            remove_hits.add( hit )
    for hit in remove_keys:
        lib_hit2seqid.pop(hit)
    remove_seqids = set()

def count(lib_seqid2hit, lib_hit2seqid):
	count0_1 = [ len(lib_seqid2hit),len(lib_hit2seqid) ]; count1 = [0,0]
	for seqid in lib_seqid2hit:
		if len(lib_seqid2hit[seqid]) == 0:
			count0_1[0]-=1
		elif len(lib_seqid2hit[seqid]) == 1:
			count1[0]+=1
			count0_1[0]-=1
	for hit in lib_hit2seqid:
		if len(lib_hit2seqid[hit]) == 0:
			count0_1[1]-=1
		elif len(lib_hit2seqid[hit]) == 1:
			count1[1]+=1
			count0_1[1]-=1
	print(count0_1)
	print(count1)
	return(count0_1, count1)
	

print("Разрешение конфликтов")
lock_pair = set()
iteration = 0
count0_1 = [1,1]
count1 = [1,2]
while count0_1 != [0,0] or count1[0] != count1[1]:
    print("\niteration: %s" % iteration)
    remove_hits = set(); remove_seqids = set()
    main(lib_hit2seqid, lib_seqid2hit, remove_hits, remove_seqids, lock_pair)
    clean(lib_hit2seqid, lib_seqid2hit, remove_hits, remove_seqids, lock_pair)
    count0_1, count1 = count(lib_seqid2hit, lib_hit2seqid)
    iteration+=1
remove_hits = set(); remove_seqids = set()
main(lib_hit2seqid, lib_seqid2hit, remove_hits, remove_seqids, lock_pair)
clean(lib_hit2seqid, lib_seqid2hit, remove_hits, remove_seqids, lock_pair)
#print(lib_seqid2hit['Efra.gene904_i0_length1302.p1'])
#print(lib_hit2seqid["evm.model.scaffold23.15"])
#print(lib_hit2seqid["evm.model.scaffold2885.1"])

fileout = open("%s.reconciled" % argv[1], "w")
for hit in lib_hit2seqid:
    if len( lib_hit2seqid[hit] ) == 0: continue
    for seqid in lib_full[hit]:
        if lib_hit2seqid[hit][0] == seqid[0]:
            fileout.write(seqid[1])
            break
fileout.close()
