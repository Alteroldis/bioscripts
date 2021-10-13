#!/usr/bin/python3
from math import isclose,fabs,log
from sys import argv

if len(argv) < 7:
    print("\n" + "##"*38 + "\n" + \
          "# USAGE:\n#  qpcr.py --cycles-file cycles.txt --samples condA,condB,condC\n" + \
          "#            --ref-genes Tubulin\n" + \
          "# Required:\n" + \
          "#  --cycles-file <string>  tabulatural file with cycles, genes and samples\n" + \
          "#  --samples     <string>  comma separated list of samples. First sample used\n" + \
          "#                as control\n" + \
          "#  --ref-genes   <string>  comma separated list of reference genes\n" + \
          "# Optional:\n" + \
          "#  --neg-contr   <string>  sample name of negative control\n" + \
          "#  --efficiency  <string>  comma separeted pairs genes:eff\n" + \
          "#                (ex. geneA:90,geneB:100). For unlisted genes efficiency is 100.\n" + \
          "#                For unlisted genes efficiency is 100%.\n" + \
          "#  --melt-temp   <string>  comma separeted pairs genes:melt-temperature\n" + \
          "#                (ex. geneA:80,geneB:66). In default mode product melting\n" + \
          "#                temperature estimated from data with range +/-1.5 Celsius\n" + \
          "#                degree. With manual mode for given genes estimated\n" + \
          "#                temperature replaced by the specified (with range +/-1.5)\n" + \
          "#                in command line\n" + \
          "\n" + "##"*38)
    exit(0)
    
file_cycles = argv[argv.index("--cycles-file") + 1]
houskeepings = argv[ argv.index("--ref-genes") + 1 ].split(",")
samples = argv[ argv.index("--samples") + 1 ].split(",")
if "--neg-contr" in argv:
    neg_contr = argv[ argv.index("--neg-contr") + 1 ]
    samples.append( neg_contr )
else:
    neg_contr = "WTF?!"
    

lib_gene2temp = {}
field = "no"
file = open(file_cycles, "r")
for line in file:
    if "Melt Temperature" in line:
        line = line.split("\n")[0]
        field = "yes"
        index_gene = line.split("\t").index( "Target" )
        index_sample = line.split("\t").index( "Sample" )
        index_cycle = line.split("\t").index( "Cq" )
        index_temp = line.split("\t").index( "Melt Temperature" )
        continue
    if field == "no": continue
    if line.split("\t")[index_sample] == neg_contr:
        continue
    gene = line.split("\t")[index_gene]
    try:
        temp = float( line.split("\t")[index_temp] )
    except: continue
    if gene not in lib_gene2temp:
        lib_gene2temp[ gene ] = [ temp ]
    else:
        lib_gene2temp[ gene ].append( temp )

for gene in lib_gene2temp:
    avg = sum( lib_gene2temp[gene] )/ len( lib_gene2temp[gene] )
    tmp = []
    for temp in lib_gene2temp[gene]:
        tmp.append([ fabs(temp - avg), temp] )
    tmp.sort()
    lib_gene2temp[ gene ] = tmp[0][1]
if "--melt-temp" in argv:
    given_temp = argv[ argv.index("--melt-temp") + 1 ].split(",")
    for gene2temp in given_temp:
        temp = float(gene2temp.split(":")[1])
        gene = gene2temp.split(":")[0]
        if isclose( temp, lib_gene2temp[ gene ], abs_tol=1 ):
            continue
        lib_gene2temp[ gene ] = temp

lib_gene2cycle = {}
lib_gene2eff = {}
for gene in lib_gene2temp:
    lib_gene2cycle[ gene ] = []
    for sample in samples:
        lib_gene2cycle[ gene ].append([])
        lib_gene2eff[ gene ] = 2
if "--efficiency" in argv:
    given_eff = argv[ argv.index("--efficiency") + 1 ].split(",")
    for gene2eff in given_eff:
        gene = gene2eff.slpit(":")[0]
        if gene in lib_gene2eff:
            eff = 2*float(gene2eff.split(":")[1])
            lib_gene2eff[ gene ] = eff

file.seek(0)
field = "no"
for line in file:
    if "Melt Temperature" in line:
        field = "yes"
        continue
    if field == "no": continue
    gene = line.split("\t")[index_gene]
    sample = line.split("\t")[index_sample]
    if sample not in samples: continue
    #if sample == neg_contr: continue
    try:
        temp = float( line.split("\t")[index_temp] )
        if isclose( temp, lib_gene2temp[ gene ], abs_tol=1 ):
            cycle = float( line.split("\t")[index_cycle] )
        else:
            cycle = 0
    except ValueError:
        cycle = 0
    index = samples.index( sample )
    lib_gene2cycle[ gene ][index].append( cycle )
file.close()

lib_gene2avg = {}
for gene in lib_gene2cycle:
    lib_gene2avg[ gene ] = []
    eff = lib_gene2eff[gene]
    for sample in samples:
        index = samples.index( sample )
        cycles = lib_gene2cycle[ gene ][index]
        if sample == neg_contr:
            #Если хоть одна из реплик с циклом=0, то всем проставляется 50
            #Это избавляет от случая, когда случайно капнул в одну из лунок контроля матрицу
            if 0 in cycles:
                for cycle in cycles:
                    cycles[ cycles.index(cycle) ] = 50
        for cycle in cycles:
            if cycle == 0:
                #Когда одна из реплик не вышла, ей назначается максимальный цикл среди реплик
                cycles[ cycles.index(0) ] = max(cycles)
        if 0 in cycles:
            #Если все реплики не вышли, то им присуждается цикл 50
            cycles = [50]
            print("WARNING: No detected specific product for %s:%s" % (gene, sample) )
        for index in range(0, len(cycles)):
            cycle = cycles[ index ]
            cycles[ index ] = log( eff**cycle, 2 )
        avg = sum(cycles)/len(cycles)
        lib_gene2avg[ gene ].append( avg )

if neg_contr in samples:
    samples = samples[:-1]
lib_gene2delta = {}
for gene in lib_gene2avg:
    if gene in houskeepings: continue
    lib_gene2delta[ gene ] = []
    for houskeeping in houskeepings:
        tmp = []
        for sample in samples:
            index = samples.index( sample )
            delta = lib_gene2avg[gene][index] - lib_gene2avg[houskeeping][index]
            tmp.append( delta )
            print(houskeeping, sample, round(lib_gene2avg[houskeeping][index], 2))
        lib_gene2delta[ gene ].append(tmp)

lib_gene2exprs = {}
for gene in lib_gene2delta:
    lib_gene2exprs[gene] = []
    for houskeeping in houskeepings:
        tmp = []
        index = houskeepings.index( houskeeping )
        refdelta = lib_gene2delta[gene][index][0]
        for delta in lib_gene2delta[gene][index][1:]:
            exprs = refdelta-delta
            tmp.append(exprs)
        lib_gene2exprs[ gene ].append(tmp)

fileout = open(file_cycles+".out", "w")
fileout.write("gene\tmelt-temp\n")
for gene in lib_gene2temp:
    lineout = "%s\t%s\n" % (gene, lib_gene2temp[gene])
    fileout.write( lineout )
fileout.write("\n\n")

fileout.write("Logfc values\n")
lineout = "gene\thouskeeping\t" + "\t".join(samples[1:]) + "\n"
fileout.write(lineout)
for gene in lib_gene2exprs:
    for houskeeping in houskeepings:
        lineout = gene + "\t" + houskeeping
        index = houskeepings.index( houskeeping )
        for exprs in lib_gene2exprs[gene][index]:
            lineout += "\t" + str(round(exprs,3))
        fileout.write(lineout+"\n")

if len(houskeeping)>1:
    fileout.write("\n\nGeometric mean or mean logfc values\n")
    lineout = "gene\t" + "\t".join(samples[1:]) + "\n"
    fileout.write(lineout)
    for gene in lib_gene2exprs:
        lineout = gene
        for sample in samples[1:]:
            index = samples[1:].index(sample)
            tmp = []
            for exprs in lib_gene2exprs[gene]:
                tmp.append(exprs[index])
            avg = sum(tmp)/len(tmp)
            lineout += "\t" + str(round(avg,3))
        fileout.write( lineout+"\n" )

if "--neg-contr" not in argv:
    fileout.close()
    exit()

fileout.write("\n\nNegative control\n")
fileout.write("gene\tavg-cycle\tneg-cycle\tFoldChange(neg/avg)\n")
for gene in lib_gene2avg:
    lineout = gene
    avgs = []
    for sample in samples:
        index = samples.index(sample)
        avgs.append( lib_gene2avg[gene][index] )
    avg = sum(avgs)/len(avgs)
    lineout += "\t" + str( round(avg,3) )
    neg = lib_gene2avg[gene][-1]
    lineout += "\t" + str( round(neg,3) )
    lineout += "\t" + str( round(neg/avg, 3) )
    fileout.write(lineout+"\n")
fileout.close()
