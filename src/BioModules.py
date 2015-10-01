#!/usr/env/bin python

"""
These scripts are for blank... needs to be filled in...

__author__ = "Andrew Garst"
__copyright__ = "Copyright 2015, The CREATE Project"
__credits__ = ["Andrew Garst", "Andrea Halweg-Edwards", "others"]
__license__ = "BSD"
__version__ = "0.1.0-dev"
__maintainer__ = "Andrew Garst"
__email__ = "ag.theseasquirt@gmail.com"
__status__ = "Development"
"""

def revcomp(s):
	import string
	complement = string.maketrans('ATCGNatcgnKRSBDMYWVHkm', 'TAGCNtagcnMYWVHKRSBDmk')
	return s.translate(complement)[::-1]

def reverse(s):
	return s[::-1]

def translate(s):
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    coding_dna = Seq(s, generic_dna) 

def digits_only(mystring):
    """
    digits_only(mystring):
    extract the numerical values from a mixed string
    
    example usage:
    digits_only('NonSynonymous K153T')
    
    output:
    153
    """
    import string
    all=string.maketrans('','')
    nodigs=all.translate(all, string.digits)
    return int(mystring.translate(all, nodigs))
def complement(s):
	"""
	Return complement of string
	"""
	import string
	complement = string.maketrans('ATCGNatcgn', 'TAGCNtagcn')
	return s.translate(complement)        

def genbank2fasta(infile, outfile): 
	""" infile and outfile are the full paths of the specified files"""
	from Bio import SeqIO
	infile=open(infile,'r')
	outfile=open(outfile,'w')
	count = SeqIO.convert(infile, "gb", outfile, "fasta")
	print "Converted %i records" % count
    
def fasta_upper(inputfile,outputfile):
    """
    Convert all records in a fasta file to uppercase
    """
    records = (rec.upper() for rec in SeqIO.parse(inputfile, "fasta"))
    SeqIO.write(records, outputfile, "fasta")
    
def fastq2fasta(inputfile,outputfile):
    """
    fastq2fasta(inputfile,outputfile). Convert fastq file into a fasta formatted file 
    using the Biopython SeqIO module"""
    from Bio import SeqIO
    from collections import Counter
    input_handle = open(inputfile,"rU")
    output_handle = open(outputfile, "w")
    sequences = SeqIO.parse(input_handle, "fastq")
    count = SeqIO.write(sequences, output_handle, "fasta")
    output_handle.close()
    input_handle.close()
    
def Biofasta2list(fastafile):
    """
    fasta2list(fastafile). Imports a fasta file with or without line wrapped sequences and 
    generates a python list of 'id', 'sequence' tuples for further manipulation. 
    """
    from Bio import SeqIO
    myseqs=SeqIO.parse(fastafile,"fasta")
    mylist=[]
    for record in myseqs:
    	mylist.append((record.id, record.seq))
    return mylist
    
def csv2list(infile,delim):
    """
    csv2list(infile,delim)
    example:
    csv2list(myfile,',')
    general function for importing tab delimited files to a python list. Will add as many rows as are in the file
    """
    import csv
    db=[]
    with open(infile, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            if row.startswith('#'):
                pass
            else:
                db.append(row)
    return db

def fasta_filecounter(fname):
    """
    Counts the number of sequences present in a fasta file
    """
    counts=0
    with open(fname) as f:
        for line in f:
            if line.startswith('>'):
                pass
            else:
                counts+=1
    print fname, counts

def csv_iter(infile):

    """
    csv_iter(infile)
    mport csv

    with open(infile, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            yield tuple(row)
    general function for importing tab delimited files to a python list. Will add as many rows as are in the file
    """
    import csv

    with open(infile, 'rU') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            yield tuple(row)

def fasta2list(infile):
    """
    fasta2list(infile)
    For entry in fasta file append name, sequence tuples to list and return the 
    list object for further information
    """
    def read_fasta(fp):
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))
    mylist=[]
    with open(infile,"r") as fp:
        for name, seq in read_fasta(fp):
            mylist.append((name, seq))
    return mylist

def tab2fasta(infile,namecol,seqcol,outfile):
    """
    tab2fasta(infile,outfile,namecol,seqcol):
    takes a tab delimited file with name and sequence columns as input and outputs a fasta file
    """
    import BioModules
    f=open(outfile,'w')
    for row in BioModules.csv_iter(infile):
        name,seq=row[namecol],row[seqcol]
        print >>f, '>'+name+'\n'+seq
    f.close()    
            
def fasta_iter(infile):
    """
    fasta_iter(infile)
    For entry in fasta file iterate over it and return the header and sequence 
    information. infile is the path of the file of interest
    """
    def read_fasta(fp):
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name.strip('>'), ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))
    mylist=[]
    
    with open(infile,"r") as fp:
        for name, seq in read_fasta(fp):
            yield name.strip('>'), seq
    
def fastq_iter(fastqfile):
    """
    fastq_iter(fastqfile)
    input: path to fastq file
    outputs ident,seq, qual from each entry iterated over in the file
    """
    f=open(fastqfile,'r')
    for line in f:
        ident,seq,space,qual=line.rstrip(),f.next().rstrip(), f.next().rstrip(), f.next().rstrip()
        yield ident, seq, qual

def fastq_gzip_iter(fastqfile):
    """
    fastq_gzip_iter(fastqfile)
    input: path to fastq.gz file
    outputs ident,seq, qual from each entry iterated over in the a .gz file
    """
    import gzip
    f=gzip.open(fastqfile,'r')
    for line in f:
        ident,seq,space,qual=line.rstrip(),f.next().rstrip(), f.next().rstrip(), f.next().rstrip()
        yield ident, seq, qual
def fastq_array(inputfile):
    """
    fastq_array(inputfile): Imports a fastq file with each tab seperated column as a 
    """
    f=open(inputfile)
    array=[]
    for line in f.readlines():
        array.append([x for x in line.strip("\n").split('\t')])
    return array

def lsdir(path,suffix):
    """
    lsdir(path,suffix):
    return a list of all files with suffix in the indicated path
    
    """
    import os
    files=[]
    for filename in os.listdir(path):
        if filename.endswith(suffix): 
            files.append(os.path.join(path, filename))
    return files

    
def N_seqfilter(mylist, Nallowed):
	"""
	def N_seqfilter(mylist):
    newlist=[]
    for id,seq in mylist:
        if seq.count('N') <=Nallowed:
            newlist.append((id, seq))# inner parentheses allow me to append both items back to tuples in the newlist
	return newlist
    """
	newlist=[]
	for id,seq in mylist:
		if seq.count('N') <=Nallowed:
			newlist.append((id, seq))# inner parentheses allow me to append both items back to tuples in the newlist
	return newlist
    
def fastafileprinter(list1,file):
    """
    Takes input list and prints to an open file container and writes a fasta file to it
    def fastafileprinter(list1,f):
    for id, seq in list1:
        print>> f, id,"\n"+seq#creates fasta file format output of all sequences for input into RNAfold server
    """
    f=open(file,'w')
    
    for id, seq in list1:
        print >> f, id,"\n"+seq#creates fasta file format output of all sequences for input into RNAfold server
    f.close()
    
    
def NCBI_genome2string(ncbi_id):
    """
    NCBI_genome2stringncbi_id). For example if I excecute the command 
    thermotogagenome=NCBI_genome2string("NC_000853.1") I will load the thermotoga
    maritima genome as an oject that I can slice as desired
    Useful genome sequences

    E. coli MG1655 'U00096.3'
    """
    from cogent.db.ncbi import EFetch
    from cogent.db.ncbi import ESearch
    ef = EFetch(id=ncbi_id)##read in E. coli genome 
    seq = ef.read().split('\n',1)[1].strip()
    genome = seq.replace('\n','')
    return genome

def usearch_counter(database,globalalignfile,cutoff):
    """
    usearch_counter(database,globalalignfile,cutoff):
    database=path to database of reference sequence
    globalalignfile= usearch results from global_alignment algorithm in blast6 format 
    (or userout with alignment included)
    cutoff=minimum percent identity to count from the matches in the globalalignfile
    will return a list of the sequence ids and counts for each
    """
    outsplit=usearchcountfile.split('/')
    outcolname=outsplit[-1].rstrip('.hits')
    
    readcounter=dict()
    outlist=[]
    
    
    d=open(database,'rU')
    for line in d:
        if line.startswith('#'):
            pass
        else:
            line=line.split('\t')
            name=str(line[0])
            readcounter[name]=0
    f=open(globalalignfile,'r')
    for line in f:
        if line.startswith('#'):
            pass
        else:
            line=line.rstrip().split('\t')
            ID,matchpercent=str(line[1]),line[2]
            if matchpercent>=cutoff:
                readcounter[ID]+=1
    
    idcountlist=sorted([(key,value) for key,value in readcounter.iteritems()])
    counts=[]
    counts.append(outcolname)
    for key, value in mylist:
        counts.append(value)
    return idcountlist

def usearch_hit_aligner(hitfile,sortcolumn,columns):
    """
    takes usearch hit files with userfields query+target+id+alnlen+mism+opens+trow+qrowdots
    and sorts by identity
    """
    import BioModules
    import os
    unique=set()
    mylist=[]
    #hitfile columns: query+target+id+alnlen+mism+opens+qrow+trow+evalue+bits+diffs
    for i in BioModules.tab_iter(hitfile):
        if i[0] not in unique:
            mylist.append(i)
            unique.add(i[0])
        else:
            pass
    #-userfields query+target+id+alnlen+mism+opens+qrow+trow+evalue+bit
    sortlist=sorted(mylist,key = lambda x : x[sortcolumn])#sort by % id
    widths = [max(map(len, col)) for col in zip(*sortlist)]
    for i in sortlist:
        #print "\t".join((val.ljust(width) for val, width in zip(i, widths)))
        print '\t'.join(['query','target','%id','aln_length','mismatches','gaps'])
        print '\t'.join(i[:5])
        print i[6]
        print i[8]
        print i[7]

def myBLAST(queryfile, dbfile, outfile, eval, fmtout, wordsize): 
    """
    myBLAST(queryfile, dbfile, outfile, eval, fmtout, wordsize): 
     	import os
     	from Bio.Blast.Applications import NcbiblastnCommandline
     	comline = NcbiblastnCommandline(query=queryfile, db=dbfile, strand="both",
                                  evalue=eval, out=outfile, outfmt=fmtout word_size=wordsize)
    	os.system(str(comline))

    """
    import os
    from Bio.Blast.Applications import NcbiblastnCommandline
    comline = NcbiblastnCommandline(query=queryfile, db=dbfile, strand="both",
                                  evalue=eval, out=outfile, outfmt=fmtout, word_size=wordsize)
    os.system(str(comline))

def PAMfinder_cas9(PAM,seq,n,d):
    """
    PAMfinder(PAM,seq,n,d). Input PAM of interest and target sequence and return tuples 
    consisting of the site/strand info and sequence with n nucletides upstream or d 
    nucleotides downstream of PAM sequence. Can accept degenerate
    nucleotides in the PAM using the Bio.SeqUtils.nt_search function
    """
    from BioModules import revcomp
    import re
    import string
    from Bio import SeqUtils
    seq=seq.upper()
    seqRev=revcomp(seq)
    positionsF=SeqUtils.nt_search(seq,PAM)
    positionsR=SeqUtils.nt_search(seqRev,PAM)
    sitelist=[]
    for item in positionsF[1:]:
        site=int(item)
        Fseq=seq[site-n:site+3+d]
        if seq!="":
        	yield str(item)+"Reverse",Fseq
    for item in positionsR[1:]:
        site=int(item)
        Fseq=seqRev[site-n:site+3+d]
        if seq!="":
        	yield str(item)+"Forward",Fseq

   
    
def hamdist(str1, str2):
	"""
	Counts the number of character differences in a string and returns this as an integer value
	"""
	from itertools import izip
	return float(sum(unicode(c1) != unicode(c2) for c1, c2 in izip(str1, str2)))


def exp_reads_generator(refseq,slicesize):
    """
    for input(refseq,slicesize)
    Returns dictionary of key:value pairs where the key is the start index and F or R for
    forward or reverse strand and the value is the slice of refseq of size indicated in slicesize
    exp_reads_generator(refseq,slicesize):
		mydict=dict()
		def revcomp(s):
			import string
			complement = string.maketrans('ATCGNatcgn', 'TAGCNtagcn')
			return s.translate(complement)[::-1]
		start=0
		while start<=len(refseq)-slicesize:
			seqslice=refseq[start:slicesize+start]
			mydict[str(start+1)+"\t"+ "F"]=seqslice
			mydict[str(start+1)+"\t"+"R"]=revcomp(seqslice)
			start+=1
		return mydict
    """
    mydict=dict()
    def revcomp(s):
        import string
        complement = string.maketrans('ATCGNatcgn', 'TAGCNtagcn')
        return s.translate(complement)[::-1]
    start=0
    while start<=len(refseq)-slicesize:
        seqslice=refseq[start:slicesize+start]
        mydict[str(start+1)+"\t"+ "F"]=seqslice
        mydict[str(start+1)+"\t"+"R"]=revcomp(seqslice)
        start+=1
    return mydict
    
def BLAST_string2file(inputstring, database,word, ecutoff, formattype, outfile):
    """
    BLAST_string(inputstring, db, outfile). Takes a string as input for a blast against the desired database and 
    outputs the results as a tabular BLAST to the outfile path.
    Example usage:
    BLAST_string("GGTTTACGCTTTACGTATAG","/Users/andrewgarst/blastdb_files/MG1655", "/Users/andrewgarst/Desktop/folAdonorhits.tab")
    """
    from Bio.Blast.Applications import NcbiblastnCommandline
    query = inputstring #your string from some external source
    blastn_cline = NcbiblastnCommandline(db=database,evalue=ecutoff, word_size=word, outfmt=formattype, out=outfile) #Blast command
    out = blastn_cline(stdin=query)
    
def BLAST_string_printer(inputstring, database,word, ecutoff, formattype):
    """
    BLAST_string_printer(inputstring, database, word, ecutoff). Takes a string and evalue word size criteria for a BLAST
    against the designated database and prints format type (5=xml, 7=tab)
    output on screen
    Example usage:
    BLAST_string("GGTTTACGCTTTACGTATAG","/Users/andrewgarst/blastdb_files/MG1655")
    """
    from Bio.Blast.Applications import NcbiblastnCommandline
    query = inputstring #your string from some external source
    blastn_cline = NcbiblastnCommandline(db=database,evalue=ecutoff, word_size=word, outfmt=formattype) #Blast command
    out, err = blastn_cline(stdin=query)
    print out

def NCBI_genomedownload(refseq_id,path2db):
    from Bio import Entrez
    from Bio import SeqIO 
    Entrez.email='ag.theseasquirt@gmail.com'
    import os
    genomehandle=path2db+refseq_id+'.gb'
    if os.path.isfile(genomehandle)==True:
        pass
    else:
        handle = Entrez.efetch(db="nucleotide", id=refseq_id ,rettype="gbwithparts",validate=True)
        genome = SeqIO.read(handle, 'genbank')
        SeqIO.write(genome,genomehandle,'genbank')
    
def BLAST_tabreader(filename):
	"""
	Takes tabular BLAST file (outfmt7 option) and reads to list composed of:
	column1:query
	column2: alnlength
	column3: qstart
	column4: mism
	column5: gap
	column6: qend
	"""
	f=open(filename,"r")
	myhitlist=[]
	for line in f.readlines():
		if "#" in line:
			pass
		else:
			line=line.split()
			query=line[0]
			target=line[1]
			percentid=float(line[2])
			alnlength=int(line[3])
			mism=int(line[4])
			gap=int(line[5])
			qstart=int(line[6])
			qend=int(line[7])
			tstart=int(line[8])
			tend=int(line[9])
			eval=float(line[10])
			bit=float(line[11])
			myhitlist.append((query,alnlength, qstart, mism, gap, qend))
	return myhitlist

def mistarget_calc(input_list):
    offtargetcounter=[]
    for query, alnlength, qstart, mism, gap, qend in input_list:
        if int(qend)>=22 and int(alnlength)>=12 and (int(mism)<=1 or int(gap)<=1):
            offtargetcounter.append(int(alnlength))
            break
        else:
            pass
    return sum(offtargetcounter)
    
def NCBI_genomedownload(ncbi_id, filepath, genomeid):
    """
    NCBI_genomedownload(ncbi_id, filepath, genomeid). Download a genome to the specified filepath as a fasta
    file with a header specified by genomeid. For example
   "NC_000913.3" will give me the E.coli genome 
    """
    from cogent.db.ncbi import EFetch
    from cogent.db.ncbi import ESearch
    ef = EFetch(id=ncbi_id)##read in E. coli genome 
    seq = ef.read().split('\n',1)[1].strip()
    genome = seq.replace('\n','')
    f=open(filepath,"w")
    f.write(">"+genomeid+"\n"+genome)
    f.close()
    
def dna2protein(DNA):
    DNA=DNA.upper()

    gencode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}           
    protein=''
    for i in range(0,len(DNA),3):
    	if i+3<=len(DNA):
    		codon=DNA[i]+DNA[i+1]+DNA[i+2]
    		protein=protein+gencode[codon]
    return protein

def list2csv(mylist,s,outfile):
    import csv
    readsfile=open(outfile,'w')
    csvwriter=csv.writer(readsfile,delimiter=s)
    for row in mylist:
        csvwriter.writerow(row)
    readsfile.close()


def list2fasta(list,filepath,namecol,seqcol):
	import codecs
	"""

	list2fastafile(list,filepath,namecol,seqcol): takes a list of (name, sequence) tuples and writes to a fasta file. Will assume that the name has '>' already included
    such as is the case from the fasta2list command
		f=open(filepath, "w")
		for x in list:
			print>> f, ">"+str(x[0]),"\n"+str(x[1]);#creates fasta file format output of all sequences for input into RNAfold server
		f.close()
	"""
	f=codecs.open(filepath, "w","utf-8")
	for x in list:
		print>> f,str(x[namecol]),"\n"+str(x[seqcol]);#creates fasta file format output of all sequences for input into RNAfold server
	f.close()
	
def gene_lookup(genename,db):
    """
    Usage: gene_lookup(genename,db):
    Given a whole genome fasta file as a database (db) lookup the coding sequence 
    for a gene (string input)
    """
    import BioModules
    for name, seq in BioModules.fasta_iter(db):
        if genename in name:
            return name, seq


def extended_gene_lookup(mylist,genedb,genomedb, distance):
    """
    Usage:
    extended_gene_lookup(mylist,genedb,genomedb, distance)
    given a fasta file of orfs (made with BioModules.genbank2fasta) and a fasta file
    with the whole genome sequence as one contiguous string (made with BioModules.NCBI_genome2string)
    lookup the gene by name and return it's sequence with additional sequence up and 
    downstream of the orf
	Default files to use for my system
	genedb='/Users/andrewgarst/blastdb_files/MG1655_genefeatures.fasta'
	genomedb='/Users/andrewgarst/blastdb_files/MG1655_genome.fasta'
	"""
    
    def gene_lookup(genename,genedb):
        """
        given a whole genome fasta file as a database (db) lookup the coding sequence for a given gene name
        """
        for name, seq in fasta_iter(genedb):
            if genename in name:
                return name, seq
    
    
    orflist=[]
    
    for gene in mylist:
        gene=gene_lookup(gene,genedb)
        orflist.append(gene)
    extendedlist=[]
    
    for name, seq in orflist:
        for ref, refseq in fasta_iter(genomedb):
            if seq in refseq:
                start=refseq.index(seq)
                stop=start+len(seq)
                orf=refseq[start:stop].lower()
                upstream=refseq[start-distance:start]
                downstream=refseq[stop:stop+distance]
                extendedlist.append((name,upstream+orf+downstream))
            else:
                seq=revcomp(seq)
                start=refseq.index(seq)
                stop=start+len(seq)
                orf=refseq[start:stop].lower()
                upstream=refseq[start-distance:start]
                downstream=refseq[stop:stop+distance]
                extendedlist.append((name,upstream+orf+downstream))
                
    
    return extendedlist
    
def Tm_calc(seq, oligo_conc, salt_conc, RNA):
	"""
	def Tm_calc(seq, oligo_conc, salt_conc, RNA):
		import Bio
    	import Bio.SeqUtils.MeltingTemp
    	Tmcalc=Bio.SeqUtils.MeltingTemp.Tm_staluc
    	print ('%0.2f' % Tmcalc(seq, oligo_conc, salt_conc, RNA))
    	
    Utilizes the tm_staluc module from Biopython to calculate hybridization energy
    
    oligo_conc is in nM, salt_conc is mM and RNA is 0 by default
    the following
    Tm_calc(seq, 1, 110, 0) will give results similar to NEB site for Taq
    Tm_calc(seq,2,600,0) will give results similar to NEB site for Phusion at 500 nM primer_conc
    
    for more info see:
    http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-pysrc.html
	"""
	import Bio
	import Bio.SeqUtils.MeltingTemp
	Tmcalc=Bio.SeqUtils.MeltingTemp.Tm_staluc
	freeenergy=float(Tmcalc(seq, oligo_conc, salt_conc, RNA))
	return round(freeenergy,2)


def RNA_folder(seq):
    """
    Takes a sequence input (as string) and outputs the dot bracket notation structure. Uses the cogent 
    application controller to execute RNAfold
    
    RNA_folder(seq):
        from cogent.app.vienna_package import RNAfold
        r = RNAfold()
        r.InputHandler ="_input_as_lines"
        seq=[seq]
        res=r(seq)
        ref_fold=res['StdOut'].read().split()
        return ref_fold[1]
    """
    from cogent.app.vienna_package import RNAfold
    r = RNAfold()
    r.InputHandler ="_input_as_lines"
    seq=[seq]
    res=r(seq)
    ref_fold=res['StdOut'].read().split()
    return ref_fold[1]

def RNAfold_hairpin(seq):
    """
    Takes a sequence input (as string) and outputs the free energy of the structure. Uses the cogent 
    application controller to execute RNAfold
    
    RNA_folder(seq):
        from cogent.app.vienna_package import RNAfold
        r = RNAfold()
        r.InputHandler ="_input_as_lines"
        seq=[seq]
        res=r(seq)
        ref_fold=res['StdOut'].read().split()
        return float(ref_fold[2].strip('()'))
    """
    from cogent.app.vienna_package import RNAfold
    r = RNAfold()
    r.InputHandler = "_input_as_lines"
    seq=[seq]
    res=r(seq)
    out=res['StdOut'].read().split()
    for i in out:
        try:
            energy=float(i.strip('()'))
            return energy
        except: 
            ValueError
            
def RNAfold_dist(seq1,seq2):
    """
    RNAfold_dist(seq1,seq2):
    Computes the distance in the minimal free energy structure of two sequences by 
    comparing their dot bracket structure and outputs a floating value of
    the percent similarity of the two dot bracket strings 
    Uses the RNA_folder and hamdist functions from BioModules
    from BioModules
    """
    
    def RNA_folder(seq):
        
        from cogent.app.vienna_package import RNAfold
        r = RNAfold()
        r.InputHandler ="_input_as_lines"
        seq=[seq]
        res=r(seq)
        ref_fold=res['StdOut'].read().split()
        return unicode(ref_fold[1])
   
    def hamdist(str1, str2):
        from itertools import izip
        assert len(str1) == len(str2)
        return float(sum(unicode(c1) != unicode(c2) for c1, c2 in izip(str1, str2)))
    
    fold1=RNA_folder(seq1)
    fold2=RNA_folder(seq2)
        
    dist=hamdist(fold1,fold2)
    percentdist=round((1.0-dist/len(seq1)),2)
    return percentdist
    

def GC_calc(seq):
	"""
	GC_calc(seq)
	Returns GC content as a float value with 2 decimal places.
	Reported as a fraction of 1 instead of percentage 
	Ex Input='ATCGGGAC'
	output=0.63
	"""
	seq=seq.upper()
	mydict=dict()
	GCtotal=int()
	Gtotal=int()
	Ctotal=int()
	for char in seq:
		if char in mydict:
			mydict[char]+=1
		else:
			mydict[char]=1
	for key, value in mydict.iteritems():
		if key =='G':
			Gtotal=value
		elif key=='C':
			Ctotal=value
	GCcontent=round(float(Gtotal+Ctotal)/len(seq),2)
	return GCcontent

def BLASTtab2list(infile):
    """

    """
    hits=open(infile,'r')
    myhits=hits.readlines()
    myhitlist=[]
    for line in myhits:
        if line.startswith("#"):
            pass
        else:
            line=line.split()
            query=line[0]
            target=line[1]
            percentid=float(line[2])
            alnlength=int(line[3])
            mism=int(line[4])
            gap=int(line[5])
            qstart=int(line[6])
            qend=int(line[7])
            tstart=int(line[8])
            tend=int(line[9])
            eval=float(line[10])
            bit=float(line[11])
            myhitlist.append((query, alnlength, qstart, mism, gap, qend))
    return myhitlist
    
def syn_codon_maxdiff(codon):
    """
    Finds the synonymous codon with the largest edit distance from the input codon. 
    """
    from BioModules import syn_codon_finder as syncodon
    from Levenshtein import distance as dist
    import operator
    codonlist=syncodon(codon)
    distlist=[]
    for i in codonlist:
        distlist.append(dist(i,codon))
    finallist=sorted(zip(codonlist,distlist),  key=operator.itemgetter(1), reverse=True)
    return finallist[0][0] 

    
def orthogonalPAMfinder(gene, PAM,spacerlen,genome):
    """
    Usage: orthogonalPAMfinder(gene, PAM,spacerlen,genome):
    for example (orthogonalPAMfinder(cas9, 'AGG', 20, '/Users/andrewgarst/blastdb_files/MG1655_genome')
    will BLAST a list of all possible AGG associated PAMS in the cas9 gene to a local MG1655 genome fasta file
    and return those that pass the stringency filters for match near the PAM.
    
    """
    import BioModules
    import os
    import Bio.SeqUtils
    genePAMS=BioModules.PAMfinder(PAM,gene,spacerlen,0)
    handle1='/tmp/foo.txt'
    handle2='/tmp/foo1.txt'
    
    #Write PAMS to temp fasta file
    BioModules.list2fastafile(genePAMS,handle1)
    #read PAMS from temporary fasta file
    BioModules.myBLAST(handle1,genome,handle2, 100, 7, 7)
    hitlist=BioModules.BLASTtab2list(handle2)
    
    #identify position of PAM in list entry
    def PAM_orientation(seq,PAM,posn):
        from Bio import SeqUtils
        PAMlocations=SeqUtils.nt_search(seq,PAM)
        if posn not in PAMlocations:
            return 'False'
        else:
            return 'True'
    #from all PAMS choose those that pass the filters in this identify_good_spacers   
    def identify_good_spacers(myhitlist,PAMS):
        from BioModules import revcomp
        badPAMS=set()
        goodPAMS=[]
        for site, seq in PAMS:
            revPAM=revcomp(PAM)
            check1=PAM_orientation(seq,PAM,20)
            check2=PAM_orientation(seq,revPAM,0)
            for query, alnlength, qstart, mism, gap, qend in myhitlist:
                if check1=='True'  and qend>=22 and alnlength>=14 and mism<=1 and gap<=1:
                    badPAMS.add(query)
                elif check2 and qstart<=0 and alnlength>=14 and mism<=1 and gap<=1:
                    badPAMS.add(query)
                else:
                    pass
               
        for site, seq in PAMS:
            if site not in badPAMS:
               goodPAMS.append((site, seq))
        return goodPAMS
    
    os.remove(handle1)
    os.remove(handle2)
    return identify_good_spacers(hitlist,genePAMS)


def string2fasta(seq,name,filepath):
    import codecs
    """
    list2fastafile(list,filepath): takes a list of (name, sequence) tuples and writes to a fasta file
        f=open(filepath, "w")
        for x in list:
            print>> f, ">"+str(x[0]),"\n"+str(x[1]);#creates fasta file format output of all sequences for input into RNAfold server
        f.close()
    """
    f=codecs.open(filepath, "w","utf-8")
    print>> f, ">"+str(name),"\n"+seq#creates fasta file format output of all sequences for input into RNAfold server
    f.close()

def oligo_check_Taq(seq):
	"""
	Usage: oligo_check(seq)
	input a string and it will return the Tm (assuming perfect complimentarity
	and the free energy (kcal/mol) of the mfe hairpin (computed using RNAfold
	so it overestimates stability a bit	
	"""
	import BioModules
	hybridenergy=BioModules.Tm_calc(seq,0.5,50,0)
	hairpin=BioModules.RNAfold_hairpin(seq)
	return hybridenergy, hairpin

def oligo_check_Phusion(seq):
	"""
	Usage: oligo_check(seq)
	input a string and it will return the Tm (assuming perfect complimentarity
	and the free energy (kcal/mol) of the mfe hairpin (computed using RNAfold
	so it overestimates stability a bit	
	"""
	import BioModules
	hybridenergy=BioModules.Tm_calc(seq,2,600,0)
	hairpin=BioModules.RNAfold_hairpin(seq)
	return hybridenergy, hairpin


def geneinforeader(genename,input):
    """
    geneinforeader(genename,input)
    input can be 'gene' for gene name, 'start' for orf start site
    'stop' for stop site, 'orientation' for orientation relative to 
    origin, or 'target_strand' for strand that you would want to target
    for recombineering
    """
    f=open('/Users/andrewgarst/databases/Recombineeringtargets.txt')
    for line in f.readlines():
        if genename in line:
            line=line.split('\t')
            gene=str(line[0])
            start=int(line[1])
            stop=int(line[2])
            orientation=str(line[3])#forward or reverse
            target_strand=str(line[4]).strip('\n')#sense or antisense
            if input=='gene':
                return gene
            elif input=='start':
                return start
            elif input=='stop':
                return stop
            elif input=='orientation':
                return orientation
            elif input=='target_strand':
                return target_strand
            else:
                pass
                

                    



def fasta2string(fastafile):
	f=open(fastafile)
	for line in f.readlines():
		seq=''
		if line.startswith('>'):
			pass
		else:
			line=line.strip('\n')
			seq=''.join(line)
			return seq
			
def gene2seq(gene,flank):
    import BioModules
    genedb='/Users/andrewgarst/blastdb_files/MG1655_genefeatures.fasta'
    genomedb='/Users/andrewgarst/blastdb_files/MG1655_genome.fasta'   
    infolist=open('/Users/andrewgarst/databases/Recombineeringtargets.txt','r')
    genomeseq=BioModules.fasta2string(genomedb)
    start=int(BioModules.geneinforeader(gene,'start'))
    stop=int(BioModules.geneinforeader(gene,'stop'))
    orientation=BioModules.geneinforeader(gene, 'orientation')
    upstream=genomeseq[start-flank:start].lower()
    orf=genomeseq[start-1:stop].upper()
    downstream=genomeseq[stop:stop+flank].lower()
    
    if orientation=='Forward':
        return downstream+orf+upstream
    else:
        return BioModules.revcomp(downstream+orf+upstream)

def spacer_printer_typeII(seq,PAM,start,stop):
    '''
    spacer_printer(PAM,seq,start,stop)
    
    Input is a PAM (i.e. 'NGG') and the sequence of interest. The start and stop inputs
    allow user choice the window in which to choose a PAM from
    '''
    from BioModules import revcomp
    from Bio import SeqUtils
    revPAM=revcomp(PAM)
    trimseq=seq[start:stop].upper()
    NGGlist=SeqUtils.nt_search(trimseq, PAM)[1:]
    CCNlist=SeqUtils.nt_search(trimseq,revPAM)[1:]
    for i in NGGlist:
        spacer=seq[i-20:i]
        PAM=seq[i:i+3].upper()
        yield str(i)+'\t'+spacer+'\t'+PAM+'\t'+'Reverse'
    for i in CCNlist:
        spacer=revcomp(seq[i+3:i+23])
        PAM=seq[i:i+3].upper()
        yield str(i)+'\t'+spacer+'\t'+PAM+'\t'+'Forward'


def spacer_printer_typeI(seq,PAM,scanstart,scanstop,offset, spacerlen):
    '''
    spacer_printer(seq,PAM,scanstart,scanstop,offset,spacerlen)
    
    Input is a PAM (i.e. 'NGG') and the sequence of interest. The start and stop inputs
    allow user choice the window in which to choose a PAM from
    '''
    from BioModules import revcomp
    from Bio import SeqUtils
    revPAM=revcomp(PAM)
    trimseq=seq[start:stop].upper()
    FPAMlist=SeqUtils.nt_search(trimseq, PAM)[1:]
    RPAMlist=SeqUtils.nt_search(trimseq,revPAM)[1:]
    myspacers=[]
    for i in FPAMlist:
        spacerF=seq[i+offset:i+offset+spacerlen]
        myPAMF=seq[i:i+len(PAM)].upper()
        yield str(i)+'\t'+myPAMF+'\t'+spacerF+'\t'+'Reverse'
    for i in  RPAMlist:
        spacerF=seq[i+offset:i+offset+spacerlen]
        myPAMF=seq[i:i+len(PAM)].upper()
        yield str(i)+'\t'+myPAMF+'\t'+spacerF+'\t'+'Forward'

def CRISPR_target_finder(seq,PAM,systemtype):
    '''
    CRISPR_target_finder(seq,PAM,systemtype)
    seq is the target sequence, PAM is the PAM to search (can use degenerate IUPAC notation)
    systemtype is either 1 for type I or 2 for type II
    returns a list of name, sequence tuples
    '''
    
    from BioModules import revcomp
    from Bio import SeqUtils
    targetlist=[]
    seq=seq.upper()
    revPAM=revcomp(PAM)
    Flist=SeqUtils.nt_search(seq, PAM)[1:]
    Rlist=SeqUtils.nt_search(seq,revPAM)[1:]
    if systemtype==1:
        for i in Flist:
            target=seq[i:i+35]
            targetlist.append((str(i)+'Reverse',target))
        for i in Rlist:
            target=revcomp(seq[i-35:i])
            targetlist.append((str(i)+'Forward',target))

    if systemtype==2:
        for i in Flist:
            target=seq[i-20:i+3]
            targetlist.append((str(i)+'Reverse',target))
        for i in Rlist:
            target=revcomp(seq[i:i+23])
            targetlist.append((str(i)+'Forward',target))
    return targetlist 
    
def framefinder(orfstart,site):
	"""
	framefinder(orfstart,site). input the start site of an ORF and
	the another short sequence (i.e. a sequencing read or oligo) and output
	which frame the use for other analysis
	for example if orfstart=100 and site=110 then the output should be
	1 meaning that to read the shorter sequence
	"""
	def framecheck(num):
		return num%3==0
	frame=int()

	if framecheck(site-orfstart)==True:
		frame=0
	elif framecheck(site+1-orfstart)==True:
		frame=1
	elif framecheck(site+2-orfstart)==True:  
		frame=2
	return frame

def first_spacer_finder(seq,PAM,start,stop):
	'''
	first_spacer_finder(seq,PAM,start,stop):

	Input is a PAM (i.e. 'NGG') and the sequence of interest. The start and stop inputs
	allow user to specify the range of the sequence in which to search for the
	PAM. Unlike spacer_printer function, this one only returns the first in the 
	sequence
	'''
	import BioModules
	from Bio import SeqUtils
	revPAM=BioModules.revcomp(PAM)
	trimseq=seq[start:stop].upper()
	NGGlist=SeqUtils.nt_search(trimseq, PAM)[1:]
	CCNlist=SeqUtils.nt_search(trimseq,revPAM)[1:]
	if firstNGG<firstCCN:
		spacer=seq[firstNGG-20:firstNGG]
		PAM=seq[firstNGG:firstNGG+3].upper()
		return firstNGG, PAM, spacer
	else:
		spacer=seq[firstCCN+3:firstCCN+23]
		PAM=seq[firstCCN:firstCCN+3].upper()
		return firstCCN, PAM, spacer

def myBLAST_NCBI(seq,db, prog,evalue,hits):
    """
     myBLAST_NCBI(seq,db, prog,evalue,hits)
     
     example ussage:
     
     target='atgccaactatccagcagctaattcgtagcgaacgctcgaaggtacagaagaaaactaaatcc'
     myBLAST_NCBI(target,'nt', 'blastn',0.001,10)
    """
    
    from Bio.Blast import NCBIWWW
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    seq=seq.upper()
    result_handle = NCBIWWW.qblast(prog, db ,seq,alignments=hits,expect=evalue, format_type='Text',hitlist_size=hits)
    print result_handle.read()
    result_handle.close()

    

def cas9_spacerfinder(seq,PAM,scanstart,scanstop,offset, spacerlen,outfile):
    """
    spacer_printer(seq,PAM,scanstart,scanstop,offset,spacerlen)
    
    Input is a PAM (i.e. 'NGG') and the sequence of interest. The start and stop inputs
    allow user choice the window in which to choose a PAM from
    """
    def spacer_printer(seq,PAM,scanstart,scanstop,offset, spacerlen):
        gRNAhandle='GTTTTAGAGCTAGAAATAGCAAGTTAAAAT'
        optimalfold='....................(((((((.((((....))))...)))))))'
        d='\t'
        def hamdist(str1, str2):
            from itertools import izip
            return float(sum(unicode(c1) != unicode(c2) for c1, c2 in izip(str1, str2)))
        
        import BioModules
        from Bio import SeqUtils
        revPAM=BioModules.revcomp(PAM)
        trimseq=seq[scanstart:scanstop]
        FPAMlist=SeqUtils.nt_search(trimseq, PAM)[1:]
        RPAMlist=SeqUtils.nt_search(trimseq, revPAM)[1:]
        try:
            for i in  FPAMlist:
                spacerF=(seq[i-spacerlen:i])
                gRNAF=spacerF+gRNAhandle
                foldF=BioModules.RNA_folder(spacerF+gRNAhandle)
                distF=1-(BioModules.hamdist(foldF,optimalfold)/len(optimalfold))
                distF=str(round(distF,2))
                gcF=str(BioModules.GC_calc(spacerF))
                myPAMF=seq[i:i+3].upper()
                yield str(i)+d+myPAMF+d+spacerF+d+'Reverse'+d+gRNAF+d+foldF+d+distF+d+gcF
            for i in  RPAMlist:
                spacerR=BioModules.revcomp(seq[i+offset:i+offset+spacerlen])
                gRNAR=spacerR+gRNAhandle
                foldR=BioModules.RNA_folder(spacerR+gRNAhandle)
                distR=1-(BioModules.hamdist(foldR,optimalfold)/len(optimalfold))
                distR=str(round(distR,2))
                gcR=str(BioModules.GC_calc(spacerR))
                myPAMR=seq[i:i+3].upper()
                yield str(i)+d+myPAMR+d+spacerR+d+'Forward'+d+gRNAR+d+foldR+d+distR+d+gcR
        except:
            TypeError
    f=open(outfile,'w')
    for i in spacer_printer(seq,PAM,scanstart,scanstop,offset, spacerlen):
        print>>f, i
    f.close()

def PAMfinder_typeI(PAM,seq,n,d):
    """
    PAMfinder(PAM,seq,n,d). Input PAM of interest and target sequence and return tuples 
    consisting of the site/strand info and sequence with n nucletides upstream or d 
    nucleotides downstream of PAM sequence. Can accept degenerate
    nucleotides in the PAM using the Bio.SeqUtils.nt_search function
    """
    from BioModules import revcomp
    import re
    import string
    from Bio import SeqUtils
    seq=seq.upper()
    seqRev=revcomp(seq)
    positionsF=SeqUtils.nt_search(seq,PAM)
    positionsR=SeqUtils.nt_search(seqRev,PAM)
    sitelist=[]
    for item in positionsF[1:]:
        site=int(item)
        Fseq=seq[site-n:site+3+d]
        if seq!="":
        	yield str(item)+'\t'+"Reverse"+'\t'+Fseq
    for item in positionsR[1:]:
        site=int(item)
        Fseq=seqRev[site-n:site+3+d]
        if seq!="":
        	yield str(item)+'\t'+"Forward"+'\t'+Fseq

def tab_iter(infile,sep,linestart):
    """
    tab_iter(infile):
    f=open(infile,'r')
    for line in f:
        line=line.split('\t')
        yield line[i]
        
    """
    f=open(infile,'rU')
    for i,line in enumerate(f):
        if i>=linestart:
            line=line.strip().rstrip().split(sep)
            yield line


def CPECdesigner(insert,backbone,tm):
    """
    CPECdesigner(insert, backbone, tm)

    insert and backbone are the sequences you want to stich together, tm is the target
    tm of the overlap between these two

    The program will out put a list primer names and sequences to construct the final vector as well
    as the final vector sequence

    This should also work for designing assemblies for gibson or homologous recombination based 
    construction approaches
    """
    import BioModules
   
    fbbprimer=''
    rbbprimer=''
    iFoverlap=''
    iRvoerlap=''
    iF=''
    iR=''
    outputlist=[]
    count=1
    while count<50:
        if BioModules.Tm_calc(backbone[:count],1, 110, 0)<tm:
            pass
        else:
            fbbprimer=backbone[:count]
            break
        count+=1
    
    count=1
    while count<50:
        if BioModules.Tm_calc(backbone[-count:],1, 110, 0)<tm:
            pass
        else:
            rbbprimer=BioModules.revcomp(backbone[-count:])
            break
        count+=1  
    
    count=1
    while count<50:
        if BioModules.Tm_calc(insert[:count],1, 110, 0)<60:
            pass
        else:
            iFoverlap=insert[:count]
            break
        count+=1        
    count=1
    while count<50:
        if BioModules.Tm_calc(insert[-count:],1, 110, 0)<60:
            pass
        else:
            iRoverlap=BioModules.revcomp(insert[-count:])
            break
        count+=1  
    outputlist=[]
    outputlist.append(("bbF",fbbprimer))
    outputlist.append(("bbR",rbbprimer))    
    outputlist.append(("insertF",BioModules.revcomp(fbbprimer)+iRoverlap))
    outputlist.append(("insertR",BioModules.revcomp(rbbprimer)+iFoverlap))      
    outputlist.append(("final_vector",insert+backbone))
    return outputlist

def CRISPR_GGassembler(repeat,offset,step, *spacers):
    """
    CRISPR_GGassembler(repeat,offset,step, *spacers)
    repeat= string of repeat sequence
    offset=position to start splicing of repeat sequences
    step= stepsize fore each move through the repeat sequence
    *spacers=list of the spacer sequences of interest
    """
    import BioModules
    import re
    Bsa_F='ggtctca'
    x='AGCC'
    Bsa_R=BioModules.revcomp(Bsa_F)
    pieces=[]
    array=''
    spacerlen=len(spacers[0])
    overlap=''
    
    def oligo_overlap(seq,step):
        while step<len(seq):
            if BioModules.Tm_calc(seq[len(seq)/2-step:len(seq)/2+step],1, 110, 0)<60:
                pass
            else:
                overlap=seq[len(seq)/2-step:len(seq)/2+step]
                start=re.search(overlap,seq).start()
                return seq[:start+len(overlap)], BioModules.revcomp(seq[start:])
                break
            step+=1
    pieces.append(('Name', 'Seq','sticky end'))
    anchor2=-len(repeat)+offset+step*len(spacers)
    pieces.append(('FgenPrimer',BioModules.revcomp(repeat[:offset+4]+Bsa_R+x),repeat[offset:offset+4]))
    count=1
    for name, seq in spacers:
        array+=repeat+seq.lower()
        l=repeat[offset:]
        r=repeat[:-len(l)+4+step]
        piece=x+Bsa_F+l+seq.lower()+r+Bsa_R+x
        oligos=oligo_overlap(piece,1)
        pieces.append((name+str(count)+'F',oligos[0],l[:4]))
        pieces.append((name+str(count)+'R',oligos[1],r[-4:]))
       
        
        offset+=step
        count+=1
    array+=repeat
    
    pieces.append(('RgenPrimer',x+Bsa_F+repeat[anchor2:], repeat[anchor2:anchor2+4]))
    pieces.append(('Target Array','',''))
    pieces.append((array,'',''))
    
    return pieces

def simple_NNKoligo_generator(seq,genename,upstream,downstream,offset,*posns):
    """
    NNKoligo_generator(seq,genename,upstream,downstream,*posns):
    seq is the target sequence that is to be mutated
    
    genename is used to annotate the oligo
    
    upstream and downstream indicate the number of nucleotides to be taken from either side of the
    NNK mutation
    
    *posns is the list of amino acid positions(this can take an input list with *listname as the arguments)
    
        ex. 
        folA='TTCGAAGGTGCGCTGAAAACCGGGCGTCTGGCACTGGAAAGTTTAGGTCTGGGGCCGTATGAAGCGCGAGAACGTGCCGATGTGTTCCGCCGCTTTAATATTCAGATGGTGGAAGAGATGGCAATGGTTGAGAACGACACCAAAGCCCGCGCGGCGGTCTATAAACGCACCAGCGCGATGTTAAGTGAGATCATTACCGAGGACCGCGAACATCTGTCATTAATTCAACGACATGGCTGGCAGGGAACCGAAGAAGGTAAACATACCGGCAACATGGCGGATGAACCGGAAACGAAACCCTCATCCTAATAAAGAGTGACGTAAATCACACTTTACAGCTAACTGTTTGTTTTTGTTTCATTGTAATGCGGCGAGTCCAGGGAGAGAGCGTGGACTCGCCAGCAGAATATAAAATTTTCCTCAACATCATCCTCGCACCAGTCGACGACggtttacgctttacgtatagTGGCGACAATTTTTTTTATCGGGAAATCTCAatgATCAGTCTGATTGCGGCGTTAGCGGTAGATCGCGTTATCGGCATGGAAAACGCCATGCCGTGGAACCTGCCTGCCGATCTCGCCTGGTTTAAACGCAACACCTTAAATAAACCCGTGATTATGGGCCGCCATACCTGGGAATCAATCGGTCGTCCGTTGCCAGGACGCAAAAATATTATCCTCAGCAGTCAACCGGGTACGGACGATCGCGTAACGTGGGTGAAGTCGGTGGATGAAGCCATCGCGGCGTGTGGTGACGTACCAGAAATCATGGTGATTGGCGGCGGTCGCGTTTATGAACAGTTCTTGCCAAAAGCGCAAAAACTGTATCTGACGCATATCGACGCAGAAGTGGAAGGCGACACCCATTTCCCGGATTACGAGCCGGATGACTGGGAATCGGTATTCAGCGAATTCCACGATGCTGATGCGCAGAACTCTCACAGCTATTGCTTTGAGATTCTGGAGCGGCGGTAATTTTGTATAGAATTTACGGCTAGCGCCGGATGCGACGCCGGTCGCGTCTTATCCGGCCTTCCTATATCAGGCTGTGTTTAAGACGCCGCCGCTTCGCCCAAATCCTTATGCCGGTTCGACGGCTGGACAAAATACTGTTTATCTTCCCAGCGCAGGCAGGTTAATGTACCACCCCAGCAGCAGCCGGTATCCAGCGCGTATATACCTTCCGGCGTACCTTTGCCCTCCAGCGATGCCCAGTGACCAAAGGCGATGCTGTATTCTTCAGCGACAGGGCCAGGAATCGCAAACCACGGTTTCAGTGGGGCAGGGGCCTCTTCCGGCGATTCTTTGCTGTACATATCCAGTTGACCGTTCGGGAAGCAAAAACGCATACGGGTAAAAGCGTTGGTGATAAAACGCAGTCTTCCCAGCCCCCGCAATTCCGGTGACCAGTTATTTGGCATATCGCCGTACATGGCATCAAGAAAGAAGGGATAGGAGTCACTCGATAGCACCGC'

        folAsites=[10,20,21,26,30,45,94,115,153,158]
        folAoffset=40

        for name, oligo, site,wtcodon in BioModules.simple_NNKoligo_generator(folA.upper(),'folA',40,20,500,*folAsites):
            print name, oligo,wtcodon

    gives: 

        folA_10V TCGGGAAATCTCAATGATCAGTCTGATTGCGGCGTTAGCGNNKGATCGCGTTATCGGCAT GTA
        folA_20M GGCGTTAGCGGTAGATCGCGTTATCGGCATGGAAAACGCCNNKCCGTGGAACCTGCCTGC ATG
        folA_21P GTTAGCGGTAGATCGCGTTATCGGCATGGAAAACGCCATGNNKTGGAACCTGCCTGCCGA CCG
        folA_26A CGTTATCGGCATGGAAAACGCCATGCCGTGGAACCTGCCTNNKGATCTCGCCTGGTTTAA GCC
        folA_30W GGAAAACGCCATGCCGTGGAACCTGCCTGCCGATCTCGCCNNKTTTAAACGCAACACCTT TGG
        folA_45H TAAACGCAACACCTTAAATAAACCCGTGATTATGGGCCGCNNKACCTGGGAATCAATCGG CAT
        folA_94I AGCCATCGCGGCGTGTGGTGACGTACCAGAAATCATGGTGNNKGGCGGCGGTCGCGTTTA ATT
        folA_115I ACAGTTCTTGCCAAAAGCGCAAAAACTGTATCTGACGCATNNKGACGCAGAAGTGGAAGG ATC
        folA_153F ATTCCACGATGCTGATGCGCAGAACTCTCACAGCTATTGCNNKGAGATTCTGGAGCGGCG TTT
        folA_158R TGCGCAGAACTCTCACAGCTATTGCTTTGAGATTCTGGAGNNKCGGTAATTTTGTATAGA CGG
    """
    import BioModules

    mutlist=[]
    ntlist=[i*3-3+offset for i in posns]
 
    for i in ntlist:
        mutlist.append((genename+'_'+str((i-offset)/3+1)+BioModules.dna2protein(seq[i:i+3]),seq[i-upstream:i]+'NNK'+seq[i+3:i+downstream], (i-offset)/3+1, seq[i:i+3]))
    return mutlist


def Gibsonassembler(tm,minprime,minoverlap,maxprime,constructname, pieces):
    """
    This is for design of primers for gibson assembly or CPEC like constructions
    Gibsonassembler(tm,minlen1,minoverlap,maxlen,constructname, pieces)
    inputs are the 
    1) tm desired for each overlap (calculated with BioModules.Tm_calc(seq[:count],1, 110, 0))
    2) the minimal length of the priming arm for each piece
    3) the maximal overlap lengths
    4) the pieces needed for construction. This is a list of name, seq tuples(i.e. [(name1,seq1), (name2,seq2), etc]
    
    
    this will output a list of primers for each piece in the order they're entered in the list and a sequence of the final vector
    
    To print oligos as a tab delimited list for ordering call this function as so
    
    for name, seq, anneal,tm, full_fold in Gibsonassembler(60,18,20,50, 'NADPH_low_red',partslist3):
        print name, seq, anneal,tm,full_fold
    for the example this would output a list like so
    [('pRR1bb_noGFPV2F', 'ttcccccgacgGGCGGGTGTCGGGGCGCAG', 72.72, 86.14, ['...((((((((......)))))))).....', -15.1])
    ('pRR1bb_noGFPV2R', 'agtgtagatcgcTGAATTCCAACTGAGCGCCGGTCGC', 72.72, 80.1, ['(((((((((.....)))))...)).))..........', -6.3])
    ('RFP_BBa_E1010V2F', 'agttggaattcaGCGATCTACACTAGCACTATCAG', 60.81, 71.42, ['.((((.....)))).....................', -2.3])
    ('RFP_BBa_E1010V2R', 'agaggcagatttATGGCTTCCTCCGAAGACG', 61.13, 73.09, ['.......(((((.((........)).)))))', -5.4])
    ('soxRmoduleV2F', 'gaggaagccatAAATCTGCCTCTTTTCAGTGTTC', 60.13, 71.42, ['((((.((........)).))))............', -4.3])
    ('soxRmoduleV2R', 'cgacacccgccCGTCGGGGGAAACCCTCCT', 65.27, 82.52, ['.((.(((...)))))((((......)))).', -7.5])]
     as well as the full construct at the end of the list

    """
    import BioModules
    junction_sites=[]
    primers=[]
    endseq=''
    junctions=[]
    junctions.append(maxprime)
    #make the full length construct and make a pseudo circular model for primer design
    for name,seq in pieces:
        endseq=endseq+seq.upper()
        junctions.append(maxprime+len(endseq))
    endseq=pieces[-1][1][-maxprime:]+endseq+pieces[0][1][:maxprime]
    #read in each part and map to full construct to make primers
    starts=[]
    counter1=0
    for name, seq in pieces:
        partstart=junctions[counter1]
        partend=junctions[counter1+1]
        Fstart=0
        Rstart=0
        #Fend=0
        #Rend=0
        
        #generate F priming site indices in endseq for each piece
        count=1
        while count<=maxprime:
            tmF=BioModules.Tm_calc(endseq[partstart:partstart+count],1,110,0)
            if count>minprime and count<maxprime and tmF>tm:
                Fend=partstart+count
                break
            else:
                Fend=partstart+count
            count+=1 
        
        #generate R priming site indices in endseq for each piece
        count=1
        while count<=maxprime:
            tmR=BioModules.Tm_calc(endseq[partend-count:partend],1,110,0)
            if count>minprime and count<maxprime and tmR>=tmF-5:
                Rstart=partend-count
                break
            else:
                Rstart=partend-count
            count+=1 
            
        #generate F overlap site using indices in endseq for each piece
        count=1
        while count<=maxprime:
            if count>minoverlap/2 and count<maxprime and BioModules.Tm_calc(endseq[partstart-count:partstart+count],1,110,0)>tm:
                #overlaps.append(endseq[i-count:i+count])
                Fstart=partstart-count
                break
            else:
                Fstart=partstart-count
            count+=1 
        #generate R overlap site using indices in endseq for each piece
        while count<=maxprime:
            if count>minoverlap/2 and count<maxprime and BioModules.Tm_calc(endseq[partend-count:partend+count],1,110,0)>tm:
                #overlaps.append(endseq[i-count:i+count])
                Rend=partend+count
                break
            else:
                Rend=partend+count
            count+=1
        counter1+=1
        #using numbers calculated from above, make primers and annealing sites
        Fprimer=endseq[Fstart:partstart].lower()+endseq[partstart:Fend]
        Fpartial=endseq[partstart:Fend]
        Rprimer=endseq[Rstart:partend]+endseq[partend:Rend].lower()
        Rpartial=endseq[Rstart:partend]
        #calculate tm for partial sites (inital primer annealing) and full length primers
        tmFp=BioModules.Tm_calc(Fprimer,1,110,0)
        tm_fpartial=BioModules.Tm_calc(Fpartial,1,110,0)
        tmRp=BioModules.Tm_calc(Rprimer,1,110,0)
        tm_rpartial=BioModules.Tm_calc(Rpartial,1,110,0)
        f_fold=[BioModules.RNA_folder(Fprimer), BioModules.RNAfold_hairpin(Fprimer)]
        r_fold=[BioModules.RNA_folder(Rprimer), BioModules.RNAfold_hairpin(Rprimer)]
        #append these to a list for each part
    
        primers.append((name+'F',Fprimer,tm_fpartial,tmFp, f_fold))  
        primers.append((name+'R',BioModules.revcomp(Rprimer),tm_rpartial,tmRp, r_fold))  
    if len(primers)< 2*len(pieces):
        return "Error: Tm to high or maxprimerlength to low"
    else:
        primers.append((constructname,endseq[maxprime:len(endseq)-maxprime],'','',''))
        return primers

def degenerate_re_convert(seq):
        
    """
    degenerate_re_convert(seq)
    make degenerate DNA strings to use with the re module
    conversiondict={'N':'(A|T|G|C)','R':'(A|G)','A':'A','C':'C','G':'G','T':'T'}
    
    example usage:
    
    degenerate_re_convert('NRG')
    
    output:
    '(A|T|G|C)(A|G)G'
    
    """
    conversiondict={'N':'(A|T|G|C)','S':'(C|G)','W':'(A|T)','R':'(A|G)','A':'A','C':'C','G':'G','T':'T'}
    newlist=list(seq)
    newseq=list()
    for i in newlist:
        for key, value in conversiondict.iteritems():
            if i==key:
                newseq.append(value)
    return ''.join(newseq)

def PAM_spacer_re_finder(seq,PAM,systemtype):
    '''
    PAM_spacer_re_finder(seq,PAM,systemtype)
    
    example usage:
    acrBtest='cattgatacaacgtgtaatcactaaggccgcgtaagcggccttttttatgcataacctacgaacattaaggagtaatt'
    PAM_spacer_re_finder(acrBtest,'NRG',2)
    [('23Reverse', 'TGATACAACGTGTAATCACTAAG'),
     ('33Reverse', 'TGTAATCACTAAGGCCGCGTAAG'),
     ('36Reverse', 'AATCACTAAGGCCGCGTAAGCGG'),
     ('67Reverse', 'ATGCATAACCTACGAACATTAAG'),
     ('70Reverse', 'CATAACCTACGAACATTAAGGAG'),
     ('0Forward', 'AGTGATTACACGTTGTATCAATG'),
     ('8Forward', 'GCGGCCTTAGTGATTACACGTTG'),
     ('11Forward', 'TACGCGGCCTTAGTGATTACACG'),
     ('19Forward', 'AGGCCGCTTACGCGGCCTTAGTG'),
     ('28Forward', 'GCATAAAAAAGGCCGCTTACGCG'),
     ('36Forward', 'TAGGTTATGCATAAAAAAGGCCG'),
     ('50Forward', 'CTCCTTAATGTTCGTAGGTTATG')]
    '''
    
    from BioModules import revcomp
    from Bio import SeqUtils
    import re
    seq=seq.upper()
    def degenerate_re_convert(seq):
        conversiondict={'N':'(A|T|G|C)','S':'(C|G)','W':'(A|T)','R':'(A|G)','A':'A','C':'C','G':'G','T':'T'}
        newlist=list(seq)
        newseq=list()
        for i in newlist:
            for key, value in conversiondict.iteritems():
                if i==key:
                    newseq.append(value)
        return ''.join(newseq)
    
    targetlist=[]
    revPAM=revcomp(PAM)
    PAM=degenerate_re_convert(PAM.upper())
    revPAM=degenerate_re_convert(revPAM.upper())
    Flist=[match.start() for match in re.finditer(PAM,seq)]
    Rlist=[match.start() for match in re.finditer(revPAM,seq)]
    if systemtype==1:
        for i in Flist:
            P=seq[i:i+3]
            target=seq[i+3:i+35]
            if len(target)==32:
                targetlist.append((str(i)+'_Reverse',target,P))
            else:
                pass
        for i in Rlist:
            target=revcomp(seq[i-32:i])
            P=revcomp(seq[i-3:i])
            targetlist.append((str(i)+'_Forward',target,P))
            if len(target)==32:
                targetlist.append((str(i)+'_Forward',target))
            else:
                pass
    if systemtype==2:
        for i in Flist:
            P=seq[i:i+3]
            target=seq[i-20:i]
            if len(target)<20:
                pass
            else:
                targetlist.append((str(i)+'_Reverse',target,P))
        for i in Rlist:
            P=revcomp(seq[i-3:i])
            target=revcomp(seq[i:i+20])
            if len(target)<20:
                pass
            else:
                targetlist.append((str(i)+'_Forward',target,P))
    return targetlist

def orthogonal_strings(genomeID,querys,degenerateposns,seedlen):
    """
    orthogonal_spacercheck2(genomeID,querys,degenerateposns,seedlen)
    
    example usage for finding orthogonal spacer/PAM sites in a genomic sequence:
    #get sequence
    acrB3prime='attgatacaacgtgtaatcactaaggccgcgtaagcggccttttttatgcataacctacgaacattaaggagtaatt'  
    # create list of one or more degenerate sites
    posns=[16,20]
    spacerlist=PAM_spacer_re_finder(acrB3,'NGG',2)
    hitlist=orthogonal_spacercheck2('U00096.3',spacerlist,posns,13)              
    for i in hitlist:
        print i
        
    out: list of spacer posn/strand, sequence, search sequence, number of sites, positions found in genome
    ('22Reverse', 'TGATACAACGTGTAATCACTAAG', 'TGTANTCACTNAG', 1, [481233])
    ('32Reverse', 'TGTAATCACTAAGGCCGCGTAAG', 'AAGGNCGCGTNAG', 2, [508547, 481223])
    ('35Reverse', 'AATCACTAAGGCCGCGTAAGCGG', 'GCCGNGTAAGNGG', 1, [481220])
    ('66Reverse', 'ATGCATAACCTACGAACATTAAG', 'TACGNACATTNAG', 2, [2168431, 481189])
    ('69Reverse', 'CATAACCTACGAACATTAAGGAG', 'GAACNTTAAGNAG', 1, [481186])
        
    
    """
    import BioModules
    from Bio import SeqUtils
    genome=BioModules.NCBI_genome2string(genomeID)
    targets=[]
    
    def degenerator(seq,posns):
        newseq=list(seq)
        for i in posns:
            newseq[i]='N'
        out=''.join(newseq)
        return out
    degenerates=[]
    for name,seq in querys:
        degenerates.append((name,seq,degenerator(seq,degenerateposns)))
    for name,seq,degen in degenerates:
        degen=degen[-seedlen:]
        F=SeqUtils.nt_search(genome,degen)[1:]
        R=SeqUtils.nt_search(genome,BioModules.revcomp(degen))[1:]
        matches=F+R
        targets.append((name,seq,degen,len(matches),matches))
    return targets


def mylocal_BLAST(dbhandle,outfile,evalue,fmtout,wordsize,queries):
    """
    for blasting sequences against local database ( can be a fasta file). Takes a python list of input queries
    
    mylocal_BLAST(dbhandle,outfile,evalue,fmtout,wordsize,queries):
    genomeID=NCBI ref number for genome to run BLAST against
    evalue=evalue cutoff
    querys is a list of name,sequence tuples that is to be blasted
    
    example usage:
    1) generate querieslist from fasta file
    queries=BioModules.fasta2list('/Users/andrewgarst/Desktop/gxl_TRACE_082314/seqsegments.fasta') 
    2) blast the files you want seperately against the file(s)
    
    ex:
    for myfile in BioModules.lsdir('/Users/andrewgarst/Desktop/gxl_TRACE_082314/', '.seq'):
        mylocal_BLAST(myfile,myfile.replace('.seq','_hits.txt'),10,0,7,queries)

    """
    import BioModules
    import os.path
    import os
    import subprocess
    handle1='/tmp/foo.txt'

    ##generate Blast database
    cmline='makeblastdb -in '+dbhandle+' -dbtype nucl -out '+dbhandle.replace('.fasta','.db')
    p = subprocess.Popen(cmline, stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    ##generate fasta file of queries for blasting
    f=open(handle1,'w')
    for name, seq in queries:
        print>>f, '>'+name+'\n+'+seq
    f.close()
    def myBLAST(queryfile, dbfile, outfile, eval, fmtout, wordsize): 
        """
        myBLAST(queryfile, dbfile, outfile, eval, fmtout, wordsize): 
            import os
            from Bio.Blast.Applications import NcbiblastnCommandline
            comline = NcbiblastnCommandline(query=queryfile, db=dbfile, strand="both",
                                      evalue=eval, out=outfile, outfmt=fmtout word_size=wordsize)
            os.system(str(comline))
    
        """
        import os
        from Bio.Blast.Applications import NcbiblastnCommandline
        comline = NcbiblastnCommandline(query=queryfile, db=dbfile, strand="both",
                                      evalue=eval, out=outfile, outfmt=fmtout, word_size=wordsize)
        os.system(str(comline))
    myBLAST(handle1,dbhandle.replace('.fasta','.db'),outfile,evalue,fmtout,wordsize)
    os.remove(handle1)
    os.remove(dbhandle.replace('.fasta','.db')+'.nin')
    os.remove(dbhandle.replace('.fasta','.db')+'.nhr')
    os.remove(dbhandle.replace('.fasta','.db')+'.nsq')

def read_fasta(fasta):
        DNAs = {}
        for line in fasta:
                if line[0] == '>':
                        name = line.strip()
                else:
                        DNA = ''.join(line.split())
  
                        if name in DNAs:
                                DNAs[name] += DNA
                        else:
                                DNAs[name] = DNA
        return DNAs

def excel2list():
    """
    excel2list()
    
    Will allow copy and paste of excel/google spreadsheets from the clipboard into a python list
    """
    import pandas as pd
    mylist=pd.read_clipboard(header=None,sep='\t')
    return mylist.values.tolist()

def gRNA_CPECinsert_designer(spacer):
    QiovlapF='CAGCTAGCTCAGTCCTAGGTATAATACTAGT'
    QiovlapR='GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC'
    spacer='ATATTCAGGGAGACCACAA'
    Fprimer=QiovlapF+spacer.lower()
    ovlapR=''
    count=1
    maxlen=30
    while count<=maxlen:
                if BioModules.Tm_calc(Fprimer[-count:],2,600,0)>60:
                    ovlapR=Fprimer[-count:]
                    break
                else:
                    pass
                count+=1 
    Rprimer=BioModules.revcomp(ovlapR+QiovlapR)
    return ['gfp_gRNAF',Fprimer,'gfp_gRNAR',Rprimer]    


def annotate_gb(seq,featurelist,outfile):
    
    """
    annotate_gb(seq,featurelist)
    input feature list with format of 
    posn,length,type,Fcolor,Rcolor,orientation ('F'is forward, 'R' is reverse' as so:
    
    [(27, 3, 'codon', 'green', 'blue', 'F'),
     (57, 3, 'codon', 'green', 'blue', 'F'), 
     (60, 3, 'codon', 'green', 'blue', 'F'),
     (75, 3, 'codon', 'green', 'blue', 'F'),
     (87, 3, 'codon', 'green', 'blue', 'F'),
     (132, 3, 'codon', 'green', 'blue', 'F'),
     (279, 3, 'codon', 'green', 'blue', 'F'),
     (342, 3, 'codon', 'green', 'blue', 'F'),
     (456, 3, 'codon', 'green', 'blue', 'F'),
     (471, 3, 'codon', 'green', 'blue', 'F')]
    
    """
    ################ A: Make a SeqRecord ################
    
    # 1. Create a sequence
    
    from Bio.Seq import Seq
    my_sequence = Seq(seq)
    
    # 2. Create a SeqRecord and assign the sequence to it
    
    from Bio.SeqRecord import SeqRecord
    my_sequence_record = SeqRecord(my_sequence)
    
    # 3. Assign an alphabet to the sequence (in this case DNA)
    
    from Bio.Alphabet import generic_dna
    my_sequence_record.seq.alphabet = generic_dna
    
    # This is the minimum required info for BioPython to be able to output 
    # the SeqRecord in Genbank format.
    # You probably would want to add other info (e.g. locus, organism, date etc)
    

    ################ B: Make a SeqFeature ################
    
    # 1. Create a start location and end location for the feature
    #    Obviously this can be AfterPosition, BeforePosition etc.,
    #    to handle ambiguous or unknown positions
    
    from Bio import SeqFeature
    from Bio.SeqFeature import FeatureLocation
    features=[]
    for posn,length,feattype,Fcolor,Rcolor,orientation in sorted(featurelist):
        mydict={'ApEinfo_fwdcolor':Fcolor,'ApEinfo_revcolor':Rcolor}
        # 2. Use the locations do define a FeatureLocation
        if orientation=='F':
            my_feature_location = FeatureLocation(posn,posn+length,strand=+1)
        
            # 3. Define a feature type as a text string 
            #     (you can also just add the type when creating the SeqFeature)
            my_feature_type = feattype
            # 4. Create a SeqFeature
            from Bio.SeqFeature import SeqFeature
            my_feature = SeqFeature(my_feature_location,type=my_feature_type,qualifiers=mydict)
            
            # 5. Append your newly created SeqFeature to your SeqRecord
            
            my_sequence_record.features.append(my_feature)
        #optional: print the SeqRecord to STDOUT in genbank format, with your new feature added.
        elif orientation=='R':
            my_feature_location = FeatureLocation(posn,posn+length,strand=-1)
        
            # 3. Define a feature type as a text string 
            #     (you can also just add the type when creating the SeqFeature)
            my_feature_type = feattype
            # 4. Create a SeqFeature
            from Bio.SeqFeature import SeqFeature
            my_feature = SeqFeature(my_feature_location,type=my_feature_type,qualifiers=mydict)
            
            # 5. Append your newly created SeqFeature to your SeqRecord
            
            my_sequence_record.features.append(my_feature)
    f=open(outfile,'w')
    f.write(my_sequence_record.format("gb"))
    f.close()


def PAMsite_finder(seq,PAM,systemtype):
    from BioModules import revcomp
    import re
    seq=seq.upper()
    def degenerate_re_convert(seq):
        conversiondict={'N':'(A|T|G|C)','R':'(A|G)','A':'A','C':'C','G':'G','T':'T'}
        newlist=list(seq)
        newseq=list()
        for i in newlist:
            for key, value in conversiondict.iteritems():
                if i==key:
                    newseq.append(value)
        return ''.join(newseq)
    
    targetlist=[]
    revPAM=revcomp(PAM)
    PAM=degenerate_re_convert(PAM.upper())
    revPAM=degenerate_re_convert(revPAM.upper())
    Flist=[match.start() for match in re.finditer(PAM,seq)]
    Rlist=[match.start() for match in re.finditer(revPAM,seq)]
    if systemtype==1:
        for i in Flist:
            target=seq[i:i+32]
            targetlist.append([i,seq[i-3:i],'Reverse',target])
        for i in Rlist:
            target=revcomp(seq[i-32:i])
            targetlist.append([i,seq[i+3],'Forward',target])

    if systemtype==2:
        for i in Flist:
            target=seq[i-20:i]
            if len(target)<20:
                pass
            else:
                targetlist.append([i,seq[i:i+3],'Reverse',target])
        for i in Rlist:
            target=revcomp(seq[i+3:i+23])
            if len(target)<20:
                pass
            else:
                targetlist.append([i,seq[i:i+3],'Forward',target])
    return targetlist    

def spacer_codon_finder(gene,PAM,aasites,start,maxdist,type):
    """
    makes a list of aa target, 
    spacer_codon_finder(gene,PAM,aasites,start,dist)
    ntsites=[(i*3-3+start,i) for i in aasites]
        edits=[]
        import BioModules
        for codonsite,aa_id in ntsites:
            for site, PAM,strand, spacer in PAMsite_finder(gene,PAM,2):
                for i in range(3,dist):
                    if 3<abs(site-codonsite)<=i:# check that the site is within a certain distance of the PAM
                        edits.append([aa_id,i, site,codonsite,PAM,strand,spacer])
                        break
                    else:
                        pass
        return sorted(edits) 
    """
    ntsites=[(i*3-3+start,i) for i in aasites]
    edits=[]
    import BioModules
    for codonsite,aa_id in ntsites:
        for site, PAM,strand, spacer in PAMsite_finder(gene,PAM,type):
            for i in range(0,maxdist):
                if abs(site-codonsite)<=maxdist:# check that the site is within a certain distance of the PAM
                    edits.append([aa_id,abs(site-codonsite), site,codonsite,PAM,strand,spacer])
                    break
                else:
                    pass
    return sorted(edits)         
    
def framefinder(orfstart,site):
    """
    framefinder(orfstart,site). input the start site of an ORF and
    the another short sequence (i.e. a sequencing read or oligo) and output
    which frame the use for other analysis
    for example if orfstart=100 and site=110 then the output should be
    1 meaning that to read the shorter sequence
    """
    def framecheck(num):
        return num%3==0
    frame=int()

    if framecheck(site-orfstart)==True:
        frame=0
    elif framecheck(site+1-orfstart)==True:
        frame=1
    elif framecheck(site-1-orfstart)==True:
        frame=-1
    elif framecheck(site+2-orfstart)==True:
        frame=2
    elif framecheck(site-2-orfstart)==True:
        frame=-2
    return frame



def syn_codon_finder(codon):
    from BioModules import dna2protein as D2P
    aa_wt=D2P(codon)
    # Ecoli_genetic code with in order of usage frequency
    gencode=  [['ATG', 'M', 1.0, 26.4, 96695.0],['TGG', 'W', 1.0, 13.9, 50991.0],['AAA', 'K', 0.74, 35.3, 129137.0],['GAA', 'E', 0.68, 39.1, 143353.0],['CAG', 'Q', 0.66, 28.4, 104171.0],
                ['GAT', 'D', 0.63, 32.7, 119939.0],['TAA', '*', 0.61, 2.0, 7356.0],['TAT', 'Y', 0.59, 17.5, 63937.0],['TTT', 'F', 0.58, 22.1, 80995.0],['CAT', 'H', 0.57, 12.5, 45879.0],
                ['TGC', 'C', 0.54, 6.1, 22188.0],['AAC', 'N', 0.51, 21.4, 78443.0],['ATT', 'I', 0.49, 29.8, 109072.0],['CCG', 'P', 0.49, 20.9, 76644.0],['AAT', 'N', 0.49, 20.6, 75436.0],
                ['CTG', 'L', 0.47, 48.4, 177210.0],['TGT', 'C', 0.46, 5.2, 19138.0],['CAC', 'H', 0.43, 9.3, 34078.0],['TTC', 'F', 0.42, 16.0, 58774.0],['TAC', 'Y', 0.41, 12.2, 44631.0],
                ['ACC', 'T', 0.4, 22.0, 80547.0],['ATC', 'I', 0.39, 23.7, 86796.0],['GAC', 'D', 0.37, 19.2, 70394.0],['GGC', 'G', 0.37, 27.1, 99390.0],['CGT', 'R', 0.36, 20.0, 73197.0],
                ['CGC', 'R', 0.36, 19.7, 72212.0],['GTG', 'V', 0.35, 24.4, 89265.0],['GGT', 'G', 0.35, 25.5, 93325.0],['CAA', 'Q', 0.34, 14.6, 53394.0],['GCG', 'A', 0.33, 30.1, 110308.0],
                ['GAG', 'E', 0.32, 18.7, 68609.0],['TGA', '*', 0.3, 1.0, 3623.0],['GTT', 'V', 0.28, 19.8, 72584.0],['GCC', 'A', 0.26, 24.2, 88721.0],['AAG', 'K', 0.26, 12.4, 45459.0],
                ['ACG', 'T', 0.25, 13.7, 50269.0],['AGC', 'S', 0.25, 15.2, 55551.0],['GCA', 'A', 0.23, 21.2, 77547.0],['GTC', 'V', 0.2, 14.3, 52439.0],['CCA', 'P', 0.2, 8.6, 31534.0],
                ['ACT', 'T', 0.19, 10.3, 37842.0],['CCT', 'P', 0.18, 7.5, 27340.0],['GCT', 'A', 0.18, 17.1, 62479.0],['GTA', 'V', 0.17, 11.6, 42420.0],['TCT', 'S', 0.17, 10.4, 38027.0],
                ['ACA', 'T', 0.17, 9.3, 33910.0],['AGT', 'S', 0.16, 9.9, 36097.0],['GGG', 'G', 0.15, 11.3, 41277.0],['TTA', 'L', 0.14, 14.3, 52382.0],['TCA', 'S', 0.14, 8.9, 32715.0],
                ['TCG', 'S', 0.14, 8.5, 31146.0],['TTG', 'L', 0.13, 13.0, 47500.0],['CCC', 'P', 0.13, 5.4, 19666.0],['GGA', 'G', 0.13, 9.5, 34799.0],['CTT', 'L', 0.12, 11.9, 43449.0],
                ['ATA', 'I', 0.11, 6.8, 24984.0],['CGG', 'R', 0.11, 5.9, 21552.0],['CTC', 'L', 0.1, 10.2, 37347.0],['TCC', 'S', 0.1, 59.1, 33430.0],['TAG', '*', 0.09, 0.3, 989.0],
                ['CGA', 'R', 0.07, 3.8, 13844.0],['AGA', 'R', 0.07, 3.6, 13152.0],['CTA', 'L', 0.04, 4.2, 15409.0],['AGG', 'R', 0.04, 2.1, 7607.0]]
    syncodons=list()
    for codonmut, aa, fraction,frequency,count in gencode:
        if aa==aa_wt and codonmut.upper()!=codon.upper():
            syncodons.append(codonmut)
    if len(syncodons)<1:
        syncodons.append(codon)
    return syncodons
    
def GenTRACR_ORF(organism,name,aasites,start,minlen,size,mutations,P1,P2):
    
    
    """
    GenTRACR_ORF(name,gene,aasites,start,minlen,size,mutations,P1,P2):
    example usage:

    proSARcodons={
              'ACG':'T','GCG':'A','TGC':'C','GAT':'D','GAA':'E',
              'TTT':'F','GGC':'G','CAT':'H','ATT':'I','AAA':'K',
              'CTG':'L','ATG':'M','AAC':'N','CCG':'P','CAG':'Q',
              'GTG':'V','TGG':'W','TAT':'Y','CGC':'R','AGC':'S'}
    mutlist=[codon for codon,aa in proSARcodons.iteritems()]
    GenTRACR_ORF('galK',galK,[145],100,3,100,['taa'],QiovlapF,QiovlapR)

    This will yield a list with the structure
    ['name', 'PAM', 'PAMmutation', 'PAMcodon', 'codon', 'distance', 'editingoligo', 'edit_size']

    
    """
    import BioModules
    from BioModules import revcomp as revc
    from BioModules import framefinder as framefinder
    from BioModules import spacer_codon_finder as spacer_codon_finder
    from BioModules import syn_codon_finder as syncodon
    gene=BioModules.NCBI_gene_fetcher(organism,name,100)[0][1]
    gencode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':u'*', 'TAG':u'*',
                'TGC':'C', 'TGT':'C', 'TGA':u'*', 'TGG':'W'}
    
    oligoslist1=[]
    oligoslist1.append('name,PAM,PAMmutation,PAMcodon,codon,distance,oligo,constructsize,edit_site_size'.split(','))
    uniqueset=set()
    #i,aa_id, PAMsite,codonsite,PAM,strand,spacer
    
    for aa_id,distance, PAMsite,codonsite,PAM,strand,spacer in spacer_codon_finder(gene,aasites,start,size-20):
        
        #find frame of PAM relative to coding sequence
        frame=PAMsite+framefinder(start,PAMsite)
        
        #assign wt target codon nt sequence
        targetcodon=gene[codonsite+1:codonsite+4]
        
        #assign target codon aa sequence
        targetaa=BioModules.dna2protein(targetcodon)
        
        #assign codon nearest PAM or just upstream (U) or downstream (D) of the immediate PAM
        PAMcodon=gene[frame:frame+3]
        PAMcodonD=gene[frame+3:frame+6]
        PAMcodonU=gene[frame-3:frame]
        
        #wt codon at PAM codon sites
        PAMaa=BioModules.dna2protein(PAMcodon)
        PAMaaD=BioModules.dna2protein(PAMcodonD)
        PAMaaU=BioModules.dna2protein(PAMcodonU)
        
        
        ### for each aa/codon in the genetic code see if there's a synonymous mutation that can simulatenously
        ### delete the PAM and make a silent mutation to the ORF of interest.
        newseq=gene
        if distance<minlen:
            pass
        else:
            for mutcodon in syncodon(PAMcodon):
                #print mutcodon
                newseq=gene[:frame]+mutcodon.lower()+gene[frame+3:]
                PAMmutation=newseq[PAMsite:PAMsite+3].upper()
                
                # check PAM orientation to ensure
                # and make sure PAM mutation is synonymous and that the PAM is disrupted
                if PAM.endswith('GG')==True and mutcodon!=PAMcodon and PAMmutation[1:]!='GG' and PAMmutation[1:]!='AG':
                    for mutation in mutations:
                        newseq2=newseq[:codonsite]+mutation.lower()+newseq[codonsite+len(mutation):]
                        targetmut=BioModules.dna2protein(mutation)
                        residue=((codonsite-start+3)/3)
                        description=targetaa+str(residue)+targetmut
                        #if codon
                        if PAMsite<codonsite:
                            for i in range(0,size/2):
                                oligo=newseq2[PAMsite-i:codonsite+i-1]
                                if len(oligo)>=size:
                                    #editingoligo=QiovlapF+oligo+'actcgag'+J23119+spacer.lower()+sgRNA
                                    editingoligo=P1+oligo+'actcgag'+J23119+spacer.lower()+P2
                                    #oligoname='-'+PAM+'/'+newseq2[PAMsite:PAMsite+3]+'-'+PAMcodon+'-'+codon+'-'+str(len(editingoligo))
                                    if description not in uniqueset:
                                        oligoslist1.append([name+description,PAM,PAMmutation,PAMcodon,mutcodon,distance,editingoligo,len(editingoligo),len(oligo)])
                                        uniqueset.add(description)
                                    else:
                                        pass
                        if PAMsite>codonsite:
                            for i in range(0,size/2):
                                oligo=newseq2[PAMsite-i:codonsite+i-1]
                                if len(oligo)>=size:
                                    #editingoligo=QiovlapF+oligo+'actcgag'+J23119+spacer.lower()+sgRNA
                                    editingoligo=P1+oligo+'actcgag'+J23119+spacer.lower()+P2
                                    #oligoname='-'+PAM+'/'+newseq2[PAMsite:PAMsite+3]+'-'+PAMcodon+'/'+codon+'-'+str(len(editingoligo))
                                    if description not in uniqueset:
                                        oligoslist1.append([name+description,PAM,PAMmutation,PAMcodon,mutcodon,distance,editingoligo,len(editingoligo),len(oligo)])
                                        uniqueset.add(description)
                                    else:
                                        pass
           
                ###make sure PAM mutation is synonymous and that the PAM is disrupted
                elif PAM.startswith('CC')==True and mutcodon!=PAMcodon and PAMmutation[:2]!='CT' and PAMmutation[:2]!='CC':
                    ### make codon substitutions at each site   
                    for mutation in mutations:
                        newseq2=newseq[:codonsite]+mutation.lower()+newseq[codonsite+len(mutation):]
                        targetmut=BioModules.dna2protein(mutation)
                        residue=((codonsite-start+3)/3)
                        description=targetaa+str(residue)+targetmut
                        if PAMsite<codonsite:
                            for i in range(0,size/2):
                                oligo=newseq2[PAMsite-i:codonsite+i-1]
                                if len(oligo)>=size:
                                    #editingoligo=QiovlapF+oligo+'actcgag'+J23119+spacer.lower()+sgRNA
                                    editingoligo=P1+oligo+'actcgag'+J23119+spacer.lower()+P2
                                    #oligoname='-'+PAM+'/'+newseq2[PAMsite:PAMsite+3]+'-'+PAMcodon+'/'+codon+'-'+str(len(editingoligo))
                                    if description not in uniqueset:
                                        oligoslist1.append([name+description,PAM,PAMmutation,PAMcodon,mutcodon,distance,editingoligo,len(editingoligo),len(oligo)])
                                        uniqueset.add(description)
                                    else:
                                        pass
                        elif PAMsite>codonsite:
                            for i in range(0,size/2):
                                oligo=newseq2[PAMsite-i:codonsite+i-1]
                                if len(oligo)>=size:
                                    #editingoligo=QiovlapF+oligo+'actcgag'+J23119+spacer.lower()+sgRNA
                                    editingoligo=P1+oligo+'actcgag'+J23119+spacer.lower()+P2
                                    #oligoname='-'+PAM++'/'+newseq2[PAMsite:PAMsite+3]+'-'+PAMcodon+'/'+codon+'-'+str(len(editingoligo))
                                    if description not in uniqueset:
                                        oligoslist1.append([name+description,PAM,PAMmutation,PAMcodon,mutcodon,distance,editingoligo,len(editingoligo),len(oligo)])
                                        uniqueset.add(description)
                                    else:
                                        pass 

    return oligoslist1


def NCBI_gene_fetcher(refseq_id,gene,flank):
    """
    NCBI_gene_fetcher(refseq_id,genenames,flank)
    
    input the a list of gene names and the reference sequence you'd like to pull them from as well as
    how much of the up and downstream flanking genomic sequence you'd like to include
    returns a list of the [gene name, seq, orf_start,orf_end]
    
    
    """
    from Bio import Entrez
    from Bio import SeqIO 
    import BioModules
    Entrez.email='ag.theseasquirt@gmail.com'
    import os
    genomehandle='/Users/andrea/repositories/design_pipeline/main/database/gb/'+refseq_id+'.gb'
    if os.path.isfile(genomehandle)==True:
        genome=SeqIO.read('/Users/andrea/repositories/design_pipeline/main/database/gb/'+refseq_id+'.gb','genbank')
    else:
        handle = Entrez.efetch(db="nucleotide", id=refseq_id ,rettype="gbwithparts",validate=True)
        genome = SeqIO.read(handle, 'genbank')
        SeqIO.write(genome,genomehandle,'genbank')
    
    geneseqslist=[]
    uniqueset=set()
    refgene=''
    for feature in genome.features:
        for key, value in feature.qualifiers.iteritems():
            if key=='gene':
                refgene=value[0]
            if key=='gene' and gene==value[0].upper():
                location=feature.location
                #position gives the actual python slice location of the start and end. need to reverse if on negative strand
                #could use to get upstream or downstream regions
                start=location.start.position
                end=location.end.position
                if feature.strand==1:
                    start=start-flank
                    end=end+flank
                    seq=str(genome.seq[start:end])
                if feature.strand==-1:
                    start=start-flank
                    end=end+flank
                    seq=str(genome.seq[start:end])
                    seq=BioModules.revcomp(seq)
                if gene not in uniqueset:
                    geneseqslist.append([gene,seq,start,end])
                    uniqueset.add(gene)
                    
                if len(geneseqslist)==0:
                    print key,'NCBI_gene_fetcher did not find gene in this organism'
            if key=='gene_synonym'and gene in value:
                location=feature.location
                #position gives the actual python slice location of the start and end. need to reverse if on negative strand
                #could use to get upstream or downstream regions
                start=location.start.position
                end=location.end.position
                if feature.strand==1:
                    start=start-flank
                    end=end+flank
                    seq=str(genome.seq[start:end])
                if feature.strand==-1:
                    start=start-flank
                    end=end+flank
                    seq=str(genome.seq[start:end])
                    seq=BioModules.revcomp(seq)
                if refgene+'/'+gene not in uniqueset:
                    geneseqslist.append([refgene+'/'+gene,seq,start,end])
                    uniqueset.add(refgene+'/'+gene)
                    
    if len(geneseqslist)==0:
        print key,'NCBI_gene_fetcher did not find gene in this organism'
    else:
        return geneseqslist
def NCBI_genebank_downloader(refseq_id,path,outfmt):
    """
     NCBI_genebank_downloader(refseq_id,outfmt):
        
    Given a refseq id and outfmt this module will download the sequence and write to fasta or genbank formats
    in the desired folder indicated by the path variable
    
    
    """
    from Bio import Entrez
    from Bio import SeqIO 
    import BioModules
    Entrez.email='ag.theseasquirt@gmail.com'
    import os
    filetype=''
    if outfmt=='fasta':
        filetype='.fasta'
        genomehandle='/Users/andrea/repositories/design_pipeline/main/database/gb/'+refseq_id+filetype
    elif outfmt=='gb':
        filetype='.gb'
        genomehandle='/Users/andrea/repositories/design_pipeline/main/database/gb/'+refseq_id+filetype
    if os.path.isfile(genomehandle)==True:
        genome=SeqIO.read('/Users/andrea/repositories/design_pipeline/main/database/gb/'+refseq_id+filtype,outfmt)
    else:
        handle = Entrez.efetch(db="nucleotide", id=refseq_id ,rettype="gb")
        genome = SeqIO.read(handle, outfmt)
        SeqIO.write(genome,genomehandle,outfmt)

def revtrans(aa_in):
    """
    revtrans(aa_in):

    will take an input amino acid (1 letter abbreviation), and return a lyst of the codons that would produce this amino acid in order of
    most frequently used in E. coli to least frequently used.
    """
    codonlist=[]
    gencode =[['ATG', 'M', 1.0, 26.4, 96695.0],['TGG', 'W', 1.0, 13.9, 50991.0],['AAA', 'K', 0.74, 35.3, 129137.0],['GAA', 'E', 0.68, 39.1, 143353.0],['CAG', 'Q', 0.66, 28.4, 104171.0],
                ['GAT', 'D', 0.63, 32.7, 119939.0],['TAA', '*', 0.61, 2.0, 7356.0],['TAT', 'Y', 0.59, 17.5, 63937.0],['TTT', 'F', 0.58, 22.1, 80995.0],['CAT', 'H', 0.57, 12.5, 45879.0],
                ['TGC', 'C', 0.54, 6.1, 22188.0],['AAC', 'N', 0.51, 21.4, 78443.0],['ATT', 'I', 0.49, 29.8, 109072.0],['CCG', 'P', 0.49, 20.9, 76644.0],['AAT', 'N', 0.49, 20.6, 75436.0],
                ['CTG', 'L', 0.47, 48.4, 177210.0],['TGT', 'C', 0.46, 5.2, 19138.0],['CAC', 'H', 0.43, 9.3, 34078.0],['TTC', 'F', 0.42, 16.0, 58774.0],['TAC', 'Y', 0.41, 12.2, 44631.0],
                ['ACC', 'T', 0.4, 22.0, 80547.0],['ATC', 'I', 0.39, 23.7, 86796.0],['GAC', 'D', 0.37, 19.2, 70394.0],['GGC', 'G', 0.37, 27.1, 99390.0],['CGT', 'R', 0.36, 20.0, 73197.0],
                ['CGC', 'R', 0.36, 19.7, 72212.0],['GTG', 'V', 0.35, 24.4, 89265.0],['GGT', 'G', 0.35, 25.5, 93325.0],['CAA', 'Q', 0.34, 14.6, 53394.0],['GCG', 'A', 0.33, 30.1, 110308.0],
                ['GAG', 'E', 0.32, 18.7, 68609.0],['TGA', '*', 0.3, 1.0, 3623.0],['GTT', 'V', 0.28, 19.8, 72584.0],['GCC', 'A', 0.26, 24.2, 88721.0],['AAG', 'K', 0.26, 12.4, 45459.0],
                ['ACG', 'T', 0.25, 13.7, 50269.0],['AGC', 'S', 0.25, 15.2, 55551.0],['GCA', 'A', 0.23, 21.2, 77547.0],['GTC', 'V', 0.2, 14.3, 52439.0],['CCA', 'P', 0.2, 8.6, 31534.0],
                ['ACT', 'T', 0.19, 10.3, 37842.0],['CCT', 'P', 0.18, 7.5, 27340.0],['GCT', 'A', 0.18, 17.1, 62479.0],['GTA', 'V', 0.17, 11.6, 42420.0],['TCT', 'S', 0.17, 10.4, 38027.0],
                ['ACA', 'T', 0.17, 9.3, 33910.0],['AGT', 'S', 0.16, 9.9, 36097.0],['GGG', 'G', 0.15, 11.3, 41277.0],['TTA', 'L', 0.14, 14.3, 52382.0],['TCA', 'S', 0.14, 8.9, 32715.0],
                ['TCG', 'S', 0.14, 8.5, 31146.0],['TTG', 'L', 0.13, 13.0, 47500.0],['CCC', 'P', 0.13, 5.4, 19666.0],['GGA', 'G', 0.13, 9.5, 34799.0],['CTT', 'L', 0.12, 11.9, 43449.0],
                ['ATA', 'I', 0.11, 6.8, 24984.0],['CGG', 'R', 0.11, 5.9, 21552.0],['CTC', 'L', 0.1, 10.2, 37347.0],['TCC', 'S', 0.1, 59.1, 33430.0],['TAG', '*', 0.09, 0.3, 989.0],
                ['CGA', 'R', 0.07, 3.8, 13844.0],['AGA', 'R', 0.07, 3.6, 13152.0],['CTA', 'L', 0.04, 4.2, 15409.0],['AGG', 'R', 0.04, 2.1, 7607.0]]
    for codon,aaref,fraction,freq_1K,count in gencode:
        if aa_in==aaref:
            codonlist.append(codon)
    return codonlist

def nt2aa_aligner_basic(ntseq,start,stop):
    """
    nt2aa_aligner_basic(ntseq,start,stop)
    input a nucleotide sequence and the positions to start and stop translation 

    """
    ntseq=ntseq[start:stop]
    trans=BioModules.dna2protein(ntseq)
    print ' '.join([ntseq[i:i+3] for i in range(0,len(ntseq),3)])
    print ' '+'   '.join(trans)


def csv2sortedlist(csv_file, delimiter,column,r):
    import operator
    import csv
    """ 
    Reads in a CSV file and returns the contents as a list sorted by the user chosen
    column,
    where every row is stored as a sublist, and each element
    in the sublist represents 1 cell in the table.
    adapted from
    http://nbviewer.ipython.org/github/rasbt/python_reference/blob/master/tutorials/sorting_csvs.ipynb
    
    """
    with open(csv_file, 'r') as csv_con:
        reader = csv.reader(csv_con, delimiter=delimiter)
        mylist=list(reader)
        for row in range(len(mylist)):
            for cell in range(len(mylist[row])):
                try:
                    mylist[row][cell] = float(mylist[row][cell])
                except ValueError:
                    pass 
        
    header = mylist[0]
    print header
    column=input()
    body = mylist[1:]
    if isinstance(column, str):  
        col_index = header.index(column)
    else:
        col_index = col
    body = sorted(body, 
           key=operator.itemgetter(col_index), 
           reverse=r)
    body.insert(0, header)
    return body
def list2df(mylist,mycolumns):
    """
    list2df(mylist,mycolumns)
    mylist is a list of input values 
    mycolumns is a list of column names to use as headers
    can handle any size 2D array
    """
    import pandas as pd
    mydf=pd.DataFrame(mylist,columns=mycolumns)
    return mydf

def csv2df(filename,s,columns):
    """
    sv2df(filename,s,columns)
    filename is the target file
    s= seperator (usually ',' or '\t')
    columns is a list of column names to use as headers i.e. ['A','B','C','D']
    can handle any size 2D array
    """
    import pandas as pd
    mydf=pd.read_csv(filename,sep=s,names=columns)
    return mydf

def primer_designer(targetname,seq,start,stop,tm,fracfold,GC,maxlen):
   
    """
    primer_designer(targetname,seq,start,stop,tm,fracfold,GC,maxlen):
    
    design primers to amplify sequence from start to stop positions.The 
    Sequence search will begin at these start/stop positions and move outward until suitable primers are found
    
    targetname is used to give each primer a name associated with what is being amplified
    
    seq is the target sequence
    
    start indicates the 3' end of where the forward primer anneals
    
    stop indicates the 3' end of where the reverse primer anneals
    
    tm is the minimal tm
    
    fracfold is the fraction of '(' or ')' characters in the RNA.fold string, used as a proxy to determine how much
    structure is in the oligo should be a number between 0-1
    
    GC is the maximal GC content desired for the primer
    
    maxlen is the maximal allowed primer length
    """
    import RNA
    import BioModules
    from Bio import SeqUtils
    from BioModules import Tm_calc as Tm
    import Bio
    for i in range(1, maxlen):
        Fprimer=seq[start-i:start]
        TmF=Tm(Fprimer,1,110,0)
        Ffold=RNA.fold(Fprimer).count('.')/float(len(Fprimer))
        GC_f=Bio.SeqUtils.GC(Fprimer)
        if TmF>=tm and Ffold>fracfold and 30<GC_f<GC:#site not in Fprimingsite and
            break
        else:
            continue
    for i in range(1, maxlen):
        Rprimer=BioModules.revcomp(seq[stop:stop+i])
        TmR=Tm(Rprimer,1,110,0)
        Rfold=RNA.fold(Rprimer)[0].count('.')/float(len(Rprimer))
        GC_r=Bio.SeqUtils.GC(Rprimer)
        if TmR>tm  and abs(TmR-TmF)<3 and Rfold>fracfold and 40<GC_r<GC:#site not in Fprimingsite and
           
            break
        else:
            continue
    mylist=[]
    mylist.append([targetname+'_F',Fprimer,TmF,Ffold,GC_f,len(Fprimer)])
    mylist.append([targetname+'_R',Rprimer,TmR,Rfold,GC_r,len(Rprimer)])
    mylist.append([targetname,seq,'','','',''])
    return mylist

    def take(n, iterable):
        """
        take(n,iterable)
        
        n=number of iterations to sample
        iterable= iterable is the iterator tosample from
        returns entries in list of tuples
        """
        from itertools import islice
        "Return first n items of the iterable as a list"
        return list(islice(iterable, n))


def BCmatcher(seq,BClist):
    """
    BCmatcher(seq,BClist)
    input a string (seq) and a barcode list that is comprised of two columns, BC:seq
    return the match with the closest distance to the barcode as an (ID,distance) tuple
    """
    import Levenshtein as Lev
    d=dict()
    for i in BClist:
         d[id]=Lev.distance(seq,BC)
    return min(d.items(), key=lambda x: x[1])

def gz2list(gzfile):
    """
    read generic .gz file to python list
    
    """
    import gzip
    mylist=[]
    f=gzip.open(gzfile,'r')
    for line in f:
        line=line.strip().rstrip()
        mylist.append(line)
    return mylist

def fastq_gz2list(gzfile):
    """
    read fastq.gz file to python list of name,seq,quality tuples
    
    """
    import gzip
    from itertools import izip
    f=gzip.open(gzfile,'rb')
    mylist1=enumerate(izip(f))
    mylist2=[]
    for i, line in mylist1:
        try:
            seq=mylist1.next()
            spacer=mylist1.next()
            qual=mylist1.next()

            ID=line[0]
            read=seq[1][0]
            quality=qual[1][0]

            mylist2.append([ID,read,quality])
        except:
            StopIteration
    return mylist2

def fastq_gz_iter(gzfile):
    """
    read fastq.gz file and yield name,seq,quality lists
    
    """
    import gzip
    from itertools import izip
    f=gzip.open(gzfile,'rb')
    mylist1=enumerate(izip(f))
    mylist2=[]
    for i, line in mylist1:
        try:
            seq=mylist1.next()
            spacer=mylist1.next()
            qual=mylist1.next()

            ID=line[0]
            read=seq[1][0]
            quality=qual[1][0]

            yield [ID,read,seq]
        except:
            StopIteration

def pyls(suffix):
    import os
    mylist=[]
    for subdir, dirs, files in os.walk('./'):
            for file in files:
                if file.endswith(suffix):
                    mylist.append(file)
    return mylist

def tab_gz_iter(infile,sep):
    """
    tab_gz_iter(infile,sep)
        import gzip
        f=gzip.open(infile,'r')
        for line in f:
            line=line.rstrip().strip().split(sep)
            yield line
    """
    import gzip
    f=gzip.open(infile,'r')
    for line in f:
        line=line.rstrip().strip().split(sep)
        yield line

def codon2aa_dist(codonwt,aamut):
    """
    Will tell you the minimal mutation distance required to convert from a given codon to a target amino acid.
    def codon2aa_dist(codonwt,aamut)
        import Levenshtein
        distlist=[]
        for i in BioModules.revtrans(aamut):
            distlist.append(Levenshtein.distance(codonwt,aamut))
        return min(distlist

    """
    
    import Levenshtein
    distlist=[]
    for i in BioModules.revtrans(aamut):
        distlist.append(Levenshtein.distance(codonwt,aamut))
    return min(distlist)   
def pyUSEARCH(queryfile,db,percentid,threads,strand,maxaccepts,maxrejects,out):
    '''
    pyUSEARCH(queryfile,db,percentid,threads,strand,maxaccepts,maxrejects,out)
    calls USEARCH8 from commandline. 
    Input files must be in fasta format
    -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits '
    '''
    import subprocess
    program = 'usearch8 '
    module=' -usearch_global '
    query= queryfile
    db=' -db '+db
    pid=' -id '+percentid
    strand=' -strand  ' +strand
    maxaccept=' -maxaccepts '+ maxaccepts 
    maxreject=' -maxrejects '+ maxrejects
    out=' -userout  '+ out
    fields=' -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits '
    cmd = [program+module+query+db+threads+pid+strand+maxaccept+maxreject+out+fields]
    subprocess.call(cmd, shell=True)

def pyUBLAST(queryfile,db,evalue,strand,outfile,outtype,maxhits):
    '''
     pyUBLAST(queryfile,db,evalue,strand,outfile,outtype,maxhits)
    
    calls USEARCH8 from commandline. Input files must be in fasta or blast db format
    can use 'alnout'(human readable) 
    or 'userout' (-userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits ')
    as outtype
    '''
    import subprocess

    program = 'usearch8 '
    module=' -ublast '
    query= queryfile
    db=' -db '+db
    evalue=' -evalue '+str(evalue)
    strand=' -strand  ' +strand
    hits=' -maxaccepts' + str(maxhits)
    if outtype=='userout':
        out=' -userout ' +outfile
        fields=' -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits '
        cmd = [program+module+query+db+evalue+strand+out+fields]
        subprocess.call(cmd, shell=True)
    if outtype=='alnout':
        out=' -alnout ' +outfile
        cmd = [program+module+query+db+evalue+strand+out]
        subprocess.call(cmd, shell=True)


def gb_writer(seq,featurelist,outfile):
    
    """
    annotate_gb(seq,featurelist,outfile)
    input feature list with format of 
    feature start,feature stop,type,Fcolor,Rcolor,orientation ('F'is forward, 'R' is reverse' as so:
    
   [[0, 140, 'oligo', 'blue', 'orange', 'F'], 
   [120, 260, 'oligo', 'blue', 'orange', 'F'], 
   [237, 377, 'oligo', 'blue', 'orange', 'F'], 
   [354, 494, 'oligo', 'blue', 'orange', 'F'], 
   [474, 614, 'oligo', 'blue', 'orange', 'F'], 
   [593, 733, 'oligo', 'blue', 'orange', 'F'], 
   [711, 851, 'oligo', 'blue', 'orange', 'F'], 
   [772, 912, 'oligo', 'blue', 'orange', 'F']]
    
    """
    ################ A: Make a SeqRecord ################
    
    # 1. Create a sequence
    
    from Bio.Seq import Seq
    my_sequence = Seq(seq)
    
    # 2. Create a SeqRecord and assign the sequence to it
    
    from Bio.SeqRecord import SeqRecord
    my_sequence_record = SeqRecord(my_sequence)
    
    # 3. Assign an alphabet to the sequence (in this case DNA)
    
    from Bio.Alphabet import generic_dna
    my_sequence_record.seq.alphabet = generic_dna
    
    # This is the minimum required info for BioPython to be able to output 
    # the SeqRecord in Genbank format.
    # You probably would want to add other info (e.g. locus, organism, date etc)
    

    ################ B: Make a SeqFeature ################
    
    # 1. Create a start location and end location for the feature
    #    Obviously this can be AfterPosition, BeforePosition etc.,
    #    to handle ambiguous or unknown positions
    
    from Bio import SeqFeature
    from Bio.SeqFeature import FeatureLocation
    features=[]
    for start,stop,feattype,Fcolor,Rcolor,orientation in sorted(featurelist):
        mydict={'ApEinfo_fwdcolor':Fcolor,'ApEinfo_revcolor':Rcolor}
        # 2. Use the locations do define a FeatureLocation
        if orientation=='F':
            my_feature_location = FeatureLocation(start,stop,strand=+1)
        
            # 3. Define a feature type as a text string 
            #     (you can also just add the type when creating the SeqFeature)
            my_feature_type = feattype
            # 4. Create a SeqFeature
            from Bio.SeqFeature import SeqFeature
            my_feature = SeqFeature(my_feature_location,type=my_feature_type,qualifiers=mydict)
            
            # 5. Append your newly created SeqFeature to your SeqRecord
            
            my_sequence_record.features.append(my_feature)
        #optional: print the SeqRecord to STDOUT in genbank format, with your new feature added.
        elif orientation=='R':
            my_feature_location = FeatureLocation(start,stop,strand=-1)
        
            # 3. Define a feature type as a text string 
            #     (you can also just add the type when creating the SeqFeature)
            my_feature_type = feattype
            # 4. Create a SeqFeature
            
            from Bio.SeqFeature import SeqFeature
            my_feature = SeqFeature(my_feature_location,type=my_feature_type,qualifiers=mydict)
            
            # 5. Append your newly created SeqFeature to your SeqRecord
            
            my_sequence_record.features.append(my_feature)
    f=open(outfile,'w')
    f.write(my_sequence_record.format("gb"))
    f.close()

def Rsitefinder(Rsite,seq):
    """
    Rsitefinder(Rsite,seq)
    ex. 
    BsaI='GGTCTC'
    Rsitefinder(BsaI,seq)
    will return a list of sites where this restriction site (Rsite) occurs
    in the longer sequence (both forward and reverse complements)

    """
    import re
    import BioModules
    
    sites=[match.start() for match in re.finditer(Rsite.upper(),seq.upper())]
    rsites=[match.start() for match in re.finditer(BioModules.revcomp(Rsite.upper()),seq.upper())]
    finallist=sites+rsites
    return sites

def syncodon_replace(seq,sites,ORFstart):
    """
    syncodon_replace(seq,sites,ORFstart)    
    seq=input sequence
    sites=sites where synonymous codon mutations are to be made
    ORFstart=position of the ATG codon in the string

    This module can be used generally or in combination with the Rsitefinder 
    function to remove cut sites from an ORF for gene construction purposes

    """
    seq=seq.lower()
    import BioModules
    wtcodon1=''
    wtcodon2=''
    mutcodon1=''
    mutcodon2=''
    
    for i in sites:
        posn=BioModules.framefinder(ORFstart,i)+i
        try:
            wtcodon1=seq[posn:posn+3]
            wtcodon2=seq[posn+3:posn+6]
            mutcodon1=BioModules.syn_codon_finder(wtcodon1)[0]
            mutcodon2=BioModules.syn_codon_finder(wtcodon2)[0]
        except:
            NameError
        
        if wtcodon1!=mutcodon1:
            outseq=seq[:posn]+mutcodon1.upper()+seq[posn+3:]
            break
        if wtcodon2!=mutcodon2:
            outseq=seq[:posn+3]+mutcodon2.upper()+seq[posn+6:]
            break
        if wtcodon1==mutcodon1 and wtcodon2==mutcodon2:
            print "No mutation found"
    #print mutcodon1,wtcodon1, mutcodon2,wtcodon2
    #print BioModules.dna2protein(seq[ORFstart:])
    #print BioModules.dna2protein(outseq[ORFstart:])
    return outseq


