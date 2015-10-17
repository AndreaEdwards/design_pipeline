import os
import urllib2
import subprocess
import shutil
from util import FileHandlers
import numpy
import pandas
from operator import itemgetter
import matplotlib.pyplot as plt
import itertools
import amino_acids
import settings
import BioModules

class MutationListGenerator:
	def __init__(self, sequence_annotations, dna_sequences, aln_diff, pdb_code='', psiblast=''):
		self.sequence_annotations = sequence_annotations
		self.pdb_code = pdb_code
		self.psiblast = psiblast
		self.dna_sequences = dna_sequences
		self.aln_diff = aln_diff
		self.ddG_results_filepath = ''
		self.ddG_data = []
		self.ddG_data_map = {}
		self.llikelihood_filepath = ''
		self.llikelihood_data = []
		self.mutation_map = {}

	def _get_outfile(self):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		out_files = file_handlers.find_files(file_paths, 'out')
		if ( self.pdb_code != '' and self.psiblast != '' ):
			for out_file in out_files:
				if ( self.pdb_code + '_mutants' ) == file_handlers.get_file_name(out_file).split('.')[0]:
					self.ddG_results_filepath = out_file
				elif ( self.psiblast + '_mutants' ) == file_handlers.get_file_name(out_file).split('.')[0]:
					self.llikelihood_filepath = out_file
		elif ( self.pdb_code != '' and self.psiblast == '' ):
			for out_file in out_files:
				if ( self.pdb_code + '_mutants' ) == file_handlers.get_file_name(out_file).split('.')[0]:
					self.ddG_results_filepath = out_file
					print "Fetching data from %s ....." % file_handlers.get_file_name(out_file)
		elif ( self.pdb_code == '' and self.psiblast != '' ):
			for out_file in out_files:
				if ( self.psiblast + '_mutants' ) == file_handlers.get_file_name(out_file).split('.')[0]:
					self.llikelihood_filepath = out_file
		else:
			print "You have not specified any results data to parse."
			exit(1)


	def _capture_ddG_data(self):
		ddG_Data = open(self.ddG_results_filepath, 'rU')
		self.ddG_data = ddG_Data.readlines()
		ddG_Data.close()

	def _capture_llikelihood_data(self):
		llikelihood_Data = open(self.llikelihood_filepath, 'rU')
		self.llikelihood_data = llikelihood_Data.readlines()
		llikelihood_Data.close()

	def _get_data(self):
		self._get_outfile()
		if ( self.ddG_results_filepath != '' and self.llikelihood_filepath != '' ):
			self._capture_ddG_data()
			self._capture_llikelihood_data()
		elif ( self.ddG_results_filepath != '' and self.llikelihood_filepath == '' ):
			self._capture_ddG_data()
		elif ( self.ddG_results_filepath == '' and self.llikelihood_filepath != '' ):
			self._capture_llikelihood_data()
		else:
			print "Program should have quit with exit code 1 in the event no data files were specified.. something weird happended."


	def _parse_ddG_data(self):
		file_handlers = FileHandlers()
		self._get_data()
		ddG_data_map = {}
		for i in range(len(self.ddG_data)):
			fields = self.ddG_data[i].split(' ')
			cleaned = file_handlers.clean(fields)
			while cleaned.count('') > 0:
				cleaned.remove('')
			if len(cleaned[1].split("-")) < 2:	## ignore first line
				pass
			else:
				chain, mutation = cleaned[2].split("-")
				wt_res, position, mut_res = mutation[0], mutation[1:-1], mutation[-1]
				ddG_data_map[(chain, wt_res, position, mut_res)] = cleaned[3]
		self.ddG_data_map = ddG_data_map

	def write_csv_output(self, aln_diff):
		self._parse_ddG_data()
		self._build_mut_seq_list()
		out_filename = self.pdb_code + '_mutation_list.csv'
		output = open(out_filename, 'w')
		output_dictionary = {}
		for sequence in self.sequence_annotations:
			output_dictionary[sequence] = {}
		for key in self.mutation_map:
			start_site = '100'
			chain = key[0]
			wild_type_residue = key[1]
			#position = str( int(key[2]) ) # for pdb numbering that is same as genomic numbering # worked for 4JDR, 3NBU, 1S7C, 2CMD, 
			#position = str( int(key[2]) + int(self.aln_diff[0][1]) )
			#position = str( int(key[2]) - 1000 ) # for 1X15 where numbers begin with 1001
			#position = str( int(key[2]) - 999 ) # for 1SRU also begins with 1001
			#position = str( int(key[2]) + 1 ) # for 4N72 where numbers begin with 383, and for 4TWZ begins with 6; worked for 1H16
			position = str( int(key[2]) + 26 ) # 3TCH
			codons = ",".join(list(self.mutation_map[key]))
			for sequence in output_dictionary:
				if chain == sequence[0]:
					output_dictionary[sequence][(start_site, wild_type_residue, position, codons)] = []
		for key in output_dictionary:
			chain = key[0]
			gene_name = key[1]
			organism = key[2]
			aa_sequence = key[3]
			for i in range(len(self.dna_sequences)):
				for sequence in self.dna_sequences[i]:
					if gene_name == sequence[0]:
						dna_sequence = sequence[1]
			for item in output_dictionary[key]:
				start_site = item[0]
				wild_type_residue = item[1]
				position = item[2]
				codons = item[3]
				output.write(gene_name + '	' + dna_sequence + '	' + start_site + '	' + wild_type_residue + '	' + position + '	' + codons + '\n')
				self._check_sequence_accuracy(wild_type_residue, position, dna_sequence, start_site)
		self._print_total_seq_count(out_filename)
		output.close

	def _print_total_seq_count(self, out_filename):
		read_output = open(out_filename, 'rU')
		read_output_lines = read_output.readlines()
		read_output.close()
		print "length of mutation_map is:", len(self.mutation_map)
		print "length of output_file is: ", len(read_output_lines)
		print "total length of library is: ", len(self.ddG_data_map)


	def _check_sequence_accuracy(self, wild_type_residue, position, dna_sequence, start_site):
		print "position is:", position
		start = ((int(position) * 3) + int(start_site)) - 3
		stop = start + 3
		residue = BioModules.dna2protein(dna_sequence[start:stop])
		print "residue is: ", residue
		print "what i think the residue is: ", wild_type_residue

	def _build_mut_seq_list(self):
		#self._parse_ddG_data()
		mutation_map = {}
		for key in self.ddG_data_map:
			if (key[0], key[1], key[2]) not in mutation_map:
				mutation_map[(key[0], key[1], key[2])] = set([settings.E_COLI_CODON_USAGE[key[3]][0][0]])
			else:
				#mutation_map[key[1]].add(key[2])
				mutation_map[(key[0], key[1], key[2])].add(settings.E_COLI_CODON_USAGE[key[3]][0][0])
		self.mutation_map = mutation_map
				
	def _test_codon_usage(self):
		for key in settings.E_COLI_CODON_USAGE:
			total = 0
			for item in settings.E_COLI_CODON_USAGE[key]:
				total += item[1]
			print key, total

class SequenceGetter:
	def __init__(self, genes, organism):
		self.genes = genes # list
		self.dna_sequences = []
		self.protein_sequences = []
		self.file_path = ''
		self.organism = settings.ORGANISM_MAP[''.join(set(organism)).upper()]

	def _get_pdb_file_path(self):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		pdb_files = file_handlers.find_files(file_paths, 'pdb')
		for pdb_file in pdb_files:
			if (self.filename + '_0001') == file_handlers.get_file_name(pdb_file).split('.')[0]:
				self.file_path = pdb_file

	def _write_fasta(self, seq_list):
		for i in range(len(seq_list)):
			seq_id = seq_list[i][0] + '_' + self.organism
			outfile = os.getcwd() + '/' + seq_id + '.fasta'
			output = open(outfile, 'w')
			output.write('>' + seq_list[i][0] + '\n' + seq_list[i][1] + '\n')
			output.close()

	def get_DNA_sequence(self):
		for gene in self.genes:
			dna_sequence = BioModules.NCBI_gene_fetcher(self.organism, gene, 100)
			self.dna_sequences.append(dna_sequence)
		return self.dna_sequences

	def get_protein_sequence(self):
		for i in range(len(self.dna_sequences)):
			for sequence in self.dna_sequences[i]:
				self.protein_sequences.append(BioModules.dna2protein(sequence[1][100:-100]))
			seq_list = []
			for gene, sequence in zip(self.genes, self.protein_sequences):
				seq_list.append((gene, sequence))
			self._write_fasta(seq_list)

	def get_pdb_sequence(self):
		pass












