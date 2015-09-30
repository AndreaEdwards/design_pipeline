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

class SeqGenerator:
	def __init__(self, pdb_code='', psiblast=''):
		self.pdb_code = pdb_code
		self.psiblast = psiblast
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
			if len(cleaned[1].split("-")) < 2:
				pass
			else:
				chain, mutation = cleaned[1].split("-")
				wt_res, position, mut_res = mutation[0], mutation[1:-1], mutation[-1]
				ddG_data_map[(wt_res, position, mut_res)] = cleaned[3]
		self.ddG_data_map = ddG_data_map

	def write_csv_output(self):
		self._parse_ddG_data()
		self._build_mut_seq_list()
		out_filename = self.pdb_code + '_mutation_list.csv'
		output = open(out_filename, 'w')
		for key in self.mutation_map:
			output.write(key[0] + ',' + key[1] + ',' + '[' + ' '.join(list(self.mutation_map[key])) + ']' + '\n')
		output.close

	def _build_mut_seq_list(self):
		self._parse_ddG_data()
		mutation_map = {}
		for key in self.ddG_data_map:
			if (key[0], key[1]) not in mutation_map:
				mutation_map[(key[0], key[1])] = set([settings.e_coli_codon_usage[key[2]][0][0]])
			else:
				#mutation_map[key[1]].add(key[2])
				mutation_map[(key[0], key[1])].add(settings.e_coli_codon_usage[key[2]][0][0])
		self.mutation_map = mutation_map
				
	def _test_codon_usage(self):
		for key in settings.e_coli_codon_usage:
			total = 0
			for item in settings.e_coli_codon_usage[key]:
				total += item[1]
			print key, total










