import os
from Bio import Entrez, SeqIO
from src.util import FileHandlers

class GenBank:
	def __init__(self):
		Entrez.email = "andrea.edwards@colorado.edu"
		self.dir_path = os.getcwd() + '/genbank_records'
		self._mkdir()
		self.genbank_record_number = ''
		self.handle = ''

	def _mkdir(self):
		file_handlers = FileHandlers()
		file_handlers.make_results_folder(self.dir_path.split('/')[-1])

	def _save_genbank_record(self):
		print "Saving .gb for record number %s" % self.genbank_record_number
		file_name = str(self.genbank_record_number) + ".gb"
		file_path = self.dir_path + '/' + file_name
		save_file = open(file_path, "w")
		save_file.write(self.handle)
		save_file.close()

	def _output_fasta(self):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		gb_files = file_handlers.find_files(file_paths, 'gb')
		for gb_file in gb_files:
			if self.genbank_record_number == file_handlers.get_file_name(gb_file).split('.')[0]:
				SeqIO.convert(gb_file, 'genbank', (file_handlers.get_file_name(gb_file) + '.fasta'), 'fasta')

	def fetch_record(self, genbank_id):
		self.genbank_record_number = str(genbank_id)
		handle = Entrez.efetch(db="protein", id=self.genbank_record_number, rettype="gb", retmode="text")
		self.handle = handle.read()
		handle.close()
		self._save_genbank_record()
		self._output_fasta()


