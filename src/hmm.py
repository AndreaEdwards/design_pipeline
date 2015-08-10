import subprocess
from util import FileHandlers

class HMM:
	def __init__(self):
		self.genbank_id = ''
		self.gb_file_path = ''

	def _get_gb_record(self):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		gb_files = file_handlers.find_files(file_paths, 'gb')
		print gb_files
		print self.genbank_id
		for gb_file in gb_files:
			if self.genbank_id == file_handlers.get_file_name(gb_file).split('.')[0]:
				self.gb_file_path = gb_file
				print self.gb_file_path

	def _hmm_to_temp_file(self):
		cmd = ['hmmscan ../src/database/Pfam-A.hmm ' + self.gb_file_path]
		subprocess.call(cmd, shell=True)

	def fetch_hmm(self, genbank_id):
		self.genbank_id = str(genbank_id)
		self._get_gb_record()
		self._hmm_to_temp_file()


