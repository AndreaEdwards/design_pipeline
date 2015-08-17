import os
import re
import subprocess

class PfamFromUniprot:
	def __init__(self, uniprot_id):
		self.uniprot_id = uniprot_id
		self.data = []

	def _pfam_to_temp_file(self, uniprot_id):
		cmd = ['python wsdbfetch_suds.py fetchData UniProtKB:' + uniprot_id + ' | grep pfam | grep -v "SUPFAM" > pfam_temp.out']
		subprocess.call(cmd, shell=True)
		
	def _parse_temp_file(self):
		temp_file = os.getcwd() + '/pfam_temp.out'
		Data = open(temp_file)
		self.data = Data.readlines()
		Data.close()
		os.remove(temp_file)

	def get_pfam_id(self):
		self._pfam_to_temp_file(self.uniprot_id)
		self._parse_temp_file()
		for line in self.data:
			pfam_id = re.search('>(.*)<', line)
			print pfam_id.group(1)

def main():
	uniprot_id = 'P0A9Q1'
	pfam_getter = PfamFromUniprot(uniprot_id)
	id_list = pfam_getter.get_pfam_id()

#main()
