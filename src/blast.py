import os
from Bio.Blast import NCBIWWW, NCBIXML
from src.util import FileHandlers

class BLAST:
	def __init__(self, genbank_record_number):
		self.genbank_record_number = genbank_record_number
		self.dir_path = os.getcwd() + '/blast_results'
		self.xml_result = ''
		self._mkdir() 

	def blast_record(self):
		print "BLASTing record number %d ..." % int(self.genbank_record_number)
		result_handle = NCBIWWW.qblast("blastp", "nr", self.genbank_record_number)
		print "extracting result..."
		self.xml_result = result_handle.read()
		result_handle.close()
		return self.xml_result
	
	def save_blast_xml(self):
		print "Saving BLAST result (.xml) for record number %d" % int(self.genbank_record_number)
		file_name = str(self.genbank_record_number) + "_blast.xml"
		file_path = self.dir_path + '/' + file_name
		save_file = open(file_path, "w")
		save_file.write(self.xml_result)
		save_file.close()

	def parse_xml(self, file_path):
		result_handle = open(file_path)
		blast_record = NCBIXML.read(result_handle)
		result_handle.close()
		return blast_record

	def _mkdir(self):
		file_handlers = FileHandlers()
		file_handlers.make_results_folder(self.dir_path.split('/')[-1])


#def main():
##	Example of how to use the xml output:
#	file_paths = file_handlers.search_directory()
#	xml_files = file_handlers.find_files(file_paths, 'xml')
#	for xml_file in xml_files:
#		blast_record = blast.parse_xml(xml_file)
#		for description in blast_record.description:
#			print description.title
#		for alignment in blast_record.alignments:
#			#print alignment
#			for hsp in alignment.hsps:
#				#print('**** Alignment ****')
#				#print('sequence:', alignment.title)
#				#print('length:', alignment.length)
#				#print('e value', hsp.expect)
#				#print(hsp.query[0:75] + '...')
#				#print(hsp.match[0:75] + '...')
#				#print(hsp.sbjct[0:75] + '...')
#				print(hsp.identities)
#

#main()
