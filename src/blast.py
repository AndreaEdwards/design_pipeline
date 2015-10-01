import os
from Bio.Blast import NCBIWWW, NCBIXML
from src.util import FileHandlers
import subprocess
import settings

class wwwBLAST:
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

class localBLAST:
	def __init__(self, pdb_code, sequence_annotations):
		self.filename = pdb_code
		self.sequence_annotations = sequence_annotations

	def _pairwise_blast(self, query_files, subject_files, protein=False, nucleotide=False):
		file_handlers = FileHandlers()
		blast_command = 'blastp' if protein == True else 'blastn'
		print query_files
		print subject_files
		for query, subject in zip(query_files, subject_files):
			outfile = ( file_handlers.get_file_name(query).split('.')[0] + file_handlers.get_file_name(subject).split('.')[0] + '_blastp.out')
			cmd = [blast_command + ' -query ' + query + ' -subject ' + subject + ' > ' + outfile]
			subprocess.call(cmd, shell=True)

	def _blast_pairwise(self, query_files, subject_files, protein=False, nucleotide=False):
		if ( protein == True ) and ( nucleotide == False ):
			self._pairwise_blast(query_files, subject_files, protein=True)
		elif ( protein == False ) and ( nucleotide == True ):
			self._pairwise_blast(query_files, subject_files, nucleotide=True)
		elif ( protein == False ) and ( nucleotide == False ):
			print "Please indicate either use of blastp or blastn by providing the keyword argument protein=True or nucleotide=True"
		else:
			print "You may only specify protein=True or nucleotide=True, not both"

	def _get_fasta_file_paths(self):
		query_files = []
		subject_files = []
		for i in range(len(self.sequence_annotations)):
			organism_id = settings.ORGANISM_MAP[self.sequence_annotations[i][2].upper()]
			print organism_id
			print self.sequence_annotations[i][1] + '_' + organism_id
			file_handlers = FileHandlers()
			file_paths = file_handlers.search_directory()
			fasta_files = file_handlers.find_files(file_paths, 'fasta')
			for fasta_file in fasta_files:
				if ( self.sequence_annotations[i][1] + '_' + organism_id ) == file_handlers.get_file_name(fasta_file).split('.fasta')[0]:
					query_files.append(fasta_file)
				elif ( self.sequence_annotations[i][1] + '_' + self.filename + '_chain-' + self.sequence_annotations[i][0] ) == file_handlers.get_file_name(fasta_file).split('.')[0]:
					subject_files.append(fasta_file)
		return query_files, subject_files

	def align_pdb_seqs(self):
		query_files, subject_files = self._get_fasta_file_paths()
		self._blast_pairwise(query_files, subject_files, protein=True)


class localPSIBLAST:
	def __init__(self):
		pass

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

#main()
