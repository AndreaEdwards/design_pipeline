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
		blast_result_files = []
		for query, subject in zip(query_files, subject_files):
			blast_result = ( file_handlers.get_file_name(query).split('.')[0] + file_handlers.get_file_name(subject).split('.')[0] + '_blastp.out')
			cmd = [blast_command + ' -query ' + query + ' -subject ' + subject + ' > ' + blast_result]
			subprocess.call(cmd, shell=True)
			blast_result_files.append(blast_result)
		return blast_result_files

	#def _blast_pairwise(self, query_files, subject_files, protein=False, nucleotide=False):
	#	if ( protein == True ) and ( nucleotide == False ):
	#		self._pairwise_blast(query_files, subject_files, protein=True)
	#	elif ( protein == False ) and ( nucleotide == True ):
	#		self._pairwise_blast(query_files, subject_files, nucleotide=True)
	#	elif ( protein == False ) and ( nucleotide == False ):
	#		print "Please indicate either use of blastp or blastn by providing the keyword argument protein=True or nucleotide=True"
	#	else:
	#		print "You may only specify protein=True or nucleotide=True, not both"

	def _get_fasta_file_paths(self):
		query_files = []
		subject_files = []
		for i in range(len(self.sequence_annotations)):
			organism_id = settings.ORGANISM_MAP[self.sequence_annotations[i][2].upper()]
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
		blast_result_files = self._pairwise_blast(query_files, subject_files, protein=True)
		aln_diff = self._get_aln_diff(blast_result_files)
		return aln_diff

	def _open_results_file(self, result_file):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		out_files = file_handlers.find_files(file_paths, 'out')
		for out_file in out_files:
			if result_file == file_handlers.get_file_name(out_file):
				Data = open(out_file, 'r')
				data = Data.readlines()
				Data.close()
		return data

	def _get_aln_diff(self, blast_result_files):
		file_handlers = FileHandlers()
		aln_diff = []
		for result_file in blast_result_files:
			data = self._open_results_file(result_file)
			cleaned_data = []
			for i in range(len(data)):
				fields = data[i].split(' ')
				cleaned = file_handlers.clean(fields)
				while cleaned.count('') > 0:
					cleaned.remove('')
				cleaned_data.append(cleaned)
			#for line in cleaned_data:
			gene_start = ''
			pdb_start = ''
			while gene_start == '':
				for line in cleaned_data:
					if ( len(line) == 4 and line[0] == 'Query' ):
						gene_start = int(line[1])
						break
					else:
						pass
			while pdb_start == '':
				for line in cleaned_data:
					if ( len(line) == 4 and line[0] == 'Sbjct' ):
						pdb_start = int(line[1])
						break
					else:
						pass
			seq_name = result_file.split('_')[0] + '_' + result_file.split('_')[3] + '_' + result_file.split('_')[4]
			start_site_difference = gene_start - pdb_start
			aln_diff.append((seq_name, start_site_difference))
		#print aln_diff
		return aln_diff
			

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
