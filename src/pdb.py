import os
import urllib2
import subprocess
import shutil
from Bio.PDB import PDBList, PDBParser
from util import FileHandlers
from prody import fetchPDBviaFTP, parsePDB, writePDB
import numpy
import pandas
from operator import itemgetter
import matplotlib.pyplot as plt
import itertools
import amino_acids

class PDBFromUniprot:
	def __init__(self):
		self.pdb_dir = os.getcwd() + '/database/pdbs/pdb'

	def fetch_pdb(self, pdb_code):
		pdb_list = PDBList()
		pdb_list.retrieve_pdb_file(pdb_code, pdir=self.pdb_dir)
		self._rename_pdb_file(pdb_code)
	
	def get_pdb_id(self, queryText):
		url = 'http://www.rcsb.org/pdb/rest/search'
		print "querying PDB..."
		req = urllib2.Request(url, data=queryText)
		f = urllib2.urlopen(req)
		result = f.readlines()
		if result:
			PDB_codes = []
			print "Found number of PDB entries: %d" % len(result)
			for pdb_code in result:
				PDB_codes.append(pdb_code.strip())
			return PDB_codes
		else:
			print "Failed to retrieve results"
			PDB_codes = []
			return PDB_codes

	def _get_downloaded_file_path(self, pdb_code):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		ent_files = file_handlers.find_files(file_paths, 'ent')
		for ent_file in ent_files:
			if pdb_code == file_handlers.get_file_name(ent_file).split('.')[0].lstrip('pdb').upper():
				return ent_file
	
	def _rename_pdb_file(self, pdb_code):
		pdb_file = self._get_downloaded_file_path(pdb_code)
		split_path = pdb_file.split('/')
		name = split_path[-1].split('.')
		new_name = name[0].lstrip('pdb').upper() + '.pdb'
		split_path[-1] = new_name
		new_path = '/'.join(split_path)
		shutil.copyfile(pdb_file, new_path)
		os.remove(pdb_file)	


class CIFFromUniprot:
	def __init__(self):
		pass

	def fetch_mmCIF(self, pdb_code, pdb_dir):
		fetchPDBviaFTP(pdb_code, format='cif', folder=pdb_dir)

	def unzip_mmCIF(self, gz_file):
		cmd = ['gunzip -d ' + gz_file]
		subprocess.call(cmd, shell=True)

	def rename_mmCIF_file(self, cif_file):
		split_path = cif_file.split('/')
		name = split_path[-1].split('.')
		new_name = name[0].upper() + '.cif'
		split_path[-1] = new_name
		new_path = '/'.join(split_path)
		shutil.copyfile(cif_file, new_path)
		os.remove(cif_file)	


class LigandBindingSite:
	def __init__(self, pdb_code):
		print "Searching for ligand binding sites...."
		self.filename = pdb_code
		self.file_path = ''
		self.out_filename = ''
		self.ligand_centroid = []
		self.ligands = []
		self.ligand_binding_pocket =[]
		self.ignore = ['MN ', 'DOD', 'SO4', 'MES']
		self.ligand_chain = []
		self.all_ligand_atoms = []

	def _get_file_path(self, ligand=False, pdb=False):
		#os.chdir("./database/pdbs/pdb")
		self.file_path = ''
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		pdb_files = file_handlers.find_files(file_paths, 'pdb')
		if ligand == True:
			for pdb_file in pdb_files:
				if self.filename == file_handlers.get_file_name(pdb_file).split('.')[0]:
					self.file_path = pdb_file
		elif pdb == True:
			file_names = []
			for pdb_file in pdb_files:
				found = file_handlers.get_file_name(pdb_file).split('.')[0]
				file_names.append(found)
				if self.filename == file_handlers.get_file_name(pdb_file).split('.')[0]:
					current_pdb = pdb_file
					#print "current_pdb", current_pdb
			if (self.filename + '_0001') not in file_names:
				# Run Rosetta Minimizer
				cmd = ['~/rosetta/main/source/bin/minimize.static.macosclangrelease -jd2 -s ' + current_pdb + ' -ignore_unrecognized_res']
				subprocess.call(cmd, shell=True)
				file_handlers = FileHandlers()
				file_paths = file_handlers.search_directory()
				pdb_files = file_handlers.find_files(file_paths, 'pdb')
				for pdb_file in pdb_files:
					#print "pdb_file:", pdb_file
					if (self.filename + '_0001') == file_handlers.get_file_name(pdb_file).split('.')[0]:
						self.file_path = pdb_file
			else:
				print "Found Rosetta minimized structure file %s" % (self.filename + '_0001')
				for pdb_file in pdb_files:
					if (self.filename + '_0001') == file_handlers.get_file_name(pdb_file).split('.')[0]:
						self.file_path = pdb_file

		#os.chdir("../../../")

	def _find_ligand(self):
		self._get_file_path(ligand=True)
		protein = parsePDB(self.file_path)
		try:
			seq = protein['A'].getSequence()
		except:
			pass
		else:
			ligand = protein.select('not protein and not water')
			repr(ligand)
			if ligand:
				self.out_filename = self.file_path.split('.')[0] + '_ligand.pdb'
				writePDB(self.out_filename, ligand)

	def _get_ligand_coordinates(self):
		self._find_ligand()
		p = PDBParser(QUIET=True)
		if self.out_filename == '':
			print "No ligand found for %s" % self.filename
		else:
			structure = p.get_structure('ligand', self.out_filename)
			chain_ids = ['A', 'B', 'C', 'D', 'X']
			for entry in chain_ids:
				try:
					chain = structure[0][entry]
					if chain:
						print "Found ligand on chain %s" % entry
						for residue in chain.get_residues():
							if residue.resname not in self.ignore:
								atom_pos_array = []
								min_array = []
								for atom in residue:
									min_array = [atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2]]
									atom_pos_array.append(min_array)
								self.ligand_centroid.append(numpy.mean(atom_pos_array,axis = 0))
								self.all_ligand_atoms = atom_pos_array
				except KeyError:
					pass
					#print 'Cannot calculate ligand centroid for chain %s of %s' % (entry, self.filename)


	def _get_ligand_name(self):
		p = PDBParser(QUIET=True)
		ligand = p.get_structure('ligand', self.out_filename)
		chain_ids = ['A', 'B', 'C', 'D', 'X']
		for entry in chain_ids:
			try:
				chain = ligand[0][entry]
				for residue in chain.get_residues():
					if residue.resname in self.ignore:
						pass
					else:
						self.ligands.append(residue.resname)
						self.ligand_chain.append(entry)
				print "Ligands found: ", self.ligands
				return self.ligands
			except KeyError:
				pass
				#print 'No ligand found in chain %s of %s' % (entry, self.filename)

	def get_residues_within_5A(self, centroid='False'):
		self._get_ligand_coordinates()
		if self.ligand_centroid == []:
			print 'No ligand bound to %s' % self.filename
		else:
			self._get_ligand_name()
			self._get_file_path(pdb=True)
			p = PDBParser(QUIET=True)
			structure = p.get_structure('protein', self.file_path)
			#chain_ids = ['A', 'B', 'C', 'D', 'X']
			for entry in self.ligand_chain:
				try:
					chain = structure[0][entry]
					ligand_binding_pocket = []
					for residue in chain.get_residues():
						for atom in residue:
							if residue.resname in self.ligands or (' ' + residue.resname.rstrip()) in self.ligands or (residue.resname.lstrip() + ' ') in self.ligands:
								pass
							elif residue.resname in self.ignore or (' ' + residue.resname.rstrip()) in self.ignore or (residue.resname.lstrip() + ' ') in self.ignore:
								pass
							else:
								if centroid == 'True':
									for centroid in self.ligand_centroid:
										if abs(float(atom.get_coord()[0]) - float(centroid[0])) <= 5 and abs(float(atom.get_coord()[1]) - float(centroid[1])) <= 5 and abs(float(atom.get_coord()[2]) - float(centroid[2])) <= 5:
											ligand_binding_pocket.append(str(residue.resname) + '\t' + str(residue.id[1]))
								else:
									for coordinate in self.all_ligand_atoms:
										if abs(float(atom.get_coord()[0]) - float(coordinate[0])) <= 5 and abs(float(atom.get_coord()[1]) - float(coordinate[1])) <= 5 and abs(float(atom.get_coord()[2]) - float(coordinate[2])) <= 5:
											ligand_binding_pocket.append(str(residue.resname) + '\t' + str(residue.id[1]))
					self.ligand_binding_pocket = set(ligand_binding_pocket)
					print "Residues in ligand binding pocket: ", [line.split('\t') for line in self.ligand_binding_pocket]
					return self.ligand_binding_pocket
				except KeyError:
					print 'No ligand found in chain %s of %s' % (entry, self.filename)


	def write_residue_output(self):
		#os.chdir('./database/pdbs/pdb')
		active_site_filename = os.getcwd() + '/' + self.filename + '_lpocket.txt'
		output_residues = open(active_site_filename, 'w')
		for line in self.ligand_binding_pocket:
			output_residues.write(line.split('\t')[1] + '\t' + str(100.00) + '\n')
		output_residues.close()
		#os.chdir('../../../')


class Rosetta:
	def __init__(self, pdb_code):
		self.filename = pdb_code
		self.pdb = []
		self.standard_res = [
								'GLY', 'PRO', 'ALA', 'VAL', 'LEU', 
								'ILE', 'MET', 'CYS', 'PHE', 'TYR',
								'TRP', 'HIS', 'LYS', 'ARG', 'GLN',
								'ASN', 'GLU', 'ASP', 'SER', 'THR'
							]


	def _get_pdb(self, rosetta_min=False, refined_pocket=False):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		pdb_files = file_handlers.find_files(file_paths, 'pdb')
		for pdb_file in pdb_files:
			if rosetta_min == True and refined_pocket == True:
				print "Invalid input"
			elif rosetta_min == True:
				if (self.filename + '_0001') == file_handlers.get_file_name(pdb_file).split('.')[0]:
					print "Found Rosetta minimized structure file ", (self.filename + '_0001.pdb')
					filepath = pdb_file
			elif refined_pocket == True:
				if ('pocket0') == file_handlers.get_file_name(pdb_file).split('.')[0]:
					print "Found pocket file pocket0.pdb"
					filepath = pdb_file
			else:
				if self.filename == file_handlers.get_file_name(pdb_file).split('.')[0]:
					print "Found target structure file ", (self.filename + '.pdb')
					filepath = pdb_file
		return filepath

	def _open_file(self, filepath):
		Data = open(filepath, 'r')
		data = Data.readlines()
		Data.close()
		os.remove(filepath)
		return data

	def _get_data(self, file_tag):
		data = []
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		txt_files = file_handlers.find_files(file_paths, 'txt')
		for txt_file in txt_files:
			if self.filename == file_handlers.get_file_name(txt_file).split('_')[0]:
				if (file_tag.split('_')[1] + '.txt') == file_handlers.get_file_name(txt_file).split('_')[1]:
					TXT = open(txt_file)
					data = TXT.readlines()
					TXT.close()
		return data

	def _get_surface_residues(self, file_tag):
		data = self._get_data(file_tag)
		residues = []
		for line in data:
			residue_number = line.split('\t')[0].strip()
			residues.append(residue_number)
		return residues

	def _get_output(self):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		out_files = file_handlers.find_files(file_paths, 'out')
		for out_file in out_files:
			if 'pocket_score' == file_handlers.get_file_name(out_file).split('.')[0]:
				TXT = open(out_file)
				data = TXT.readlines()
				TXT.close()
				#os.remove(out_file)
		return data

	def _get_score(self):
		data = self._get_output()
		for line in data:
			score = float(line.split(' ')[-1].strip())
		return score

	def _get_resnums(self, filtered_data):
		resnums = []
		for data in filtered_data:
			resnums.append(data[0])
		return resnums

	def _write_scores(self, data, handle):
		out_file = os.getcwd() + '/' + self.filename + handle
		out = open(out_file, 'w')
		for tuple_element in data:
			out.write(str(tuple_element[0]) + '\t' + str(tuple_element[1]) + '\n')
		out.close()

	def _write_pocket_residues(self, handle):
		pocket_residues = self._get_pocket_residues()
		out_file = os.getcwd() + '/' + self.filename + handle
		out = open(out_file, 'w')
		for resnum in pocket_residues:
			out.write(str(resnum) + '\t' + str(pocket_residues[resnum]) + '\n')
		out.close()

	def _plot_data(self, data):
		x = []
		y = []
		for d in data:
			x.append(d[0])
			y.append(d[1])
		plt.scatter(x, y)
		plt.show()

	def _sort_scores(self, data, plt=False):
		sorted_data = sorted(data, key=itemgetter(1))
		if plt == True:
			self._plot_data(sorted_data)
		return sorted_data

	def _pocket_finder(self, filepath, residue_number, output=False):
		if len(residue_number) == 1:
			chain = self._get_chain_id(filepath, residue_number[0])
			cmd = ['~/rosetta/main/source/bin/pocket_measure.static.macosclangrelease -s ' + filepath + ' -central_relax_pdb_num ' + str(residue_number[0]) + ':' + chain + ' -pocket_num_angles 100 | grep Largest >> pocket_score.out']
			subprocess.call(cmd, shell=True)
		elif len(residue_number) == 2 and output == False:
			chain1 = self._get_chain_id(filepath, residue_number[0])
			chain2 = self._get_chain_id(filepath, residue_number[1])
			cmd = ['~/rosetta/main/source/bin/pocket_measure.static.macosclangrelease -s ' + filepath + ' -central_relax_pdb_num ' + str(residue_number[0]) + ':' + chain1 + ',' + str(residue_number[1]) + ':' + chain2 + ' -pocket_num_angles 100 | grep Largest >> pocket_score.out']
			subprocess.call(cmd, shell=True)
		elif len(residue_number) == 2 and output == True:
			chain1 = self._get_chain_id(filepath, residue_number[0])
			chain2 = self._get_chain_id(filepath, residue_number[1])
			cmd = ['~/rosetta/main/source/bin/pocket_measure.static.macosclangrelease -s ' + filepath + ' -central_relax_pdb_num ' + str(residue_number[0]) + ':' + chain1 + ',' + str(residue_number[1]) + ':' + chain2 + ' -pocket_dump_pdbs -pocket_num_angles 100 | grep Largest >> pocket_score.out']
			subprocess.call(cmd, shell=True)
		else:
			print "invalid input"

	def _get_chain_id(self, filepath, resnum):
		p = PDBParser(QUIET=True)
		structure = p.get_structure('protein', filepath)
		for chain in structure.get_chains():
			chain_id = chain.id
			chain = structure[0][chain_id]
			for residue in chain.get_residues():
				if str(residue.id[1]) == str(resnum):
					return chain_id

	def _is_standard_res(self, filepath, resnum):
		p = PDBParser(QUIET=True)
		structure = p.get_structure('protein', filepath)
		for chain in structure.get_chains():
			chain_id = chain.id
			chain = structure[0][chain_id]
			for residue in chain.get_residues():
				if str(residue.id[1]) == resnum:
					if residue.resname in self.standard_res:
						return True
					else:
						return False

	def _scan_surface(self):
		data = []
		filepath = self._get_pdb(rosetta_min=True)
		residues = self._get_surface_residues("_SurfRes")
		print 'Scanning surface for pockets...'
		for residue in residues:
			if self._is_standard_res(filepath, residue):
				self._pocket_finder(filepath, [residue])
				score = self._get_score()
				print 'Residue number: %s Score: %s' % (residue, str(score))
				data.append((residue, score))
		self._write_scores(data, '_pocket_scores.out')
		return data

	def _refine_pockets(self, filtered_data):
		print "Refining pocket location..."
		filepath = self._get_pdb(rosetta_min=True)
		resnums = self._get_resnums(filtered_data)
		combos = itertools.combinations(resnums, 2)
		results = []
		for combo in combos:
			print "combo is:", combo
			res1, res2 = combo[0], combo[1]
			self._pocket_finder(filepath, combo)
			score = self._get_score()
			results.append((combo, score))
		sorted_results = self._sort_scores(results, plt=False)
		self._pocket_finder(filepath, list(sorted_results[-1][0]), output=True)

	def _get_pocket_coordinates(self):
		pockets_filepath = self._get_pdb(refined_pocket=True)
		pockets = self._open_file(pockets_filepath)
		pocket_coordinates = []
		pocket_data = []
		for line in pockets:
			if line[13:16] == 'TPB':
				pocket_data.append(line)
				coordinates = (line[31:38], line[39:46], line[47:54])
				pocket_coordinates.append(coordinates)
		self._extract_pocket(pocket_data)
		return pocket_coordinates

	def _extract_pocket(self, pocket_data):
		out_file = os.getcwd() + '/' + self.filename + 'main_pocket.pdb'
		out = open(out_file, 'w')
		for line in pocket_data:
			out.write(line + '\n')
		out.close() 

	def _get_pocket_residues(self):
		pocket_residues = {}
		pocket_coordinates = self._get_pocket_coordinates()
		protein_filepath = self._get_pdb(rosetta_min=True)
		p = PDBParser(QUIET=True)
		structure = p.get_structure('protein', protein_filepath)
		for chain in structure.get_chains():
			chain = structure[0][chain.id]
			for residue in chain.get_residues():
				if residue.id[0] == ' ':
					for atom in residue:
						for coordinate in pocket_coordinates:
							if abs(float(atom.get_coord()[0]) - float(coordinate[0])) <= 2.5 and abs(float(atom.get_coord()[1]) - float(coordinate[1])) <= 2.5 and abs(float(atom.get_coord()[2]) - float(coordinate[2])) <= 2.5:
								pocket_residues[str(residue.id[1])] = str(100.00)
				else:
					pass
		return pocket_residues

	def find_pockets(self):
		data = self._scan_surface()
		sorted_data = self._sort_scores(data, plt=False)
		highest_score = float(sorted_data[-1][1])
		lowest_score = float(sorted_data[0][1])
		filtered_data = []
		for data in sorted_data:
			frac = abs((lowest_score - float(data[1])) / (highest_score - lowest_score))
			if frac >= 0.50:
				filtered_data.append(data)
		self._refine_pockets(filtered_data)
		self._write_scores(filtered_data, '_pockets.txt')
		self._write_pocket_residues('_pocketres.txt')


class MutantListMaker:
	def __init__(self, pdb_code):
		self.handles = []
		self.filename = pdb_code
		self.resnums = set([])
		self.codes = {	'GLY' : 'G', 'PRO' : 'P', 'ALA' : 'A', 'VAL' : 'V', 'LEU' : 'L', 
						'ILE' : 'I', 'MET' : 'M', 'CYS' : 'C', 'PHE' : 'F', 'TYR' : 'Y',
						'TRP' : 'W', 'HIS' : 'H', 'LYS' : 'K', 'ARG' : 'R', 'GLN' : 'Q',
						'ASN' : 'N', 'GLU' : 'E', 'ASP' : 'D', 'SER' : 'S', 'THR' : 'T'
					}
		self.inverse_codes = {}
		for key in self.codes:
			self.inverse_codes[self.codes[key]] = key
		self.amino_acids = set([])
		for key in self.codes:
			self.amino_acids.add(self.codes[key]) 

	def _open_file(self, filepath):
		Data = open(filepath, 'r')
		data = Data.readlines()
		Data.close()
		return data

	def _get_filepath(self, handle, data_file=False, pdb_file=False):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		if data_file == True:
			files = file_handlers.find_files(file_paths, 'txt')
			for path in files:
				if (self.filename + handle) == file_handlers.get_file_name(path).split('.')[0]:
					return path
		elif pdb_file == True:
			files = file_handlers.find_files(file_paths, 'pdb')
			for path in files:
				if (self.filename + '_0001') == file_handlers.get_file_name(path).split('.')[0]:
					return path
		else:
			print "Specify file type"

	def _get_resnums(self):
		resnums = set([])
		for handle in self.handles:
			filepath = self._get_filepath(handle, data_file=True)
			data = self._open_file(filepath)
			for line in data:
				resnum = line.split('\t')[0]
				resnums.add(resnum)
		self.resnums = resnums
		return resnums

	def _get_resmapping(self):
		res_mapping = []
		filepath = self._get_filepath('', pdb_file=True)
		p = PDBParser(QUIET=True)
		structure = p.get_structure('protein', filepath)
		chain = structure[0]['A']
		for residue in chain.get_residues():
			if str(residue.id[1]) in self.resnums:
				res_mapping.append((self.codes[residue.resname], residue.id[1]))
		return res_mapping

	def _get_mutants(self, res_mapping):
		mutants = {}
		for mapping in res_mapping:
			resid = set(mapping[0])
			aminoacids = set(self.inverse_codes.keys())
			mutant_list = list(aminoacids.difference(resid))
			mutants[mapping] = mutant_list
		return mutants

	def _write_mutant_list(self, mutant_dict, handle):
		out_file = os.getcwd() + '/' + self.filename + handle
		out = open(out_file, 'w')
		for key_tuple in mutant_dict:
			for aa in mutant_dict[key_tuple]:
				out.write('A ' + str(key_tuple[0]) + ' ' + str(key_tuple[1]) + ' ' + aa + '\n')
		out.close()

	def generate_mutant_list(self, pocketres=True, lpocket=True, SurfRes=True):
		if pocketres == True and lpocket == True and SurfRes == True:
			#self.handles = ['_pocketres', '_lpocket', '_SurfRes']
			self.handles = ['_pocketres', '_lpocket']
			self._get_resnums()
		res_mapping = self._get_resmapping()
		mutants = self._get_mutants(res_mapping)
		self._write_mutant_list(mutants, '_mutant_list.txt')


class DDGMonomer:
	def __init__(self, pdb):
		self.filename = pdb

	def _get_filepath(self, data_file=False, pdb_file=False):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		if data_file == True:
			files = file_handlers.find_files(file_paths, 'txt')
			for path in files:
				if (self.filename + '_mutant_list') == file_handlers.get_file_name(path).split('.')[0]:
					return path
		elif pdb_file == True:
			files = file_handlers.find_files(file_paths, 'pdb')
			for path in files:
				if (self.filename + '_0001') == file_handlers.get_file_name(path).split('.')[0]:
					return path
		else:
			print "Specify file type"

	def _run_ddg_monomer(self, pdb_filepath, mutant_list, threshold):
		print "Running ddg_monomer..."
		cmd = ['~/rosetta/main/source/bin/pmut_scan_parallel.static.macosclangrelease -jd2 -s ' + pdb_filepath + ' -ex1 -ex2 -extrachi_cutoff 1 -use_input_sc -ignore_unrecognized_res -no_his_his_pairE -multi_cool_annealer 10 -mute basic core -mutants_list ' + mutant_list + ' -DDG_cutoff ' + threshold + ' | grep PointMut|grep -v "go()"|grep -v "main()" > ' + self.filename + '_mutants.out']
		subprocess.call(cmd, shell=True)

	def get_targets(self, threshold):
		mutant_list = self._get_filepath(data_file=True)
		pdb_filepath = self._get_filepath(pdb_file=True)
		self._run_ddg_monomer(pdb_filepath, mutant_list, str(threshold))


class SurfaceResidues:
	def __init__(self, pdb_code):
		print "Calculating solvent accessible surface area using POPS...."
		self.dir_path = os.getcwd() + '/pops_results'
		self._mkdir()
		self.filename = pdb_code
		self.file_path = ''
		self.out_filename = ''
		self.out_file_path = ''
		self.data = []
		self.SASA_dict = {}
		self.found_modres = False

	def _get_file_path(self, modres=False):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		pdb_files = file_handlers.find_files(file_paths, 'pdb')
		if modres == True:
			for pdb_file in pdb_files:
				if (self.filename + '_0001_edit_modres') == file_handlers.get_file_name(pdb_file).split('.')[0]:
					self.file_path = pdb_file
					self.out_file = file_handlers.get_file_name(pdb_file).split('.')[0] + '_pops.out'
					self.out_file_path = self.dir_path + '/' + self.out_file
		else:
			for pdb_file in pdb_files:
				if (self.filename + '_0001') == file_handlers.get_file_name(pdb_file).split('.')[0]:
					self.file_path = pdb_file
					self.out_file = file_handlers.get_file_name(pdb_file).split('.')[0] + '_pops.out'
					self.out_file_path = self.dir_path + '/' + self.out_file

	def _mkdir(self):
		file_handlers = FileHandlers()
		file_handlers.make_results_folder(self.dir_path.split('/')[-1])

	def _check_for_modres(self, protein_filepath):
		p = PDBParser(QUIET=True)
		structure = p.get_structure('protein', protein_filepath)
		modres_list = []
		swapped_list = []
		for chain in structure.get_chains():
			chain = structure[0][chain.id]
			for residue in chain.get_residues():
				if residue.resname not in amino_acids.longer_names.keys():
					if residue.resname in amino_acids.modres.keys():
						self.found_modres = True
						modres_list.append((str(residue.id[1]), residue.resname))
						swapped_list.append((str(residue.id[1]), amino_acids.modres[residue.resname]))
		if self.found_modres:
			print 'Modified residues found in structures:', modres_list
			print 'Swapping modified residues with canonical residues in PDB for SASA calculation.'
			print 'Swapped residues are: ', swapped_list
			pdb_editor = EditPDB(self.filename)
			pdb_editor.edit_resname(modres_list, '_0001_edit_modres.pdb')

	def _run_POPS(self):
		self._get_file_path()
		self._check_for_modres(self.file_path)
		if self.found_modres:
			self.file_path = ''
			self._get_file_path(modres=True)
		cmd = ['pops --pdb ' + self.file_path + ' --coarse --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut --popsOut ' + self.out_file_path]
		with open(os.devnull, 'w') as fnull:
			subprocess.call(cmd, stdout=fnull, stderr=fnull, shell=True)

	def _get_data(self):
		Data = open(self.out_file_path)
		self.data = Data.readlines()
		Data.close()
		#os.remove(self.out_file_path)

	def _build_SASA_dict(self):
		file_handlers = FileHandlers()
		self.SASA_dict[self.filename] = {}
		self._run_POPS()
		self._get_data()
		for line in self.data:
			fields = line.split('\t')
			cleaned = file_handlers.clean(fields)
			if len(cleaned) == 9: 
				(position, 
				aa, 
				tot_SA, 
				SASA, 
				frac_SASA, 
				phob, 
				phil) = (cleaned[2],
						cleaned[0],
						cleaned[8],
						cleaned[5],
						cleaned[6],
						cleaned[3],
						cleaned[4])
				self.SASA_dict[self.filename][position] = [aa, 
													tot_SA, 
													SASA, 
													frac_SASA, 
													phob, 
													phil]
		#return self.SASA_dict

	def write_resi_sasa_output(self):
		self._build_SASA_dict()
		sasa_b_factor_filename = os.getcwd() + '/' + self.filename + '_sasa.txt'
		output_b_factors = open(sasa_b_factor_filename, 'w')
		for residue_position in self.SASA_dict[self.filename]:
			#print "residue id:", self.SASA_dict[self.filename][residue_position][0]
			output_b_factors.write(residue_position + '\t' + self.SASA_dict[self.filename][residue_position][2] + '\n')
		output_b_factors.close()

	def write_frac_sasa_output(self):
		self._build_SASA_dict()
		sasa_b_factor_filename = os.getcwd() + '/' + self.filename + '_fracsasa.txt'
		output_b_factors = open(sasa_b_factor_filename, 'w')
		for residue_position in self.SASA_dict[self.filename]:
			output_b_factors.write(residue_position + '\t' + self.SASA_dict[self.filename][residue_position][3] + '\n')
		output_b_factors.close()

	def write_surface_resi_output(self, threshold):
		surface_residues = self._get_surface_residues(threshold)
		out_filename = self.filename + '_SurfRes.txt'
		output = open(out_filename, 'w')
		for line in surface_residues:
			resi_position = line.split('\t')[1]
			output.write(resi_position + '\t100.00\n')
		output.close

	def _get_surface_residues(self, threshold):
		surface_residues = []
		if self.SASA_dict == {}:
			self._build_SASA_dict()
			for pdb_code in self.SASA_dict:
				for resi_position, SASA_info_list in self.SASA_dict[pdb_code].iteritems():
					if float(SASA_info_list[3]) > threshold:
						surface_residues.append(SASA_info_list[0] + '\t' + resi_position)
		else:
			for pdb_code in self.SASA_dict:
				for resi_position, SASA_info_list in self.SASA_dict[pdb_code].iteritems():
					if float(SASA_info_list[3]) > threshold:
						surface_residues.append(SASA_info_list[0] + '\t' + resi_position)
		print "Surface residues are: ", [line.split('\t') for line in surface_residues]
		return surface_residues


class EditPDB:
	def __init__(self, pdb_code):
		self.filename = pdb_code
		self.pdb = []
		self.data = []
		self.data_dict = {}

	def _get_pdb(self, modres=False):
		#os.chdir("./database/pdbs/pdb")
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		pdb_files = file_handlers.find_files(file_paths, 'pdb')
		for pdb_file in pdb_files:
			if modres == True:
				if (self.filename + '_0001') == file_handlers.get_file_name(pdb_file).split('.')[0]:
					PDB = open(pdb_file)
					self.pdb = PDB.readlines()
					PDB.close()
			else:
				if self.filename == file_handlers.get_file_name(pdb_file).split('.')[0]:
					PDB = open(pdb_file)
					self.pdb = PDB.readlines()
					PDB.close()

	def _get_data(self, file_tag):
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		txt_files = file_handlers.find_files(file_paths, 'txt')
		for txt_file in txt_files:
			if self.filename == file_handlers.get_file_name(txt_file).split('_')[0]:
				if (file_tag.split('_')[1] + '.txt') == file_handlers.get_file_name(txt_file).split('_')[1]:
					TXT = open(txt_file)
					self.data = TXT.readlines()
					TXT.close()

	def _build_data_dict(self, file_tag):
		self.data_dict = {}
		self._get_data(file_tag)
		file_handlers = FileHandlers()
		for line in self.data:
			fields = line.split('\t')
			cleaned = file_handlers.clean(fields)
			self.data_dict[int(cleaned[0])] = float(cleaned[1])

	def _edit_bfactor(self, file_tag):
		self._get_pdb()
		self._build_data_dict(file_tag)
		lines_out = []
		for line in self.pdb:
			if line[0:6] == "ATOM  ":
				resnum = int(line[22:26].strip())
				if resnum in self.data_dict.keys():
					lines_out.append("%s%6.2F%s" % (line[:60], self.data_dict[resnum], line[66:]))
				else:
					lines_out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
			elif line[0:6] == "HETATM":
				lines_out.append("%s%6s%s" % (line[:60],"NA",line[66:]))   
			else:
				lines_out.append(line)
		return lines_out

	def edit_resname(self, reslist, file_tag):
		self._get_pdb(modres=True)
		lines_out = []
		for line in self.pdb:
			if line[0:6] == "ATOM  ":
				for resnum in reslist[0]:
					res_num = int(line[22:26].strip())
					if res_num == resnum:
						lines_out.append("%s%s%s" % (line[:17], amino_acids.modres[reslist[1]], line[20:]))
					else:
						lines_out.append(line)
		self._write_output(lines_out, file_tag)

	def _write_output(self, lines, file_tag):
		#os.chdir('./database/pdbs/pdb') 
		output_file = os.getcwd() + "/" + self.filename + file_tag
		outfile = open(output_file, 'w')
		lines = self._edit_bfactor(file_tag)
		for line in lines:
			outfile.write(line)
		outfile.close()

	def _pull_bfactor(self):
		self._get_pdb()
		lines_out = {}
		for line in self.pdb:
			if line[0:6] == "ATOM  ":
				resnum = line[23:26].strip()
				b_factor = line[60:66].strip()
				lines_out[resnum] = resnum + '\t' + b_factor + '\n'
		return lines_out

	def _sort_keys(self, input_dict):
		temp_list = []
		for key in input_dict:
			temp_list.append(int(key))
		temp_list.sort()
		return temp_list

	def edit_bfactor_sasa(self):
		lines = self._edit_bfactor("_sasa")
		self._write_output(lines, "_sasa.pdb")

	def edit_bfactor_percent_sasa(self):
		lines = self._edit_bfactor("_fracsasa")
		self._write_output(lines, "_fracsasa.pdb")

	def edit_bfactor_ligand_binding_pocket(self):
		lines = self._edit_bfactor("_lpocket")
		self._write_output(lines, "_lpocket.pdb")

	def edit_bfactor_surface_residues(self):
		lines = self._edit_bfactor("_SurfRes")
		self._write_output(lines, "_SurfRes.pdb")

	def edit_bfactor_pockets(self):
		lines = self._edit_bfactor("_pockets")
		self._write_output(lines, "_pockets.pdb")

	def edit_bfactor_pocket_residues(self):
		lines = self._edit_bfactor("_pocketres")
		self._write_output(lines, "_pocketres.pdb")

	def get_bfactor(self):
		lines = self._pull_bfactor()
		sorted_list = self._sort_keys(lines)
		#os.chdir('./database/pdbs/pdb')
		output_file = os.getcwd() + "/" + self.filename + '_pulled_bfactors.txt'
		outfile = open(output_file, 'w')
		for resnum in sorted_list:
			outfile.write(lines[str(resnum)])
		outfile.close()
		#os.chdir('../../../')
















