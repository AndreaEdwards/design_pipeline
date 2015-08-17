import os
import urllib2
import subprocess
import shutil
from Bio.PDB import PDBList, PDBParser
from util import FileHandlers
from prody import fetchPDBviaFTP, parsePDB, writePDB
import numpy
import pandas

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
		self.filename = pdb_code
		self.file_path = ''
		self.out_filename = ''
		self.ligand_centroid = []
		self.ligands = []
		self.ligand_binding_pocket =[]
		self.ignore = ['MN ', 'DOD']

	def _get_file_path(self):
		os.chdir("./database/pdbs/pdb")
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		pdb_files = file_handlers.find_files(file_paths, 'pdb')
		for pdb_file in pdb_files:
			if self.filename == file_handlers.get_file_name(pdb_file).split('.')[0]:
				self.file_path = pdb_file
		os.chdir("../../../")

	def _find_ligand(self):
		self._get_file_path()
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

	def _get_ligand_centroid(self):
		self._find_ligand()
		p = PDBParser(QUIET=True)
		structure = p.get_structure('ligand', self.out_filename)
		chain = structure[0]['A']
		for residue in chain.get_residues():
		    if residue.resname not in self.ignore:
		    	atom_pos_array = []
    			min_array = []
    			for atom in residue:
		    		min_array = [atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2]]
		    		atom_pos_array.append(min_array)
				self.ligand_centroid.append(numpy.mean(atom_pos_array,axis = 0))

	def _get_ligand_name(self):
		p = PDBParser(QUIET=True)
		ligand = p.get_structure('ligand', self.out_filename)
		chain = ligand[0]['A']
		for residue in chain.get_residues():
			if residue.resname in self.ignore:
				pass
			else:
				self.ligands.append(residue.resname)
		print "Ligands found: ", self.ligands

	def get_residues_within_5A(self):
		self._get_ligand_centroid()
		self._get_ligand_name()
		p = PDBParser(QUIET=True)
		structure = p.get_structure('protein', self.file_path)
		chain = structure[0]['A']
		ligand_binding_pocket = []
		for residue in chain.get_residues():
			for atom in residue:
				if residue.resname in self.ligands:
					pass
				elif residue.resname in self.ignore:
					pass
				else:
					for centroid in self.ligand_centroid:
						if abs(float(atom.get_coord()[0]) - float(centroid[0])) <= 5 and abs(float(atom.get_coord()[1]) - float(centroid[1])) <= 5 and abs(float(atom.get_coord()[2]) - float(centroid[2])) <= 5:
							ligand_binding_pocket.append(str(residue.resname) + '\t' + str(residue.id[1]))
		self.ligand_binding_pocket = ligand_binding_pocket

	def write_residue_output(self):
		os.chdir('./database/pdbs/pdb')
		active_site_filename = os.getcwd() + '/' + self.filename + '_lpocket.txt'
		output_residues = open(active_site_filename, 'w')
		for line in self.ligand_binding_pocket:
			output_residues.write(line.split('\t')[1] + '\t' + str(100.00) + '\n')
		output_residues.close()
		os.chdir('../../../')


class Rosetta:
	def __init__(self):
		pass

	def minimize_pdb(pdb_file):
		cmd = ['~/rosetta/main/source/bin/minimize.static.macosclangrelease -s ' + pdb_file + ' -ignore_unrecognized_res']
		subprocess.call(cmd, shell=True)

	def run_ddg_monomer(pdb, mutation_list):
		# Not tested yet.
		cmd = ['~/rosetta/main/source/bin/pmut_scan_parallel.static.macosclangrelease  -s 4OY6_0001.pdb -ex1 -ex2 -extrachi_cutoff 1 -use_input_sc -ignore_unrecognized_res -no_his_his_pairE -multi_cool_annealer 10 -mute basic core -mutants_list mutant_list -DDG_cutoff 0 | grep PointMut|grep -v "go()"|grep -v "main()" > mutants.out']
		subprocess.call(cmd, shell=True)

	def pocket_finder():
		pass


class SurfaceResidues:
	def __init__(self, pdb_code):
		self.dir_path = os.getcwd() + '/pops_results'
		self._mkdir()
		self.filename = pdb_code
		self.file_path = ''
		self.out_filename = ''
		self.out_file_path = ''
		self.data = []
		self.SASA_dict = {}

	def _get_file_path(self):
		os.chdir("./database/pdbs/pdb")
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		pdb_files = file_handlers.find_files(file_paths, 'pdb')
		for pdb_file in pdb_files:
			if self.filename == file_handlers.get_file_name(pdb_file).split('.')[0]:
				self.file_path = pdb_file
				self.out_file = file_handlers.get_file_name(pdb_file).split('.')[0] + '_pops.out'
				self.out_file_path = self.dir_path + '/' + self.out_file

	def _mkdir(self):
		file_handlers = FileHandlers()
		file_handlers.make_results_folder(self.dir_path.split('/')[-1])

	def _run_POPS(self):
		self._get_file_path()
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
			output_b_factors.write(residue_position + '\t' + self.SASA_dict[self.filename][residue_position][2] + '\n')
		output_b_factors.close()

	def write_surface_resi_output(self, threshold):
		surface_residues = self._get_surface_residues(threshold)
		out_filename = self.filename + '_SurfRes.txt'
		output = open(out_filename, 'w')
		for resi_position in surface_residues:
			output.write(resi_position + '\t100.00\n')
		output.close

	def _get_surface_residues(self, threshold):
		surface_residues = []
		if self.SASA_dict == {}:
			self._build_SASA_dict()
			for pdb_code in self.SASA_dict:
				for resi_position, SASA_info_list in self.SASA_dict[pdb_code].iteritems():
					if float(SASA_info_list[3]) > threshold:
						surface_residues.append(resi_position)
		else:
			for pdb_code in self.SASA_dict:
				for resi_position, SASA_info_list in self.SASA_dict[pdb_code].iteritems():
					if float(SASA_info_list[3]) > threshold:
						surface_residues.append(resi_position)
		return surface_residues


class EditPDB:
	def __init__(self, pdb_code):
		self.filename = pdb_code
		self.pdb = []
		self.data = []
		self.data_dict = {}

	def _get_pdb(self):
		#os.chdir("./database/pdbs/pdb")
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		pdb_files = file_handlers.find_files(file_paths, 'pdb')
		for pdb_file in pdb_files:
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
			self.data_dict[cleaned[0]] = float(cleaned[1])

	def _edit_bfactor(self, file_tag):
		self._get_pdb()
		self._build_data_dict(file_tag)
		lines_out = []
		for line in self.pdb:
			if line[0:6] == "ATOM  ":
				resnum = line[23:26].strip()
				if resnum in self.data_dict.keys():
					lines_out.append("%s%6.2F%s" % (line[:60], self.data_dict[resnum], line[66:]))
				else:
					lines_out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
			elif line[0:6] == "HETATM":
				lines_out.append("%s%6s%s" % (line[:60],"NA",line[66:]))   
			else:
				lines_out.append(line)
		return lines_out

	def _write_output(self, lines, file_tag):
		#os.chdir('./database/pdbs/pdb') 
		output_file = os.getcwd() + "/" + self.filename + file_tag
		outfile = open(output_file, 'w')
		lines = self._edit_bfactor(file_tag)
		for line in lines:
			outfile.write(line)
		outfile.close()

	def edit_bfactor_sasa(self):
		lines = self._edit_bfactor("_sasa")
		self._write_output(lines, "_sasa.pdb")

	def edit_bfactor_ligand_binding_pocket(self):
		lines = self._edit_bfactor("_lpocket")
		self._write_output(lines, "_lpocket.pdb")

	def edit_bfactor_surface_residues(self):
		lines = self._edit_bfactor("_SurfRes")
		self._write_output(lines, "_SurfRes.pdb")












