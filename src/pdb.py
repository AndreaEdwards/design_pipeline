import os
import urllib2
from Bio.PDB import PDBList
from util import FileHandlers
from prody import fetchPDBviaFTP

class PDBFromUniprot:
	def __init__(self):
		pass

	def fetch_pdb(self, pdb_code, pdb_dir):
		pdb_list = PDBList()
		pdb_list.retrieve_pdb_file(pdb_code, pdir=pdb_dir)
	
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
	
	def rename_pdb_file(self, pdb_file):
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
	def __init__(self):
		pass

	def find_ligand(pdb_file):
		protein = parsePDB(pdb_file)
		try:
			seq = protein['A'].getSequence()
		except:
			pass
		else:
			ligand = protein.select('not protein and not water')
			repr(ligand)
			if ligand:
				writePDB(pdb_file.split('.')[0] + '_ligand.pdb', ligand)

	def get_ligand_centroid(ligand_pdb):
		for atom in residue:
		    min_array = [atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2]]
		    atom_pos_array.append(min_array)
		centroid_position = mean(atom_pos_array,axis = 0)


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
	def __init__(self):
		pass

	def get_surface_residues():
		pass


class EditPDB:
	def __init__(self, pdb_code):
		self.filename = pdb_code
		self.data = ''

	def _open_file(self):
		#os.chdir("../src/database/pdbs")
		file_handlers = FileHandlers()
		file_paths = file_handlers.search_directory()
		pdb_files = file_handlers.find_files(file_paths, 'pdb')
		for pdb_file in pdb_files:
			if self.filename == file_handlers.get_file_name(pdb_file).split('.')[0]:
				Data = open(pdb_file)
				self.data = Data.readlines()
				Data.close

	def b_factor(self):
		self._open_file()
		for line in self.data:
			if line[0:6] == "HETATM":
				print line[60:66]


#	import sys
#
#	def loadDataFile(data_file,col1=0,col2=1,abs_value=False):
#	    """
#	    Parses input from a two column data file, creating a dictionary of
#	    column[0]:column[1] key/value pairs.  Returns dictionary.
#	    """
#				    
#		# Read in file
#		f = open(data_file,'r')
#		data = f.readlines()
#		f.close()
#
#		# Strip out blank lines and comments
#		data = [d for d in data if d[0] != "#" and d.strip() != ""]
#
#		data_dict = {}
#		for record in data:
#		    try:
#		        field = record.split()
#		        key = field[col1]
#		        data_dict.update([(key,float(field[col2]))])
#
#		        if abs_value:
#		            data_dict[key] = abs(data_dict[key])
#				            
#				    except IndexError:
#				        sys.stderr.write("Mangled data, skipping line:\n%s" % record)
#				        continue
#				    except ValueError:
#				        sys.stderr.write("Mangled data, skipping line:\n%s" % record)
#				        continue
#				        
#				return data_dict
#
#
#	def pdbBfactor(pdb,data_dict):
#	    """
#	    Goes through pdb line by line.  If residue is in dictionary data_dict,
#	    the b-factor is replaced in output_file by the dictionary value.  If the
#	    residue is not in the dictionary, the b-factor is given value 0.0.
#	    Returns void.
#	    """
#
#		out = []
#		for line in pdb:
#		    if line[0:6] == "ATOM  ":
#		        resnum = line[23:26].strip()
#		        if resnum in data_dict.keys():
#		            out.append("%s%6.2F%s" % (line[:60],data_dict[resnum],
#		                                      line[66:]))
#		        else:
#		            out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
#		    
#		    elif line[0:6] == "HETATM":
#		        out.append("%s%6s%s" % (line[:60],"NA",line[66:]))
#		    
#		    else:
#		        out.append(line)
#		
#		return out












