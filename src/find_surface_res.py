"""
Parse the output files for POPS* run using the --coarse flag.
Output a dictionary: Position: [aa, tot_SA, SASA, Phob, Phil]

key questions:
1. What is the composition of amino acids on surface?
2. Draw a heat map with amino acids in rows and pdbs in columns. 
heat indicates frequency on surface. use frac_SASA as cutoff for definition as
'on surface' 

"""

from util import FileHandlers
import numpy
import matplotlib.pyplot as plt

def load_out_files():
	file_handlers = FileHandlers()
	file_paths = file_handlers.search_directory()
	out_files = file_handlers.find_files(file_paths, 'out')
	return out_files


def build_SASA_dict(out_files):
	SASA_dict = {}
	for path in out_files:
		file_handlers = FileHandlers()
		file_name = file_handlers.get_file_name(path)
		SASA_dict[file_name] = {}
		for line in open(path):
			file_handlers2 = FileHandlers()
			fields = line.split('\t')
			cleaned = file_handlers2.clean(fields)
			if len(cleaned) == 9: #and int(cleaned[2]) >= 1:
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
				SASA_dict[file_name][position] = [aa, 
													tot_SA, 
													SASA, 
													frac_SASA, 
													phob, 
													phil]
	return SASA_dict

def calculate_aa_comp(pdb_dict):
	aa_comp_dict = {}
	for resi_position, SASA_info_list in pdb_dict.iteritems():
		if float(SASA_info_list[3]) > 0.3:
			if SASA_info_list[0] in aa_comp_dict:
				aa_comp_dict[SASA_info_list[0]] += 1
			else:
				aa_comp_dict[SASA_info_list[0]] = 1
	return aa_comp_dict


def calculate_frac_surface_resi(pdb_dict, aa_comp_dict):
	total = len(pdb_dict)
	surface_aa = 0
	for aa in aa_comp_dict:
		surface_aa += aa_comp_dict[aa]
	frac_surface_resi = float(surface_aa) / total
	return frac_surface_resi 

def build_lists(SASA_dict):
	total_length_list = []
	frac_surface_resi_list = []
	for file_name, pdb_dict in SASA_dict.iteritems():
		total_length_list.append(len(pdb_dict))
		aa_comp_dict = calculate_aa_comp(pdb_dict)
		frac_surface_resi = calculate_frac_surface_resi(pdb_dict, aa_comp_dict)
		frac_surface_resi_list.append(frac_surface_resi)
	return total_length_list, frac_surface_resi_list


def plot_FracSurfaceResi_LengthProtein(total_length_list, frac_surface_resi_list):
	plt.scatter(total_length_list, frac_surface_resi_list)
	plt.draw()
	plt.savefig('output.png')

def draw_heat_map():
	pass

def main():
	out_files = load_out_files()
	SASA_dict = build_SASA_dict(out_files)
	#for key, value in SASA_dict.iteritems():
	#	aa_comp_dict = calculate_aa_comp(value)
	#	print "aa composition is: ", aa_comp_dict
	#	frac_surface_resi = calculate_frac_surface_resi(value, aa_comp_dict)
	#	print "fraction of aa on surface is: ", frac_surface_resi

	total_length_list, frac_surface_resi_list = build_lists(SASA_dict)
	#print "shortest sequence: ", min(total_length_list) 
	#print "longest sequence: ", max(total_length_list) 
	#print "lowest fraction surface residues: ", min(frac_surface_resi_list)
	#print "highest fraction surface residues: ", max(frac_surface_resi_list)
	#print "average +/- stddev", numpy.mean(frac_surface_resi_list), numpy.std(frac_surface_resi_list)
	plot_FracSurfaceResi_LengthProtein(total_length_list, frac_surface_resi_list)

	#for file_name, len_dict in output_dict:
	#	print min(output_dict), max(output_dict)
	
	#for file_name, pdb_dict in SASA_dict.iteritems():
	#	aa_comp_dict = calculate_aa_comp(pdb_dict)
	#	frac_surface_resi = calculate_frac_surface_resi(pdb_dict, aa_comp_dict)
	#	print aa_comp_dict
	#	print frac_surface_resi		
	

main()
