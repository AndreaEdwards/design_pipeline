#Iterate over all atoms in a structure:
import numpy
from Bio.PDB import *
p = PDBParser()
structure = p.get_structure('CBM','1FAT.pdb')
def inRange(value,limit_list):
    """"
    You input a value and check is it lies withing a range
    Input:
    value: The number you are checking
    limit_list: list of 2 values within which the number should lie.
    Output:
    Boolean
    """
    if value<max(limit_list[0],limit_list[1]) and value > min(limit_list[0],limit_list[1]):
        return True
    else:
        return False

count_vector = ()
chain = structure[0]['A']
residue_dict = {}
residue_list = [100,101,102,103,105,110]
for residue in residue_list:
    residue_info = chain[residue]
    res_no = residue
    res_list = []
    for atom in residue_info:
        x_limit = (atom.coord[0]+5,atom.coord[0]-5)
        y_limit = (atom.coord[1]+5,atom.coord[1]-5)
        z_limit = (atom.coord[2]+5,atom.coord[2]-5)
        count = 0
        for atom1 in chain.get_atoms():
            residue_1 = atom1.get_parent()
            res1_id = residue_1.id[1]
            if inRange(atom1.coord[0],x_limit) and inRange(atom1.coord[1],y_limit) and inRange(atom1.coord[2],z_limit):
                if res1_id not in res_list:
                    res_list.append(res1_id)
    residue_dict[res_no] = res_list
            
            
            
"""          
            #distance_vector = atom.coord - atom1.coord
            #distance = numpy.linalg.norm(distance_vector)
            #if distance <= 5:
                #count += 1
    #count_vector = count_vector + (count,) 
#print count_vector
                
"""
             

            