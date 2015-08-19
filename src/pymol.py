import __main__
__main__.pymol_argv = ['pymol', '-qei']
import pymol
pymol.finish_launching()

from pymol import cmd


pdb_file = '4PDJ.pdb'
pdb_colors = '4PDJ_lpocket.pdb'
pdb_name = '4PDJ'
cmd.load(pdb_colors, pdb_name)
cmd.spectrum("b", "blue_white_red", "%s and %s" % (pdb_name, pdb_name))
cmd.show("cartoon", pdb_name)
cmd.show("surface", pdb_name)

