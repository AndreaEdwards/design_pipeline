from rosetta import *
rosetta.init()

scorefxn = get_fa_scorefxn()
print scorefxn

scorefxn2 = ScoreFunction()
scorefxn2.set_weight(fa_atr,1)
scorefxn2.set_weight(fa_rep,1)

from toolbax import pose_from_rcsb
ras = pose_from_rcsb("6Q21")
