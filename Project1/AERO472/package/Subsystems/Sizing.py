"""
Do size calculations that can be used in the top level design.
"""

import package.Subsystems as sub
from .. import given


def range_relation(R = None,vel = None,LoD = None,spcf= None,wiOverwf= None):
	"""
	Use the	Breguet Range equation for this relation.
	if the plane is a piston prop, pass prop efficincy instead of velocity
	**Ensure that values units properelly cancel**
	"""
		
		
	import sys
	from math import log,exp
	#count of the number of None Types
	cNone = [R,vel,LoD,spcf,wiOverwf].count(None)
		
	if cNone != 1:
		sys.exit("Incorrect number of None Variables in AERO472.Sizing.range_relation. Make sure to pass only one None for variable of interest")

	if R == None:
		return vel / spcf * LoD * log(wiOverwf)
	elif vel == None:
		return pow(1 / R / spcf * LoD * log(wiOverwf),-1)
	elif LoD == None:
		return pow(vel / R/ spcf * log(wiOverwf),-1)
	elif spcf == None:
		return vel / R * LoD * log(wiOverwf)
	elif wiOverwf == None:
		x = R * spcf / LoD / vel
		return 1 - exp(-x)


def endurance_relation(E = None, spcf = None, eff = None, vel = None, LoD = None, wiOverwf = None):
	"""
	Using the loiter endurance equaiton
	this is true for a piston prop engine. 
	"""

	import sys
	from math import log,exp
	#count of the number of None Types
	cNone = [E,spcf,LoD,wiOverwf].count(None)
	
	if cNone != 1:
		sys.exit("Incorrect number of None Variables in AERO472.Sizing.endurance_relations. Make sure to pass only one None for variable of interest")

	if E == None:
		return 1 / spcf * LoD * log(wiOverwf)
	elif spcf == None:
		return 1 / E * LoD * log(wiOverwf)
	elif LoD == None:
		return E * spcf / log(wiOverwf)
	elif wiOverwf == None:
		x = E * spcf * vel / LoD / eff
		print(x)
		return 1 - exp(-x)



#find initial weight guess woi
class weight:
    def __init__(self):
        self.wpp = given["Wpp"]
        print(self.wpp)