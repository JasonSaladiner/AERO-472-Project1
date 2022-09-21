"""
This package will consist of given variables and general functions
"""



##givens##
alt = 12000     #ft                 #altitude of cruise
Wpp = 1000      #lbf                #Weight of passengars and payload
wash = 2        #degree             #Linear washout of the wings
RC = 400        #ft/min             #Rate of climb
spfc = 0.45     #lbf/hp/hr          #specific fuel consumption
Eprop = .85     #"%"                #prop efficiency
cf = 0.00045    #coef.              #skin friction coef
Etail = .9      #"%"                #efficinecy (eta) tail
kappa = .5      #non-dim            #downwash effect

#givens dict
#pass this using **given for ease
given = {"alt" : alt, "Wpp" : Wpp, "wash" : wash, "RC" : RC, "spfc" : spfc, "Eprop" : Eprop, "cf" : cf, "Etail" : Etail, "kappa": kappa}


#Given relations
VLength = lambda wo : 4.47*pow(wo,0.23)         #length [ft], wo[lbf]           #Vehicle length
Sfuse = lambda l : 3.1*pow(l,1.3)               #wetted fuse [ft^2], length fuse [ft]       #Wetted fuselage surface area
from math import pi
Power = lambda w,v, cdo, q, s, rc, eta,e,ar : w*v/eta * (cdo*q*s/w + w/s/pi/e/ar/q + rc/v)

##%%can be passed but less likely to be constly passed so im not going to as a default%%
givenRelations = {"VLength": VLength, "Sfuse" : Sfuse, "Power" : Power}


##Other Givens##
#Cruise : {300nm @ 150kts, 2hr}
#SM >= 15%
#stall speed <= 75kts (@ sea level, 60F, W=wo)
#Tail is NACA 0012
#Tail AR = ARwing - 2
#6% fuel reserves


