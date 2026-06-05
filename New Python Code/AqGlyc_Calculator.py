"""Glycerol_calculator.py"""

## This python code is based on the MatLab code orginaly provided by Chris Westbrook
## http://www.met.reading.ac.uk/~sws04cdw/viscosity_calc.html

__author__  = "Matthew Partridge"
__license__ = "GPL"
__version__ = "1.0"
__credits__ = "Chris Westbrook"


#Required packages ----------------

import numpy
import math

def AqGly_rho_eta(T,glycerolMass,waterMass):
    glycerolDen = (1273.3-0.6121*T)/1000 
    waterDen = (1-math.pow(((abs(T-4))/622),1.7))

    totalMass=glycerolMass+waterMass
    mass_fraction=glycerolMass/totalMass

    ##Andreas Volk polynomial method
    contraction_av = 1-math.pow(3.520E-8*((mass_fraction*100)),3)+math.pow(1.027E-6*((mass_fraction*100)),2)+2.5E-4*(mass_fraction*100)-1.691E-4
    contraction = 1+contraction_av/100

    waterVol = waterMass/waterDen
    glycerolVol= glycerolMass/glycerolDen
    vol_fraction= glycerolVol/(glycerolVol+waterVol)
    ## Distorted sine approximation method
    #contraction_pc = 1.1*math.pow(math.sin(numpy.radians(math.pow(mass_fraction,1.3)*180)),0.85)
    #contraction = 1 + contraction_pc/100

    rho = (glycerolDen*vol_fraction+waterDen*(1-vol_fraction))*contraction * 1000

    #Viscosity calcualtor ----------------

    glycerolVisc=0.001*12100*numpy.exp((-1233+T)*T/(9900+70*T))
    waterVisc=0.001*1.790*numpy.exp((-1230-T)*T/(36100+360*T))

    a=0.705-0.0017*T
    b=(4.9+0.036*T)*numpy.power(a,2.5)
    alpha=1-mass_fraction+(a*b*mass_fraction*(1-mass_fraction))/(a*mass_fraction+b*(1-mass_fraction))
    A=numpy.log(waterVisc/glycerolVisc)

    eta = glycerolVisc*numpy.exp(A*alpha) * 1000
    return rho, eta
"""
#Variables ----------------

T = 35	#temperature (degrees Celcius)

#Densities ----------------

glycerolDen = (1273.3-0.6121*T)/1000 			#Density of Glycerol (g/cm3)
waterDen = (1-math.pow(((abs(T-4))/622),1.7)) 	#Density of water (g/cm3)


#Fraction cacluator ----------------

#glycerolMass=glycerolDen*glycerolVol
#waterMass=waterDen*waterVol
glycerolMass = 0
waterMass = 1
totalMass=glycerolMass+waterMass
mass_fraction=glycerolMass/totalMass

waterVol = waterMass/waterDen
glycerolVol= glycerolMass/glycerolDen
#print(waterVol)
#print(glycerolVol)
vol_fraction= glycerolVol/(glycerolVol+waterVol)

print ("Mass fraction of mixture =", mass_fraction)
print ("Volume fraction of mixture =", vol_fraction)

#Density calculator ----------------

##Andreas Volk polynomial method
contraction_av = 1-math.pow(3.520E-8*((mass_fraction*100)),3)+math.pow(1.027E-6*((mass_fraction*100)),2)+2.5E-4*(mass_fraction*100)-1.691E-4
contraction = 1+contraction_av/100

## Distorted sine approximation method
#contraction_pc = 1.1*math.pow(math.sin(numpy.radians(math.pow(mass_fraction,1.3)*180)),0.85)
#contraction = 1 + contraction_pc/100

density_mix=(glycerolDen*vol_fraction+waterDen*(1-vol_fraction))*contraction

print ("Density of mixture =",density_mix*1000,"kg/m3")


#Viscosity calcualtor ----------------

glycerolVisc=0.001*12100*numpy.exp((-1233+T)*T/(9900+70*T))
waterVisc=0.001*1.790*numpy.exp((-1230-T)*T/(36100+360*T))

a=0.705-0.0017*T
b=(4.9+0.036*T)*numpy.power(a,2.5)
alpha=1-mass_fraction+(a*b*mass_fraction*(1-mass_fraction))/(a*mass_fraction+b*(1-mass_fraction))
A=numpy.log(waterVisc/glycerolVisc)

viscosity_mix=glycerolVisc*numpy.exp(A*alpha)

print ("Viscosity of mxiture =",viscosity_mix*1000, "cP")"""