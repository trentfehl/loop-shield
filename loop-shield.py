# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# description
print "============================================================================================="
print "                    Active Spacecraft Radiation Shielding Calculations                       "
print "============================================================================================="
print "This program is for calculating the displacement an incoming energetic particle will experie-"
print "nce due to a local magnetic field generated from current through a loop of wire. The particle"
print "is assumed to be coming from the positive X direction and is deflected in one of the Y direc-"
print "tions. " 

print "                                         /@\                                                  " 
print "                                     Z  //@\\\                                                 "
print "                                          @                                                   " 
print "                                          @                                                   " 
print "                                          @                                                   " 
print "                                          @                                                   " 
print "                                          @                                                   " 
print "                                          @                                                   " 
print "                      ,,,-----**````````` @ `````````**---,,,,                                " 
print "              __--````                    @       R   \#0|    `````--__                       "
print "           .@/                            @        __-/*\|             `*&_                   "  
print "        %%                                @    __-/*                       `@             X   " 
print "      *(                                  @__-/*                             #        |,      "  
print "      @                                 _/@+++++++++++++++++++++++++++++++++++++++++++@@>     "                             
print "      |@                              _/* @                                  #        |*      "
print "        *#_                         _/*   @                               ,@                  " 
print "           ``-_                   _/*     @                         __--`                     "    
print "               ``**--,,,_       _/*       @                 _,,,--**``` _____                 "     
print "                       ```****||*<<<<<<<<<@>>>>>>>>*******```          __\$||                 "
print "                            _/*           @                   _,,,--**```  \`                 "
print "                          _/*             @            *****```         I                     "     
print "                       |v/,               @                                                   "      
print "                   Y   7/**`              @                                                   "
print "                                          @                                                   "
print "                                          @                                                   "
print "                                          @                                                   "
print "                                          @                                                   "
print "                                          @                                                   "
print "                                                                                              "

#importing packages
import math
import sys
import random
import time
import numpy
import scipy
import scipy.integrate
from tqdm import *
# from plotly.graph_objs import *
# import plotly.plotly as py
# from decimal import Decimal

while True:
        value = raw_input('Do you want to continue? (y/n) ')
        try:
                answer = str(value)
        except ValueError:
                print 'Valid answer, please'
                continue
        if answer == "y":
                break
        else:
                quit()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# variables

# magnetic field variables
print "========================================"
print "Magnetic Field Variable Inputs"
y = 0 # coordinates
z = 0 # coordinates
phi = 0 # phi angle?
mu_0 = math.pi*4*(10 ** -7)

# loop radius
while True:
	value = raw_input('Input a value for the radius of the current loop (meters): ')
	try:
		global a
        	a = float(value)
	except ValueError:
		print 'Valid number, please'
		continue
	if 0 < a:
		print "Loop Radius = " + str(a)
        	break
    	else:
       		print 'Valid range, please: 0 < value'
# loop radius
while True:
        value = raw_input('Input a value for the current magnitude (amps): ')
        try:
		global i
                i = float(value)
        except ValueError:
                print 'Valid number, please'
                continue
        if 0 < i:
                print "Loop Current = " + str(i)
                break
        else:
                print 'Valid range, please: 0 < value'
'''
# x minimum
while True:
        value = raw_input('Input the closest distance to the ring to stop integrating (meters): ')
        try:
		global xmin
                xmin = float(value)
        except ValueError:
                print 'Valid number, please'
                continue
        if a < xmin:
                print "X Min = " + str(xmin)
                break
        else:
                print 'Valid range, please: Loop Radius < value'
'''
# removed user input for this for simplicity
xmin = a*1.001

# x maximum
while True:
        value = raw_input('Input the distance from the loop to start integrating B field (meters): ')
        try:
		global xmax
                xmax = (float(value) - xmin)
        except ValueError:
                print 'Valid number, please'
                continue
        if xmin < xmax:
                print "Distance = " + str(xmax)
                break
        else:
                print 'Valid range, please: X Min < value'


# incoming particle variables
print
print "========================================"
print "Incoming Particle Variable Inputs"
q = 1.60218*(10 ** -19) 
c = 2.998*(10 ** 8)
mass_H = 1.6726*(10 ** -27)
mass_e = 0.9109*(10 ** -30)

# particle species
while True:
        value = raw_input('Input incoming particle value (Z # for ions and -1 for electrons): ')
        try:
		global Z
                Z = int(value)
        except ValueError:
                print 'Valid number, please'
                continue
        if Z == 1 or Z == -1:
                print "Z Value = " + str(Z)
                break
        else:
                print 'Sorry only Hydrogen and electrons are supported right now'
# incoming particle energy
while True:
        value = raw_input('Input kinetic energy of incoming particle (MeV): ')
        try:
		global E
                E = float(value)
        except ValueError:
                print 'Valid number, please'
                continue
        if 0 < E:
                print "E = " + str(E)
                break
        else:
                print 'Valid range, please: 0 < value'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to integrate
def f(x):	
	global y
	global z
	global alpha_sq
	global beta_sq
	global k_sq

	C = mu_0*i/math.pi	
	rho_sq = (x ** 2) + (y ** 2)
	r_sq = (x ** 2) + (y ** 2) + (z ** 2)
	alpha_sq =  (a ** 2) + r_sq - 2*a*math.sqrt(rho_sq)
	beta_sq = (a ** 2) + r_sq + 2*a*math.sqrt(rho_sq)
	k_sq = 1 - alpha_sq/beta_sq
	s = math.sin(phi)

	# elliptical integrals
	def E(k_sq, phi):	
		return math.sqrt(1 - k_sq*(s ** 2))
	e_k, err = scipy.integrate.quad(E, 0, math.pi/2, args=(k_sq)) 

	def K(k_sq, phi):
        	return 1/math.sqrt(1 - k_sq*(s ** 2))
	k_k, err = scipy.integrate.quad(K, 0, math.pi/2, args=(k_sq))  
	
	return C*(((a ** 2)-r_sq)*e_k - alpha_sq*k_k)/(2*alpha_sq*math.sqrt(beta_sq))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# x range
xrange = xmax - xmin
# print "X Range = " + str(xrange)

# b range determination
print
print "========================================"
print "Determining B(x) Range in Tesla"

# number of steps for B(x) determination
while True:
        value = raw_input('Input the number of steps used in B(x) determination: ')
        try:
                global numSteps
                numSteps = int(value)
        except ValueError:
                print 'Valid number, please'
                continue
        if 0 < numSteps:
                print "Number of Steps = " + str(numSteps)
                break
        else:
                print 'Valid range, please: 0 < value'

bmin = f(xmin)
bmax = bmin
brange = bmax - bmin
for i in tqdm(range(numSteps)):
	#	global x
	#	global b

	x = xmin + (xrange*float(i)/numSteps)
	b = f(x)
	if b < bmin: bmin = b
	if b > bmax: bmax = b
brange = bmax - bmin
print "B Min =		%.4E" % bmin
print "B Max =		%.4E" % bmax
print "B Range =	%.4E" % brange

print
print "========================================"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# monte carlo
print "Monte Carlo Runs for Integrated B Field in Tesla"

# number of points per run 
while True:
        value = raw_input('Input the number of points per run: ')
        try:
                global numPoints
                numPoints = int(value)
        except ValueError:
                print 'Valid number, please'
                continue
        if 0 < numPoints:
                print "Number of Points = " + str(numPoints)
                break
        else:
                print 'Valid range, please: 0 < value'

# number of monte carlo runs
while True:
        value = raw_input('Input the number of Monte Carlo runs to perform: ')
        try:
                global numRuns
                numRuns = int(value)
        except ValueError:
                print 'Valid number, please'
                continue
        if 0 < numRuns:
                print "Number of Runs = " + str(numRuns)
                break
        else:
                print 'Valid range, please: 0 < value'

if bmax >= brange:
	rectArea = xrange * math.fabs(bmax)
	# sys.stdout.write("Integration Guess:	%.4E" % rectArea)
	# print "      NOTE: Guess is result of B Max multiplied by X range"
elif math.fabs(bmin) >= brange:
	rectArea = xrange * math.fabs(bmin)
	# sys.stdout.write("Integration Guess:	%.4E" % rectArea)
	# print "	     NOTE: Guess is result of B Min multiplied by X range"
else:
	rectArea = xrange * brange #initial guess at integral
	# sys.stdout.write("Integration Guess:	%.4E" % rectArea)
	# print "	     NOTE: Guess is result of B Range multiplied by X range"

b_tot = 0
b_tot_sq = 0
for k in range(numRuns):
	ctr = 0
	b_int = 0
	# sys.stdout.write("Run {:d} of {:d}".format((k+1),numRuns))
	for j in tqdm(range(numPoints)):
		#	global x
		#	global b
	
		x_ran = random.SystemRandom()
		b_ran = random.SystemRandom()
		x = xmin + (xrange*x_ran.random())
		b = bmin + (brange*b_ran.random())
		if math.fabs(b) <= math.fabs(f(x)):
			if f(x) > 0 and b > 0:
				ctr += 1
			if f(x) < 0 and b < 0:
				ctr -= 1
		# print "f(x): " + str(f(x))
		# print "B: " + str(b)
	b_int = rectArea * float(ctr) / numPoints
	print "Run {:d} of {:d} complete. Value: {:.4E}".format((k+1),numRuns,b_int)
	b_tot += b_int
	b_tot_sq += (b_int ** 2)
b_avg = b_tot/numRuns

print "MC Integration Result:	%.4E" % b_avg

I = b_avg 
I_sq = b_tot_sq/numRuns

sig = I_sq - (I ** 2)
print "1 Sigma Uncertainty:	{:.4E}" .format(sig) 
print "Note: uncertainty proportional to 1/sqrt(N) where N is the number of points per run" 

print
print "========================================"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# incoming particle relativistic velocity and mass
if Z == 1:
	m_rest = mass_H
elif Z < 0:
	m_rest = mass_e

E_joules = 1.60218*(10 ** -13)*E

v = c*math.sqrt(1 - ((m_rest*c*c/(E_joules + m_rest*c*c)) ** 2))

gamma = 1/math.sqrt(1 - ((v/c) ** 2))

m_rel = gamma * m_rest

print "Relativistic Incoming Particle Velocity and Mass"
print "Kinetic Energy (MeV) = 		%.2E" % E
print "Velocity (m/s) = 		%.4E" % v
print "Velocity (v/c) =		{0:.0f}%".format(v/c * 100)
print "Rest Mass (kg) =		%.5E" % m_rest
print "Relativistic Mass (kg) =	%.5E" % m_rel
print "Gamma = 			%.7E" % gamma

print
print "========================================"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# V x B force calculations
t = xrange / v
# print "T = " + str(t)

F = q*v*b_int
accel = F/m_rel
# print "Acceleration = " + str(accel)

y_delta = 0.5*accel*(t ** 2)

print "Incoming Particle Displacement due to qVxB Active Shield"
print "Y Axis Displacement (meters):	%.5E" % math.fabs(y_delta)
print "========================================"
