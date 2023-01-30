using Plots, SpecialFunctions, Random, DifferentialEquations
include("functions.jl")

#|||---VARIABLES----||||#
XLENGTH = 20
YLENGTH = 20
initialDensity = 0.2

times = [5,6,7,8]
numSimulations = 10
probMovement = 1
probProliferation = 0.1
BIAS = 0

STEPSIZE = 1
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)+1

#Produce averaged lattices with proliferation for a number of time steps
#Produce numerical solutions for a number of time steps using the logistic equation
#Graph these together

#Now set up an initial condition independent of the vertical position
#Run simulation with proliferation
#Solve the Fisher-Kolmogorov model numerically
#Plot together



#||||----TO DO----||||#
#  - Start to fill in this CODE
#  - Write function for the Fisher-Kolmogorov model