using Plots, SpecialFunctions, Random

#|||---VARIABLES----||||#
XLENGTH = 500
YLENGTH = 500
initialDensity = 0.5


#Initialise a grid with randomly placed agents

simGrid = zeros(XLENGTH,YLENGTH)
for i in simGrid
    occupiedCheck = rand(1)
    if
end#for


#Run proliferation simulations. Record averages for a number of timesteps

#Graph number of agents vs the logistic equation