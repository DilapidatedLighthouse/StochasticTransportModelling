using Plots, SpecialFunctions, Random
include("functions.jl")

#|||---VARIABLES----||||#
XLENGTH = 10
YLENGTH = 10
initialDensity = 0.5

totalTime = 10
numSimulations = 1
probMovement = 1
BIAS = 0

#Initialise a grid with randomly placed agents

simGrid = zeros(XLENGTH,YLENGTH)
for i in eachindex(simGrid)
    occupiedCheck = rand(1)[1]
    if(occupiedCheck < initialDensity)
        simGrid[i] = 1.0
    end#if
end#for

a = StochasticExclusionWalkAverage([XLENGTH,YLENGTH], totalTime, simGrid, numSimulations, probMovement, BIAS)
grid = PrepGridForDisplay(a)
p1 = scatter(grid)
display(p1)
#Run proliferation simulations. Record averages for a number of timesteps

#Graph number of agents vs the logistic equation