using Plots, SpecialFunctions, Random, DifferentialEquations
include("functions.jl")

#|||---VARIABLES----||||#
XLENGTH = 20
YLENGTH = 20
initialDensity = 0.2

times = 1:200
numSimulations = 1
probMovement = 1
probProliferation = 0.1
BIAS = 0

STEPSIZE = 0.50
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)+1
#Initialise a grid with randomly placed agents

simGrid = zeros(XLENGTH,YLENGTH)
for i in eachindex(simGrid)
    occupiedCheck = rand(1)[1]
    if(occupiedCheck < initialDensity)
        simGrid[i] = 1.0
    end#if
end#for

#Run simulations and record for given times

a = StochasticExclusionWalkAverageWithProliferationMultTimes([XLENGTH,YLENGTH], times, simGrid, numSimulations, probMovement, probProliferation)

# #Graph solutions on a grid
# for i in times
#     grid = PrepGridForDisplay(a[i])
#     p1 = scatter(grid[1],grid[2])
#     display(p1)
# end#for

densities = zeros(length(times))
for i in eachindex(times)
    densities[i] = sum(a[i])
end#for
C0 = fill(initialDensity, length(times))

numericSolutions = zeros(length(Times),NUMBEROFSTEPS) 
theoreticSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS, probMovement, probProliferation], C0, times, Logistic!)
p1 = scatter(times, densities)
p1 = plot!(times, theoreticSolutions)
display(p1)
