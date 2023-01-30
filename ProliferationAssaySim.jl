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

STEPSIZE = 0.50
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)+1
#Initialise a grid with randomly placed agents


#Randomly distributing elements on array
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
# for i in eachindex(times)
#     gridmy = PrepGridForDisplay(a[1])
#     p1 = scatter(gridmy[1],gridmy[2])
#     display(p1)
# end#for

# densities = zeros(length(times))
# for i in eachindex(times)
#     densities[i] = sum(a[i])
# end#for

densities = map(x -> calculateDensities(x), a)
C0 = fill(initialDensity, NUMBEROFSTEPS)

numericSolutions = zeros(length(times),NUMBEROFSTEPS) 
theoreticSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS, probMovement, probProliferation], C0, times, FisherKolmogorov!)
p1 = scatter(-XLENGTH/2:XLENGTH/2 , densities[4])
p1 = plot!(-XLENGTH/2:STEPSIZE:XLENGTH/2, theoreticSolutions[4,:])
display(p1)
