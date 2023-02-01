using Plots, SpecialFunctions, Random, DifferentialEquations
include("functions.jl")

#|||---VARIABLES----||||#
XLENGTH = 50
YLENGTH = 50
initialDensity = 0.2
times = 0:1000
numSimulations = 20
probMovement = 1
probProliferation = 0.01
BIAS = 0

STEPSIZE = 1
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)+1


#Setup evenly distributed agents on a grid

#Randomly distributing elements on array
# randCount = 0
# simGrid = zeros(XLENGTH,YLENGTH)
# for i in eachindex(simGrid)
#     occupiedCheck = rand(1)[1]
#     if(occupiedCheck < initialDensity)
#         randCount += 1
#         simGrid[i] = 1.0
#     end#if
#     if(randCount >= initialDensity*XLENGTH*YLENGTH)
#         break
#     end#if
# end#for


#Produce averaged lattices with proliferation for a number of time steps
averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimesWithRandomIC([XLENGTH,YLENGTH], times, numSimulations, probMovement, probProliferation, initialDensity)
densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
densities_full = map(x -> sum(x)/length(x), densities_y) #Average results over x direction
initialDensity = densities_full[1]
#Produce numerical solutions for a number of time steps using the logistic equation
theoreticalDensities = map(t -> initialDensity/((1-initialDensity)*exp(-probProliferation*t) + initialDensity), times)
#Graph these together
timeChangePlot = plot(times, densities_full)
timeChangePlot = plot!(times, theoreticalDensities)
display(timeChangePlot)

#Now set up an initial condition independent of the vertical position
XLENGTH = 20
YLENGTH = 20
#Run simulation with proliferation
#Solve the Fisher-Kolmogorov model numerically
#Plot together



#||||----TO DO----||||#
#  - Start to fill in this CODE
#  - Write function for the Fisher-Kolmogorov model