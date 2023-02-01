using Plots, SpecialFunctions, Random, DifferentialEquations
include("functions.jl")

#|||---VARIABLES----||||#
XLENGTH = 50
YLENGTH = 50
initialDensity = 0.2
times = 0:1000
numSimulations = 2
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


#||||----QUESTION 1----||||#


# #Produce averaged lattices with proliferation for a number of time steps
# averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimesWithRandomIC([XLENGTH,YLENGTH], times, numSimulations, probMovement, probProliferation, initialDensity)
# densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
# densities_full = map(x -> sum(x)/length(x), densities_y) #Average results over x direction
# initialDensity = densities_full[1]
# #Produce numerical solutions for a number of time steps using the logistic equation
# theoreticalDensities = map(t -> initialDensity/((1-initialDensity)*exp(-probProliferation*t) + initialDensity), times)
# #Graph these together
# timeChangePlot = plot(times, densities_full)
# timeChangePlot = plot!(times, theoreticalDensities)
# display(timeChangePlot)



#||||----QUESTION 2----||||#

#Now set up an initial condition independent of the vertical position
XLENGTH = 250
YLENGTH = 50
STEPSIZE = 1
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)
times = [100]
H = 50

xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...]

simGrid = zeros(XLENGTH,YLENGTH)

simGrid = createBlock(0,2*H, YLENGTH, simGrid,xAxisValues)
#Run simulation with proliferation
averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimes([XLENGTH,YLENGTH], times, simGrid, numSimulations, probMovement, probProliferation)
densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction

#Solve the Fisher-Kolmogorov model numerically
C0 = calculateDensities(simGrid) 
numericSolutions = zeros(length(times),NUMBEROFSTEPS) 
theoreticSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS, probMovement, probProliferation], C0, times, FisherKolmogorov!)
#Plot together



#||||----TO DO----||||#
#  - Start to fill in this CODE
#  - Write function for the Fisher-Kolmogorov model