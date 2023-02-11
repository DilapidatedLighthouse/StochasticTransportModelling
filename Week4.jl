using Plots, SpecialFunctions, Random, DifferentialEquations
include("functions.jl")
include("functionsWeek4.jl")




#Grid parameters
XLENGTH = 10
YLENGTH = 10
lengths = [XLENGTH, YLENGTH]
xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...]
NUMBEROFSIMULATIONS = 1
STEPSIZE = 1 

times = [0]
NUMBEROFSIMULATIONS = 1
biases = 0 #not properly implemented

probMovements = [1,1]
#if one of the grids contain a cell at the specified coordinates, the index of that grid will be returned. Else the function will return 0.

#Set up initial conditions
function createBlock(center, width, height, simGrid, xAxisValues)

    #ensure the height of the block is at most the height of the grid
    if(height > size(simGrid,2))
        height = size(simGrid,2)
    end#if
    
    #ensure the width of the block is at most the width of the grid. (That'd be a borring simulation, but you do you.)
    if(width > size(simGrid,1))
        width = size(simGrid,1)
    end#if

    #place the block
    for i in axes(simGrid,1)
        if(abs(xAxisValues[i]-center) <= width/2)
            for j in 1:height
                simGrid[i,j]=1.0
            end#for
        end#if
    end#for
    
    return simGrid
end#function
#population 1
H = 2
simGrid1 = zeros(XLENGTH,YLENGTH)

simGrid1 = createBlock(0,2*H, 100000, simGrid1,xAxisValues)
#population 2
H = 25
simGrid2 = zeros(XLENGTH,YLENGTH)

simGrid2 = createBlock(100,2*H, YLENGTH, simGrid2,xAxisValues)

simGridsMaster = [copy(simGrid1), copy(simGrid2)]







simulations = StochasticExclusionWalkAverageMultTimesMultiplePopulations(lengths, times, copy(simGridsMaster), NUMBEROFSIMULATIONS, probMovements, biases)

densities = map(sim -> map(pop -> calculateDensities(pop), sim), simulations)

scatter(xAxisValues, densities[1][1])
#Numerical solution for two populations

#Plot the two together