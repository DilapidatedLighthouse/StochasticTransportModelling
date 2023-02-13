using Plots, SpecialFunctions, Random, DifferentialEquations
include("functions.jl")
include("FunctionsWeek4.jl")




#||||----Grid parameters----||||#
XLENGTH = 500
YLENGTH = 50
lengths = [XLENGTH, YLENGTH]
xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...]
STEPSIZE = 1 
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)

times = [0, 500]
NUMBEROFSIMULATIONS = 3
biases = 0 #not properly implemented

probMovements = [1,0.5]


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

#   population 1
H = 20
simGrid1 = zeros(XLENGTH,YLENGTH)

simGrid1 = createBlock(100,2*H, YLENGTH, simGrid1,xAxisValues)
#   population 2
H = 20
simGrid2 = zeros(XLENGTH,YLENGTH)

simGrid2 = createBlock(0,2*H, YLENGTH, simGrid2,xAxisValues)

#Had to make some name changes here and pass in copies because of some odd behaviour. Arrays were being changed when they shouldn't have been.
simGridsMaster = [copy(simGrid1), copy(simGrid2)]





#||||----Run Simulation----||||#

#Simulations is a nested array as follows [population][time][x,y]
simulations = StochasticExclusionWalkAverageMultTimesMultiplePopulations(lengths, times, deepcopy(simGridsMaster), NUMBEROFSIMULATIONS, probMovements, biases)
#densities is a nested array as follows [population][time][x]
densities = map(sim -> map(pop -> calculateDensities(pop), sim), simulations)


#||||----Numerical Solutions----||||#
nInitialConditions = zeros(2,NUMBEROFSTEPS)
nInitialConditions[1,:] = calculateDensities(simGrid1)
nInitialConditions[2,:] = calculateDensities(simGrid2)

#Access these solutions with nSolutions[population, x-value, time].
nSolutions = pdesolverWeek4(XLENGTH, STEPSIZE, NUMBEROFSTEPS, times, nInitialConditions, probMovements[1]/4,probMovements[2]/4)



#||||----Plot Solutions----||||#


# myGraph = scatter(xAxisValues, densities[1][2])
# myGraph = scatter!(xAxisValues, densities[2][2])
# myGraph = scatter!(xAxisValues, densities[1][2]+densities[2][2])
# display(myGraph)

# myheatmap = heatmap(transpose(simulations[1][2] + simulations[2][2]*2))
# display(myheatmap)
#Numerical solution for two populations

#Plot the two together
timeIndex = 2
myPlot = plot(xAxisValues, nSolutions[1,:,timeIndex])
myplot = plot!(xAxisValues, nSolutions[2,:,timeIndex])
myplot = scatter!(xAxisValues, densities[1][timeIndex])
myplot = scatter!(xAxisValues, densities[2][timeIndex])