using Plots, SpecialFunctions, Random, DifferentialEquations
include("functions.jl")
include("FunctionsWeek4.jl")




#||||----Parameters----||||#
XLENGTH = 500
YLENGTH = 50
lengths = [XLENGTH, YLENGTH]
xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...] #For Graphing against later

STEPSIZE = 1
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)

biases = 0 #Not properly implemented, insufficient time to remove


TIMES = [0, 500, 1000, 1500] #Times for which the calculations will be saved
#TIMES = [0,10,15,20]
NUMBEROFSIMULATIONS = 10

probMovements = [0.8,0.8] #Probabilty that an agent will move for each population


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
H1 = 20
simGrid1 = zeros(XLENGTH,YLENGTH)
center1 = -100
simGrid1 = createBlock(center1,2*H1, YLENGTH, simGrid1,xAxisValues)
#   population 2
H2 = 20
center2 = 0
simGrid2 = zeros(XLENGTH,YLENGTH)

simGrid2 = createBlock(center2,2*H2, YLENGTH, simGrid2,xAxisValues)

#Had to make some name changes here and pass in copies because of some odd behaviour. Arrays were being changed when they shouldn't have been.
simGridsMaster = [copy(simGrid1), copy(simGrid2)]





#||||----Run Simulation----||||#

#Simulations is a nested array as follows [population][time][x,y]
simulations = StochasticExclusionWalkAverageMultTimesMultiplePopulations(lengths, TIMES, deepcopy(simGridsMaster), NUMBEROFSIMULATIONS, probMovements, biases)
#densities is a nested array as follows [population][time][x]
densities = map(sim -> map(pop -> calculateDensities(pop), sim), simulations)


#||||----Numerical Solutions----||||#
nInitialConditions = zeros(2,NUMBEROFSTEPS)
nInitialConditions[1,:] = calculateDensities(simGrid1)
nInitialConditions[2,:] = calculateDensities(simGrid2)

#Access these solutions with nSolutions[population, x-value, time].
nSolutions = pdesolverWeek4(XLENGTH, STEPSIZE, NUMBEROFSTEPS, TIMES, nInitialConditions, probMovements[1]/4,probMovements[2]/4)



#||||----Plot Solutions----||||#




function prepareGraphNoLabel(timeIndex, times, probMovements)
    
    time = times[timeIndex]
    timeString = string(times[timeIndex])

    #total stochastic
    plot1 = scatter(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex],markerstrokewidth = 0, markershape = :rect, markersize = 4)

    #Stochastic values
    for populationIndex in eachindex(probMovements)
        plot1 = scatter!(xAxisValues,  densities[populationIndex][timeIndex], markerstrokewidth = 0, markersize = 2.6)
    end#for
    #Numeric Values
    for populationIndex in eachindex(probMovements)
        plot1 = plot!(xAxisValues,nSolutions[populationIndex, :, timeIndex], lw = 2.5, xlabel = "x", ylabel="Density", framestyle = :box)
    end#for
    #total numeric
    plot1 = plot!(xAxisValues,  nSolutions[1, :, timeIndex] + nSolutions[2, :, timeIndex],xlabel = "x", title = string("t = ", time), lc = :black, ls=:dash, ylabel="Density", framestyle = :box, xlims = (-XLENGTH/2, XLENGTH/2), ylims = (0,1) )

    return plot1

end#function


#plots for 
plot1 = prepareGraphNoLabel(1, TIMES, probMovements)
plot2 = prepareGraphNoLabel(2, TIMES, probMovements)
plot3 = prepareGraphNoLabel(3, TIMES, probMovements)
plot4 = prepareGraphNoLabel(4, TIMES, probMovements)

combinedPlot = plot(plot1,plot2,plot3,plot4, legend = false, xlims = (-XLENGTH/2, XLENGTH/2), ylims = (0,1) )
display(combinedPlot)



#||||-----Special Conditions Graph----||||#
#When probabilities of movement are identical for all populations, the PDE can be solved exactly. 

function prepareGraphSpecial(timeIndex, times, probMovements)
    
    time = times[timeIndex]
    timeString = string(times[timeIndex])
    
    D = probMovements[1]/4
    C(x)=0.5*(erf((H1-(x-center1))/sqrt(4*D*time))+erf((H1+(x-center1))/sqrt(4*D*time)) + erf((H2-(x-center2))/sqrt(4*D*time))+erf((H2+(x-center2))/sqrt(4*D*time)));

    plot1 = scatter(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex],markerstrokewidth = 0, markershape = :rect, markersize = 4)
    plot1=plot!(C, -XLENGTH/2:XLENGTH/2-1,lw=4,lc=:red,ls=:dash,label="Exact",xlabel="x",ylabel="C(x, 0)",xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1),framestyle=:box)
    plot1 = plot!(xAxisValues,  nSolutions[1, :, timeIndex] + nSolutions[2, :, timeIndex],xlabel = "x", title = string("t = ", time), lc = :black, ls=:dash, ylabel="Density", framestyle = :box, xlims = (-XLENGTH/2, XLENGTH/2), ylims = (0,1) )

    return plot1

end#function


plot1 = prepareGraphSpecial(1, TIMES, probMovements)
plot2 = prepareGraphSpecial(2, TIMES, probMovements)
plot3 = prepareGraphSpecial(3, TIMES, probMovements)
plot4 = prepareGraphSpecial(4, TIMES, probMovements)

combinedPlot2 = plot(plot1,plot2,plot3,plot4, legend = false, xlims = (-XLENGTH/2, XLENGTH/2), ylims = (0,1) )
display(combinedPlot2)






