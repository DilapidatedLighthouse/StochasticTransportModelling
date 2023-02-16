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

times = [0, 500, 1000, 1500]
#times = [0,10,200,300]
NUMBEROFSIMULATIONS = 10
biases = 0 #not properly implemented

probMovements = [0.8,0.8]


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




# function prepareGraphNoLabel(timeIndex, times, probMovements)
    
#     time = times[timeIndex]
#     timeString = string(times[timeIndex])

#     plot1 = scatter(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex],markerstrokewidth = 0, markershape = :rect, markersize = 4, legend = false)

#     #Stochastic values
#     for populationIndex in eachindex(probMovements)
#         plot1 = scatter!(xAxisValues,  densities[populationIndex][timeIndex], markerstrokewidth = 0, markersize = 2.6, legend = false)
#     end#for
#     #Numeric Values
#     for populationIndex in eachindex(probMovements)
#         plot1 = plot!(xAxisValues,nSolutions[populationIndex, :, timeIndex], lw = 2.5, xlabel = "x", ylabel="Density", framestyle = :box, legend = false)
#     end#for

#     #plot1 = scatter!(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex],markerstrokewidth = 0, markershape = :rect, label="Stochastic total")
#     plot1 = plot!(xAxisValues,  nSolutions[1, :, timeIndex] + nSolutions[2, :, timeIndex],xlabel = "x", title = string("t = ", time), lc = :black, ls=:dash, ylabel="Density", framestyle = :box, legend = false, xlims = (-XLENGTH/2, XLENGTH/2), ylims = (0,1) )

#     return plot1

# end#function


# #plots for 
# plot1 = prepareGraphNoLabel(1, times, probMovements)
# plot2 = prepareGraphNoLabel(2, times, probMovements)
# plot3 = prepareGraphNoLabel(3, times, probMovements)
# plot4 = prepareGraphNoLabel(4, times, probMovements)

# combinedPlot = plot(plot1,plot2,plot3,plot4, legend = false, xlims = (-XLENGTH/2, XLENGTH/2), ylims = (0,1) )
# display(combinedPlot)



#||||-----Special Conditions Graph----||||#


function prepareGraphSpecial(timeIndex, times, probMovements)
    
    time = times[timeIndex]
    timeString = string(times[timeIndex])
    
    D = probMovements[1]/4
    C(x)=0.5*(erf((H1-(x-center1))/sqrt(4*D*time))+erf((H1+(x-center1))/sqrt(4*D*time)) + erf((H2-(x-center2))/sqrt(4*D*time))+erf((H2+(x-center2))/sqrt(4*D*time)));

    plot1 = scatter(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex],markerstrokewidth = 0, markershape = :rect, markersize = 4, legend = false)
    plot1=plot!(C, -XLENGTH/2:XLENGTH/2-1,lw=4,lc=:red,ls=:dash,label="Exact",xlabel="x",ylabel="C(x, 0)",xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1),framestyle=:box)
    #plot1 = scatter!(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex],markerstrokewidth = 0, markershape = :rect, label="Stochastic total")
    plot1 = plot!(xAxisValues,  nSolutions[1, :, timeIndex] + nSolutions[2, :, timeIndex],xlabel = "x", title = string("t = ", time), lc = :black, ls=:dash, ylabel="Density", framestyle = :box, legend = false, xlims = (-XLENGTH/2, XLENGTH/2), ylims = (0,1) )

    return plot1

end#function


plot1 = prepareGraphSpecial(1, times, probMovements)
plot2 = prepareGraphSpecial(2, times, probMovements)
plot3 = prepareGraphSpecial(3, times, probMovements)
plot4 = prepareGraphSpecial(4, times, probMovements)

combinedPlot = plot(plot1,plot2,plot3,plot4, legend = false, xlims = (-XLENGTH/2, XLENGTH/2), ylims = (0,1) )
display(combinedPlot)


#display(prepareGraph(3,times, probMovements))







function prepareGraph(timeIndex, times, probMovements)
    
    time = times[timeIndex]
    timeString = string(times[timeIndex])

    plot1 = scatter(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex],markerstrokewidth = 0, markershape = :rect, markersize = 4, label="Stochastic total")

    #Stochastic values
    for populationIndex in eachindex(probMovements)
        plot1 = scatter!(xAxisValues,  densities[populationIndex][timeIndex], markerstrokewidth = 0, markersize = 2.6, label=string("Stochastic ", populationIndex))
    end#for
    #Numeric Values
    for populationIndex in eachindex(probMovements)
        plot1 = plot!(xAxisValues,nSolutions[populationIndex, :, timeIndex], lw = 2.5, xlabel = "x", ylabel="Density", framestyle = :box, label=string("Numeric ", populationIndex))
    end#for

    #plot1 = scatter!(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex],markerstrokewidth = 0, markershape = :rect, label="Stochastic total")
    plot1 = plot!(xAxisValues,  nSolutions[1, :, timeIndex] + nSolutions[2, :, timeIndex],xlabel = "x", title = string("t = ", time), lc = :black, ls=:dash, ylabel="Density", framestyle = :box, label=string("Numeric total"))

    return plot1

end#function










# # myGraph = scatter(xAxisValues, densities[1][2])
# # myGraph = scatter!(xAxisValues, densities[2][2])
# # myGraph = scatter!(xAxisValues, densities[1][2]+densities[2][2])
# # display(myGraph)

# # myheatmap = heatmap(transpose(simulations[1][2] + simulations[2][2]*2))
# # display(myheatmap)
# #Numerical solution for two populations

# #Plot the two together
# timeIndex = 2
# myPlot = plot(xAxisValues, nSolutions[1,:,timeIndex])
# myplot = plot!(xAxisValues, nSolutions[2,:,timeIndex])
# myplot = scatter!(xAxisValues, densities[1][timeIndex])
# myplot = scatter!(xAxisValues, densities[2][timeIndex])

# #||||----Plot standard----||||#


# #plot1 = scatter(xAxisValues, nSolutions[1, :, timeIndex] + nSolutions[2, :, timeIndex], markerstrokewidth = 0, markersize = 0.5, markershape = :rect, label="Stochastic total")
# #plot1 = plot(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex], lw=2, ls=:dash, xlabel = "x", ylabel="Density", framestyle = :box, label=string("Numeric total"))


# timeIndex = 3
# timeString = string(times[timeIndex])

# plot1 = scatter(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex],markerstrokewidth = 0, markershape = :rect, markersize = 2.8, label="Stochastic total")

# populationIndex = 1
# plot1 = scatter!(xAxisValues,  densities[populationIndex][timeIndex], markerstrokewidth = 0, markersize = 2, label=string("Stochastic ", populationIndex))
# #plot1 = plot!(xAxisValues,nSolutions[populationIndex, :, timeIndex],  xlabel = "x", ylabel="Density", framestyle = :box, label=string("Numeric ", populationIndex))


# populationIndex = 2
# plot1 = scatter!(xAxisValues, densities[populationIndex][timeIndex], markercolor = :lightgreen, markerstrokewidth = 0, markersize = 2, label=string("Stochastic ", populationIndex))
# #plot1 = plot!(xAxisValues, nSolutions[populationIndex, :, timeIndex], xlabel = "x", ylabel="Density", framestyle = :box, label=string("Numeric ", populationIndex))


# populationIndex = 1
# plot1 = plot!(xAxisValues,nSolutions[populationIndex, :, timeIndex], ls=:dash, lc=:red, xlabel = "x", ylabel="Density", framestyle = :box, label=string("Numeric ", populationIndex))


# populationIndex = 2
# plot1 = plot!(xAxisValues, nSolutions[populationIndex, :, timeIndex], ls=:dash, lc = :green, xlabel = "x", ylabel="Density", framestyle = :box, label=string("Numeric ", populationIndex))



# #plot1 = scatter!(xAxisValues, densities[1][timeIndex]+ densities[2][timeIndex],markerstrokewidth = 0, markershape = :rect, label="Stochastic total")
# plot1 = plot!(xAxisValues,  nSolutions[1, :, timeIndex] + nSolutions[2, :, timeIndex],xlabel = "x", lc = :black, ls=:dash, ylabel="Density", framestyle = :box, label=string("Numeric total"))


# display(plot1)


