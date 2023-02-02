using Plots, SpecialFunctions, Random, DifferentialEquations, Measures
include("functions.jl")

#|||---VARIABLES----||||#
XLENGTH = 100
YLENGTH = 100
initialDensity = 0.2
times = 1:200
numSimulations = 10
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


# #||||----QUESTION 1----||||#


# #Produce averaged lattices with proliferation for a number of time steps
# probProliferation = 0.1
# averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimesWithRandomIC([XLENGTH,YLENGTH], times, numSimulations, probMovement, probProliferation, initialDensity)
# densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
# densities_full = map(x -> sum(x)/length(x), densities_y) #Average results over x direction
# initialDensity = densities_full[1]
# #Produce numerical solutions for a number of time steps using the logistic equation
# theoreticalDensities = map(t -> initialDensity/((1-initialDensity)*exp(-probProliferation*t) + initialDensity), times)
# #Graph these together
# timeChangePlot1 = scatter(times, densities_full, xlabel = "t", ylabel = "N", title = "P_m = 1, P_p = 0.1", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
# timeChangePlot1 = plot!(times, theoreticalDensities, lc = :black, label = "Exact", lw = 3)

# probProliferation = 0.01
# averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimesWithRandomIC([XLENGTH,YLENGTH], times, numSimulations, probMovement, probProliferation, initialDensity)
# densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
# densities_full = map(x -> sum(x)/length(x), densities_y) #Average results over x direction
# initialDensity = densities_full[1]
# #Produce numerical solutions for a number of time steps using the logistic equation
# theoreticalDensities = map(t -> initialDensity/((1-initialDensity)*exp(-probProliferation*t) + initialDensity), times)
# #Graph these together
# timeChangePlot2 = scatter(times, densities_full, xlabel = "t", ylabel = "N", title = "P_m = 1, P_p = 0.01", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
# timeChangePlot2 = plot!(times, theoreticalDensities, lc = :black, label = "Exact", lw = 3)

# probProliferation = 0.001
# averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimesWithRandomIC([XLENGTH,YLENGTH], times, numSimulations, probMovement, probProliferation, initialDensity)
# densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
# densities_full = map(x -> sum(x)/length(x), densities_y) #Average results over x direction
# initialDensity = densities_full[1]
# #Produce numerical solutions for a number of time steps using the logistic equation
# theoreticalDensities = map(t -> initialDensity/((1-initialDensity)*exp(-probProliferation*t) + initialDensity), times)
# #Graph these together
# timeChangePlot3 = scatter(times, densities_full, xlabel = "t", ylabel = "N", title = "P_m = 1, P_p = 0.001", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
# timeChangePlot3 = plot!(times, theoreticalDensities, lc = :black, label = "Exact", lw = 3)

# probMovement = 0.01
# probProliferation = 0.01
# averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimesWithRandomIC([XLENGTH,YLENGTH], times, numSimulations, probMovement, probProliferation, initialDensity)
# densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
# densities_full = map(x -> sum(x)/length(x), densities_y) #Average results over x direction
# initialDensity = densities_full[1]
# #Produce numerical solutions for a number of time steps using the logistic equation
# theoreticalDensities = map(t -> initialDensity/((1-initialDensity)*exp(-probProliferation*t) + initialDensity), times)
# #Graph these together
# timeChangePlot4 = scatter(times, densities_full, xlabel = "t", ylabel = "N", title = "P_m = 0.01, P_p = 0.01", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
# timeChangePlot4 = plot!(times, theoreticalDensities, lc = :black, label = "Exact", lw = 3)

# totalPlot = plot(timeChangePlot1,timeChangePlot2, timeChangePlot3, timeChangePlot4, layout = (1,4), size = (2000,500))
# display(totalPlot)


#||||----QUESTION 2----||||#

#Now set up an initial condition independent of the vertical position
numSimulations = 10
XLENGTH = 250
YLENGTH = 100
println("YLENGTH = ", YLENGTH)
STEPSIZE = 1
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)
times = [100]
H = 25
probMovement = 1

xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...]



simGrid = zeros(XLENGTH,YLENGTH)
simGrid = createBlock(0,2*H, YLENGTH, simGrid,xAxisValues)


probProliferation = 0.1
#Run simulation with proliferation
averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimes([XLENGTH,YLENGTH], times, simGrid, numSimulations, probMovement, probProliferation)
densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
#Solve the Fisher-Kolmogorov model numerically
C0 = calculateDensities(simGrid) 
numericSolutions = zeros(length(times),NUMBEROFSTEPS) 
theoreticSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS, probMovement, probProliferation], C0, times, FisherKolmogorov!)
#Plot together
plot1 = scatter(xAxisValues, densities_y, xlabel = "x", ylabel = "N", title = "P_m = 1, P_p = 0.1", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
plot1 = plot!(xAxisValues, theoreticSolutions[1,:], lc = :black, label = "Exact", lw = 3)



probProliferation = 0.1
#Run simulation with proliferation
averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimes([XLENGTH,YLENGTH], times, simGrid, numSimulations, probMovement, probProliferation)
densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
#Solve the Fisher-Kolmogorov model numerically
C0 = calculateDensities(simGrid) 
numericSolutions = zeros(length(times),NUMBEROFSTEPS) 
theoreticSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS, probMovement, probProliferation], C0, times, FisherKolmogorov!)
#Plot together
plot2 = scatter(xAxisValues, densities_y, xlabel = "x", ylabel = "N", title = "P_m = 1, P_p = 0.01", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
plot2 = plot!(xAxisValues, theoreticSolutions[1,:], lc = :black, label = "Exact", lw = 3)



probProliferation = 0.1
#Run simulation with proliferation
averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimes([XLENGTH,YLENGTH], times, simGrid, numSimulations, probMovement, probProliferation)
densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
#Solve the Fisher-Kolmogorov model numerically
C0 = calculateDensities(simGrid) 
numericSolutions = zeros(length(times),NUMBEROFSTEPS) 
theoreticSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS, probMovement, probProliferation], C0, times, FisherKolmogorov!)
#Plot together
plot3 = scatter(xAxisValues, densities_y, xlabel = "x", ylabel = "N", title = "P_m = 1, P_p = 0.001", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 1cm, label = "Stochastic")
plot3 = plot!(xAxisValues, theoreticSolutions[1,:], lc = :black, label = "Exact", lw = 3)


probMovement = 1
probProliferation = 1
#Run simulation with proliferation
averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimes([XLENGTH,YLENGTH], times, simGrid, numSimulations, probMovement, probProliferation)
densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
#Solve the Fisher-Kolmogorov model numerically
C0 = calculateDensities(simGrid) 
numericSolutions = zeros(length(times),NUMBEROFSTEPS) 
theoreticSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS, probMovement, probProliferation], C0, times, FisherKolmogorov!)
#Plot together
plot4 = scatter(xAxisValues, densities_y, xlabel = "x", ylabel = "N", title = "P_m = 0.01, P_p = 0.01", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
plot4 = plot!(xAxisValues, theoreticSolutions[1,:], lc = :black, label = "Exact", lw = 3)

totalPlot = plot(plot1,plot2, plot3, plot4, layout = (1,4), size = (2000,500))
display(totalPlot)


# #||||----QUESTION 2.b----||||#

# numSimulations = 10
# XLENGTH = 250
# YLENGTH = 100
# println("YLENGTH = ", YLENGTH)
# STEPSIZE = 1
# NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)
# times = [100, 500, 1000, 1500]
# H = 25
# probMovement = 1

# xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...]


# probMovement = 0.01
# probProliferation = 0.01
# #Run simulation with proliferation
# averagedSimulationGrids = StochasticExclusionWalkAverageWithProliferationMultTimes([XLENGTH,YLENGTH], times, simGrid, numSimulations, probMovement, probProliferation)
# densities_y = map(x -> calculateDensities(x), averagedSimulationGrids) #Average results over the y direction
# #Solve the Fisher-Kolmogorov model numerically
# C0 = calculateDensities(simGrid) 
# numericSolutions = zeros(length(times),NUMBEROFSTEPS) 
# theoreticSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS, probMovement, probProliferation], C0, times, FisherKolmogorov!)
# #Plot together
# plot1 = scatter(xAxisValues, densities_y[1], xlabel = "x", ylabel = "N", title = "t = 100", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
# plot1 = plot!(xAxisValues, theoreticSolutions[1,:], lc = :black, label = "Exact", lw = 3)

# plot2 = scatter(xAxisValues, densities_y[2], xlabel = "x", ylabel = "N", title = "t = 500", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
# plot2 = plot!(xAxisValues, theoreticSolutions[1,:], lc = :black, label = "Exact", lw = 3)

# plot3 = scatter(xAxisValues, densities_y[3], xlabel = "x", ylabel = "N", title = "t = 1000", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
# plot3 = plot!(xAxisValues, theoreticSolutions[1,:], lc = :black, label = "Exact", lw = 3)

# plot4 = scatter(xAxisValues, densities_y[4], xlabel = "x", ylabel = "N", title = "t = 1500", ylims = (0.0,1.0), size = (400,400), left_margin=1cm, right_margin=1cm, top_margin=1cm, bottom_margin = 0.5cm, label = "Stochastic")
# plot4 = plot!(xAxisValues, theoreticSolutions[1,:], lc = :black, label = "Exact", lw = 3)

# totalPlot = plot(plot1,plot2, plot3, plot4, layout = (1,4), size = (2000,500))
# display(totalPlot)