#TO DO:
#   *Add option for bias to movement


using Plots, SpecialFunctions, Random
include("functions.jl")
gr() #A graphical thing? I need to look it up



#||||----VARIABLES----||||#

#Declaring variables for simulation
XLENGTH = 500
YLENGTH = 100
PROBABILTYOFMOVE = 1.0
PROBABILTYOFPROLIFERATION = 1.0
TIMEOFSIMULATION = 0
NUMBEROFSIMULATIONS = 1

BIAS = -0.125

#Constants for initial conditions
H = 50

#Constants for heat equation
D = PROBABILTYOFMOVE/4

#||||----CODE----||||#

#Setting up x-axis for graphing and grid initialisation
xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...]

#Initialise grid with agents
simGrid = zeros(XLENGTH,YLENGTH)

simGrid = createBlock(0,2*H, YLENGTH, simGrid,xAxisValues)
center=-150
#simGrid = createBlockEven(center,2*H, YLENGTH, simGrid,xAxisValues)
#Heat equation. Will have to be manually changed to match the initial conditions and simulation behaviour
T=TIMEOFSIMULATION
C(x)=0.25*(erf((H-x)/sqrt(4*D*T))+erf((H+x)/sqrt(4*D*T)));
#C(x)= 0.5*(erf((H-(x-center))/sqrt(4*D*T))+erf((H+(x-center))/sqrt(4*D*T)));
#C(x)=0.5*(erf((H-x)/sqrt(4*D*T))+erf((H+x)/sqrt(4*D*T))) + 0.25*(erf((H-(x-center))/sqrt(4*D*T))+erf((H+(x-center))/sqrt(4*D*T)));

#Calculate average densities as function of space.
densities = calculateDensities(StochasticExclusionWalkAverage([XLENGTH, YLENGTH], TIMEOFSIMULATION, simGrid, NUMBEROFSIMULATIONS, PROBABILTYOFMOVE, BIAS))
println("Calculations Completed")
#||||----PLOTS----||||#
myPlot = scatter(xAxisValues,densities,mc=:black,msc=:match,label="Stochastic")
myPlot=plot!(C, -XLENGTH/2:XLENGTH/2-1,lw=4,lc=:green,ls=:dash,label="Exact",xlabel="x",ylabel="C(x, 0)",xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1),framestyle=:box)
display(myPlot)
