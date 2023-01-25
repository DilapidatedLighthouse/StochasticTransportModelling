using DifferentialEquations
include("functions.jl")

#||||----Variables----||||#
XLENGTH = 500
YLENGTH = 50
Times = [0,1000,2000,3000]
STEPSIZE = 0.50
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)+1
x=LinRange(-XLENGTH/2,XLENGTH/2,NUMBEROFSTEPS)
NUMSIMULATIONS = 10
VELOCITY = 0.05


#Associated with exact solutions
DIFFUSIVITY=0.25
D=DIFFUSIVITY

#Associated with initial conditions
H=20
C0=zeros(NUMBEROFSTEPS)



numericSolutions = zeros(length(Times),NUMBEROFSTEPS) 


for i in 1:NUMBEROFSTEPS
    if abs(x[i]) <= H
        C0[i]=1
    end#if
end#for

xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...]
simGrid = zeros(XLENGTH,YLENGTH)
simGrid = createBlock(0,2*H, YLENGTH, simGrid,xAxisValues)


#Numeric solution of the nonlinear advection diffusion PDE
numericSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS,D, VELOCITY],C0, Times, AdvectionDiffusionFunction!)

#Biased migration with diffusion
densities = calculateDensities(StochasticExclusionWalkAverage([XLENGTH,YLENGTH],Times[2],simGrid,NUMSIMULATIONS,1,VELOCITY))

p1 = scatter(xAxisValues,densities,mc=:black,msc=:match,label="Stochastic")

p1=plot!(x,numericSolutions[2,:],lw=2,xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1.2), ylabel="C(x,t)",legend=false)
# p1!=plot(x,densities,lw=2,xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1.2), ylabel="C(x,t)",legend=false)
display(p1)
