using DifferentialEquations
include("functions.jl")

#||||----Variables----||||#
XLENGTH = 500
YLENGTH = 50
Times = [1, 150, 200, 250]
STEPSIZE = 1

NUMSIMULATIONS = 10
VELOCITY = 0.4


#Associated with exact solutions
DIFFUSIVITY=0.25
D=DIFFUSIVITY

#Associated with initial conditions
H=20

xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...]

simGrid = zeros(XLENGTH,YLENGTH)
simGrid = createBlock(0,2*H, YLENGTH, simGrid,xAxisValues)




#Biased migration with diffusion

#p1 = scatter(xAxisValues,densities,mc=:black,msc=:match,label="Stochastic")

#Numeric solution of the nonlinear advection diffusion PDE

# # p1!=plot(x,densities,lw=2,xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1.2), ylabel="C(x,t)",legend=false)



function CalcNumericSolution(STEPSIZE, Times, D, VELOCITY,XLENGTH)
    NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)+1
    println("NUMBER OF STEPS = ", NUMBEROFSTEPS)
    x=LinRange(-XLENGTH/2,XLENGTH/2,NUMBEROFSTEPS)

    
    C0=zeros(NUMBEROFSTEPS)



     
    for i in 1:NUMBEROFSTEPS
        if abs(x[i]) <= H
            C0[i]=1
        end#if
    end#for

    

    
    numericSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS,D, VELOCITY],C0, Times, AdvectionDiffusionFunction!)

    
    return [x,numericSolutions]
    
end#function
#p1=plot(x,numericSolutions[2,:],lw=2,xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1.2), ylabel="C(x,t)",legend=false)

# x, numericSolutions = CalcNumericSolution(STEPSIZE/2,Times, D, VELOCITY*1.25, XLENGTH)

# plotArray = []
# for i in 1:size(numericSolutions)[1]
#     append!(plotArray, plot(x,numericSolutions[i,:],lw=2,xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1.2), ylabel="C(x,t)",legend=false))
# end#for




# #||--Calculating Differences--||#
# stepSizes = [0.1, 0.2, 0.5]
# refNumericSolutions = CalcNumericSolution(1, [100], D, VELOCITY, XLENGTH)[2]
# errors = []
# for i in eachindex(stepSizes)
#     sols = CalcNumericSolution(stepSizes[i], [100], D, VELOCITY, XLENGTH)[2]
#     tempErrors = zeros(length(refNumericSolutions))
#     for j in eachindex(refNumericSolutions)
#         tempErrors[j] = refNumericSolutions[j] - sols[1+(j-1)*Int(1/stepSizes[i])] 
#     end#for
#     append!(errors, [tempErrors])
# end#for

# p2 = plot(-XLENGTH/2:1:XLENGTH/2, errors[2], lc=:black, ylims=(-0.1,0.1),framestyle=:box, legend = false, title = "Step-size 1 vs 0.2", xlabel = "x", ylabel = "N(x,100)")
# p3 = plot(-XLENGTH/2:1:XLENGTH/2, errors[3], lc=:black, ylims=(-0.1,0.1),framestyle=:box, legend = false, title = "Step-size 1 vs 0.5", xlabel = "x", ylabel = "N(x,100)")
# p1 = plot(-XLENGTH/2:1:XLENGTH/2, errors[1], lc=:black, ylims=(-0.1,0.1),framestyle=:box, legend = false, title = "Step-size 1 vs 0.1", xlabel = "x", ylabel = "N(x,100)")

# fullPlot = plot(p3,p2,p1, layout = (3,1), titlefontsize=10)

#||--Plotting with Stochastic--||#
stochasticSolutions = StochasticExclusionWalkAverageMultTimes([XLENGTH,YLENGTH],Times,simGrid,NUMSIMULATIONS,1,VELOCITY)
# densities = map(x -> calculateDensities(x), stochasticSolutions)
densities = []
for i in stochasticSolutions
    append!(densities, [calculateDensities(i)])
end#for


numericSolutions = CalcNumericSolution(STEPSIZE, Times, D, VELOCITY/2, XLENGTH)

p1 = scatter(-XLENGTH/2:XLENGTH/2,densities[1], label="Stochastic")
p1 = plot!(numericSolutions[1], numericSolutions[2][1,:],lw=2,lc=:black,ls=:dash,label="Numeric",xlabel="x",ylabel="C(x, 1)",xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1),framestyle=:box)
p2 = scatter(-XLENGTH/2:XLENGTH/2,densities[2], label="Stochastic")
p2 = plot!(numericSolutions[1], numericSolutions[2][2,:],lw=2,lc=:black,ls=:dash,label="Numeric",xlabel="x",ylabel="C(x, 50)",xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1),framestyle=:box)
p3 = scatter(-XLENGTH/2:XLENGTH/2,densities[3], label="Stochastic")
p3 = plot!(numericSolutions[1], numericSolutions[2][3,:],lw=2,lc=:black,ls=:dash,label="Numeric",xlabel="x",ylabel="C(x, 100)",xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1),framestyle=:box)
p4 = scatter(-XLENGTH/2:XLENGTH/2,densities[4], label="Stochastic")
p4 = plot!(numericSolutions[1], numericSolutions[2][4,:],lw=2,lc=:black,ls=:dash,label="Numeric",xlabel="x",ylabel="C(x, 150)",xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1),framestyle=:box)

fullPlot = plot(p1,p2,p3,p4)

#scatter(calculateDensities(tempstochasticSolutions[1]))