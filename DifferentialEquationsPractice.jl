using Plots, DifferentialEquations, SpecialFunctions

include("functions.jl")
#=
#||||----PDE Solver----||||#


#Solves a DE problem and generates a matrix of solution values with the time on one access and the x position on the other
#Requires the existence of a list of zeros "numericSolutions" in the scope in which the function is called.
#Paramters is an array of variables according to the needs of the ODEFunction
#C0 is the initial conditions for the DE
#ODEFunction is a function that defines the DE to be solved.
#Note: ODEFunction will not be called directly but ny ODEProblem from the DifferentialEquations package
function PDESolver(parameters, C0, Times,ODEFunction)
    timeSpan = (0.0, maximum(Times))
    problem = ODEProblem(ODEFunction,C0,timeSpan,parameters) #create a 'problem' for use with DifferentialEquations package
    
    solution = solve(problem,saveat=Times)
    
    for i in 1:length(solution[:,])
        numericSolutions[i,:]=solution[:,i]
    end#for

    return numericSolutions
    
end#function


#A possible ODEFunction for PDESolver. Diffusion equation
#Should have parameters = [step-size, Number of steps(along x direction), Diffusion constant] 
function DiffusionFunction!(du, u, parameters, t)
    dx,N,D = parameters #unpacking into variables
    #Setting up numeric scheme for other derivatives
    N = Int(N) #For some reason it will implicitly convert N to a float64... not sure why
    for i in 2:N-1
        du[i] = D*(u[i+1]-2*u[i]+u[i-1])/dx^2
    end#for

    du[1]=0.0
    du[N]=0.0
end#function


#A possible ODEFunction for PDESolver. Diffusion and advection equation
#Should have parameters = [step-size, Number of steps(along x direction), Diffusion constant, velocity] 
function AdvectionDiffusionFunction!(du, u, parameters, t)
    dx,N,D,v = parameters #unpacking into variables
    N = Int(N)
    #Setting up numeric scheme for other derivatives
    for i in 2:N-1
        du[i] = D*(u[i+1]-2*u[i]+u[i-1])/dx^2 - (v/(2*dx))*(u[i+1]-u[i-1])
    end#for

    du[1]=0.0
    du[N]=0.0
end#function
=#

#||||----Variables----||||#
XLENGTH = 500
Times = [0,1000,2000,3000]
STEPSIZE = 0.50
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)+1
x=LinRange(-XLENGTH/2,XLENGTH/2,NUMBEROFSTEPS)

VELOCITY = 0.05


#Associated with exact solutions
DIFFUSIVITY=0.25

#Associated with initial conditions
h=20
C0=zeros(NUMBEROFSTEPS)



numericSolutions = zeros(length(Times),NUMBEROFSTEPS) 


for i in 1:NUMBEROFSTEPS
    if abs(x[i]) <= h
        C0[i]=1
    end#if
end#for



#Diffusion
#=
C1(x)=0.5*(erf((h-x)/sqrt(4*D*Times[1]))+erf((h+x)/sqrt(4*D*Times[1])));
C2(x)=0.5*(erf((h-x)/sqrt(4*D*Times[2]))+erf((h+x)/sqrt(4*D*Times[2])));
C3(x)=0.5*(erf((h-x)/sqrt(4*D*Times[3]))+erf((h+x)/sqrt(4*D*Times[3])));
C4(x)=0.5*(erf((h-x)/sqrt(4*D*Times[4]))+erf((h+x)/sqrt(4*D*Times[4])));

numericSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS,DIFFUSIVITY], C0, Times,DiffusionFunction!)
=#


#Advection

C1(x)=0.5*(erf((h-(x-VELOCITY*Times[1]))/sqrt(4*D*Times[1]))+erf((h+x-VELOCITY*Times[1])/sqrt(4*D*Times[1])));
C2(x)=0.5*(erf((h-x+VELOCITY*Times[2])/sqrt(4*D*Times[2]))+erf((h+x-VELOCITY*Times[2])/sqrt(4*D*Times[2])));
C3(x)=0.5*(erf((h-x+VELOCITY*Times[3])/sqrt(4*D*Times[3]))+erf((h+x-VELOCITY*Times[3])/sqrt(4*D*Times[3])));
C4(x)=0.5*(erf((h-x+VELOCITY*Times[4])/sqrt(4*D*Times[4]))+erf((h+x-VELOCITY*Times[4])/sqrt(4*D*Times[4])));

numericSolutions = PDESolver([STEPSIZE, NUMBEROFSTEPS,DIFFUSIVITY,VELOCITY], C0, Times,AdvectionDiffusionFunction!)


p1=plot(x,numericSolutions[1,:],lw=2,xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1.2), ylabel="C(x,t)",legend=false)
p1=plot!(C1,-XLENGTH/2,XLENGTH/2,lw=2,ls=:dash,lc=:red)
p2=plot(x,numericSolutions[2,:],lw=2,xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1.2), ylabel="C(x,t)",legend=false)
p2=plot!(C2,-XLENGTH/2,XLENGTH/2,lw=2,ls=:dash,lc=:red)
p3=plot(x,numericSolutions[3,:],lw=2,xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1.2), ylabel="C(x,t)",legend=false)
p3=plot!(C3,-XLENGTH/2,XLENGTH/2,lw=2,ls=:dash,lc=:red)
p4=plot(x,numericSolutions[4,:],lw=2,xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1.2), ylabel="C(x,t)",legend=false)
p4=plot!(C4,-XLENGTH/2,XLENGTH/2,lw=2,ls=:dash,lc=:red,xlabel="x")
p5=plot(p1,p2,p3,p4,layout=(4,1))
display(p5)