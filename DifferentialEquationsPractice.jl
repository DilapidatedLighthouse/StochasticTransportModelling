using Plots, DifferentialEquations, SpecialFunctions

include("functions.jl")

#||||----Variables----||||#
XLENGTH = 500
Times = [0,1000,2000,3000]
STEPSIZE = 0.50
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)+1
x=LinRange(-XLENGTH/2,XLENGTH/2,NUMBEROFSTEPS)

VELOCITY = 0.05


#Associated with exact solutions
DIFFUSIVITY=0.25
D=DIFFUSIVITY

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