#Created for a later course called 'Medical Mathematics'

using Plots, SpecialFunctions, Random, DifferentialEquations
include("functions.jl")

#|||---VARIABLES----||||#
XLENGTH = 20 #Width of grid
YLENGTH = 20 #Height of grid
lengths = [XLENGTH, YLENGTH]
initialDensity = 0.5 #Density of 

numSimulations = 2
probMovement = 1
probProliferation = 1
BIAS = 0
probDeath = 0.7

foodDensity = 0.8

MAXTIME = 100

STEPSIZE = 1
NUMBEROFSTEPS = Int(XLENGTH/STEPSIZE)+1

#||||----Functions----||||#

#In this model we have the  normal model for movement, but, instead of each cell having a certain probability of proliferating, 
#we have that each node on the graph has a certain probability of causing any present cell to proliferate
#(modelling a restricted availability of food) and each cell also has a certain probability of dying each time
function ResourceAndSpaceLimitedProliferationWithDeath(lengths, MAXTIME, numSimulations, probMovement, probProliferation, initialDensity, probDeath, foodDensity)
    
    sumGrids = fill(zeros(lengths...),1:MAXTIME) #The result of each simulation will be added to this variable so it can be averaged later
    maxTime = MAXTIME #Remnant of the code I changed to make this



    for sim in 1:numSimulations
        simGrid = zeros(lengths...)
        randCount = 0
        simGrid = zeros(XLENGTH,YLENGTH)

        #Prepare Grid
        for i in eachindex(simGrid)
            occupiedCheck = rand(1)[1]
            if(occupiedCheck < initialDensity && simGrid[i] != 1.0)
                randCount += 1
                simGrid[i] = 1.0
                

            end#if
        end#for
        totalAgents = sum(simGrid)


        println("Simulation: ",sim)#Print the current number of simulations

        tempGrid = copy(simGrid)

        for t in 1:maxTime
            
            local count = 0
            
            println("Moving:")

            #Move the agents
            while count < totalAgents #choose random cell to try to find an agent to move
               randCoordinates = []
               for i in eachindex(lengths)
                    append!(randCoordinates,rand(1:lengths[i]))
               end#for

               #Check if the cell contains an agent
               if(tempGrid[randCoordinates...]==1)
                    count+=1
                    #Does it move?
                    moveQuery = rand(1)[1]
                    if(moveQuery <= probMovement)
                        #decide on direction. 1: up 2: right 3:down 4:left
                        moveDirection = rand(1:4)
                        #Move the agent
                        tempGrid = attemptActionWithDirection(moveDirection,tempGrid,randCoordinates, moveAgent)
                    end#if
                end#if
            end#while
            
            println("Proliferating")
            println(totalAgents)

            #proliferate
            tempTotalAgents = totalAgents
            count = 0
            while count < round(foodDensity*lengths[1]*lengths[2]) #choose random cell to try to find an agent to proliferate
                randCoordinates = []
                for i in eachindex(lengths)
                     append!(randCoordinates,rand(1:lengths[i]))
                end#for
 
                #Check if the cell contains an agent
                if(tempGrid[randCoordinates...]==1)
                     
                     
                     #Does it proliferate?
                     proliferateQuery = rand(1)[1]
                     if(proliferateQuery <= probProliferation)
                         #decide on direction. 1: up 2: right 3:down 4:left
                         moveDirection = rand(1:4)
                     
                         #Proliferate
                         tempGrid = attemptActionWithDirection(moveDirection,tempGrid,randCoordinates, proliferate)
                         tempTotalAgents += 1
                     end#if
                 end#if
                 count += 1
             end#while
             totalAgents = tempTotalAgents
             println(totalAgents)



             println("Dying")
             count=0
             #death
             while count < totalAgents #choose random cell to try to find an agent to move
                randCoordinates = []
                for i in eachindex(lengths)
                     append!(randCoordinates,rand(1:lengths[i]))
                end#for
 
                #Check if the cell contains an agent
                if(tempGrid[randCoordinates...]==1)
                    count += 1
                     
                     #Does it die?
                     deathQuery = rand(1)[1]
                     if(deathQuery <= probDeath)
                         tempGrid[randCoordinates...] = 0;
                         tempTotalAgents -= 1
                     end#if
                 end#if
                 
             end#while
            totalAgents = tempTotalAgents
            println(totalAgents)
            
            sumGrids[t] +=tempGrid
        end#for
    end#for
    
    #Average placement of agents across the grid throughout all simulations
    averageGrid = sumGrids/numSimulations
    

    return averageGrid #averageGrid[time] gives the averaged grid for the time 'time' for all simulations. Note, it returns a grid

end#function

function calculateFullDensities(simGrid)
    lengths=[size(simGrid)...]#make array with the dimensions of the grid
    #----Relevant to a particular problem. Returning just averageGrid may be preferable
    #calculate density along vertical lines
    density=zeros(lengths[1])
    for i in 1:lengths[1]
        for j in 1:size(simGrid[1])[1]*size(simGrid[1])[2]
            #for k in 1:size(simGrid[1])[2]
                density[i]+=simGrid[i][j]
            #end#for
        end#for
    end#for

    #Average number of agents along each vertical line on the grid, across all simulations
    averageDensity = density/(size(simGrid[1])[1]*size(simGrid[1])[2])

    return averageDensity
end#function




#||||----Calculations----||||#
simGrids = ResourceAndSpaceLimitedProliferationWithDeath(lengths, MAXTIME, numSimulations, probMovement, probProliferation, initialDensity, probDeath, foodDensity)
densities = calculateFullDensities(simGrids)
plot(densities)