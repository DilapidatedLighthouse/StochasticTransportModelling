using Plots, SpecialFunctions, Random, DifferentialEquations
include("functions.jl")





#Grid parameters
XLENGTH = 500
YLENGTH = 50
lengths = [XLENGTH, YLENGTH]
xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...]
NUMBEROFSIMULATIONS = 10
STEPSIZE = 1 

times = [0,2]
NUMBEROFSIMULATIONS = 1
biases = 0 #not properly implemented

probMovements = [1,1]

#Set up initial conditions

#population 1
H = 25
simGrid1 = zeros(XLENGTH,YLENGTH)

simGrid1 = createBlockEven(-100,2*H, 100000, simGrid1,xAxisValues)
#population 2
H = 25
simGrid2 = zeros(XLENGTH,YLENGTH)

simGrid2 = createBlockEven(100,2*H, YLENGTH, simGrid2,xAxisValues)

simGrids = [simGrid1, simGrid2]

#Run simulation for two populations. Here we run each population on its own grid. The function just makes sense on all simgGrids. 
function StochasticExclusionWalkAverageMultTimesMultiplePopulations(lengths, times, simGrids, numSimultaions, probMovements, biases)
    totalAgents = sum(sum(simGrids)) #We sum over each simGrid and then over all simgrids
    sumGrids = map(x -> fill(zeros(lengths...),length(times)), simGrids)#The result of each simulation will be added to this variable so it can be averaged later
    maxTime = maximum(times)

    for sim in 1:numSimultaions
        println("Simulation: ",sim)#Print the current number of simulations

        tempGrids = copy(simGrids)

        for t in 0:maxTime
            
            count=0
            
            #Move the agents
            while count < totalAgents #choose random cell to try to find an agent to move
               randCoordinates = []
               for i in eachindex(lengths)
                    append!(randCoordinates,rand(1:lengths[i]))
               end#for

               #Check if the cell contains an agentc
               occupiedIndex = CheckCellEmpty(randCoordinates, simGrids)
               if(occupiedIndex != 0)
                    count+=1
                    #Does it move?
                    moveQuery = rand(1)[1]
                    if(moveQuery <= probMovements[occupiedIndex])
                        moveDirection = chooseDirectionWithXBias(biases)

                        #Move the agent
                        tempGrids = attemptActionWithDirectionMultiPopulation(moveDirection,tempGrids, occupiedIndex, randCoordinates, moveAgent)
                    end#if
                end#if
            end#while
            if(t in times)
                for i in eachindex(probMovements) #dumb way of getting number of populations
                    # println("sumGrids, size = ", size(sumGrids),": ", typeof(sumGrids))
                    # println("tempGrids, size = ", size(tempGrids),": ", typeof(tempGrids))

                    sumGrids[i][findlast(time -> time==t, times)]+=tempGrids[i]
                end
            end#if

        end#for
        
        # #Sum the results of all simulations
        
        # sumGrid+=tempGrid
        
    end#for
    
    #Average placement of agents across the grid throughout all simulations
    averageGrid = sumGrids/numSimultaions
    

    return averageGrid

end#function

function attemptActionWithDirectionMultiPopulation(moveDirection, simGrids, actingIndex, coordinates,action)
    moveVector = zeros(length(coordinates)) #Used to move in the right direction

    #set move vector based on input direction
    moveVector = [Int(round(cos(pi/2*(2-moveDirection)))),Int(round(sin(pi/2*(2-moveDirection))))]



    #Point to move to. Converts each element into an integer so it doesn't break when accessing the index.
    moveCoord = [Int((coordinates+moveVector)[i]) for i in eachindex(coordinates+moveVector)]
    #Point to move to if movement reflected. Used for boundary conditions. Converts each element into an integer so it doesn't break when accessing the index.
    oppositeMoveCoord = [Int((coordinates-moveVector)[i]) for i in eachindex(coordinates-moveVector)]
    #converts coordinates to integers
    intCoordinates = [Int(coordinates[i]) for i in eachindex(coordinates)]

    #Check to see if move is within bounds
    if(checkInBoundaries(simGrids[actingIndex],moveCoord))

        #Check not moving to an occupied space
        if(CheckCellEmpty(moveCoord, simGrids)!=0)
            simGrids[actingIndex] =action(simGrids[actingIndex],intCoordinates,moveCoord)
        end#if
    
    #Otherwise bounce off the boundary (naive. just a move in the other direction)
    elseif(checkInBoundaries(simGrids[actingIndex],oppositeMoveCoord))
        #Check not moving to an occupied space
        if(CheckCellEmpty(oppositeMoveCoord, simGrids)!=0)
            #=
            simGrid[(oppositeMoveCoord)...]=1.0
            simGrid[intCoordinates...] = 0.0
            =#
            simGrids[actingIndex] = action(simGrids[actingIndex],intCoordinates,oppositeMoveCoord)
        end#if
    end#if elseif
    return simGrids
end#Function

#if one of the grids contain a cell at the specified coordinates, the index of that grid will be returned. Else the function will return 0.
function CheckCellEmpty(coordinates, simGrids)
    for gridIndex in eachindex(simGrids)
        if(simGrids[gridIndex][coordinates...] == 1)
            return gridIndex
        end#if
    end#for
    return 0
end#Function

simulations = StochasticExclusionWalkAverageMultTimesMultiplePopulations(lengths, times, simGrids, NUMBEROFSIMULATIONS, probMovements, biases)
densities = map(sim -> map(pop -> calculateDensities(pop), sim), simulations)

scatter(xAxisValues, densities[1][1])
#Numerical solution for two populations

#Plot the two together