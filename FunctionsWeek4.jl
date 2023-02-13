#using Base.DeepCopy
#if one of the grids contain a cell at the specified coordinates, the index of that grid will be returned. Else the function will return 0.
function CheckCellEmpty(coordinates, simGridsTemp3)
    for gridIndex in eachindex(simGridsTemp3)
        if(Int(simGridsTemp3[gridIndex][coordinates...]) == 1)
            
            return gridIndex
        end#if
    end#for
    return 0
end#Function

function moveAgentUpdate(simGridTemp5, intCoordinates, moveCoord)
    
    simGridTemp5[(moveCoord)...]=1.0
    simGridTemp5[intCoordinates...] = 0.0
    
    return simGridTemp5
end#function


function attemptActionWithDirectionMultiPopulation(moveDirection, simGridsTemp2, actingIndex, coordinates,action)
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
    if(checkInBoundaries(simGridsTemp2[actingIndex],moveCoord))

        #Check not moving to an occupied space
        if(CheckCellEmpty(moveCoord, simGridsTemp2)==0)
            
            simGridsTemp2[actingIndex] =action(simGridsTemp2[actingIndex],intCoordinates,moveCoord)
            
        end#if
    
    #Otherwise bounce off the boundary (naive. just a move in the other direction)
    elseif(checkInBoundaries(simGridsTemp2[actingIndex],oppositeMoveCoord))
        #Check not moving to an occupied space
        if(CheckCellEmpty(oppositeMoveCoord, simGridsTemp2)==0)
            #=
            simGrid[(oppositeMoveCoord)...]=1.0
            simGrid[intCoordinates...] = 0.0
            =#
            simGridsTemp2[actingIndex] = action(simGridsTemp2[actingIndex],intCoordinates,oppositeMoveCoord)
        end#if
    end#if elseif
    return simGridsTemp2
end#Function



#Run simulation for two populations. Here we run each population on its own grid. The function just makes sense on all simgGrids. 
function StochasticExclusionWalkAverageMultTimesMultiplePopulations(lengths, times, simGrids, numSimultaions, probMovements, biases)
    
    simGridsTemp = deepcopy(simGrids)

    totalAgents = sum(sum(simGridsTemp)) #We sum over each simGrid and then over all simgrids
    sumGrids = map(x -> fill(zeros(lengths...),length(times)), simGridsTemp)#The result of each simulation will be added to this variable so it can be averaged later
    maxTime = maximum(times)
    
    for sim in 1:numSimultaions
        println("Simulation: ",sim)#Print the current number of simulations

        tempGrids = deepcopy(simGridsTemp)
        display(heatmap(tempGrids[1]))

        for t in 0:maxTime
            
            count=0
            if(t in times)
                
                for i in eachindex(probMovements) #dumb way of getting number of populations
                    # println("sumGrids, size = ", size(sumGrids),": ", typeof(sumGrids))
                    # println("tempGrids, size = ", size(tempGrids),": ", typeof(tempGrids))

                    sumGrids[i][findlast(time -> time==t, times)]+=deepcopy(tempGrids[i])
                end
            end#if

            #Move the agents
            while count < totalAgents #choose random cell to try to find an agent to move
               randCoordinates = []
               for i in eachindex(lengths)
                    append!(randCoordinates,rand(1:lengths[i]))
               end#for

               #Check if the cell contains an agent
               occupiedIndex = CheckCellEmpty(randCoordinates, tempGrids)
               if(occupiedIndex != 0)
                
                    count+=1
                    #Does it move?
                    moveQuery = rand(1)[1]
                    if(moveQuery <= probMovements[occupiedIndex])
                        moveDirection = chooseDirectionWithXBias(biases)
                        
                        #Move the agent
                        
                        tempGrids = attemptActionWithDirectionMultiPopulation(moveDirection,tempGrids, occupiedIndex, randCoordinates, moveAgentUpdate)
                        
                    end#if
                end#if
            end#while
            
            
        
        end#for
        display(heatmap(tempGrids[1]))
        # #Sum the results of all simulations
        
        # sumGrid+=tempGrid
        
    end#for
    
    #Average placement of agents across the grid throughout all simulations
    averageGrid = sumGrids/numSimultaions
    

    return averageGrid

end#function

function diffusionTwoPopulations!(du,u,p,t)
    dx,N,D1,D2=p
    S = u[1,:] + u[2,:]
    for i in 2:N-1
        du[1,i]=D1*((1-S[i])*(u[1,i+1]-2*u[1,i]+u[1,i-1])/dx^2 + (u[1,i]*(S[i+1]-2*S[i]+S[i-1])/dx^2))
        du[2,i]=D2*((1-S[i])*(u[2,i+1]-2*u[2,i]+u[2,i-1])/dx^2 + (u[2,i]*(S[i+1]-2*S[i]+S[i-1])/dx^2))
    end#for
    du[1,1]=0.0
    du[1,N]=0.0
    du[2,1]=0.0
    du[2,N]=0.0
end#function


#Access these solutions with storageVariable[population, x-value, time].
#L is the x-length, dx is the stepsize, N is the number of steps, t is an array of times to store the solution at
#C0 is an array of arrays storing the intial conditions, D1 is the diffusivity for the 1st population, D2 is the diffusivity for the 2nd population
function pdesolverWeek4(L,dx,N,T,C0,D1,D2)
    p=(dx,N,D1,D2)
    tspan=(0.0,maximum(T))
    prob=ODEProblem(diffusionTwoPopulations!,C0,tspan,p)
    sol=solve(prob,saveat=T);
    return sol
end#function