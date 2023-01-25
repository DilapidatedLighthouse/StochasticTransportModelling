
using Plots, SpecialFunctions, Random

#||||----GENERAL----||||#

#Can be plotted with scatter(gridPoints[1],gridPoints[2])
function PrepGridForDisplay(simGrid)
    gridPoints = [[],[]]
    for i in axes(simGrid,1)
            for j in axes(simGrid, 2)
                if(simGrid[i,j] == 1.0)
                    append!(gridPoints[1],i)
                    append!(gridPoints[2],j)
                end#if
            end#for
    end#for
    return gridPoints
end#Function


#||||----STOCHASTIC PROCESSES----||||#

#Meat and bones. Runs simulations of cells moving on a grid and averages the results.
#Lengths is an array of length n specifying the dimensions of the grid.
#totalTime is the number of time-steps for each simulation
#numSimulations is the number of simulations
#probMovement is the probabilty that an agent will try to move if given the chance 
#the nullVariable is so the proliferation probability doesnt have to be removed when switching between functions
function StochasticExclusionWalkAverage(lengths, totalTime, simGrid, numSimultaions, probMovement, bias = 0, nullVariable = "")
    totalAgents = sum(simGrid)
    integerTime = Int(totalTime)
    sumGrid = zeros(lengths...)#The result of each simulation will be added to this variable so it can be averaged later
    

    for sim in 1:numSimultaions
        println("Simulation: ",sim)#Print the current number of simulations

        tempGrid = copy(simGrid)

        for t in 1:integerTime
            count=0
            
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
                        # #decide on direction. 1: up 2: right 3:down 4:left
                        # randMoveIndex = rand(1:length(biases))
                        # moveDirection = biases[randMoveIndex]
                        moveDirection = chooseDirectionWithXBias(bias)

                        #Move the agent
                        tempGrid = attemptActionWithDirection(moveDirection,tempGrid,randCoordinates, moveAgent)
                    end#if
                end#if
            end#while
        end#for
        
        #Sum the results of all simulations
        sumGrid+=tempGrid
        
    end#for
    
    #Average placement of agents across the grid throughout all simulations
    averageGrid = sumGrid/numSimultaions

    return averageGrid

    # #----Relevant to a particular problem. Returning just averageGrid may be preferable
    # #calculate density along vertical lines
    # for i in 1:lengths[1]
    #     for j in 1:lengths[2]
    #         density[i]+=averageGrid[i,j]
    #     end#for
    # end#for

    # #Average number of agents along each vertical line on the grid, across all simulations
    # averageDensity = density/lengths[2]
    # return averageDensity

end#function


function chooseDirectionWithXBias(Bias)
    randomVar = rand(1)[1]
    direction = (randomVar < 0.25) ? 1 : (randomVar >= 0.25 && randomVar < 0.5 + Bias/4) ? 2 : (randomVar >= 0.5 + Bias/4 && randomVar < 0.75) ? 4 : 3
    return direction
end#function


#ensure all times in 'times' array are integers
function StochasticExclusionWalkAverageWithProliferationMultTimes(lengths, times, simGrid, numSimultaions, probMovement, probProliferation)
    totalAgents = sum(simGrid)
    sumGrids = fill(zeros(lengths...),length(times))#The result of each simulation will be added to this variable so it can be averaged later
    maxTime = maximum(times)

    for sim in 1:numSimultaions
        println("Simulation: ",sim)#Print the current number of simulations

        tempGrid = copy(simGrid)

        for t in 1:maxTime
            
            count=0
            
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

            #proliferate
            tempTotalAgents = totalAgents
            count = 0
            while count < totalAgents #choose random cell to try to find an agent to move
                randCoordinates = []
                for i in eachindex(lengths)
                     append!(randCoordinates,rand(1:lengths[i]))
                end#for
 
                #Check if the cell contains an agent
                if(tempGrid[randCoordinates...]==1)
                     count += 1
                     tempTotalAgents += 1
                     #Does it proliferate?
                     proliferateQuery = rand(1)[1]
                     if(proliferateQuery <= probProliferation)
                         #decide on direction. 1: up 2: right 3:down 4:left
                         moveDirection = rand(1:4)
                     
                         #Proliferate
                         tempGrid = attemptActionWithDirection(moveDirection,tempGrid,randCoordinates, proliferate)
                     end#if
                 end#if
             end#while

            #Making sumGrid
            if(t in times)
                sumGrids[findlast(time -> time==t, times)]+=tempGrid
            end#if

        end#for
        
        # #Sum the results of all simulations
        
        # sumGrid+=tempGrid
        
    end#for
    
    #Average placement of agents across the grid throughout all simulations
    averageGrid = sumGrids/numSimultaions
    

    return averageGrid

end#function

function StochasticExclusionWalkAverageWithProliferation(lengths, totalTime, simGrid, numSimultaions, probMovement, probProliferation)
    totalAgents = sum(simGrid)
    integerTime = Int(totalTime)
    sumGrid = zeros(lengths...)#The result of each simulation will be added to this variable so it can be averaged later
    

    for sim in 1:numSimultaions
        println("Simulation: ",sim)#Print the current number of simulations

        tempGrid = copy(simGrid)

        for t in 1:integerTime
            count=0
            
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

            #proliferate
            tempTotalAgents = totalAgents
            count = 0
            while count < totalAgents #choose random cell to try to find an agent to move
                randCoordinates = []
                for i in eachindex(lengths)
                     append!(randCoordinates,rand(1:lengths[i]))
                end#for
 
                #Check if the cell contains an agent
                if(tempGrid[randCoordinates...]==1)
                     count += 1
                     tempTotalAgents += 1
                     #Does it proliferate?
                     proliferateQuery = rand(1)[1]
                     if(proliferateQuery <= probProliferation)
                         #decide on direction. 1: up 2: right 3:down 4:left
                         moveDirection = rand(1:4)
                     
                         #Proliferate
                         tempGrid = attemptActionWithDirection(moveDirection,tempGrid,randCoordinates, proliferate)
                     end#if
                 end#if
             end#while
        end#for
        
        #Sum the results of all simulations
        
        sumGrid+=tempGrid
        
    end#for
    
    #Average placement of agents across the grid throughout all simulations
    averageGrid = sumGrid/numSimultaions
    

    return averageGrid

end#function



#----Handles the movement of an agent---- 
#Will need to be modified for higher dimensions
#moveDirection expects an integer between 1 and 4 (inclusive) as follows => 1: up 2: right 3: down 4: left 
#simGrid is an nxn array expecting values of either 0 or 1
#coordinates is an array with the coordinates of the point you want to move
function attemptActionWithDirection(moveDirection, simGrid, coordinates,action)
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
    if(checkInBoundaries(simGrid,moveCoord))
        #Check not moving to an occupied space
        
        if(simGrid[(moveCoord)...]!=1)
            #=
            simGrid[(moveCoord)...]=1.0
            simGrid[intCoordinates...] = 0.0
            =#
            simGrid =action(simGrid,intCoordinates,moveCoord)
        end#if
    
    #Otherwise bounce off the boundary (naive. just a move in the other direction)
    elseif(checkInBoundaries(simGrid,oppositeMoveCoord))
        #Check not moving to an occupied space
        if(simGrid[(oppositeMoveCoord)...]!=1)
            #=
            simGrid[(oppositeMoveCoord)...]=1.0
            simGrid[intCoordinates...] = 0.0
            =#
            simGrid = action(simGrid,intCoordinates,oppositeMoveCoord)
        end#if
    end#if elseif
    return simGrid
end#Function

function moveAgent(simGrid, intCoordinates, moveCoord)
    simGrid[(moveCoord)...]=1.0
    simGrid[intCoordinates...] = 0.0
    return simGrid
end#function

function proliferate(simGrid, intCoordinates, moveCoord)
    simGrid[(moveCoord)...]=1.0
    return simGrid
end#function



#----Check if a given point lies with the bounds of a matrix----
#Takes in an nxn array and an array of indexes of length n.
function checkInBoundaries(simGrid,position)
    for i in eachindex(position)
        if(position[i] < 1 || position[i] > size(simGrid,i))
            return false
        end#if
    end#for
    return true
end

function calculateDensities(simGrid)
    lengths=[size(simGrid)...]#make array with the dimensions of the grid
    #----Relevant to a particular problem. Returning just averageGrid may be preferable
    #calculate density along vertical lines
    density=zeros(lengths[1])
    for i in 1:lengths[1]
        for j in 1:lengths[2]
            density[i]+=simGrid[i,j]
        end#for
    end#for

    #Average number of agents along each vertical line on the grid, across all simulations
    averageDensity = density/lengths[2]
    return averageDensity

end#function


#Used for setting up initial conditions. Will create a solid block of agents on grid 'simgGrid', centered at 'center', with width 'width' and height 'height'.
#xAxisValues is an array containing the position of each column of the grid in cartesian space.
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

function createBlockEven(center, width, height, simGrid, xAxisValues)

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
            for j in axes(simGrid,2)
                
                if(j%(round(size(simGrid,2)/height))==0)
                    simGrid[i,j]=1.0
                end#if
            end#for
            
        end#if
    end#for
    
    return simGrid
end#function



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
    println(Times)
    solution = solve(problem,saveat=[Times...])
    
    for i in 1:length(solution[:,])
        numericSolutions[i,:]=solution[:,i]
    end#for

    return numericSolutions
    
end#function

    #||||----Differential Equations for numerical solving----||||#
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
        du[i] = D*(u[i+1]-2*u[i]+u[i-1])/dx^2 - (v*(1-2*u[i])/(2*dx))*(u[i+1]-u[i-1])
    end#for

    du[1]=0.0
    du[N]=0.0
end#function


#Should have parameters = [step-size, Number of steps(along x direction), Diffusion constant, Logistic constant] Logistic constant is a rate parameter for the logisitic growth equation
function Logistic!(du, u, parameters, t)
    dx,N, D, lmbda = parameters #unpacking into variables
    N = Int(N)
    #Setting up numeric scheme for other derivatives
    for i in 2:N-1
        du[i] = D*(u[i+1]-2*u[i]+u[i-1])/dx^2 + lmbda*u[i]*(1-u[i])
    end#for

    du[1]=0.0
    du[N]=0.0
end#function