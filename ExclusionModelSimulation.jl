#TO DO:
#   *Add option for bias to movement
#   *Make the plot prettier


using Plots, SpecialFunctions, Random
gr() #A graphical thing? I need to look it up


#||||----FUNCTIONS----||||#

#Meat and bones. Runs simulations of cells moving on a grid and averages the results.
#Lengths is an array of length n specifying the dimensions of the grid.
#totalTime is the number of time-steps for each simulation
#numSimulations is the number of simulations
#probMovement is the probabilty that an agent will try to move if given the chance 
#the nullVariable is so the proliferation probability doesnt have to be removed when switching between functions
function StochasticExclusionWalkAverage(lengths, totalTime, simGrid, numSimultaions, probMovement, biases = [1,2,3,4], nullVariable = "")
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
                        #decide on direction. 1: up 2: right 3:down 4:left
                        randMoveIndex = rand(1:length(biases))
                        moveDirection = biases[randMoveIndex]
                    
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

    #----Relevant to a particular problem. Returning just averageGrid may be preferable
    #calculate density along vertical lines
    for i in 1:lengths[1]
        for j in 1:lengths[2]
            density[i]+=averageGrid[i,j]
        end#for
    end#for

    #Average number of agents along each vertical line on the grid, across all simulations
    averageDensity = density/lengths[2]
    return averageDensity

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

#||||----VARIABLES----||||#

#Declaring variables for simulation
XLENGTH = 500
YLENGTH = 50
PROBABILTYOFMOVE = 1.0
PROBABILTYOFPROLIFERATION = 1.0
TIMEOFSIMULATION = 400
NUMBEROFSIMULATIONS = 4

#Constants for initial conditions
H = 50

#Constants for heat equation
D = PROBABILTYOFMOVE/4




#||||----CODE----||||#

#Setting up x-axis for graphing and grid initialisation
xAxisValues = [((-XLENGTH/2):((XLENGTH/2)-1))...]

#Initialise grid with agents
simGrid = zeros(XLENGTH,YLENGTH)
#=Place vertical block of agents in center of grid with width 2*H
for i in 1:XLENGTH
    if(abs(xAxisValues[i])<=H)
        for j in 1:YLENGTH
            simGrid[i,j]=1.0
        end#for
    end#if
end#for
=#
simGrid = createBlock(0,2*H, 1000, simGrid,xAxisValues)
simGrid = createBlock(-100,H, Int(round(3*YLENGTH/4)), simGrid,xAxisValues)
#Heat equation
T=TIMEOFSIMULATION
C(x)=0.5*(erf((H-x)/sqrt(4*D*T))+erf((H+x)/sqrt(4*D*T)));

#Calculate average densities as function of space.
densities = calculateDensities(StochasticExclusionWalkAverage([XLENGTH, YLENGTH], TIMEOFSIMULATION, simGrid, NUMBEROFSIMULATIONS, PROBABILTYOFMOVE))

#||||----PLOTS----||||#
myPlot = scatter(xAxisValues,densities,mc=:blue,msc=:match,label="Stochastic")
myPlot=plot!(C, -XLENGTH/2:XLENGTH/2-1,lw=4,lc=:red,ls=:dash,label="Exact",xlabel="x",ylabel="C(x,t), ⟨ C(x,t) ⟩",xlims=(-XLENGTH/2,XLENGTH/2),ylims=(0,1))
display(myPlot)