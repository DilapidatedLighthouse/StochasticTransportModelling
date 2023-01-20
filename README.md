# Stochastic Transport Modelling
This repository contains code written for the class Stochastic Transport Modelling at the 2023 AMSI Summer School. 

## ExclusionModelSimulation.jl
This script runs simulations for random walks on a grid where individual elements cannot occupy the same node. There is also the option to allow the agents to multiply with a certain frequency. The columns of the grid are averaged over several simulations and the result is graphed against the solution to the diffusion equation:

<p align="center">
  <img src="https://user-images.githubusercontent.com/122573155/212587057-dbd8f928-4178-4ece-8ddd-a342c0b82fed.svg"/>
</p>
##functions.jl
This file contains all major functions called by other files in the project.
