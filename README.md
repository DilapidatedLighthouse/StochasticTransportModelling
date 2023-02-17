# Stochastic Transport Modelling
This repository contains code written for the class Stochastic Transport Modelling at the 2023 AMSI Summer School. 

Most of the scripts in this repository run simulations for random walks on a grid where individual elements cannot occupy the same node. The columns of the grid are averaged over several simulations and the result is graphed. These were then compared to the solutions of appropriate partial differential equations. For most of the files, the output is something like this: 

<p align="center">
  <img src="https://user-images.githubusercontent.com/122573155/212587057-dbd8f928-4178-4ece-8ddd-a342c0b82fed.svg"/>
</p>

## Files of interest:
* ExclusionModelSimulation.jl: Created for the first week, it has the simplest model and most concise code.
* functions.jl: Contains most of the functions created for the course. 
* Week4.jl: The final piece written for the course. Default model but adapted for multiple populations of agents.
* FunctionWeek4.jl: Contains some additional functions for use in Week4.jl.
