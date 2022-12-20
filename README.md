# Robust Binary Optimization
Released: December 20, 2023

Written by Timo Gersing

## Introduction
This project provides algorithms for solving robust binary optimization problems with budgeted uncertainty in the objective function.
The algorithms are implemented as described in the paper "A Branch & Bound Algorithm for Robust Binary Optimization with Budget Uncertainty", published in Mathematical Programming Computation by Christina Büsing, Timo Gersing and Arie Koster.

## Dependencies
The algorithms use Gurobi as a MILP solver. Make sure to first install Gurobi on your machine. Also specify the path to the gurobi.jar when running the code.

## Specifying Robust Optimization Problems
Robust optimization problems are handed to the algorithms via two separate files. The first file states the nominal problem without uncertainties. This can for example be an .mps or .lp file, as long as it can be read by Gurobi.

The second file specifies the uncertainty in the problem, that is the robustness Budget Γ and deviations of the objective coefficients. The first line in this file specifies the robustness budget. For example `Gamma:25.0` states that the robustness budget is 25. The following lines specify the deviations of objective coefficients. For example, `x15:81` states that the deviation of variable x15 is 81. Here, the name of the variable corresponds to the name in the first file.

An Example for a robust knapsack problem can be found in the testInstances/ folder.

## Running Algorithms
To solve a problem, run the AlgorithmExecuter class. You can specify the problem to solve and the algorithm to use either by responding in the console or by directly giving the arguments.

In the latter case, the first argument is the path to the nominal problem file. The second argument is the path to the robustness components. The third argument is the algorithm to use. Here you can choose from `bnb`, `dnc`, `ref`, `cut`, `bss`, `rp1`, `rp2`, `rp3`, `rp4`, where bnb is the branch & bound algorithm in the paper.

You can optionally specify a time limit in seconds, a destination for a log file and a destination for a solution file. For this, optionally append further arguments of the form `timelimit=*`, `logpath=*`, and `solutionpath=*` in any order.

Arguments for solving the knapsack example with the branch & bound algorithm and specifying optionals would for example be:
```testInstances/knapsack.mps.gz testInstances/knapsack_robust.txt bnb timelimit=3600 logpath=log.txt solutionpath=solution.txt```
