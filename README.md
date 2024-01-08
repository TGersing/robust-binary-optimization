# Robust Binary Optimization
Released: December 20, 2022

Last Updated: May 16, 2023

Written by Timo Gersing

## Introduction
This project provides algorithms for solving robust binary optimization problems with budgeted uncertainty in the objective function.
The algorithms are implemented as described in the thesis "Algorithms for Robust Combinatorial Optimization with Budgeted Uncertainty and Fair Planning of the Out-of-Hours Service for Pharmacies" by Timo Gersing. They have been published before in slightly different versions in the papers "A Branch and Bound Algorithm for Robust Binary Optimization with Budget Uncertainty", published in Mathematical Programming Computation by Christina Büsing, Timo Gersing and Arie Koster, and "Recycling Inequalities for Robust Combinatorial Optimization with Budget Uncertainty" presented by the same authors at IPCO 2023.

## Dependencies
The algorithms use Gurobi as a MILP solver. Make sure to first install Gurobi on your machine. Also specify the path to the gurobi.jar when running the code.

## Specifying Robust Optimization Problems
Robust optimization problems are handed to the algorithms via two separate files. The first file states the nominal problem without uncertainties. This can for example be an .mps or .lp file, as long as it can be read by Gurobi.

The second file specifies the uncertainty in the problem, that is the robustness Budget Γ and deviations of the objective coefficients. The first line in this file specifies the robustness budget. For example `Gamma:25.0` states that the robustness budget is 25. The following lines specify the deviations of objective coefficients. For example, `x15:81` states that the deviation of variable x15 is 81. Here, the name of the variable corresponds to the name in the first file.

An Example for a robust knapsack problem can be found in the `testInstances/` folder.

## Running Algorithms
To solve a problem, run the `AlgorithmExecuter` class. You can specify the problem to solve, the algorithm to use, and optional strategies in two ways. First, you may add arguments to the execution of the `AlgorithmExecuter` class (see example below). Second, you may add no arguments and simply respond to the pop-up windows after starting the program.

In the former case, the first argument is the path to the nominal problem file. The second argument is the path to the robustness components. The third argument is the algorithm to use. Here you can choose from `bnb`, `rec`, `dnc`, `ref`, `cut`, `bss`, `rp1`, `rp2`, `rp3`, `rp4`, `sub`, where bnb is the branch and bound algorithm from "A Branch and Bound Algorithm for Robust Binary Optimization with Budget Uncertainty" and rec uses recycled inequalities from "Recycling Inequalities for Robust Combinatorial Optimization with Budget Uncertainty".

You can optionally specify a time limit in seconds, a destination for a log file, a destination for a results file (with objective values/gaps/computation time...), and a destination for a solution file. For this, optionally append further arguments of the form `timelimit=*`, `logpath=*`, `resultspath=*`, and `solutionpath=*` in any order.

Optional strategies for the algorithm are chosen by appending them to the arguments above. For example, the branch and bound algorithm has a parameter `BoundingLPOptimalityCutsStrategy`, in which we choose whether to use optimality cuts when solving LPs or not. The corresponding choices to specify in the arguments are `LPOPTCUTS_ENABLE` or `LPOPTCUTS_DISABLE`. We refer to the strategy classes of each algorithm for the different available options.

Arguments for solving the knapsack example with the branch and bound algorithm and specifying optionals would for example be:
```testInstances/knapsack.mps.gz testInstances/knapsack_robust.txt bnb timelimit=3600 logpath=log.txt solutionpath=solution.txt LPOPTCUTS_DISABLE ESTIMATORS_IMPROVED```
