package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Optional;

import alg.AbstractAlgorithm.AlgorithmParameters;
import alg.RobustAlgorithm.RobustAlgorithmStrategies;
import alg.SubproblemNominalGurobi.NOSStrategies.NOSImprovingZStrategy;
import alg.SubproblemNominalGurobi.NOSStrategies.NOSTerminationStrategy;
import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBConstr;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import util.PossibleZ;
import util.Variable;

/**
 * This class implements the nominal subproblems solved during the divide and conquer algorithm.
 * 
 * @author Timo Gersing
 */
class SubproblemNominalGurobi extends SubproblemGurobi {
	/**
	 * The current value of z.
	 */
	private PossibleZ chosenZ;
		
	/**
	 * Specifies strategies.
	 */
	private NOSStrategies nosStrategies;
	
	/**
	 * The value of z for the improved solution.
	 */
	private Double improvedZ = null;

	/**
	 * The global primal bound obtained from a prior subproblem or an improved solution.
	 * Is used to store the solution values of improved solutions and for the termination.
	 */
	private double globalPrimalBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
	/**
	 * Values of variables in a new incumbent. Needs to be stored because an improved sub-optimal solution
	 * might be better than the optimal solution of the subproblem. 
	 */
	private double[] incumbentValues;
	
	/**
	 * Constructor receiving paths to problem files as well as parameters and strategies.
	 */
	SubproblemNominalGurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, NOSStrategies nosStrategies) throws IOException, GRBException {
		super(problemPath, robustPath, algorithmParameters);
		this.nosStrategies = nosStrategies;
	}
	
	/**
	 * Resets the problem to an unsolved state.
	 */
	@Override
	protected void resetProblem() throws GRBException {
		super.resetProblem();
		improvedZ = null;
	}
		
	/**
	 * Updates the current z, alters the objective function, and removes constraints
	 * added for the previous subproblem.
	 */
	void updateZ(PossibleZ chosenZ) throws GRBException {
		this.chosenZ = chosenZ;
		
		//Removes all temporary constraints added for the previous subproblem.
		for (GRBConstr grbConstr : temporaryConstraints) {
			model.remove(grbConstr);
		}
		temporaryConstraints = new ArrayList<GRBConstr>();

		//Sets the callback for the termination and improving.
		model.setCallback(new Callback());
		
		//Alters the objective function.
		GRBLinExpr objLinExpr = new GRBLinExpr();
		for (Variable var : nominalVariables) {
			objLinExpr.addTerm(var.getObjectiveCoefficient(), var.getModelVariable());
		}
		for (Variable var : uncertainVariables) {
			objLinExpr.addTerm(Math.max(var.getDeviation()-chosenZ.getValue(), 0), var.getModelVariable());
		}
		objLinExpr.addConstant(Gamma*chosenZ.getValue());
		model.setObjective(objLinExpr);
	}
	
	/**
	 * Obtains an incumbent solution, computes the optimal z and improves the primal bound if
	 * the new solution is the new best.
	 */
	private void improveZ (double incumbentObjectiveValue, double[] incumbentNominalVariablesValues, double[] incumbentUncertainVariablesValues) {
		//Computes optimal value for z.
		double optimalZ = computeOptimalBilinearZ(incumbentUncertainVariablesValues);
		//Computes the optimal objective value.
		double objectiveValue = computeBilinearSolutionValue(incumbentNominalVariablesValues, incumbentUncertainVariablesValues, optimalZ);
		//Updates the global primal bound if the new incumbent is better.
		if (objectiveValue < globalPrimalBound) {
			globalPrimalBound = objectiveValue;
			improvedZ = optimalZ;
			incumbentValues = incumbentNominalVariablesValues;
			if (objectiveValue < incumbentObjectiveValue - algorithmParameters.getAbsoluteGapTolerance()) {
				String output = "###Improved solution with optimal z = "+optimalZ+" to the new global primal bound = "+objectiveValue;
				try {
					AbstractAlgorithm.writeOutput(output, algorithmParameters);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			else {
				String output = "###Found new global primal bound = "+objectiveValue;
				try {
					AbstractAlgorithm.writeOutput(output, algorithmParameters);
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	
	/**
	 * Implements the callbacks for termination and improving of incumbent solutions.
	 */
	private class Callback extends GRBCallback {
		/**
		 * Saves the last time stamp at which we were able to prune a possible value of z via estimators. 
		 */
		private Long timeStampLastCutoff = null;
		
		/**
		 * Remaining values of z that are not yet pruned for which we have computed an estimator
		 * from the current z.
		 */
		private ArrayList<PossibleZ> comparedPossibleZs = null;
		
		@Override
		protected void callback() {
			try {
				//Tries to terminate the subproblem respecting estimators if the option is chosen.
				if (nosStrategies.getTerminationStrategy() == NOSTerminationStrategy.TERMINATE_ESTIMATORS && where == GRB.CB_MIP) {
					//Queries the current primal and dual bound from the subproblem.
					primalBound = getDoubleInfo(GRB.CB_MIP_OBJBST);
					dualBound = getDoubleInfo(GRB.CB_MIP_OBJBND);
					//Checks whether the subproblem can be terminated with respect to the global primal
					//bound and the current dual bound.
					if (AbstractAlgorithm.isOptimal(Math.min(primalBound, globalPrimalBound), dualBound, algorithmParameters)) {
						if (timeStampLastCutoff == null) {
							timeStampLastCutoff = System.nanoTime();
							comparedPossibleZs = new ArrayList<PossibleZ>(chosenZ.getEstimators().keySet());
						}
						//Checks whether we can hope for pruning additional possible z.
						boolean pruningPossible = false;
						for (int i = comparedPossibleZs.size()-1; i >= 0; i--) {
							//Computes the dual bound for the possible z using the estimator.
							PossibleZ comparedPossibleZ = comparedPossibleZs.get(i);
							double estimator = chosenZ.getEstimators().get(comparedPossibleZ);
							double dualBoundComparedPossibleZ = Math.max(comparedPossibleZ.getDualBound(), dualBound - estimator);
							//Checks whether the dual bound is good enough for pruning.
							if (AbstractAlgorithm.isOptimal(Math.min(primalBound, globalPrimalBound), dualBoundComparedPossibleZ, algorithmParameters)) {
								//Removes the possible z and sets the time stamp.
								comparedPossibleZs.remove(i);
								timeStampLastCutoff = System.nanoTime();
							}
							//If raising the dual bound up to the primal bound would be sufficient
							//to prune the possible z then we can hope for pruning it later and
							//don't terminate the subproblem if the last time stamp is not too far in the past.
							else if (AbstractAlgorithm.isOptimal(Math.min(primalBound, globalPrimalBound), primalBound - estimator, algorithmParameters)) {
								pruningPossible = true;
							}
						}
						//If there is hope for pruning further possible z and the last pruning
						//has happened in the last 10 seconds then we continue solving the subproblem.
						if (!pruningPossible || (System.nanoTime()-timeStampLastCutoff)/Math.pow(10, 9) >= 10) {
							String output = "###Terminated due to global primal bound";
							try {
								AbstractAlgorithm.writeOutput(output, algorithmParameters);
							} catch (IOException e) {}

							abort();
						}
					}
				}

				//Tries to terminate the subproblem directly using the global primal bound if the option is chosen.
				if (nosStrategies.getTerminationStrategy() == NOSTerminationStrategy.TERMINATE_DIRECT && where == GRB.CB_MIP) {
					dualBound = getDoubleInfo(GRB.CB_MIP_OBJBND);
					if (AbstractAlgorithm.isOptimal(globalPrimalBound, dualBound, algorithmParameters)) {
						String output = "###Terminated due to global primal bound";
						try {
							AbstractAlgorithm.writeOutput(output, algorithmParameters);
						} catch (IOException e) {}
						abort();
					}
				}
			} catch (GRBException e1) {
				e1.printStackTrace();
			}


			try {
				if (nosStrategies.getImprovingZStrategy() == NOSImprovingZStrategy.IMPROVE_Z && where == GRB.CB_MIPSOL) {
					double[] incumbentNominalVariablesValues = getSolution(nominalModelVariables);
					double[] incumbentUncertainVariablesValues = getSolution(uncertainModelVariables);
					improveZ(getDoubleInfo(GRB.CB_MIPSOL_OBJ), incumbentNominalVariablesValues, incumbentUncertainVariablesValues);
				}
			} catch (GRBException e) {
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Solves the model using Gurobi.
	 */
	@Override
	void solve(Optional<Double> timeLimit) throws IOException, GRBException {
		super.solve(timeLimit);
		
		//Updates global primal bound
		try {
			if (primalBound < globalPrimalBound) {
				globalPrimalBound = primalBound;
				incumbentValues = model.get(GRB.DoubleAttr.X, nominalModelVariables);
			}
		} catch (GRBException e) { }
		
		//Writes information to the node log.
		String output = "\n#####Finished MILP\n"
				+ "Primal Bound: "+primalBound+"\n"
				+ "Dual Bound: "+dualBound;
		AbstractAlgorithm.writeOutput(output, algorithmParameters);
	}
	
	/**
	 * Returns the improved value for z.
	 */
	Double getImprovedZ() {
		return this.improvedZ;
	}
	
	/**
	 * Returns the global primal bound.
	 */
	double getGlobalPrimalBound() {
		return this.globalPrimalBound;
	}
	
	/**
	 * Returns solution values of the incumbent.
	 */
	double[] getIncumbentValues() {
		return this.incumbentValues;
	}

	
	/**
	 * Sets the global primal bound.
	 */
	void setGlobalPrimalBound(double globalPrimalBound) {
		this.globalPrimalBound = globalPrimalBound;
	}
	
	/**
	 * Specifies strategies for algorithms solving nominal subproblems;
	 */
	abstract static class NOSStrategies extends RobustAlgorithmStrategies{
		/**
		 * Enum type specifying whether we terminate nominal subproblems prematurely.
		 * We can either terminate the problem respecting estimators or directly once the
		 * dual bound is strong enough compared to the global primal bound.
		 */
		public enum NOSTerminationStrategy {
			DONT_TERMINATE,
	 		TERMINATE_ESTIMATORS,
	 		TERMINATE_DIRECT;
		}
		
		/**
		 * Enum type specifying whether we improve incumbent solutions by computing an optimal choice for z.
		 */
		public enum NOSImprovingZStrategy {
			DONT_IMPROVE_Z,
			IMPROVE_Z;
		}

		protected NOSImprovingZStrategy improvingZStrategy;
		protected NOSTerminationStrategy terminationStrategy;
		
		public NOSStrategies() {
			super();
			improvingZStrategy = NOSImprovingZStrategy.IMPROVE_Z;
			terminationStrategy = NOSTerminationStrategy.TERMINATE_DIRECT;
		}

		
		public NOSImprovingZStrategy getImprovingZStrategy() {
			return improvingZStrategy;
		}
		public void setImprovingZStrategy(NOSImprovingZStrategy improvingZStrategy) {
			this.improvingZStrategy = improvingZStrategy;
		}

		public NOSTerminationStrategy getTerminationStrategy() {
			return terminationStrategy;
		}
		public void setTerminationStrategy(NOSTerminationStrategy terminationStrategy) {
			this.terminationStrategy = terminationStrategy;
		}
	}

}