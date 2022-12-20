package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Optional;

import alg.AbstractAlgorithm.AlgorithmParameters;
import alg.AlgBranchAndBound.BnBStrategies.BnBLPOptimalityCutsStrategy;
import alg.RobustAlgorithm.RobustAlgorithmStrategies;
import alg.SubproblemBoundedGurobi.BoundingSubproblemsStrategies.BoundingImprovingZStrategy;
import alg.SubproblemBoundedGurobi.BoundingSubproblemsStrategies.BoundingTerminationStrategy;
import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBConstr;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBVar;
import util.BnBNode;
import util.PossibleZ;
import util.Variable;
import util.Variable.Type;

/**
 * This class implements the integer and relaxed robust subproblems solved during the branch and bound algorithm.
 * 
 * @author Timo Gersing
 */
public class SubproblemBoundedGurobi extends SubproblemGurobi{
	/**
	 * The current lower bound on z.
	 */
	private PossibleZ lowerBoundZ;
	
	/**
	 * The current upper bound on z.
	 */
	private PossibleZ upperBoundZ;
	
	/**
	 * The substituted variables p from the reformulation.
	 */
	private GRBVar[] pPrime;
	
	/**
	 * The substituted variable z from the (clique) reformulation.
	 */
	private GRBVar zPrime;
	
	/**
	 * The value of the (resubstituted) variable z.
	 */
	private Double zSolutionValue;
	
	/**
	 * The global primal bound obtained from a prior subproblem, the solution to this subproblem or an improved solution.
	 */
	private double globalPrimalBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
	/**
	 * The value of z for the improved solution.
	 */
	private Double improvedZ = null;
	/**
	 * Values of variables in a new incumbent. Needs to be stored because an improved sub-optimal solution
	 * might be better than the optimal solution of the subproblem. 
	 */
	private double[] incumbentValues;
	
	/**
	 * Specifies strategies.
	 */
	private BoundingSubproblemsStrategies boundingStrategies;
	
	/**
	 * Constructor receiving paths to problem files as well as parameters and strategies.
	 */
	SubproblemBoundedGurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, BoundingSubproblemsStrategies boundingStrategies) throws IOException, GRBException {
		super(problemPath, robustPath, algorithmParameters);
		
		this.boundingStrategies = boundingStrategies;
		
		zPrime = model.addVar(0, Double.MAX_VALUE, Gamma, GRB.CONTINUOUS, "z");
		
		if (cliquePartitioning != null) {
			pPrime = new GRBVar[cliquePartitioning.getCliques().size()];
		}
		else {
			pPrime = new GRBVar[uncertainVariables.length];
		}
		for (int i = 0; i < pPrime.length; i++) {
			pPrime[i] = model.addVar(0, Double.MAX_VALUE, 1, GRB.CONTINUOUS, "p"+i);
		}
	}
	
	/**
	 * Resets the problem to an unsolved state.
	 */
	@Override
	protected void resetProblem() throws GRBException {
		super.resetProblem();
		zSolutionValue = null;
		improvedZ = null;
	}
	
	/**
	 * Solves the linear relaxation of the subproblem.
	 */
	void solveRelaxed(Optional<Double> timeLimit) throws IOException, GRBException  {
		//Sets alls variables to continuous
		for (GRBVar grbVar : nominalModelVariables) {
			grbVar.set(GRB.CharAttr.VType, GRB.CONTINUOUS);
		}
		//Removes callback.
		model.setCallback(null);
		
		this.solve(timeLimit);
		
		//Writes information to the node log.
		String output = "\n##### Finished LP"
				+ "\nDual bound = "+dualBound;
		AbstractAlgorithm.writeOutput(output, algorithmParameters);
	}
	
	/**
	 * Solves the integer subproblem.
	 */
	void solveInteger(Optional<Double> timeLimit) throws IOException, GRBException  {
		//Sets variables to their original type
		for (Variable var : nominalVariables) {
			if (var.getType() == Type.BINARY) {
				var.getModelVariable().set(GRB.CharAttr.VType, GRB.BINARY);
			}
			else if (var.getType() == Type.INTEGER) {
				var.getModelVariable().set(GRB.CharAttr.VType, GRB.INTEGER);
			}
			else if (var.getType() == Type.CONTINUOUS) {
				var.getModelVariable().set(GRB.CharAttr.VType, GRB.CONTINUOUS);
			}
		}
		//Adds callback
		model.setCallback(new Callback());
		
		this.solve(timeLimit);
		
		//Writes information to the node log
		String output = "\n#####Finished MILP\n"
				+ "Primal Bound: "+primalBound+"\n"
				+ "Dual Bound: "+dualBound;
		AbstractAlgorithm.writeOutput(output, algorithmParameters);
		
		//Updates global primal bound
		try {
			if (primalBound < globalPrimalBound) {
				globalPrimalBound = primalBound;
				incumbentValues = model.get(GRB.DoubleAttr.X, nominalModelVariables);
			}
		} catch (GRBException e) { }
	}
	
	/**
	 * Updates bounds on z and the reformulation for a given node of the branching tree.
	 */
	void updateBounds(BnBNode node) throws GRBException, IOException {
		//Updates bounds
		this.lowerBoundZ=node.getLowerBoundZ();
		this.upperBoundZ=node.getUpperBoundZ();
		
		//Removes all temporary constraints added for the previous subproblem
		for (GRBConstr grbConstr : temporaryConstraints) {
			model.remove(grbConstr);
		}
		temporaryConstraints = new ArrayList<GRBConstr>();
		
		//Alters the objective function.
		GRBLinExpr objLinExpr = new GRBLinExpr();
		for (Variable nominalVariable : nominalVariables) {
			objLinExpr.addTerm(nominalVariable.getObjectiveCoefficient(), nominalVariable.getModelVariable());
		}
		for (Variable uncertainVariable : uncertainVariables) {
			if (uncertainVariable.getDeviation() > upperBoundZ.getValue()) {
				objLinExpr.addTerm(uncertainVariable.getDeviation() - upperBoundZ.getValue(), uncertainVariable.getModelVariable());
			}
		}
		for (int i = 0; i < pPrime.length; i++) {
			objLinExpr.addTerm(1, pPrime[i]);
		}
		objLinExpr.addTerm(Gamma, zPrime);
		objLinExpr.addConstant(Gamma*lowerBoundZ.getValue());
		model.setObjective(objLinExpr);
		
		//Adds the constraints from the (clique) reformulation
		if (cliquePartitioning != null) {
			for (int i = 0; i < cliquePartitioning.getCliques().size(); i++) {
				//We only add a constraint if at least one of the deviations - lower bound on z is larger
				//than the feasibility tolerance. Otherwise, adding the constraint is probably insatisfied 
				//anyway and may cause numerical problems.
				boolean addConstraint = false;
				for (int varIndex : cliquePartitioning.getCliques().get(i)) {
					if (uncertainVariables[varIndex].getDeviation() - lowerBoundZ.getValue() > algorithmParameters.getFeasibilityTolerance()) {
						addConstraint = true;
						break;
					}
				}
				if (addConstraint) {
					GRBLinExpr robustnessExpr = new GRBLinExpr();
					robustnessExpr.addTerm(1, pPrime[i]);
					robustnessExpr.addTerm(1, zPrime);
					for (int varIndex : cliquePartitioning.getCliques().get(i)) {
						Variable var = uncertainVariables[varIndex];
						if (var.getDeviation() > lowerBoundZ.getValue()) {
							robustnessExpr.addTerm(-(Math.min(var.getDeviation(), upperBoundZ.getValue())-lowerBoundZ.getValue()), var.getModelVariable());
						}
					}
					temporaryConstraints.add(model.addConstr(robustnessExpr, GRB.GREATER_EQUAL, 0, ""));
				}
			}
		}
		else {
			for (int i = uncertainVariables.length-1; i >= 0; i--) {
				Variable var = uncertainVariables[i];
				//We only add a constraint if deviations - lower bound is larger than the feasibility tolerance.
				if (var.getDeviation() - lowerBoundZ.getValue() > algorithmParameters.getFeasibilityTolerance()) {
					GRBLinExpr robustnessExpr = new GRBLinExpr();
					robustnessExpr.addTerm(1, pPrime[i]);
					robustnessExpr.addTerm(1, zPrime);
					robustnessExpr.addTerm(-(Math.min(var.getDeviation(), upperBoundZ.getValue())-lowerBoundZ.getValue()), var.getModelVariable());
					temporaryConstraints.add(model.addConstr(robustnessExpr, GRB.GREATER_EQUAL, 0, ""));
				}
			}
		}
		model.update();
	}
	
	/**
	 * Obtains an incumbent solution, computes the optimal z and improves the primal bound if
	 * the new solution is the new best.
	 */
	private void improveZ (double incumbentObjectiveValue, double[] incumbentNominalVariablesValues, double[] incumbentUncertainVariablesValues) {
		//Computes optimal value for z
		double optimalZ = computeOptimalBilinearZ(incumbentUncertainVariablesValues);
		//Computes the optimal objective value
		double objectiveValue = computeBilinearSolutionValue(incumbentNominalVariablesValues, incumbentUncertainVariablesValues, optimalZ);
		//Updates the global primal bound if the new incumbent is better
		if (objectiveValue < globalPrimalBound ) {
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
		private Long timeStampLastPrune = null;
		
		/**
		 * Remaining values of z that are not yet pruned for which we have computed an estimator
		 * from the upper bound of z.
		 */
		private ArrayList<PossibleZ> remainingPossibleZsForUpperBoundZ = null;
		
		/**
		 * Remaining values of z that are not yet pruned for which we have computed an estimator
		 * from the lower bound of z.
		 */
		private ArrayList<PossibleZ> remainingPossibleZsForLowerBoundZ = null;

		@Override
		protected void callback() {
			try {
				//Tries to terminate the subproblem if the option is chosen.
				if (boundingStrategies.getTerminationStrategy() == BoundingTerminationStrategy.TERMINATE && where == GRB.CB_MIP) {
					//Queries the current primal and dual bound from the subproblem.
					primalBound = getDoubleInfo(GRB.CB_MIP_OBJBST);
					dualBound = getDoubleInfo(GRB.CB_MIP_OBJBND);
					//Checks whether the subproblem can be terminated with respect to the global primal
					//bound and the current dual bound.
					if (AbstractAlgorithm.isOptimal(Math.min(primalBound, globalPrimalBound), dualBound, algorithmParameters)) {
						//Initializes the time stamp and remaining possible z for pruning.
						if (timeStampLastPrune == null) {
							timeStampLastPrune = System.nanoTime();
							remainingPossibleZsForUpperBoundZ = new ArrayList<PossibleZ>(upperBoundZ.getEstimators().keySet());
							remainingPossibleZsForLowerBoundZ = new ArrayList<PossibleZ>(lowerBoundZ.getEstimators().keySet());
						}
						//Checks whether we can hope for pruning additional possible z.
						boolean pruningPossible = false;
						for (int i = remainingPossibleZsForUpperBoundZ.size()-1; i >= 0; i--) {
							//Computes the dual bound for the possible z using the estimator.
							PossibleZ comparedPossibleZ = remainingPossibleZsForUpperBoundZ.get(i);
							double estimator = upperBoundZ.getEstimators().get(comparedPossibleZ);
							double dualBoundComparedPossibleZ = Math.max(comparedPossibleZ.getDualBound(), dualBound - estimator);
							//Checks whether the dual bound is good enough for pruning.
							if (AbstractAlgorithm.isOptimal(Math.min(primalBound, globalPrimalBound), dualBoundComparedPossibleZ, algorithmParameters)) {
								//Removes the possible z and sets the time stamp.
								remainingPossibleZsForUpperBoundZ.remove(i);
								timeStampLastPrune = System.nanoTime();
							}
							//If raising the dual bound up to the primal bound would be sufficient
							//to prune the possible z then we can hope for pruning it later and
							//don't terminate the subproblem if the last time stamp is not too far in the past.
							else if (AbstractAlgorithm.isOptimal(Math.min(primalBound, globalPrimalBound), primalBound - estimator, algorithmParameters)) {
								pruningPossible = true;
							}
						}
						//Same as above for the possible z for which we have computed an estimator
						//from the lower bound of z.
						for (int i = remainingPossibleZsForLowerBoundZ.size()-1; i >= 0; i--) {
							PossibleZ comparedPossibleZ = remainingPossibleZsForLowerBoundZ.get(i);
							double estimator = lowerBoundZ.getEstimators().get(comparedPossibleZ);
							double dualBoundComparedPossibleZ = Math.max(comparedPossibleZ.getDualBound(), dualBound - estimator);
							if (AbstractAlgorithm.isOptimal(globalPrimalBound, dualBoundComparedPossibleZ, algorithmParameters)) {
								remainingPossibleZsForLowerBoundZ.remove(i);
								timeStampLastPrune = System.nanoTime();
							}
							else if (AbstractAlgorithm.isOptimal(globalPrimalBound, primalBound - estimator, algorithmParameters)) {
								pruningPossible = true;
							}
						}
						//If there is hope for pruning further possible z and the last pruning
						//has happened in the last 10 seconds then we continue solving the subproblem.
						if (!pruningPossible || (System.nanoTime()-timeStampLastPrune)/Math.pow(10, 9) >= 10) {
							String output = "###Terminated due to global primal bound";
							try {
								AbstractAlgorithm.writeOutput(output, algorithmParameters);
							} catch (IOException e) {}
							
							abort();
						}
					}
				}
			} catch (GRBException e1) {
				e1.printStackTrace();
			}

			//Tries to improve an incumbent solution computing an optimal value for z.
			try {
				if (boundingStrategies.getImprovingZStrategy() == BoundingImprovingZStrategy.IMPROVE_Z && where == GRB.CB_MIPSOL) {
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
	void solve(Optional<Double> timeLimit) throws IOException, GRBException  {
		super.solve(timeLimit);
		//Tries to query the solution values for z (resubstituted).
		//Value is null if this fails.
		try {
			zSolutionValue = lowerBoundZ.getValue() + zPrime.get(GRB.DoubleAttr.X);
		} catch (GRBException e) { }
	}
	
	/**
	 * Returns whether a branching is effective by checking if the current solution is cut-off
	 * in the child nodes defined by the new lower bound and the new upper bound.
	 * Takes optimality-cuts into consideration if the option is chosen.
	 */
	boolean isBranchingEffective(double Gamma, double newUpperBound, double newLowerBound, BnBLPOptimalityCutsStrategy optimalityCutsStrategy) {
		//If z obeys the new upper bound then we have to check whether the new optimality-cuts
		//or the stronger formulation cut-off the current solution.
		//If this is the case then the branching is effective, as z does not obey the new lower bound.
		if (zSolutionValue <= newUpperBound) {
			boolean isEffective = false;
			if (optimalityCutsStrategy == BnBLPOptimalityCutsStrategy.CUTS) {
				double lhs = 0;
				for (int i = 0; i < uncertainVariables.length; i++) {
					Variable var = uncertainVariables[i];
					if (var.getDeviation() > newUpperBound){
						lhs+=uncertainVariablesSolutionValues[i];
					}
				}
				if (lhs > Math.floor(Gamma)+algorithmParameters.getFeasibilityTolerance()) {
					isEffective = true;
				}
			}
			if (!isEffective) {
				isEffective = isNewUpperBoundEffective(newUpperBound);
			}
			return isEffective;
		}
		//If z obeys the new lower bound then we have to check whether the new optimality-cuts
		//or the stronger formulation cut-off the current solution.
		//If this is the case then the branching is effective, as z does not obey the new upper bound.
		if (zSolutionValue >= newLowerBound) {
			boolean isEffective = false;
			if (optimalityCutsStrategy == BnBLPOptimalityCutsStrategy.CUTS) {
				double lhs = 0;
				for (int i = 0; i < uncertainVariables.length; i++) {
					Variable var = uncertainVariables[i];
					if (var.getDeviation() >= newLowerBound){
						lhs+=uncertainVariablesSolutionValues[i];
					}
				}
				if (lhs < Math.ceil(Gamma)-algorithmParameters.getFeasibilityTolerance()) {
					isEffective = true;
				}
			}
			if (!isEffective) {
				isEffective = isNewLowerBoundEffective(newLowerBound);
			}
			return isEffective;
		}
		//Returns true of z obeys neither the new lower not the new upper bound. 
		return true;
	}
	
	/**
	 * Checks whether the current solution is cut-off for a new upper bound on z.
	 */
	private boolean isNewUpperBoundEffective(double newUpperBound) {
		//Queries the values of the substituted variable p. If this fails then there exists
		//no computed solution and we return true.
		double[] pPrimeSolutionValues;
		try {
			pPrimeSolutionValues = model.get(GRB.DoubleAttr.X, pPrime);
		} catch (GRBException e) {
			return true;
		}
		
		//Computes the left-hand and right-hand side of the new (not substituted) constraint
		//p >= (deviation - newUpperBound)x.
		//Returns true if lhs < rhs for any constraint, false otherwise.
		if (cliquePartitioning != null) {
			for (int i = 0; i < pPrime.length; i++) {
				double lhs = pPrimeSolutionValues[i];
				double rhs = 0;
				for (int varIndex : cliquePartitioning.getCliques().get(i)) {
					Variable var = uncertainVariables[varIndex];
					if (uncertainVariablesSolutionValues[varIndex] > 0) {
						//Resubstitutes the variables p.
						if (var.getDeviation() > upperBoundZ.getValue()) {
							lhs += (var.getDeviation() - upperBoundZ.getValue()) * uncertainVariablesSolutionValues[varIndex];
						}
						if (var.getDeviation() > newUpperBound) {
							rhs += (var.getDeviation() - newUpperBound) * uncertainVariablesSolutionValues[varIndex];
						}
					}
				}
				if (lhs < rhs - algorithmParameters.getFeasibilityTolerance()) {
					return true;
				}
			}
		}
		else {
			for (int i = 0; i < pPrime.length; i++) {
				if (uncertainVariablesSolutionValues[i] > 0) {
					Variable var = uncertainVariables[i];

					//Resubstitutes the variable p.
					double lhs = pPrimeSolutionValues[i];
					if (var.getDeviation() > upperBoundZ.getValue()) {
						lhs += (var.getDeviation() - upperBoundZ.getValue()) * uncertainVariablesSolutionValues[i];
					}
					
					double rhs = 0;
					if (var.getDeviation() > newUpperBound) {
						rhs += (var.getDeviation() - newUpperBound) * uncertainVariablesSolutionValues[i];
					}
					if (lhs < rhs - algorithmParameters.getFeasibilityTolerance()) {
						return true;
					}
				}
			}
		}
		return false;
	}

	/**
	 * Checks whether the current solution is cut-off for a new upper bound on z.
	 */
	private boolean isNewLowerBoundEffective(double newLowerBound) {
		//Queries the values of the substituted variable p. If this fails then there exists
		//no computed solution and we return true.
		double[] pPrimeSolutionValues;
		try {
			pPrimeSolutionValues = model.get(GRB.DoubleAttr.X, pPrime);
		} catch (GRBException e) {
			return true;
		}
		
		//Computes the left-hand and right-hand side of the new (not substituted) constraint
		//p + z >= (deviation - newLowerBound)x + newLowerBound.
		//Returns true if lhs < rhs for any constraint, false otherwise.
		if (cliquePartitioning != null) {
			for (int i = 0; i < pPrime.length; i++) {
				double lhs = zSolutionValue + pPrimeSolutionValues[i];
				double rhs = newLowerBound;
				
				for (int varIndex : cliquePartitioning.getCliques().get(i)) {
					Variable var = uncertainVariables[varIndex];
					if (uncertainVariablesSolutionValues[varIndex] > 0) {
						//Resubstitutes the variables p.
						if (var.getDeviation() > upperBoundZ.getValue()) {
							lhs += (var.getDeviation() - upperBoundZ.getValue()) * uncertainVariablesSolutionValues[varIndex];
						}
						if (var.getDeviation() > newLowerBound) {
							rhs += (var.getDeviation() - newLowerBound) * uncertainVariablesSolutionValues[varIndex];
						}
					}
				}
				
				if (lhs < rhs - algorithmParameters.getFeasibilityTolerance()) {
					return true;
				}
			}
		}
		else {
			for (int i = 0; i < pPrime.length; i++) {
				if (uncertainVariablesSolutionValues[i] > 0) {
					Variable var = uncertainVariables[i];
					double lhs = zSolutionValue + pPrimeSolutionValues[i];
					
					//Resubstitutes the variables p.
					if (var.getDeviation() > upperBoundZ.getValue()) {
						lhs += (var.getDeviation() - upperBoundZ.getValue()) * uncertainVariablesSolutionValues[i];
					}
					
					double rhs = newLowerBound;
					if (var.getDeviation() > newLowerBound) {
						rhs += (var.getDeviation() - newLowerBound) * uncertainVariablesSolutionValues[i];
					}
					
					if (lhs < rhs - algorithmParameters.getFeasibilityTolerance()) {
						return true;
					}
				}
			}
		}
		return false;
	}
	
	/**
	 * Checks whether the current solution is integer feasible.
	 */
	boolean isSolutionIntegerFeasible() {
		boolean feasible = true;
		for (int i = 0; i < nominalVariables.length; i++) {
			Variable variable = nominalVariables[i];
			if (variable.getType() == Type.BINARY || variable.getType() == Type.INTEGER) {
				if (Math.abs(nominalVariablesSolutionValues[i] - Math.round(nominalVariablesSolutionValues[i])) > algorithmParameters.getIntegerFeasibilityTolerance()) {
					feasible = false;
					break;
				}
			}
		}
		return feasible;
	}
	
	/**
	 * Returns the value for z computed by Gurobi.
	 */
	Double getComputedZ() {
		return zSolutionValue;
	}
	
	/**
	 * Returns the global primal bound.
	 */
	double getGlobalPrimalBound() {
		return this.globalPrimalBound;
	}
	
	/**
	 * Returns the improved value for z.
	 */
	Double getImprovedZ() {
		return this.improvedZ;
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
	 * Specifies strategies for algorithms solving bounded subproblems.
	 */
	abstract static class BoundingSubproblemsStrategies extends RobustAlgorithmStrategies{
		/**
		 * Enum type specifying whether we terminate robust subproblems prematurely.
		 */
		public enum BoundingTerminationStrategy {
			DONT_TERMINATE,
	 		TERMINATE;
		}
		
		/**
		 * Enum type specifying whether we improve incumbent solutions by computing an optimal choice for z.
		 */
		public enum BoundingImprovingZStrategy {
			DONT_IMPROVE_Z,
			IMPROVE_Z;
		}
		
		public BoundingSubproblemsStrategies() {
			super();
			improvingZStrategy = BoundingImprovingZStrategy.IMPROVE_Z;
			terminationStrategy = BoundingTerminationStrategy.TERMINATE;
		}


		protected BoundingImprovingZStrategy improvingZStrategy;
		protected BoundingTerminationStrategy terminationStrategy;
		
		public BoundingImprovingZStrategy getImprovingZStrategy() {
			return improvingZStrategy;
		}
		public void setImprovingZStrategy(BoundingImprovingZStrategy improvingZStrategy) {
			this.improvingZStrategy = improvingZStrategy;
		}
		public BoundingTerminationStrategy getTerminationStrategy() {
			return terminationStrategy;
		}
		public void setTerminationStrategy(BoundingTerminationStrategy terminationStrategy) {
			this.terminationStrategy = terminationStrategy;
		}
	}

}
