package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import alg.AbstractAlgorithm.AlgorithmParameters;
import alg.SubproblemGurobi.SubproblemsStrategies.ImprovingZStrategy;
import alg.SubproblemNominalGurobi.NOSStrategies.NOSTerminationStrategy;
import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBConstr;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import util.PossibleZ;
import util.Variable;

/**
 * This class implements the nominal subproblems solved during the Bertsimas Sim approach and the divide and conquer algorithm.
 * 
 * @author Timo Gersing
 */
class SubproblemNominalGurobi extends SubproblemGurobi {
	/**
	 * Specifies strategies.
	 */
	private NOSStrategies nosStrategies;
	
	/**
	 * The current value of z defining the nominal subproblem.
	 */
	private PossibleZ chosenZ;
	
	/**
	 * The improved (optimal) value of z for the incumbent solution.
	 */
	private Double improvedZ = null;

	/**
	 * The primal bound provided by a potentially improved solution.
	 * Is always at least as good as the primal bound for the current problem.
	 */
	private double improvedPrimalBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
	
	/**
	 * Values of variables in a new incumbent. Needs to be stored because an improved sub-optimal solution
	 * might be better than the optimal solution of the subproblem. 
	 */
	private double[] improvedNominalVariablesSolutionValues;
	
	/**
	 * Callback passed to the subproblem.
	 * Reports bounds from the subproblem to the master and asks for termination.
	 */
	private MasterCallback masterCallback;

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
		improvedPrimalBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
		improvedNominalVariablesSolutionValues = null;
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
	private void improveZ (double[] incumbentNominalVariablesValues, double[] incumbentUncertainVariablesValues) {
		//Computes optimal value for z.
		double optimalZ = computeOptimalBilinearZ(incumbentUncertainVariablesValues);
		//Computes the optimal objective value.
		double objectiveValue = computeBilinearSolutionValue(incumbentNominalVariablesValues, incumbentUncertainVariablesValues, optimalZ);
		//Updates the improved primal bound if the new incumbent is better.
		if (objectiveValue < improvedPrimalBound) {
			improvedPrimalBound = objectiveValue;
			improvedZ = optimalZ;
			improvedNominalVariablesSolutionValues = incumbentNominalVariablesValues;
		}
	}

	
	/**
	 * Implements the callbacks for termination and improving of incumbent solutions.
	 */
	private class Callback extends GRBCallback {	
		@Override
		protected void callback() {
			//Stores incumbent solutions and potentially tries to improve an incumbent solution computing an optimal value for z.
			try {
				if (where == GRB.CB_MIPSOL) {
					if (nosStrategies.getImprovingZStrategy() == ImprovingZStrategy.IMPROVINGZ_ENABLE) {
						double[] incumbentNominalVariablesValues = getSolution(nominalModelVariables);
						double[] incumbentUncertainVariablesValues = getSolution(uncertainModelVariables);
						improveZ(incumbentNominalVariablesValues, incumbentUncertainVariablesValues);
					}
					else {
						if (getDoubleInfo(GRB.CB_MIPSOL_OBJ) < improvedPrimalBound) {
							improvedPrimalBound = getDoubleInfo(GRB.CB_MIPSOL_OBJ);
							improvedZ = chosenZ.getValue();
							improvedNominalVariablesSolutionValues = getSolution(nominalModelVariables);
						}
					}
				}
			} catch (GRBException e) {
				e.printStackTrace();
			}
			
			//Reports current primal and dual bounds to the master and possibly asks for termination
			try {
				if (where == GRB.CB_MIP) {
					dualBound = getDoubleInfo(GRB.CB_MIP_OBJBND);
					
					if (nosStrategies.getTerminationStrategy() != NOSTerminationStrategy.TERMINATION_DISABLE) {
						if (masterCallback.updateBoundsAndDecideTermination(improvedPrimalBound, primalBound, dualBound)) {
							String output = "###Terminated due to global primal bound";
							try {
								AbstractAlgorithm.writeOutput(output, algorithmParameters);
							} catch (IOException e) {}
							abort();
						}
					}
					else {
						masterCallback.updatePrimalDualBounds(improvedPrimalBound, dualBound);
					}
				}
			} catch (GRBException e1) {
				e1.printStackTrace();
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
		if (primalBound < improvedPrimalBound) {
			improvedPrimalBound = primalBound;
			try {
				improvedNominalVariablesSolutionValues = model.get(GRB.DoubleAttr.X, nominalModelVariables);
			} catch (GRBException e) { }
		}
		masterCallback.updatePrimalDualBounds(improvedPrimalBound, dualBound);
		
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
	 * Returns the improved primal bound.
	 */
	double getImprovedPrimalBound() {
		return this.improvedPrimalBound;
	}
	
	/**
	 * Returns solution values of the potentially improved solution.
	 */
	double[] getImprovedNominalSolutionValues() {
		return this.improvedNominalVariablesSolutionValues;
	}
	
	public void setMasterCallback(MasterCallback masterCallback) {
		this.masterCallback = masterCallback;
	}
	
	/**
	 * Interface for the master callback that reports bounds from the subproblem to the master and asks for termination.
	 * Is implemented by the master algorithm (@AlgBertsimasSimSequence or @AlgDivideAndConquer).
	 */
	protected interface MasterCallback {

		boolean updateBoundsAndDecideTermination(double improvedPrimalBound, double primalBound, double dualBound);

		boolean updatePrimalDualBounds(double improvedPrimalBound, double dualBound);
		
	}
	
	/**
	 * Specifies strategies for algorithms solving nominal subproblems;
	 */
	abstract static class NOSStrategies extends SubproblemsStrategies{
		/**
		 * Enum type specifying whether we terminate nominal subproblems prematurely.
		 * We can either terminate the problem respecting estimators or directly once the
		 * dual bound is strong enough compared to the global primal bound.
		 */
		public enum NOSTerminationStrategy {
	 		TERMINATION_ESTIMATORS,
	 		TERMINATION_DIRECT,
			TERMINATION_DISABLE;
		}
		
		NOSTerminationStrategy terminationStrategy;
		
		/**
		 * Constructor obtaining arguments which are matched to the enums defining strategies.
		 */
		public NOSStrategies(List<String> argList, AlgorithmParameters algorithmParameters) throws IOException {
			super(argList, algorithmParameters);
		}

		public NOSTerminationStrategy getTerminationStrategy() {
			return terminationStrategy;
		}
		public void setTerminationStrategy(NOSTerminationStrategy terminationStrategy) {
			this.terminationStrategy = terminationStrategy;
		}
	}

}