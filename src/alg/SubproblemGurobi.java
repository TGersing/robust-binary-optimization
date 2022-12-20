package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import alg.AbstractAlgorithm.AlgorithmParameters;
import gurobi.GRB;
import gurobi.GRBConstr;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import util.PossibleZ;
import util.Variable;

/**
 * Parent class to subproblems of a robust problem.
 * 
 * @author Timo Gersing
 */
abstract class SubproblemGurobi extends RobustProblemGurobi{	
	/**
	 * Constraints temporarily added to the subproblem.
	 */
	protected List<GRBConstr> temporaryConstraints = new ArrayList<GRBConstr>();
	
	/**
	 * Optimality cut using the lower bound on z.
	 */
	protected GRBConstr lowerOptCut;
	/**
	 * Optimality cut using the upper bound on z.
	 */
	protected GRBConstr upperOptCut;
	
	SubproblemGurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters) throws GRBException, IOException {
		super(problemPath, robustPath, algorithmParameters);
	}
	
	/**
	 * Returns the optimal z for the bilinear formulation using the computed solution values
	 * for the uncertain variables.
	 */
	double computeOptimalBilinearZ() {
		return computeOptimalBilinearZ(getUncertainVariablesSolutionValues());
	}
	
	/**
	 * Returns the optimal z for the bilinear formulation using given solution values for the
	 * uncertain variables. An optimal choice is the deviation of the greatest variable (w.r.t deviation)
	 * such that the sum of the solution values of all variables with higher or equal deviation 
	 * is greater than or equal to gamma.
	 */
	protected double computeOptimalBilinearZ(double[] uncertainVariablesSolutionValues) {
		double optimalZ = 0;
		double sum = 0;
		int varIndex = getUncertainVariables().length-1;
		while (varIndex >= 0) {
			sum += uncertainVariablesSolutionValues[varIndex];
			if (sum >= getGamma()) {
				optimalZ = getUncertainVariables()[varIndex].getDeviation();
				break;
			}
			varIndex--;
		}
		return optimalZ;
	}
		
	/**
	 * Computes the optimal solution value for the bilinear formulation using the computed solution
	 * and a given value of z.
	 */
	double computeBilinearSolutionValue(double chosenZ) {
		return computeBilinearSolutionValue(getNominalVariablesSolutionValues(), getUncertainVariablesSolutionValues(), chosenZ);
	}
	
	/**
	 * Computes the optimal solution value for the bilinear formulation using a given solution
	 * and a given value of z.
	 */
	protected double computeBilinearSolutionValue(double[] nominalVariablesSolutionValues, double[] uncertainVariablesSolutionValues, double chosenZ) {
		double solutionValue = getGamma() * chosenZ;
		Variable[] nominalVariables = getNominalVariables();
		for (int i = 0; i < nominalVariables.length; i++) {
			solutionValue += nominalVariables[i].getObjectiveCoefficient() * nominalVariablesSolutionValues[i];
		}
		Variable[] uncertainVariables = getUncertainVariables();
		for (int i = 0; i < uncertainVariables.length; i++) {
			if (uncertainVariables[i].getDeviation() > chosenZ) {
				solutionValue += (uncertainVariables[i].getDeviation() - chosenZ)*uncertainVariablesSolutionValues[i];
			}
		}
		return solutionValue;
	}	
	
	/**
	 * Adds optimality-cuts to the model for a given lower and upper bound on z.
	 */
	void addOptimalityCuts(PossibleZ lowerBoundZ, PossibleZ upperBoundZ) throws GRBException, IOException {
		if (lowerBoundZ.getValue() > 0) {
			GRBLinExpr grbLinExpr = new GRBLinExpr();
			grbLinExpr = new GRBLinExpr();
			for (Variable var :uncertainVariables) {
				if (var.getDeviation() >= lowerBoundZ.getValue()){
					grbLinExpr.addTerm(1.0, var.getModelVariable());
				}
			}
			lowerOptCut = model.addConstr(grbLinExpr, GRB.GREATER_EQUAL, Math.ceil(Gamma), "LowerOptimalityCut"+lowerBoundZ.getValue());			
			temporaryConstraints.add(lowerOptCut);
		}

		GRBLinExpr grbLinExpr = new GRBLinExpr();
		for (Variable var :uncertainVariables) {
			if (var.getDeviation() > upperBoundZ.getValue()){
				grbLinExpr.addTerm(1.0, var.getModelVariable());
			}
		}
		upperOptCut = model.addConstr(grbLinExpr, GRB.LESS_EQUAL, Math.floor(Gamma), "UpperOptimalityCut"+upperBoundZ.getValue());
		temporaryConstraints.add(upperOptCut);
		model.update();
		String output = "Use Optimality Cuts ["+lowerBoundZ+","+upperBoundZ+"]";
		AbstractAlgorithm.writeOutput(output, algorithmParameters);
	}
}
