package alg;

import java.io.IOException;
import java.util.Optional;

import alg.AbstractAlgorithm.AlgorithmParameters;
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import gurobi.GRB.DoubleAttr;
import util.Variable;

/**
 * This class imports a problem from a file and solves it using Gurobi.
 * We only support minimization problems. If the problem is a maximization problem, then we convert it.
 * 
 * @author Timo Gersing
 */
public class ProblemGurobi{
	/**
	 * The Gurobi model which we import from the file.
	 */
	protected GRBModel model;
	
	/**
	 * Parameters of the algorithm, specifying logging, threads and tolerances.
	 */
	protected AlgorithmParameters algorithmParameters = new AlgorithmParameters();

	/**
	 * An array containing the variables of the problem.
	 */
	protected Variable[] nominalVariables;
	
	/**
	 * An array containing the variables of the GRBmodel.
	 */
	protected GRBVar[] nominalModelVariables;
	
	/**
	 * The variables' solution values.
	 */
	protected double[] nominalVariablesSolutionValues;
	
	/**
	 * A primal bound on the optimal objective value.
	 */
	protected double primalBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;

	/**
	 * A dual bound on the optimal objective value.
	 */
	protected double dualBound = AbstractAlgorithm.DEFAULT_DUAL_BOUND;
	
	/**
	 * Indicates whether the problem was observed to be infeasible.
	 */
	protected boolean isStatusInfeasible = false;
	
	/**
	 * Indicates whether we solved the problem to optimality.
	 */
	protected boolean isStatusOptimal = false;
	
	/**
	 * Constructor importing a problem from a file and reading the variables.
	 * If the imported problem is a maximization problem, then we convert it to a minimization problem.
	 */
	ProblemGurobi(String problemPath, AlgorithmParameters algorithmParameters) throws GRBException {
		this.algorithmParameters = algorithmParameters;
		
		//Reads the problem from file and disables output while importing.
		GRBEnv env = new GRBEnv();
		env.set(GRB.IntParam.OutputFlag, 0);
		model = new GRBModel(env, problemPath);
		model.set(GRB.IntParam.OutputFlag, 1);
		
		//Ensures minimization problem
		if (model.get(GRB.IntAttr.ModelSense) == GRB.MAXIMIZE) {
			GRBLinExpr objExpr = (GRBLinExpr) model.getObjective();
			GRBLinExpr newObjExpr = new GRBLinExpr(); 
			for (int i = 0; i < objExpr.size(); i++) {
				newObjExpr.addTerm(-objExpr.getCoeff(i), objExpr.getVar(i));
			}
			newObjExpr.addConstant(-objExpr.getConstant());
			model.setObjective(newObjExpr, GRB.MINIMIZE);
			model.update();
		}
		
		//Reads all model variables from the model and creates a Variable array with their objective coefficient.
		nominalModelVariables = model.getVars();
		nominalVariables = new Variable[nominalModelVariables.length];
		for (int i = 0; i < nominalModelVariables.length; i++) {
			GRBVar grbVar = nominalModelVariables[i];
			double objCoeff = grbVar.get(DoubleAttr.Obj);
			nominalVariables[i] = new Variable(i, grbVar, objCoeff);
		}
		
		algorithmParameters.applyParameters(this);
	}
	
	/**
	 * Resets the problem information to an unsolved state.
	 */
	protected void resetProblem() throws GRBException {
		primalBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
		dualBound = AbstractAlgorithm.DEFAULT_DUAL_BOUND;
		isStatusInfeasible = false;
		isStatusOptimal = false;
		nominalVariablesSolutionValues = null;
	}	
	
	/**
	 * Solves the model using Gurobi.
	 */
	void solve(Optional<Double> timeLimit) throws IOException, GRBException {
		//Sets the problem to an unsolved state for the case that it is reused.
		this.resetProblem();
		
		//Tries to solve the model. If an exception occurs, the problem is marked to be unsolved with
		//default dual and primal bounds. If the problem is only a subproblem then the algorithm may continue.
		try {
			if (timeLimit.isPresent()) {
				model.set(GRB.DoubleParam.TimeLimit, timeLimit.get());
			}
			model.optimize();
		} catch (GRBException e) {
			String output = "\n\nEncountered GRBException when solving model.\n"
					+"Primal and dual bound are set to default, infeasibility and optimality status to false.\n"
					+e.toString();
			AbstractAlgorithm.writeOutput(output, algorithmParameters);
			return;
		}
		
		//Tries to query the primal bound after the problem is terminated.
		//Primal bound is default if this fails.
		try {
			primalBound = model.get(GRB.DoubleAttr.ObjVal);
		} catch (GRBException e) {}

		//Tries to query the dual bound after the problem is terminated.
		//Dual bound is default if this fails.
		try {
			dualBound = model.get(GRB.DoubleAttr.ObjBound);
		} catch (GRBException e) {}
		
		//Checks whether the problem was solved to optimality after the problem is terminated.
		//Status is false if this fails.
		try {
			if (model.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL) {
				isStatusOptimal = true;
			}
		} catch (GRBException e) {}
		
		//Checks whether the problem is infeasible after the problem is terminated.
		//Status is false if this fails.
		try {
			if (model.get(GRB.IntAttr.Status) == GRB.Status.INFEASIBLE) {
				isStatusInfeasible = true;
				dualBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
			}
		} catch (GRBException e) {}
		
		//Tries to query the solution values after the problem is terminated.
		//Values are null if this fails.
		try {
			nominalVariablesSolutionValues = model.get(GRB.DoubleAttr.X, nominalModelVariables);
		} catch (GRBException e) { }
	}
	
	/**
	 * Returns the Gurobi model.
	 */
	public GRBModel getModel() {
		return model;
	}
	
	/**
	 * Returns the array of variables.
	 */
	Variable[] getNominalVariables(){
		return nominalVariables;
	}
	
	/**
	 * Returns the array of variables of the model.
	 */
	GRBVar[] getNominalModelVariables(){
		return nominalModelVariables;
	}
	
	/**
	 * Returns the solutions values.
	 */
	double[] getNominalVariablesSolutionValues() {
		return nominalVariablesSolutionValues;
	}
	
	/**
	 * Returns the primal bound computed by Gurobi.
	 */
	double getPrimalBound() {
		return primalBound;
	}
	
	/**
	 * Returns the dual bound computed by Gurobi.
	 */
	double getDualBound() {
		return dualBound;
	}
	
	/**
	 * Returns the status whether the problem is infeasible.
	 */
	boolean isStatusInfeasible() {
		return isStatusInfeasible;
	}
	
	/**
	 * Returns the status whether the problem was solved to optimality.
	 */
	boolean isStatusOptimal() {
		return isStatusOptimal;
	}	
}
