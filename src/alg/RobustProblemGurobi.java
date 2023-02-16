package alg;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import alg.AbstractAlgorithm.AlgorithmParameters;
import gurobi.GRB;
import gurobi.GRBException;
import gurobi.GRBVar;
import util.CliquePartitioning;
import util.Variable;

/**
 * This class extends a nominal problem to a robust problem.
 * 
 * @author Timo Gersing
 */
class RobustProblemGurobi extends ProblemGurobi {
	/**
	 * Path to the file stating the robustness components of the problem.
	 */
	protected String robustPath;

	/**
	 * The robustness budget.
	 */
	protected double Gamma;
	
	/**
	 * An array containing the uncertain variables of the problem.
	 */
	protected Variable[] uncertainVariables;

	/**
	 * An array containing the uncertain variables of the model.
	 */
	protected GRBVar[] uncertainModelVariables;

	/**
	 * The uncertain variables' solution values.
	 */
	protected double[] uncertainVariablesSolutionValues;
	
	/**
	 * A partitioning of the uncertain variables into cliques.
	 */
	protected CliquePartitioning cliquePartitioning;
	
	/**
	 * Constructor importing a problem from a file. Reads the value of Gamma and the deviations from a file.
	 */
	RobustProblemGurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters) throws GRBException, IOException {
		super(problemPath, algorithmParameters);
		
		this.robustPath = robustPath;
		
		Map<String, Variable> nameToVariable = new HashMap<String, Variable>();
		for (Variable variable : nominalVariables) {
			nameToVariable.put(variable.getModelVariable().get(GRB.StringAttr.VarName), variable);
		}
		
		BufferedReader robustReader = new BufferedReader(new FileReader(robustPath));
		String line = robustReader.readLine();
		String[] lineArray = line.split(":");
		Gamma = Double.parseDouble(lineArray[1]);
		
		List<Variable> uncertainVariablesList = new ArrayList<Variable>();
		while ((line = robustReader.readLine()) != null) {
			lineArray = line.split(":");
			Variable variable = nameToVariable.get(lineArray[0]);
			variable.setDeviation(Double.parseDouble(lineArray[1]));
			uncertainVariablesList.add(variable);
		}
		robustReader.close();
		
		//The list of variables with uncertain objective coefficients is sorted non-decreasing with respect to the deviations.
		Collections.sort(uncertainVariablesList);
		uncertainVariables = new Variable[uncertainVariablesList.size()];
		uncertainVariablesList.toArray(uncertainVariables);
		
		uncertainModelVariables = new GRBVar[uncertainVariables.length];
		for (int i = 0; i < uncertainVariables.length; i++) {
			uncertainModelVariables[i] = uncertainVariables[i].getModelVariable();
			uncertainVariables[i].setUncertainIndex(i);
		}
	}
	
	/**
	 * Sets the clique partitioning.
	 */
	void setCliquePartitioning(CliquePartitioning cliquePartitioning) {
		this.cliquePartitioning = cliquePartitioning;
	}
		
	/**
	 * Resets the problem to an unsolved state.
	 */
	@Override
	protected void resetProblem() throws GRBException {
		super.resetProblem();
		uncertainVariablesSolutionValues = null;
	}

	/**
	 * Solves the model using Gurobi.
	 */
	@Override
	void solve(Optional<Double> timeLimit) throws IOException, GRBException {
		super.solve(timeLimit);
		//Tries to query the solution values for the uncertain variables after the problem is terminated.
		//Values are null if this fails.
		try {
			uncertainVariablesSolutionValues = model.get(GRB.DoubleAttr.X, uncertainModelVariables);
		} catch (GRBException e) { }
	}
	
	/**
	 * Returns the path to the imported robustness components of the problem.
	 */
	String getRobustPath() {
		return robustPath;
	}
	
	/**
	 * Returns the robustness budget.
	 */
	double getGamma() {
		return Gamma;
	}
	
	/**
	 * Returns the array of uncertain variables.
	 */
	Variable[] getUncertainVariables() {
		return uncertainVariables;
	}
	
	/**
	 * Returns the array of uncertain variables in the model.
	 */
	GRBVar[] getUncertainModelVariables() {
		return uncertainModelVariables;
	}
	
	/**
	 * Returns the solution values of the uncertain variables.
	 */
	double[] getUncertainVariablesSolutionValues() {
		return uncertainVariablesSolutionValues;
	}
}
