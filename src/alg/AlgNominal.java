package alg;

import java.io.IOException;
import java.util.LinkedHashMap;

import gurobi.GRBException;
import util.Variable;

/**
 * This class solves nominal problems containing no uncertainty.
 * 
 * @author Timo Gersing
 */
public class AlgNominal extends AbstractAlgorithm{
	/**
	 * The problem to be solved.
	 */
	private ProblemGurobi problem;
	
	/**
	 * Constructor receiving problem path and parameters.
	 */
	public AlgNominal(String problemPath , AlgorithmParameters algorithmParameters) throws IOException, GRBException {
		super(algorithmParameters);
		this.problem = new ProblemGurobi(problemPath, algorithmParameters);
	}

	/**
	 * Tries to solve the problem. If present, we set the remaining time as a limit.
	 */
	@Override
	protected void executeAlgorithm() throws IOException, GRBException {
		problem.solve(getRemainingTime());
		setPrimalBound(problem.getPrimalBound());
		setDualBound(problem.getDualBound());
		
		//Stores best solution found, if available
		if (problem.getNominalVariablesSolutionValues() != null) {
			solution = new LinkedHashMap<Variable, Double>(problem.getNominalVariables().length);
			for (int i = 0; i < problem.getNominalVariables().length; i++) {
				solution.put(problem.getNominalVariables()[i], problem.getNominalVariablesSolutionValues()[i]);
			}
		}

	}
}