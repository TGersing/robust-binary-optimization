package alg;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

import alg.RobustAlgorithm.RobustAlgorithmStrategies.CliqueStrategy;
import alg.RobustAlgorithm.RobustAlgorithmStrategies.FilterStrategy;
import gurobi.GRB;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import util.CliquePartitioning;
import util.ConflictGraph;
import util.Variable;

/**
 * This class implements the standard reformulation by Bertsimas and Sim. The formulation is improved
 * as described in "A Branch & Bound Algorithm for Robust Binary Optimization with Budget Uncertainty".
 * 
 * @author Timo Gersing
 */
public class AlgReformulationGurobi extends AbstractAlgorithm implements RobustAlgorithm{
	/**
	 * The robust problem to be solved.
	 */
	private RobustProblemGurobi robustProblem;
	
	/**
	 * Specifies strategies.
	 */
	private RobustAlgorithmStrategies robustAlgorithmStrategies;

	/**
	 * Constructor receiving problem paths as well as parameters and strategies.
	 */
	public AlgReformulationGurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, RobustAlgorithmStrategies robustAlgorithmStrategies) throws GRBException, IOException {
		super(algorithmParameters);
		this.robustProblem = new RobustProblemGurobi(problemPath, robustPath, algorithmParameters);
		
		this.robustAlgorithmStrategies = robustAlgorithmStrategies;
	}

	/**
	 * Executes the algorithm by reformulating the problem and solving it within the remaining time limit.
	 */
	@Override
	protected void executeAlgorithm() throws IOException, GRBException {
		String output = "######################################################\n"
				+ "##### Solving Problem via Reformulation"
				+ "\n######################################################";
		writeOutput(output);

		reformulateModel();
		
		robustProblem.solve(getRemainingTime());
		primalBound = robustProblem.getPrimalBound();
		dualBound = robustProblem.getDualBound();
		
		//Stores best solution found
		solution = new LinkedHashMap<Variable, Double>(robustProblem.getNominalVariables().length);
		for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
			solution.put(robustProblem.getNominalVariables()[i], robustProblem.getNominalVariablesSolutionValues()[i]);
		}
	}

	/**
	 * Reformulates the given model.
	 */
	private void reformulateModel() throws GRBException, IOException {
		//Computes a conflict graph and clique partitioning if the corresponding strategy is chosen
		if (robustAlgorithmStrategies.getCliqueStrategy() == CliqueStrategy.CLIQUES || robustAlgorithmStrategies.getFilterStrategy() == FilterStrategy.FILTERED_Z) {
			String output = "\n###########################\n"
					+ "##### Start Preprocessing"
					+ "\n###########################";
			writeOutput(output);
		}
		ConflictGraph conflictGraph = null;
		CliquePartitioning cliquePartitioning = null;
		if (robustAlgorithmStrategies.getCliqueStrategy() == RobustAlgorithmStrategies.CliqueStrategy.CLIQUES) {
			conflictGraph = new ConflictGraph(robustProblem.getModel(), robustProblem.getUncertainModelVariables(), algorithmParameters);
			cliquePartitioning = new CliquePartitioning(robustProblem.getUncertainVariables(), conflictGraph, algorithmParameters);
		}
		//Computes the highest possible optimal choice for z respecting the chosen filtering and clique strategies.
		double highestPossibleZ = computeHighestPossibleZ(robustProblem.getUncertainVariables(), robustProblem.getGamma(), robustAlgorithmStrategies, algorithmParameters);

		//Queries the model and uncertain variables from the given robust problem.
		GRBModel model = robustProblem.getModel();		
		Variable[] uncertainVariables = robustProblem.getUncertainVariables();
		
		//Adds the variable z and the substituted variables p from the (clique) reformulation.
		GRBVar z = model.addVar(0, Double.MAX_VALUE, robustProblem.getGamma(), GRB.CONTINUOUS, "z");
		//The size of p may vary depending on whether we merge variables into cliques.
		GRBVar[] pPrime;
		if (cliquePartitioning != null) {
			pPrime = new GRBVar[cliquePartitioning.getCliques().size()];
		}
		else {
			pPrime = new GRBVar[uncertainVariables.length];
		}
		for (int i = 0; i < pPrime.length; i++) {
			pPrime[i] = model.addVar(0, Double.MAX_VALUE, 1, GRB.CONTINUOUS, "p"+i);
		}
		
		//Alters the objective function by adding Gamma*z and the sum over all variables p.
		GRBLinExpr objLinExpr = (GRBLinExpr) model.getObjective();
		for (int i = 0; i < pPrime.length; i++) {
			objLinExpr.addTerm(1, pPrime[i]);
		}
		for (Variable uncertainVariable : uncertainVariables) {
			if (uncertainVariable.getDeviation() > highestPossibleZ) {
				objLinExpr.addTerm(uncertainVariable.getDeviation() - highestPossibleZ, uncertainVariable.getModelVariable());
			}
		}
		objLinExpr.addTerm(robustProblem.getGamma(), z);
		model.setObjective(objLinExpr);
		
				
		//Adds the constraint of the form z+p>=d*x (using cliques if the option is chosen correspondingly), with d being the deviation.
		if (cliquePartitioning != null) {
			int pIndex = 0;
			for (List<Integer> clique : cliquePartitioning.getCliques()) {
				GRBLinExpr robustnessExpr = new GRBLinExpr();
				robustnessExpr.addTerm(1, pPrime[pIndex]);
				robustnessExpr.addTerm(1, z);
				for (int varIndex : clique) {
					Variable var = uncertainVariables[varIndex];
					robustnessExpr.addTerm(-Math.min(var.getDeviation(), highestPossibleZ), var.getModelVariable());
				}
				model.addConstr(robustnessExpr, GRB.GREATER_EQUAL, 0, "");
				pIndex++;
			}
		}
		else {
			for (int i = uncertainVariables.length-1; i >= 0; i--) {
				Variable var = uncertainVariables[i];
				GRBLinExpr robustnessExpr = new GRBLinExpr();
				robustnessExpr.addTerm(1, pPrime[i]);
				robustnessExpr.addTerm(1, z);
				robustnessExpr.addTerm(-Math.min(var.getDeviation(), highestPossibleZ), var.getModelVariable());
				model.addConstr(robustnessExpr, GRB.GREATER_EQUAL, 0, "");
			}
		}
		
		model.update();
	}
}