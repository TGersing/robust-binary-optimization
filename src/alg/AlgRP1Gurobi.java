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
import util.ConflictGraph;
import util.PossibleZ;
import util.Variable;

/**
 * This class implements the formulation RP1 from the paper "Strong formulations of robust mixed 0–1 programming" by Atamturk.
 * 
 * @author Timo Gersing
 */
public class AlgRP1Gurobi extends AbstractAlgorithm implements RobustAlgorithm{
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
	public AlgRP1Gurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, RobustAlgorithmStrategies robustAlgorithmStrategies) throws GRBException, IOException {
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
				+ "##### Solving Problem via RP1"
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
	 * Reformulates the given model into RP1.
	 */
	private void reformulateModel() throws GRBException, IOException {
		if (robustAlgorithmStrategies.getCliqueStrategy() == CliqueStrategy.CLIQUES || robustAlgorithmStrategies.getFilterStrategy() == FilterStrategy.FILTERED_Z) {
			String output = "\n###########################\n"
					+ "##### Start Preprocessing"
					+ "\n###########################";
			writeOutput(output);
		}

		//Computes a conflict graph if the corresponding strategy is chosen.
		ConflictGraph conflictGraph = null;
		if (robustAlgorithmStrategies.getCliqueStrategy() == CliqueStrategy.CLIQUES) {
			conflictGraph = new ConflictGraph(robustProblem.getModel(), robustProblem.getUncertainModelVariables(), algorithmParameters);
		}
		//Computes a list of possible optimal choices for z respecting the chosen filtering and clique strategies.
		List<PossibleZ> possibleZs = computePossibleZs(robustProblem.getUncertainVariables(), robustProblem.getGamma(), robustAlgorithmStrategies, conflictGraph, algorithmParameters);
		
		//Queries the model and uncertain variables from the given robust problem.
		GRBModel model = robustProblem.getModel();		
		Variable[] uncertainVariables = robustProblem.getUncertainVariables();

		//Adds the variables z, p, lambda, and omega to the model
		GRBVar z = model.addVar(0, Double.MAX_VALUE, robustProblem.getGamma(), GRB.CONTINUOUS, "z");
		GRBVar[] p = new GRBVar[uncertainVariables.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = model.addVar(0, Double.MAX_VALUE, 1, GRB.CONTINUOUS, "p"+i);
		}
		GRBVar[] lambda = new GRBVar[possibleZs.size()];
		for (int k = 0; k < possibleZs.size(); k++) {
			lambda[k] = model.addVar(0, Double.MAX_VALUE, 0, GRB.CONTINUOUS, "lambda"+k);
		}
		GRBVar[][] omega = new GRBVar[uncertainVariables.length][possibleZs.size()];
		for (int i = 0; i < uncertainVariables.length; i++) {
			for (int k = 0; k < possibleZs.size(); k++) {
				omega[i][k] = model.addVar(0, Double.MAX_VALUE, 0, GRB.CONTINUOUS, "omega"+i+","+k);
			}
		}
		
		//Alters the model's objective function by adding Gamma*z and the sum over all variables p.
		GRBLinExpr objLinExpr = (GRBLinExpr) model.getObjective();
		for (int i = 0; i < p.length; i++) {
			objLinExpr.addTerm(1, p[i]);
		}
		objLinExpr.addTerm(robustProblem.getGamma(), z);
		model.setObjective(objLinExpr);
		
		
		//The sum over all variables lambda has to be equal to one.
		GRBLinExpr chooseLambda = new GRBLinExpr();
		for (int k = 0; k < possibleZs.size(); k++) {
			chooseLambda.addTerm(1, lambda[k]);
		}
		model.addConstr(chooseLambda, GRB.EQUAL, 1, "chooseLambda");
		
		
		//Variable x[i] is equal to the sum over all omega[i][].
		for (int i = 0; i < uncertainVariables.length; i++) {
			GRBLinExpr chooseX = new GRBLinExpr();
			for (int k = 0; k < possibleZs.size(); k++) {
				chooseX.addTerm(1, omega[i][k]);
			}
			model.addConstr(chooseX, GRB.EQUAL, uncertainVariables[i].getModelVariable(), "chooseX"+i);
		}

		
		//We can only choose omega[][k] if lambda[k] is chosen.
		for (int i = 0; i < uncertainVariables.length; i++) {
			for (int k = 0; k < possibleZs.size(); k++) {
				model.addConstr(lambda[k], GRB.GREATER_EQUAL, omega[i][k], "omegaLambdaLink"+i+","+k);
			}
		}
		
		//z is at least the sum over all variables lambda multiplied by the corresponding deviations.
		GRBLinExpr chooseZ = new GRBLinExpr();
		for (int k = 0; k < possibleZs.size(); k++) {
			chooseZ.addTerm(possibleZs.get(k).getValue(), lambda[k]);
		}
		model.addConstr(z, GRB.GREATER_EQUAL, chooseZ, "chooseZ");
		
		//p[i] is at least the sum over all variables omega[i][k] multiplied by the maximum of zero and the deviation of i minus the corresponding deviation to k.
		for (int i = 0; i < uncertainVariables.length; i++) {
			GRBLinExpr chooseP = new GRBLinExpr();
			for (int k = 0; k < possibleZs.size(); k++) {
				if (uncertainVariables[i].getDeviation() > possibleZs.get(k).getValue()) {
					chooseP.addTerm(uncertainVariables[i].getDeviation()-possibleZs.get(k).getValue(), omega[i][k]);
				}
				else {
					break;
				}
			}
			model.addConstr(p[i], GRB.GREATER_EQUAL, chooseP, "chooseP"+i);
		}
		
		model.update();
	}
}