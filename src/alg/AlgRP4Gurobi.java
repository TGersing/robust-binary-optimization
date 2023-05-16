package alg;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import alg.RobustAlgorithm.RobustAlgorithmStrategies.CliqueStrategy;
import alg.RobustAlgorithm.RobustAlgorithmStrategies.FilterStrategy;
import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBConstr;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import util.ConflictGraph;
import util.PossibleZ;
import util.Variable;

/**
 * This class implements the formulation RP4 from the paper "Strong formulations of robust mixed 0–1 programming" by Atamturk.
 * 
 * @author Timo Gersing
 */
public class AlgRP4Gurobi extends AbstractAlgorithm implements RobustAlgorithm{
	/**
	 * The robust problem to be solved.
	 */
	private RobustProblemGurobi robustProblem;
	
	/**
	 * Specifies strategies.
	 */
	private RobustAlgorithmStrategies robustAlgorithmStrategies;
	
	/**
	 * Omega variables used for reconstructing the nominal solution values.
	 */
	private GRBVar[][] omega;
	
	/**
	 * List of possible z for which omega is defined.
	 */
	List<PossibleZ> possibleZs;

	/**
	 * Constructor storing the given robust problem and strategies.
	 */
	public AlgRP4Gurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, RobustAlgorithmStrategies robustAlgorithmStrategies) throws GRBException, IOException {
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
				+ "##### Solving Problem via RP4"
				+ "\n######################################################";
		writeOutput(output);

		reformulateModel();
		
		robustProblem.solve(getRemainingTime());
		primalBound = robustProblem.getPrimalBound();
		dualBound = robustProblem.getDualBound();
		
		//Stores best solution found
		solution = new LinkedHashMap<Variable, Double>(robustProblem.getNominalVariables().length);
		for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
			double solutionValue = robustProblem.getNominalModelVariables()[i].get(GRB.DoubleAttr.LB);
			for (int k = 0; k < possibleZs.size(); k++) {
				solutionValue += omega[i][k].get(GRB.DoubleAttr.X);
			}
			solution.put(robustProblem.getNominalVariables()[i], solutionValue);
		}
	}

	

	/**
	 * Reformulates the given model into RP4.
	 */
	private void reformulateModel() throws GRBException, IOException {
		//Computes a conflict graph if the corresponding strategy is chosen
		if (robustAlgorithmStrategies.getCliqueStrategy() == CliqueStrategy.CLIQUES_ENABLE || robustAlgorithmStrategies.getFilterStrategy() == FilterStrategy.FILTERINGZ_ENABLE) {
			String output = "\n###########################\n"
					+ "##### Start Preprocessing"
					+ "\n###########################";
			writeOutput(output);
		}
		ConflictGraph conflictGraph = null;
		if (robustAlgorithmStrategies.getCliqueStrategy() == CliqueStrategy.CLIQUES_ENABLE) {
			conflictGraph = new ConflictGraph(robustProblem.getModel(), robustProblem.getUncertainModelVariables(), algorithmParameters);
		}
		//Computes a list of possible optimal choices for z respecting the chosen filtering and clique strategies.
		possibleZs = computePossibleZs(robustProblem.getUncertainVariables(), robustProblem.getGamma(), robustAlgorithmStrategies, conflictGraph, algorithmParameters);

		//Queries the model and variables from the given robust problem.
		GRBModel model = robustProblem.getModel();
		Variable[] nominalVariables = robustProblem.getNominalVariables();
		Variable[] uncertainVariables = robustProblem.getUncertainVariables();

		//Adds the variables z, p, lambda, and omega.
		GRBVar z = model.addVar(0, Double.MAX_VALUE, robustProblem.getGamma(), GRB.CONTINUOUS, "z");
		GRBVar[] p = new GRBVar[uncertainVariables.length];
		for (int i = 0; i < uncertainVariables.length; i++) {
			p[i] = model.addVar(0, Double.MAX_VALUE, 1, GRB.CONTINUOUS, "p"+i);
		}
		GRBVar[] lambda = new GRBVar[possibleZs.size()];
		for (int k = 0; k < possibleZs.size(); k++) {
			lambda[k] = model.addVar(0, Double.MAX_VALUE, 0, GRB.CONTINUOUS, "lambda"+k);
		}
		omega = new GRBVar[nominalVariables.length][possibleZs.size()];
		
		//To construct RP4, we have to alter the whole constraint matrix of the nominal problem.
		//First, we map the variables of the imported model to an index such that we can quickly replace it in the constraints with its corresponding variable omega. 
		HashMap<GRBVar, Integer> modelVariableToOmegaIndex = new HashMap<GRBVar, Integer>();
		for (int i = 0; i < nominalVariables.length; i++) {
			modelVariableToOmegaIndex.put(nominalVariables[i].getModelVariable(), i);
		}
		
		//Formulation RP4 can only deal with variables where the lower bound is zero.
		//If the imported problem contains variables not fulfilling this property then these have to be shifted.
		for (int i = 0; i < nominalVariables.length; i++) {
			Variable variable = nominalVariables[i];
			GRBVar modelVariable = variable.getModelVariable();
			double UB = modelVariable.get(GRB.DoubleAttr.UB) - modelVariable.get(GRB.DoubleAttr.LB);

			if (variable.getType() == Variable.Type.BINARY) {
				for (int k = 0; k < possibleZs.size(); k++) {
					omega[i][k] = model.addVar(0, 1, 0, GRB.BINARY, "omega"+i+","+k);
				}
			}
			else if (variable.getType() == Variable.Type.INTEGER) {
				for (int k = 0; k < possibleZs.size(); k++) {
					omega[i][k] = model.addVar(0, UB, 0, GRB.INTEGER, "omega"+i+","+k);
				}
			}
			else {
				for (int k = 0; k < possibleZs.size(); k++) {
					omega[i][k] = model.addVar(0, UB, 0, GRB.CONTINUOUS, "omega"+i+","+k);
				}
			}
		}
		model.update();
		
		//Add the objective function to the problem.
		//If the variable is shifted, we add the corresponding value via the constant term of the objective.
		double objConstant = 0;
		GRBLinExpr objExpr = new GRBLinExpr();
		for (int i = 0; i < nominalVariables.length; i++) {
			Variable nominalVariable = nominalVariables[i];
			for (int k = 0; k < possibleZs.size(); k++) {
				objExpr.addTerm(nominalVariable.getObjectiveCoefficient(), omega[i][k]);
			}
			if (nominalVariable.getModelVariable().get(GRB.DoubleAttr.LB) != 0) {
				objConstant += nominalVariable.getObjectiveCoefficient()*nominalVariable.getModelVariable().get(GRB.DoubleAttr.LB);
			}
		}
		objExpr.addConstant(objConstant);
		//Adds Gamma*z and the sum over all variables p.
		objExpr.addTerm(robustProblem.getGamma(), z);
		for (int i = 0; i < uncertainVariables.length; i++) {
			objExpr.addTerm(1, p[i]);
		}
		model.setObjective(objExpr, GRB.MINIMIZE);
		
		//The sum over all variables lambda has to be equal to one.
		GRBLinExpr chooseLambda = new GRBLinExpr();
		for (int k = 0; k < possibleZs.size(); k++) {
			chooseLambda.addTerm(1, lambda[k]);
		}
		model.addConstr(chooseLambda, GRB.EQUAL, 1, "chooseLambda");
		
		//z is at least the sum over all variables lambda multiplied by the corresponding deviations.
		GRBLinExpr chooseZ = new GRBLinExpr();
		for (int k = 0; k < possibleZs.size(); k++) {
			chooseZ.addTerm(possibleZs.get(k).getValue(), lambda[k]);
		}
		model.addConstr(z, GRB.GREATER_EQUAL, chooseZ, "chooseZ");
		
		//We can only choose omega[][k] if lambda[k] is chosen.
		for (int i = 0; i < nominalVariables.length; i++) {
			double UB = omega[i][0].get(GRB.DoubleAttr.UB);
			for (int k = 0; k < possibleZs.size(); k++) {
				GRBLinExpr rhs = new GRBLinExpr();
				rhs.addTerm(UB, lambda[k]);
				model.addConstr(omega[i][k], GRB.LESS_EQUAL, rhs, "");
			}
		}
		
		//p[i] is at least the sum over all variables omega[i][k] multiplied by the maximum of zero and the deviation of i minus the corresponding deviation to k.
		for (int i = 0; i < uncertainVariables.length; i++) {
			Variable uncertainVariable = uncertainVariables[i];
			GRBVar modelVariable = uncertainVariable.getModelVariable();
			if (uncertainVariable.getDeviation() > 0) {
				GRBLinExpr chooseP = new GRBLinExpr();
				int omegaIndex = modelVariableToOmegaIndex.get(modelVariable);
				if (uncertainVariable.getModelVariable().get(GRB.DoubleAttr.LB) != 0) {
					double constant = 0;
					for (int k = 0; k < possibleZs.size(); k++) {
						if (uncertainVariable.getDeviation() > possibleZs.get(k).getValue()) {
							chooseP.addTerm(uncertainVariable.getDeviation()-possibleZs.get(k).getValue(), omega[omegaIndex][k]);
							constant += uncertainVariable.getModelVariable().get(GRB.DoubleAttr.LB)*(uncertainVariable.getDeviation()-possibleZs.get(k).getValue());
						}
						else {
							break;
						}
					}
					chooseP.addConstant(constant);
				}
				else {
					for (int k = 0; k < possibleZs.size(); k++) {
						if (uncertainVariable.getDeviation() > possibleZs.get(k).getValue()) {
							chooseP.addTerm(uncertainVariable.getDeviation()-possibleZs.get(k).getValue(), omega[omegaIndex][k]);
						}
						else {
							break;
						}
					}
				}
				model.addConstr(p[i], GRB.GREATER_EQUAL, chooseP, "chooseP");
			}
		}
		
		//We store the constraints of the nominal problem in an array.
		GRBConstr[] constraints = model.getConstrs();
		
		//Replaces all constraints by one counterpart for every lambda. The original variables are replaced by their corresponding omega variables.
		for (GRBConstr constraint : constraints) {
			//Count the constant term implied by the shifted variables.
			double constant = 0;
			//Initialize the left hand side expression for the current constraint for each lambda.
			GRBLinExpr[] expressions = new GRBLinExpr[possibleZs.size()];
			for (int k = 0; k < possibleZs.size(); k++) {
				expressions[k] = new GRBLinExpr();
			}
			
			GRBLinExpr lhsExpr = model.getRow(constraint);
			for (int i = 0; i< lhsExpr.size(); i++) {
				GRBVar modelVariable = lhsExpr.getVar(i);
				double coeff = lhsExpr.getCoeff(i);
				//Count the constant term for fixed variables
				if (modelVariable.get(GRB.DoubleAttr.LB) != 0) {
					constant += coeff * modelVariable.get(GRB.DoubleAttr.LB);
				}
				//Replace the original variable by omega for every k.
				int omegaIndex = modelVariableToOmegaIndex.get(modelVariable);
				for (int k = 0; k < possibleZs.size(); k++) {
					expressions[k].addTerm(coeff, omega[omegaIndex][k]);
				}
			}
			
			//Add the new constraints to the problem.
			for (int k = 0; k < possibleZs.size(); k++) {
				if (constraint.get(GRB.CharAttr.Sense) == GRB.LESS_EQUAL) {
					GRBLinExpr rhsExpr = new GRBLinExpr();
					rhsExpr.addTerm(constraint.get(GRB.DoubleAttr.RHS)-constant, lambda[k]);
					model.addConstr(expressions[k], GRB.LESS_EQUAL, rhsExpr, "");
				}
				else {
					GRBLinExpr rhsExpr = new GRBLinExpr();
					rhsExpr.addTerm(constraint.get(GRB.DoubleAttr.RHS)-constant, lambda[k]);
					model.addConstr(expressions[k], GRB.GREATER_EQUAL, rhsExpr, "");
				}
			}
		}
		
		for (GRBConstr constr : constraints) {
			model.remove(constr);
		}
		model.setCallback(new PrimalDualCallback());
		model.update();
	}
	
	/**
	 * Callback for computing the primal dual integral
	 */
	private class PrimalDualCallback extends GRBCallback {
		@Override
		protected void callback() {
			if (where == GRB.CB_MIP) {
				try {
					primalDualIntegral.update(getDoubleInfo(GRB.CB_MIP_OBJBST), getDoubleInfo(GRB.CB_MIP_OBJBND), false);
				} catch (GRBException e) {}
			}
		}
	}
}