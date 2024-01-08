package alg;

import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.TreeSet;

import alg.RobustAlgorithm.RobustAlgorithmStrategies.CliqueStrategy;
import alg.RobustAlgorithm.RobustAlgorithmStrategies.FilterStrategy;
import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import util.CliquePartitioning;
import util.ConflictGraph;
import util.Variable;

/**
 * This class implements the standard reformulation by Bertsimas and Sim
 * together with submodular cuts, as proposed by Joung and Park.
 * The underlying standard reformulation is improved as described in
 * "A Branch and Bound Algorithm for Robust Binary Optimization with Budget Uncertainty".
 * 
 * @author Timo Gersing
 */
public class AlgSubmodularCuts extends AbstractAlgorithm implements RobustAlgorithm{
	/**
	 * The robust problem to be solved.
	 */
	private RobustProblemGurobi robustProblem;
	
	/**
	 * Specifies strategies.
	 */
	private RobustAlgorithmStrategies robustAlgorithmStrategies;
	
	/**
	 * The (substituted) p variables from the Bertsimas Sim reformulation.
	 */
	private GRBVar[] pPrime;
	
	/**
	 * The z variable from the Bertsimas Sim reformulation.
	 */
	private GRBVar z;
	
	/**
	 * The highest possible optimal value of z, used for strengthening the formulation. 
	 */
	private double highestPossibleZ;

	/**
	 * Constructor receiving problem paths as well as parameters and strategies.
	 */
	public AlgSubmodularCuts(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, RobustAlgorithmStrategies robustAlgorithmStrategies) throws GRBException, IOException {
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
				+ "##### Solving Problem via Submodular Cuts"
				+ "\n######################################################";
		writeOutput(output);

		reformulateModel();
		
		robustProblem.solve(getRemainingTime());
		setPrimalBound(robustProblem.getPrimalBound());
		setDualBound(robustProblem.getDualBound());
		
		//Stores best solution found, if available
		if (robustProblem.getNominalVariablesSolutionValues() != null) {
			solution = new LinkedHashMap<Variable, Double>(robustProblem.getNominalVariables().length);
			for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
				solution.put(robustProblem.getNominalVariables()[i], robustProblem.getNominalVariablesSolutionValues()[i]);
			}
		}
	}

	/**
	 * Reformulates the given model.
	 */
	private void reformulateModel() throws GRBException, IOException {
		//Computes a conflict graph and clique partitioning if the corresponding strategy is chosen
		if (robustAlgorithmStrategies.getCliqueStrategy() == CliqueStrategy.CLIQUES_ENABLE || robustAlgorithmStrategies.getFilterStrategy() == FilterStrategy.FILTERINGZ_ENABLE) {
			String output = "\n###########################\n"
					+ "##### Start Preprocessing"
					+ "\n###########################";
			writeOutput(output);
		}
		ConflictGraph conflictGraph = null;
		CliquePartitioning cliquePartitioning = null;
		if (robustAlgorithmStrategies.getCliqueStrategy() == RobustAlgorithmStrategies.CliqueStrategy.CLIQUES_ENABLE) {
			conflictGraph = new ConflictGraph(robustProblem.getModel(), robustProblem.getUncertainModelVariables(), algorithmParameters);
			cliquePartitioning = new CliquePartitioning(robustProblem.getUncertainVariables(), conflictGraph, algorithmParameters);
		}
		
		//Computes the highest possible optimal choice for z respecting the chosen filtering and clique strategies.
		highestPossibleZ = computeHighestPossibleZ(robustProblem.getUncertainVariables(), robustProblem.getGamma(), robustAlgorithmStrategies, algorithmParameters);

		//Queries the model and uncertain variables from the given robust problem.
		GRBModel model = robustProblem.getModel();		
		Variable[] uncertainVariables = robustProblem.getUncertainVariables();
		
		//Adds the variable z and the substituted variables p from the (clique) reformulation.
		z = model.addVar(0, Double.MAX_VALUE, robustProblem.getGamma(), GRB.CONTINUOUS, "z");
		//The size of p may vary depending on whether we merge variables into cliques.
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
		objLinExpr.addTerm(robustProblem.getGamma(), z);
		//If a variable's deviation is higher than the highest possible optimal value for z,
		//then we already add the difference to the obj coeff of the variable.
		for (Variable uncertainVariable : uncertainVariables) {
			if (uncertainVariable.getDeviation() > highestPossibleZ) {
				objLinExpr.addTerm(uncertainVariable.getDeviation() - highestPossibleZ, uncertainVariable.getModelVariable());
			}
		}
		model.setObjective(objLinExpr);
		
				
		//Adds the constraint of the form z+p>=d*x (using cliques if the option is chosen correspondingly).
		//Here, d is the minimum of the deviation and the highest possible value for z
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
		//Adds the submodular cut callback
		//PreCrush has to be set to 1 for the separation.
		model.set(GRB.IntParam.PreCrush, 1);
		model.setCallback(new CutCallback());
		model.update();
	}
	
	
	/**
	 * Obtains an LP solution and returns a submodular cut.
	 */
	private GRBLinExpr computeSubmodularCut(double[] uncertainVariablesSolutionValues) {
		double Gamma = robustProblem.getGamma();
		if (Gamma <= 0) {
			return null;
		}
		
		GRBLinExpr cutExpr = new GRBLinExpr();
		
		//Sorts non-increasing with respect to solution value and breaks ties with deviations (non-decreasing) (=sorting by index).
		TreeSet<Variable> sortedUncertainVariables = new TreeSet<>((var1, var2) -> {
			if (uncertainVariablesSolutionValues[var1.getUncertainIndex()] == uncertainVariablesSolutionValues[var2.getUncertainIndex()]) {
				return Double.compare(var1.getUncertainIndex(), var2.getUncertainIndex());
			}
			else {
				return Double.compare(uncertainVariablesSolutionValues[var2.getUncertainIndex()], uncertainVariablesSolutionValues[var1.getUncertainIndex()]);
			}
		});
		
		for (Variable variable : robustProblem.getUncertainVariables()) {
			sortedUncertainVariables.add(variable);
		}
		
		//The set of variables with the highest deviations of the already considered vars (=sorting by index).
		TreeSet<Variable> highestDevVariables = new TreeSet<>((var1, var2) -> Double.compare(var1.getUncertainIndex(), var2.getUncertainIndex()));
		
		//Iterate over all uncertain variables w.r.t. the ordering above (LP solution value)
		Iterator<Variable> iterator = sortedUncertainVariables.iterator();
		while (iterator.hasNext()) {
			//Computes the coefficient for the next uncertain variable
			Variable uncertainVariable = iterator.next();
			//Adds the current variable to the variables with the highest deviation (might be removed immediately)
			highestDevVariables.add(uncertainVariable);
			double coeff = Math.min(highestPossibleZ, uncertainVariable.getDeviation());
			//If there now exist more than floor(Gamma) variables with the highest deviation, then the coefficient needs to be reduced
			if (highestDevVariables.size() > Math.floor(Gamma)) {
				//If there now exist not more than ceil(Gamma) deviating variables, then the variable with the lowest deviation will now only deviate fractionally
				//We reduce the coefficient by the difference of the full versus fractional deviation
				if (highestDevVariables.size() <= Math.ceil(Gamma)) {
					Variable fractionalVariable = highestDevVariables.first();
					coeff -= (Math.ceil(Gamma)-Gamma)*Math.min(highestPossibleZ, fractionalVariable.getDeviation());
				}
				//If there exist more than ceil(Gamma) variables with the highest deviation, then one variable has to be removed
				else {
					//The variable with the lowest deviation is removed
					Variable leavingVariable = highestDevVariables.pollFirst();
					//If Gamma is fractional, then the leaving variable deviated fractionally. This deviation is removed.
					//Moreover, the variable with the second lowest deviation now deviates fractionally. 
					if (Math.ceil(Gamma) != Gamma) {
						Variable fractionalVariable = highestDevVariables.first();
						coeff -= (Math.ceil(Gamma)-Gamma)*Math.min(highestPossibleZ, fractionalVariable.getDeviation());
						coeff -= (Gamma-Math.floor(Gamma))*Math.min(highestPossibleZ, leavingVariable.getDeviation());
					}
					//If Gamma is integer, the leaving variables deviation is subtracted.
					else {
						coeff -= Math.min(highestPossibleZ, leavingVariable.getDeviation());
					}
				}
			}
			if (coeff > 0) {
				cutExpr.addTerm(coeff, uncertainVariable.getModelVariable());
			}
		}
		
		//Adds p and z variables such that the robustness term is greater than the submodular scenario.
		cutExpr.addTerm(-Gamma, z);
		for (GRBVar pVar : pPrime) {
			cutExpr.addTerm(-1, pVar);
		}
		return cutExpr;
	}
	
	/**
	 * Callback for computing the primal dual integral and the submodular cuts
	 */
	private class CutCallback extends GRBCallback {
		@Override
		protected void callback() {
			if (where == GRB.CB_MIP) {
				try {
					setPrimalBound(getDoubleInfo(GRB.CB_MIP_OBJBST));
					setDualBound(getDoubleInfo(GRB.CB_MIP_OBJBND));
				} catch (GRBException e) {}
			}
			

			try {
				//We separate submodular cuts while we are in the root node.
				if ((where == GRB.CB_MIPNODE && getDoubleInfo(GRB.CB_MIPNODE_NODCNT) == 0 && getIntInfo(GRB.CB_MIPNODE_STATUS) == GRB.Status.OPTIMAL)) {
					//Queries the values of the uncertain variables.
					double[] uncertainVarValues = getNodeRel(robustProblem.getUncertainModelVariables());
					
					//Computes and adds the submodular cut.
					GRBLinExpr cutExpr = computeSubmodularCut(uncertainVarValues);
					addCut(cutExpr, GRB.LESS_EQUAL, 0);
				}
			} catch (GRBException e) {
				e.printStackTrace();
			}
		}
	}
}