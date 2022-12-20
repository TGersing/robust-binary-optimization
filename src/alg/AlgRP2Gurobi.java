package alg;

import java.io.IOException;
import java.util.LinkedHashMap;

import util.Variable;
import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

/**
 * This class implements the formulation RP2 from the paper "Strong formulations of robust mixed 0–1 programming" by Atamturk.
 * 
 * @author Timo Gersing
 */
public class AlgRP2Gurobi extends AbstractAlgorithm implements RobustAlgorithm{
	/**
	 * The robust problem to be solved.
	 */
	private RobustProblemGurobi robustProblem;
	
	/**
	 * Variables z,p of the reformulation, which need to be defined globally for the separation callback.
	 */
	private GRBVar z;
	private GRBVar[] p;

	
	/**
	 * Constructor storing the given robust problem.
	 */
	public AlgRP2Gurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters) throws GRBException, IOException {
		super(algorithmParameters);
		this.robustProblem = new RobustProblemGurobi(problemPath, robustPath, algorithmParameters);
	}
	
	/**
	 * Executes the algorithm by reformulating the problem and solving it within the remaining time limit.
	 */
	@Override
	protected void executeAlgorithm() throws IOException, GRBException {
		String output = "######################################################\n"
				+ "##### Solving Problem via RP2"
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
	 * Reformulates the given model into RP2.
	 */
	private void reformulateModel() throws GRBException {
		//Queries the model and uncertain variables from the given robust problem.
		GRBModel model = robustProblem.getModel();		
		Variable[] uncertainVariables = robustProblem.getUncertainVariables();
		
		//Adds the variables z, p.
		z = model.addVar(0, Double.MAX_VALUE, robustProblem.getGamma(), GRB.CONTINUOUS, "z");
		p = new GRBVar[uncertainVariables.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = model.addVar(0, Double.MAX_VALUE, 1, GRB.CONTINUOUS, "p"+i);
		}
		
		//Alters the objective function by adding Gamma*z and the sum over all variables p.
		GRBLinExpr objLinExpr = (GRBLinExpr) model.getObjective();
		for (int i = 0; i < p.length; i++) {
			objLinExpr.addTerm(1, p[i]);
		}
		objLinExpr.addTerm(robustProblem.getGamma(), z);
		model.setObjective(objLinExpr);
		
				
		//Adds the constraints from the standard reformulation of the form z+p>=d*x, with d being the deviation.
		for (int i = uncertainVariables.length-1; i >= 0; i--) {
			Variable var = uncertainVariables[i];
			GRBLinExpr robustnessExpr = new GRBLinExpr();
			robustnessExpr.addTerm(1, p[i]);
			robustnessExpr.addTerm(1, z);
			robustnessExpr.addTerm(-var.getDeviation(), var.getModelVariable());
			model.addConstr(robustnessExpr, GRB.GREATER_EQUAL, 0, "");
		}
		
		//Adds the callback for separation of valid inequalities.
		//PreCrush has to be set to 1 for the separation.
		model.set(GRB.IntParam.PreCrush, 1);
		model.setCallback(new Callback());
		
		model.update();
	}
	
	/**
	 * Implements the separation algorithm searching for violated cuts by solving a shortest path problem.
	 * We compute a shortest path from the super source to nodes corresponding to the uncertain variables.
	 * Every negative length path corresponds to a violated cut.
	 * The graph is leveled with respect to the deviations and thus acyclic.
	 */
	private class Callback extends GRBCallback {
		@Override
		protected void callback() {
			try {
				//Cuts can be added when the callback is called from a MIP node and the node's LP relaxation is solved to optimality. 
				if (where == GRB.CB_MIPNODE && getIntInfo(GRB.CB_MIPNODE_STATUS) == GRB.Status.OPTIMAL) {
					Variable[] uncertainVariables = robustProblem.getUncertainVariables();
					
					//Distance from the super source to the variables.
					double[] distances = new double[uncertainVariables.length];
					//Predecessor of each variable.
					int[] predecessors = new int[uncertainVariables.length];
					//The value of z in the solution of the current node.
					double zVal = getNodeRel(z);
					
					//Computes the shortest path to each variable, starting with the variable that has the lowest deviation.
					for (int j = 0; j < uncertainVariables.length; j++) {
						//The value of the uncertain variable in the solution of the current node.
						double xVal = getNodeRel(uncertainVariables[j].getModelVariable());
						//If the variable is zero then it is not interesting for finding cuts.
						if (xVal <= algorithmParameters.getIntegerFeasibilityTolerance()) {
							distances[j] = Double.POSITIVE_INFINITY;
						}
						else {
							//The value of p[j] in the solution of the current node.
							double pVal = getNodeRel(p[j]);
							//Initializes the distance directly from the super source.
							distances[j] = zVal + pVal - uncertainVariables[j].getDeviation() * xVal;
							//Initializes the super source as predecessor.
							predecessors[j] = -1;
							
							//We have already computed the shortest path for all variables with lower deviation.
							//Checks whether a path via these variables is shortest. 
							for (int i = 0; i < j; i++) {
								//We only consider arcs between variables of different deviations, as these paths correspond exactly to facet defining inequalities.
								if (uncertainVariables[j].getDeviation() > uncertainVariables[i].getDeviation()) {
									double newDist = distances[i] + pVal - (uncertainVariables[j].getDeviation()-uncertainVariables[i].getDeviation()) * xVal;
									if (newDist < distances[j]) {
										distances[j] = newDist;
										predecessors[j] = i;
									}
								}
								else {
									break;
								}
							}
							//If the computed shortest path has negative length, we add it to our model.
							if (distances[j] < 0) {
								GRBLinExpr cutExpr = new GRBLinExpr();
								//Add z to the cut.
								cutExpr.addTerm(1, z);
								//Add the current variable and all variables on the path from the super source (except the first) with their corresponding p.
								//The coefficient is the difference of the current vairable's and the predecessor's deviation. 
								int index = j;
								while (predecessors[index] != -1) {
									cutExpr.addTerm(1, p[index]);
									cutExpr.addTerm(-(uncertainVariables[index].getDeviation()-uncertainVariables[predecessors[index]].getDeviation()), uncertainVariables[index].getModelVariable());
									index = predecessors[index];
								}
								//The first variable that is visited from the super source has to be treated specially, as the coefficient is the variable's deviation.
								cutExpr.addTerm(1, p[index]);
								cutExpr.addTerm(-uncertainVariables[index].getDeviation(), uncertainVariables[index].getModelVariable());
								addCut(cutExpr, GRB.GREATER_EQUAL, 0);
							}
						}
					}
				}
			} catch (GRBException e) {
				e.printStackTrace();
			}
		}
	}
}