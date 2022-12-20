package alg;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.PriorityQueue;

import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import util.Variable;

/**
 * This class implements the cutting plane approach that separates scenarios from the exponential uncertainty set.
 * 
 * @author Timo Gersing
 */
public class AlgCuttingPlanesGurobi extends AbstractAlgorithm implements RobustAlgorithm{
	/**
	 * The robust problem to be solved.
	 */
	private RobustProblemGurobi robustProblem;
	/**
	 * A variable modeling the sum of the deviations in the worst case scenario.
	 */
	private GRBVar robustnessVar;
	
	/**
	 * Constructor receiving problem paths as well as parameters and strategies.
	 */
	public AlgCuttingPlanesGurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters) throws GRBException, IOException {
		super(algorithmParameters);
		this.robustProblem = new RobustProblemGurobi(problemPath, robustPath, algorithmParameters);
	}
	
	/**
	 * Executes the algorithm by reformulating the problem and solving it within the remaining time limit.
	 */
	@Override
	protected void executeAlgorithm() throws IOException, GRBException {
		String output = "######################################################\n"
				+ "##### Solving Problem via Separating Scenarios"
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
	private void reformulateModel() throws GRBException {
		GRBModel model = robustProblem.getModel();		
		
		robustnessVar = model.addVar(0, Double.MAX_VALUE, 1, GRB.CONTINUOUS, "robustnessSum");
		
		//Alters the objective function by adding the robustnessVar.
		GRBLinExpr objLinExpr = (GRBLinExpr) model.getObjective();
		objLinExpr.addTerm(1, robustnessVar);
		model.setObjective(objLinExpr);
		
		//Adds the callback for separation of lazy constraints and valid inequalities.
		//PreCrush and lazyConstraints have to be set to 1 for the separation.
		model.set(GRB.IntParam.PreCrush, 1);
		model.set(GRB.IntParam.LazyConstraints, 1);
		model.setCallback(new Callback());
		
		model.update();
	}


	
	/**
	 * Implements the separation algorithm searching for violated constraints.
	 */
	private class Callback extends GRBCallback {
		//If we cut-off an incumbent solution then we can repair it by setting the appropriate value for robustnessVar.
		//The repaired incumbent cannot be added in the same call where the lazy constraint is added.
		//We store it until we reach the next node to add the incumbent.
		private double[] incumbentNominalValues;
		private double incumbentRobustnessValue;
		private boolean hasIncumbent = false;
		
		@Override
		protected void callback() {
			try {
				//We separate cuts corresponding to scenarios while we are in the root node.
				//Also, every time we find a new incumbent, we also check whether it violates a scenario.
				if ((where == GRB.CB_MIPNODE && getDoubleInfo(GRB.CB_MIPNODE_NODCNT) == 0 && getIntInfo(GRB.CB_MIPNODE_STATUS) == GRB.Status.OPTIMAL)
						|| where == GRB.CB_MIPSOL) {
					Variable[] uncertainVariables = robustProblem.getUncertainVariables();
					GRBVar[] uncertainModelVariables = robustProblem.getUncertainModelVariables();
					double Gamma = robustProblem.getGamma();

					//Queries the values of the uncertain variables and the rubustnessVar.
					double[] varValues;
					Double robustnessVarValue;
					if (where == GRB.CB_MIPNODE) {
						varValues = getNodeRel(uncertainModelVariables);
						robustnessVarValue = getNodeRel(robustnessVar);
					}
					else {
						varValues = getSolution(uncertainModelVariables);
						robustnessVarValue = getSolution(robustnessVar);
					}
					
					//The worst case scenario includes the variables where d[i]*x[i] is as high as possible, with d being the deviation.
					//We store all variables with their corresponding value d[i]*x[i] as a tuple in a priority queue.
					//The tuples are sorted with respect to d[i]*x[i].
					//The queue contains at most ceil(Gamma) variables, since not more variables will deviate.
					PriorityQueue<IndexValueTuple> highestValuesTimesDeviations = new PriorityQueue<IndexValueTuple>();
					for (int i = 0; i < varValues.length; i++) {
						double valueTimesDeviation = varValues[i]*uncertainVariables[i].getDeviation();
						if (highestValuesTimesDeviations.size() < Math.ceil(Gamma)) {
							highestValuesTimesDeviations.add(new IndexValueTuple(i, valueTimesDeviation));
						}
						else if (valueTimesDeviation > highestValuesTimesDeviations.peek().getValue()) {
							highestValuesTimesDeviations.remove();
							highestValuesTimesDeviations.add(new IndexValueTuple(i, valueTimesDeviation));
						}
					}
					
					//Computes the sum of the deviation for the worst case scenario
					double robustnessValue = 0;
					for (IndexValueTuple tuple : highestValuesTimesDeviations) {
						robustnessValue += tuple.getValue();
					}
					//If Gamma is not integer then the variable with the smallest value d[i]*x[i] differs only by a fraction.
					robustnessValue -= (Math.ceil(Gamma) - Gamma)*highestValuesTimesDeviations.peek().getValue();
					
					//If the worst case value is greater than robustnessVarValue then we have found a violated cut. 
					if (robustnessValue > robustnessVarValue) {
						//Add the deviating variables of the worst case scenario to the cut.
						GRBLinExpr cutExpr = new GRBLinExpr();
						for (IndexValueTuple tuple : highestValuesTimesDeviations) {
							cutExpr.addTerm(uncertainVariables[tuple.getIndex()].getDeviation(), uncertainModelVariables[tuple.getIndex()]);
						}
						//If Gamma is not integer then the variable with the smallest value d[i]*x[i] differs only by a fraction.
						cutExpr.addTerm(-(Math.ceil(Gamma) - Gamma)*uncertainVariables[highestValuesTimesDeviations.peek().getIndex()].getDeviation(), uncertainModelVariables[highestValuesTimesDeviations.peek().getIndex()]);
						
						//The sum has to be smaller than the value of robustnessVar.
						cutExpr.addTerm(-1, robustnessVar);
						if (where == GRB.CB_MIPNODE) {
							addCut(cutExpr, GRB.LESS_EQUAL, 0);
						}
						else {
							addLazy(cutExpr, GRB.LESS_EQUAL, 0);
							//If we just cut-off an incumbent solution then we can repair it by setting the appropriate value for robustnessVar.
							hasIncumbent = true;
							incumbentNominalValues = getSolution(robustProblem.getNominalModelVariables());
							incumbentRobustnessValue = robustnessValue;
						}
					}
				}
				//Adds the repaired incumbent solution at the next possible MIP node.
				if (hasIncumbent && where == GRB.CB_MIPNODE) {
					setSolution(robustnessVar, incumbentRobustnessValue);
					setSolution(robustProblem.getNominalModelVariables(), incumbentNominalValues);
					hasIncumbent = false;
				}
			} catch (GRBException e1) {
				e1.printStackTrace();
			}
		}
	}
	
	/**
	 * Auxiliary class for storing indices of variables together with a value in a tuple.
	 * To allow for sorting, the tuple is comparable with respect to the value.  
	 */
	private class IndexValueTuple implements Comparable<IndexValueTuple>{
		int index;
		double value;
		
		public IndexValueTuple(int index, double value) {
			this.index = index;
			this.value = value;
		}

		public int getIndex() {
			return index;
		}
		public double getValue() {
			return value;
		}

		@Override
		public int compareTo(IndexValueTuple tuple) {
			return Double.compare(this.value, tuple.value);
		}
	}
}