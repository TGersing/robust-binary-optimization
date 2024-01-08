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
		model.setCallback(new CutCallback());
		
		model.update();
	}


	
	/**
	 * Implements the separation algorithm searching for violated constraints.
	 */
	private class CutCallback extends GRBCallback {
		//If we cut-off an incumbent solution then we can repair it by setting the appropriate value for robustnessVar.
		//The repaired incumbent cannot be added in the same call where the lazy constraint is added.
		//We store it until we reach the next node to add the incumbent.
		private double[] incumbentNominalValues;
		private double incumbentRobustnessValue;
		private boolean hasIncumbent = false;
		
		@Override
		protected void callback() {
			if (where == GRB.CB_MIP) {
				try {
					setPrimalBound(getDoubleInfo(GRB.CB_MIP_OBJBST));
					setDualBound(getDoubleInfo(GRB.CB_MIP_OBJBND));
				} catch (GRBException e) {}
			}

			try {
				//We separate cuts corresponding to scenarios while we are in the root node.
				//Also, every time we find a new incumbent, we check whether it violates a scenario.
				if ((where == GRB.CB_MIPNODE && getDoubleInfo(GRB.CB_MIPNODE_NODCNT) == 0 && getIntInfo(GRB.CB_MIPNODE_STATUS) == GRB.Status.OPTIMAL)
						|| where == GRB.CB_MIPSOL) {
					Variable[] uncertainVariables = robustProblem.getUncertainVariables();
					GRBVar[] uncertainModelVariables = robustProblem.getUncertainModelVariables();
					double Gamma = robustProblem.getGamma();

					//Queries the values of the uncertain variables and the rubustnessVar.
					double[] varValues;
					double robustnessVarValue;
					if (where == GRB.CB_MIPNODE) {
						varValues = getNodeRel(uncertainModelVariables);
						robustnessVarValue = getNodeRel(robustnessVar);
					}
					else {
						varValues = getSolution(uncertainModelVariables);
						robustnessVarValue = getSolution(robustnessVar);
					}
					
					//The worst case scenario includes the variables where d[i]*x[i] is as high as possible, with d being the deviation.
					//We store up to ceil(Gamma) indices with the highest d[i]*x[i], since not more variables will deviate.
					PriorityQueue<Integer> highestVarIndicesValuesTimesDeviations = new PriorityQueue<Integer>(
							(index1, index2) -> Double.compare(varValues[index1]*uncertainVariables[index1].getDeviation(), varValues[index2]*uncertainVariables[index2].getDeviation()));
					for (int i = 0; i < varValues.length; i++) {
						highestVarIndicesValuesTimesDeviations.add(i);
						if (highestVarIndicesValuesTimesDeviations.size() > Math.ceil(Gamma)) {
							highestVarIndicesValuesTimesDeviations.remove();
						}
					}
					
					//Computes the sum of the deviation for the worst case scenario
					double robustnessValue = 0;
					for (Integer index : highestVarIndicesValuesTimesDeviations) {
						robustnessValue += varValues[index]*uncertainVariables[index].getDeviation();
					}
					//If Gamma is not integer then the variable with the smallest value d[i]*x[i] differs only by a fraction.
					robustnessValue -= (Math.ceil(Gamma) - Gamma)*varValues[highestVarIndicesValuesTimesDeviations.peek()]*uncertainVariables[highestVarIndicesValuesTimesDeviations.peek()].getDeviation();
					
					//If the worst case value is greater than robustnessVarValue then we have found a violated cut. 
					if (robustnessValue > robustnessVarValue) {
						//Add the deviating variables of the worst case scenario to the cut.
						GRBLinExpr cutExpr = new GRBLinExpr();
						for (Integer index : highestVarIndicesValuesTimesDeviations) {
							cutExpr.addTerm(uncertainVariables[index].getDeviation(), uncertainModelVariables[index]);
						}
						//If Gamma is not integer then the variable with the smallest value d[i]*x[i] differs only by a fraction.
						cutExpr.addTerm(-(Math.ceil(Gamma) - Gamma)*uncertainVariables[highestVarIndicesValuesTimesDeviations.peek()].getDeviation(), uncertainModelVariables[highestVarIndicesValuesTimesDeviations.peek()]);
						
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
}