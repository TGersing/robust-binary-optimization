package alg;

import java.io.IOException;
import java.util.LinkedHashMap;

import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;
import util.Variable;

/**
 * This class implements the formulation RP3 from the paper "Strong formulations of robust mixed 0–1 programming" by Atamturk.
 * 
 * @author Timo Gersing
 */
public class AlgRP3Gurobi extends AbstractAlgorithm implements RobustAlgorithm{
	/**
	 * The robust problem to be solved.
	 */
	private RobustProblemGurobi robustProblem;		

	/**
	 * Constructor storing the given robust problem.
	 */
	public AlgRP3Gurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters) throws GRBException, IOException {
		super(algorithmParameters);
		this.robustProblem = new RobustProblemGurobi(problemPath, robustPath, algorithmParameters);
	}
	
	/**
	 * Executes the algorithm by reformulating the problem and solving it within the remaining time limit.
	 */
	@Override
	protected void executeAlgorithm() throws IOException, GRBException {
		String output = "######################################################\n"
				+ "##### Solving Problem via RP3"
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
	 * Reformulates the given model into RP3.
	 */
	private void reformulateModel() throws GRBException {
		//Queries the model and uncertain variables from the given robust problem.
		GRBModel model = robustProblem.getModel();		
		Variable[] uncertainVariables = robustProblem.getUncertainVariables();

		//Adds the variable z, p, ny.
		GRBVar z = model.addVar(0, Double.MAX_VALUE, robustProblem.getGamma(), GRB.CONTINUOUS, "z");
		GRBVar[] p = new GRBVar[uncertainVariables.length];
		for (int i = 0; i < uncertainVariables.length; i++) {
			p[i] = model.addVar(0, Double.MAX_VALUE, 1, GRB.CONTINUOUS, "p"+i);
		}
		GRBVar[] ny = new GRBVar[uncertainVariables.length+2];
		for (int i = 0; i < uncertainVariables.length+2; i++) {
			ny[i] = model.addVar(0, Double.MAX_VALUE, 0, GRB.CONTINUOUS, "ny"+i);
		}

		//Alters the objective function by adding Gamma*z and the sum over all variables p.
		GRBLinExpr objLinExpr = (GRBLinExpr) model.getObjective();
		for (int i = 0; i < uncertainVariables.length; i++) {
			objLinExpr.addTerm(1, p[i]);
		}
		objLinExpr.addTerm(robustProblem.getGamma(), z);
		model.setObjective(objLinExpr);
		
		//Add the constraints of the form (d[j]-d[i])x[j]+ny[j]-ny[i] <= p[j], with d being the deviation, for all 0<=i<j<=n with d[i]<d[j].
		//Note that we need to shift the indices, as ny ranges from 0 to n+1, while p and x range from 1 to n. (We have d[0]=0).
		//We start with i=0 (i=-1 with variable shift)
		for (int j = 0; j < uncertainVariables.length; j++) {
			GRBLinExpr expr = new GRBLinExpr();
			expr.addTerm(uncertainVariables[j].getDeviation()-0, uncertainVariables[j].getModelVariable());
			expr.addTerm(1, ny[j+1]);
			expr.addTerm(-1,ny[0]);
			model.addConstr(expr, GRB.LESS_EQUAL, p[j], "");
		}
		//Continue with i>=1 (i>=0 with variable shift)
		for (int i = 0; i < uncertainVariables.length; i++) {
			for (int j = i+1; j < uncertainVariables.length; j++) {
				//Only consider d[i]<d[j], as these inequalities are facet defining.
				if (uncertainVariables[j].getDeviation() > uncertainVariables[i].getDeviation()) {
					GRBLinExpr expr = new GRBLinExpr();
					expr.addTerm(uncertainVariables[j].getDeviation()-uncertainVariables[i].getDeviation(), uncertainVariables[j].getModelVariable());
					expr.addTerm(1, ny[j+1]);
					expr.addTerm(-1,ny[i+1]);
					model.addConstr(expr, GRB.LESS_EQUAL, p[j], "");
				}
			}
		}
		
		//Add the constraints of the form ny[n+1]-ny[i]<=z.
		for (int i = 0; i <= uncertainVariables.length; i++) {
			GRBLinExpr expr = new GRBLinExpr();
			expr.addTerm(1, ny[uncertainVariables.length+1]);
			expr.addTerm(-1, ny[i]);
			model.addConstr(expr, GRB.LESS_EQUAL, z, "");
		}
		
		//Add the constraint of the form ny[n+1]-ny[0]>=0.
		GRBLinExpr expr = new GRBLinExpr();
		expr.addTerm(1, ny[uncertainVariables.length+1]);
		expr.addTerm(-1, ny[0]);
		model.addConstr(expr, GRB.GREATER_EQUAL, 0, "");
		
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