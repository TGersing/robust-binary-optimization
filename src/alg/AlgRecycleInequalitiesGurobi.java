package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import alg.AlgRecycleInequalitiesGurobi.RecyclingStrategies.Relaxation;
import util.ConflictGraph;
import util.Variable;
import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBConstr;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

/**
 * This class implements the standard reformulation by Bertsimas and Sim together with the recycling of valid inequalities,
 * as described in "Recycling Inequalities for Robust Combinatorial Optimization with Budget Uncertainty" by Büsing, Gersing and Koster.
 * 
 * @author Timo Gersing
 */
public class AlgRecycleInequalitiesGurobi extends AbstractAlgorithm implements RobustAlgorithm{
	/**
	 * The robust problem to be solved.
	 */
	private RobustProblemGurobi robustProblem;

	/**
	 * Specifies strategies.
	 */
	private RecyclingStrategies recyclingStrategies;
	
	/**
	 * The p variables from the Bertsimas Sim reformulation.
	 */
	private GRBVar[] p;
	
	/**
	 * The z variable from the Bertsimas Sim reformulation.
	 */
	private GRBVar z;
	
	/**
	 * List of nominal constraints used for recycling.
	 */
	private NominalConstr[] nominalConstrs;
	
	/**
	 * Maps GRBVars to the corresponding Variable object.
	 */
	private Map<GRBVar, Variable> modelVarToVariable;
		
	/**
	 * Class for executing the separation.
	 */
	private Separator separator;
	
	/**
	 * Quotient for scaling deviations to improve numerical stability of recycled inequalities.
	 * Substitutes p and z by p*scalingQuotient and z*scalingQuotient (in objective function and constraints).
	 * Afterwards divides the robustness constraints and recycled inequalities by scalingQuotient.
	 */
	private double scalingQuotient;
	
	/**
	 * Constructor receiving problem paths as well as parameters and strategies.
	 */
	public AlgRecycleInequalitiesGurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, RecyclingStrategies recyclingStrategies) throws GRBException, IOException {
		super(algorithmParameters);
		this.robustProblem = new RobustProblemGurobi(problemPath, robustPath, algorithmParameters);
		this.recyclingStrategies = recyclingStrategies;
	}

	/**
	 * Executes the algorithm by reformulating the problem and solving it within the remaining time limit.
	 * If we consider the linear relaxation and separate cuts then we do this iteratively until the time limit is reached or no violated cuts are found.
	 */
	@Override
	protected void executeAlgorithm() throws IOException, GRBException {
		reformulateModel();
		
		GRBModel model = robustProblem.getModel();
		
		//Relaxes variables if the option is chosen
		if (recyclingStrategies.relaxation == Relaxation.SOLVE_RELAXED) {
			for (GRBVar grbVar : robustProblem.getNominalModelVariables()) {
				grbVar.set(GRB.CharAttr.VType, GRB.CONTINUOUS);
			}
		}
		//If we solve the integer problem and separate cuts then we add the callback
		//In this case, we set the precrush value to 1 such that cuts can always be used by Gurobi
		else {
			model.setCallback(new RecyclingCallback());
			if (separator != null) {
				model.set(GRB.IntParam.PreCrush, 1);
			}
		}
		
		//Disables Gurobi's cuts if the option is chosen
		if (recyclingStrategies.gurobiCutStrategy == RecyclingStrategies.GurobiCutStrategy.GCUTS_DISABLE) {
			model.set(GRB.IntParam.Cuts, 0);
		}
		
		//Solves the problem in the remaining time
		model.update();
		robustProblem.solve(getRemainingTime());
		primalBound = robustProblem.getPrimalBound();
		dualBound = robustProblem.getDualBound();
		
		//Initializes solution
		solution = new LinkedHashMap<Variable, Double>(robustProblem.getNominalVariables().length);
		
		//If we consider the linear relaxation and separate cuts then we proceed separating and reoptimizing until no more violated cuts are found
		while (recyclingStrategies.relaxation == Relaxation.SOLVE_RELAXED && separator != null) {
			if (model.get(GRB.IntAttr.Status) != GRB.OPTIMAL) {
				writeOutput("######################################################"
						+ "\n##### LP not solved to optimality. Terminate separation."
						+ "\n######################################################");
				break;
			}
			
			primalBound = robustProblem.getPrimalBound();
			dualBound = robustProblem.getDualBound();
			
			//Queries the values of the current relaxed solution.
			double[] xValues = model.get(GRB.DoubleAttr.X, robustProblem.getNominalModelVariables());
			double[] pValues = model.get(GRB.DoubleAttr.X, p);
			Double zValue = z.get(GRB.DoubleAttr.X);
			
			//Stores solution
			for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
				solution.put(robustProblem.getNominalVariables()[i], xValues[i]);
			}
			
			//Separates violated cuts
			Collection<RobustInequality> violatedRobustInequalities = separator.getViolatedRobustInequalities(xValues, pValues, zValue, algorithmParameters.getFeasibilityTolerance());
			//Adds cuts to the model
			for (RobustInequality robustInequality : violatedRobustInequalities) {
				double violation = robustInequality.computeViolation(xValues, pValues, zValue);
				robustInequality.writeSeparationMessage(violation);
				model.addConstr(robustInequality.getLHSexpr(), GRB.GREATER_EQUAL, robustInequality.getConstant(), "violatedRecycledIneq");
			}
			if (violatedRobustInequalities.isEmpty()) {
				writeOutput("######################################################"
						+ "\n##### Found no violated cut"
						+ "\n######################################################");
				break;
			}
			
			model.update();
			if (getRemainingTime().isPresent() && getRemainingTime().get() > 0) {
				robustProblem.solve(getRemainingTime());
			}
			else {
				writeOutput("######################################################"
						+ "\n##### Reached time limit. Terminate separation."
						+ "\n######################################################");
			}
		}
		
		//Stores solution
		if (robustProblem.getNominalVariablesSolutionValues() != null) {
			for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
				solution.put(robustProblem.getNominalVariables()[i], robustProblem.getNominalVariablesSolutionValues()[i]);
			}

		}
	}

	/**
	 * Reformulates the given model to the Bertsimas Sim reformulation.
	 * Also initializes the separation of recycled inequalities or adds recycled inequalities directly to the model.
	 */
	private void reformulateModel() throws GRBException, IOException {
		//Queries the model and uncertain variables from the given robust problem.
		GRBModel model = robustProblem.getModel();
		Variable[] uncertainVariables = robustProblem.getUncertainVariables();
		
		//Stores the nominal constraints for recycling before we add the robustness constraints.
		GRBConstr[] grbConstrs = model.getConstrs();
		nominalConstrs = new NominalConstr[grbConstrs.length];
		for (int j = 0; j < grbConstrs.length; j++) {
			nominalConstrs[j] = new NominalConstr(grbConstrs[j]);
		}
		
		//Maps GRBVars to the corresponding Variable object.
		modelVarToVariable = new HashMap<GRBVar, Variable>();
		for (Variable variable : robustProblem.getNominalVariables()) {
			modelVarToVariable.put(variable.getModelVariable(), variable);
		}
		
		//Computes scalingQuotient
		double minDev = Double.POSITIVE_INFINITY;
		double maxDev = Double.NEGATIVE_INFINITY;
		for (Variable var : robustProblem.getUncertainVariables()) {
			minDev = Math.min(minDev, var.getDeviation());
			maxDev = Math.max(maxDev, var.getDeviation());
		}
		scalingQuotient = Math.sqrt(minDev*maxDev);
		
		//Adds the variables z and p from the reformulation, scaled with the scalingQuotient.
		z = model.addVar(0, Double.MAX_VALUE, scalingQuotient*robustProblem.getGamma(), GRB.CONTINUOUS, "z");
		p = new GRBVar[uncertainVariables.length];
		for (int i = 0; i < p.length; i++) {
			p[i] = model.addVar(0, Double.MAX_VALUE, scalingQuotient, GRB.CONTINUOUS, "p"+i);
		}
		
		//Alters the objective function by adding z and p.
		GRBLinExpr objLinExpr = (GRBLinExpr) model.getObjective();
		for (int i = 0; i < p.length; i++) {
			objLinExpr.addTerm(scalingQuotient, p[i]);
		}
		objLinExpr.addTerm(scalingQuotient*robustProblem.getGamma(), z);
		model.setObjective(objLinExpr);
				
		//Adds the constraint of the form z+p>=d/scalingQuotient*x, with d being the deviation.
		for (int i = uncertainVariables.length-1; i >= 0; i--) {
			Variable var = uncertainVariables[i];
			GRBLinExpr robustnessExpr = new GRBLinExpr();
			robustnessExpr.addTerm(1, p[i]);
			robustnessExpr.addTerm(1, z);
			robustnessExpr.addTerm(-var.getDeviation()/scalingQuotient, var.getModelVariable());
			model.addConstr(robustnessExpr, GRB.GREATER_EQUAL, 0, "");
		}
		model.update();
		
		List<RecyclableInequality> recyclableInequalities = null;
		//Sets separators if the corresponding options are chosen 
		if (recyclingStrategies.separationStrategy == RecyclingStrategies.SeparationRecyclingStrategy.SEPRECYCLE_CONSTRAINTS) {
			recyclableInequalities = computeRecyclableInequalities();
			separator = new RecyclableSeparator(recyclableInequalities);
		}
		else if (recyclingStrategies.separationStrategy == RecyclingStrategies.SeparationRecyclingStrategy.SEPRECYCLE_LP) {
			recyclableInequalities = computeRecyclableInequalities();
			separator = new LPSeparator(recyclableInequalities);
		}
		else if (recyclingStrategies.separationStrategy == RecyclingStrategies.SeparationRecyclingStrategy.SEPRECYCLE_LIFT) {
			separator = new LiftingSeparator();
		}
		else if (recyclingStrategies.separationStrategy == RecyclingStrategies.SeparationRecyclingStrategy.SEPRECYCLE_CLIQUES) {
			separator = new CliqueSeparator();
		}
		
		//Adds recycled inequalities directly to the model if the corresponding strategy is chosen.
		if (recyclingStrategies.directRecyclingStrategy == RecyclingStrategies.DirectRecyclingStrategy.DIRRECYCLE_CONSTRAINTS) {
			if (recyclableInequalities == null) {
				recyclableInequalities = computeRecyclableInequalities();
			}
			for (RecyclableInequality recylableInequality: recyclableInequalities) {
				RobustInequality recycledInequality = recylableInequality.getRecycledInequality();
				model.addConstr(recycledInequality.getLHSexpr(), GRB.GREATER_EQUAL, recycledInequality.getConstant(), "");
			}
		}
		
		model.update();
	}
	
	/**
	 * Callback is called at the root node of the branch and bound.
	 * Passes the current solution x,p,z to a separator and adds robust inequalities if they are violated.
	 */
	private class RecyclingCallback extends GRBCallback {
		
		@Override
		protected void callback() {
			//Updates the primal dual integral
			if (where == GRB.CB_MIP) {
				try {
					primalDualIntegral.update(getDoubleInfo(GRB.CB_MIP_OBJBST), getDoubleInfo(GRB.CB_MIP_OBJBND), false);
				} catch (GRBException e) {}
			}
			
			try {
				//Checks whether we are in the root node and are allowed to add user cuts
				if (separator != null && where == GRB.CB_MIPNODE && getDoubleInfo(GRB.CB_MIPNODE_NODCNT) == 0 && getIntInfo(GRB.CB_MIPNODE_STATUS) == GRB.Status.OPTIMAL) {
					//Queries the values of the current relaxed solution.
					double[] xValues = getNodeRel(robustProblem.getNominalModelVariables());
					double[] pValues = getNodeRel(p);
					Double zValue = getNodeRel(z);
					
					//Passes the current solution to the separator 
					Collection<RobustInequality> robustInequalities = separator.getViolatedRobustInequalities(xValues, pValues, zValue, algorithmParameters.getFeasibilityTolerance());
					
					//Adds recycled inequalities if violated
					for (RobustInequality robustInequality : robustInequalities) {
						double violation = robustInequality.computeViolation(xValues, pValues, zValue);
						robustInequality.writeSeparationMessage(violation);
						addCut(robustInequality.getLHSexpr(), GRB.GREATER_EQUAL, robustInequality.getConstant());
					}
				}
			} catch (GRBException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Interface for separators, which obtain a solution x,p,z and return robust inequalities.
	 */
	private interface Separator {
		/**
		 * Returns only violated robust inequalities
		 */
		default Collection<RobustInequality> getViolatedRobustInequalities(double[] xValues, double[] pValues, double zValue, double violationThreshold) throws GRBException, IOException{
			Collection<RobustInequality> robustInequalities = this.computeRobustInequalities(xValues, pValues, zValue);
			robustInequalities.removeIf(recycledInequality -> recycledInequality.computeViolation(xValues, pValues, zValue) <= violationThreshold);
			return robustInequalities;
		}
		
		/**
		 * Computes robust inequalities according to the separation strategy.
		 * @throws IOException 
		 */
		abstract Collection<RobustInequality> computeRobustInequalities(double[] xValues, double[] pValues, double zValue) throws GRBException, IOException;
	}
	
	/**
	 * Preprocesses recyclable inequalities and provides their optimal sub-inequalities during separation.
	 */
	private class RecyclableSeparator implements Separator {
		/**
		 * A list of recyclable inequalities
		 */
		private List<RecyclableInequality> recyclableInequalities;
		
		public RecyclableSeparator(List<RecyclableInequality> recyclableInequalities) throws GRBException {
			this.recyclableInequalities = recyclableInequalities;
		}

		@Override
		public List<RobustInequality> computeRobustInequalities(double[] xValues, double[] pValues, double zValue) {
			List<RobustInequality> robustInequalities = new ArrayList<RobustInequality>();
			for (RecyclableInequality recycledInequality : recyclableInequalities) {
				robustInequalities.add(recycledInequality.getRecycledInequality(xValues, pValues));
			}
			return robustInequalities;
		}
	}
	
	/**
	 * Computes a conflict graph for the nominal problem and uses it to separate recycled clique inequalities.
	 */
	private class CliqueSeparator implements Separator{
		/**
		 * Maps variables to the corresponding neighborhood in the conflict graph.
		 */
		private Map<Variable, Set<Variable>> uncertainVariableToConflictingNeighborhood;
		
		/**
		 * Constructor computing the neighborhoods
		 */
		public CliqueSeparator() throws GRBException, IOException {
			ConflictGraph conflictGraph = new ConflictGraph(robustProblem.getModel(), robustProblem.getUncertainModelVariables(), algorithmParameters);
			uncertainVariableToConflictingNeighborhood = new HashMap<Variable, Set<Variable>>();
			for (GRBVar uncertainModelVariable : robustProblem.getUncertainModelVariables()) {
				Set<Variable> neighborhood = new HashSet<Variable>();
				for (GRBVar neighborVar : conflictGraph.getNeighbors(uncertainModelVariable)) {
					neighborhood.add(modelVarToVariable.get(neighborVar));
				}
				uncertainVariableToConflictingNeighborhood.put(modelVarToVariable.get(uncertainModelVariable), neighborhood);
			}
		}
				
		@Override
		public Set<RobustInequality> computeRobustInequalities(double[] xValues, double[] pValues, double zValue) {
			//Assigns each variable its impact on the violation in the recycled inequality.
			Function<Variable, Double> violationImpactFunc = variable -> xValues[variable.getNominalIndex()]*variable.getDeviation()/scalingQuotient - pValues[variable.getUncertainIndex()];
			//The clique inequalities to be recycled
			Set<RobustInequality> cliqueInequalities = new HashSet<RobustInequality>();
			//Set of already considered cliques to avoid duplicate recycling
			Set<Set<Variable>> cliques = new HashSet<Set<Variable>>();
			//Extends each variable with a positive violation impact to a clique
			for (Variable uncertainVariable : robustProblem.uncertainVariables) {
				if (violationImpactFunc.apply(uncertainVariable) <= 0) {
					continue;
				}
				//Initializes and extends clique
				Set<Variable> clique = new HashSet<Variable>(uncertainVariableToConflictingNeighborhood.get(uncertainVariable).size());
				clique.add(uncertainVariable);
				extendClique(clique, violationImpactFunc);
				//Adds clique inequality if not already considered
				if (!cliques.contains(clique)) {
					RobustInequality recycledCliqueInequality =  new RobustInequality();
					recycledCliqueInequality.zCoeff = 1;
					for (Variable variable : clique) {
						recycledCliqueInequality.addVariable(variable, variable.getDeviation()/scalingQuotient, 1);
					}
					cliqueInequalities.add(recycledCliqueInequality);
					cliques.add(clique);
				}
			}
			return cliqueInequalities;
		}
		
		/**
		 * Extends a clique greedily with respect to an impact function.
		 */
		private void extendClique (Set<Variable> clique, Function<Variable, Double> impactFunction) {
			Iterator<Variable> iterator = clique.iterator();
			//Potential clique members that are in the neighborhood of all current clique members.
			Set<Variable> potentialMembers = new HashSet<Variable>(uncertainVariableToConflictingNeighborhood.get(iterator.next()));
			while (iterator.hasNext()) {
				potentialMembers.retainAll(uncertainVariableToConflictingNeighborhood.get(iterator.next()));
			}
			
			//Removes all potential clique members with a negative impact 
			iterator = potentialMembers.iterator();
			while (iterator.hasNext()) {
				if (impactFunction.apply(iterator.next()) < 0) {
					iterator.remove();
				}
			}
			
			//Adds the potential member with the highest impact and removes all potential members that are not within its neighborhood.
			if (!potentialMembers.isEmpty()) {
				Comparator<Variable> comparator = (var1, var2) -> Double.compare(impactFunction.apply(var1), impactFunction.apply(var2));
				while (!potentialMembers.isEmpty()) {
					Variable newMember = Collections.max(potentialMembers, comparator);
					clique.add(newMember);
					potentialMembers.retainAll(uncertainVariableToConflictingNeighborhood.get(newMember));
				}
			}
		}
	}

	/**
	 * First tries to separate recyclable constraints.
	 * If this isn't successful, solves an LP combining nominal constraints for recycling.
	 */
	private class LPSeparator implements Separator {
		/**
		 * Separator for recyclable constraints.
		 */
		private RecyclableSeparator recyclableSeparator;

		/**
		 * The LP model
		 */
		private GRBModel sepModel;
		
		/**
		 * Variables defining the coefficient with which a constraint is used to combine the valid inequality
		 */
		private GRBVar[] coeffConstr;
		
		/**
		 * Variables defining the coefficient with which a x>=LB constraint is used to combine the valid inequality
		 */
		private GRBVar[] coeffLB;
		
		/**
		 * Variables defining the coefficient with which a x<=UB constraint is used to combine the valid inequality
		 */
		private GRBVar[] coeffUB;
		
		/**
		 * Left-hand side of the combined inequality
		 */
		private GRBVar[] pi;
		
		/**
		 * Right-hand side of the combined inequality
		 */
		private GRBVar piZero;
		
		/**
		 * Sets of constraints that can be combined to a recyclable inequality
		 */
		private List<Set<NominalConstr>> combinableConstraintsSets;
		
		public LPSeparator(List<RecyclableInequality> recyclableInequalities) throws GRBException {
			recyclableSeparator = new RecyclableSeparator(recyclableInequalities);
			//Partitions the nominal constraints into sets that can be combined
			computeCombinableConstraints();			
			//Initializes separation LP
			sepModel = new GRBModel(new GRBEnv());
			//Initializes pi variables, which have to be >= 0 so that we can recycle the resulting valid inequality
			pi = new GRBVar[robustProblem.getUncertainVariables().length];
			for (int i = 0; i < robustProblem.getUncertainVariables().length; i++) {
				pi[i] = sepModel.addVar(0, Double.POSITIVE_INFINITY, 0, GRB.CONTINUOUS, "");
			}
			//Initializes the pizero variable, which we fix to 1 in order to normalize the valid inequality
			piZero = sepModel.addVar(1, 1, 0, GRB.CONTINUOUS, "");
			
			//Initializes coeff variables of the constraints
			coeffConstr = new GRBVar[nominalConstrs.length];
			for (int i = 0; i < nominalConstrs.length; i++) {
				//ax<=b constraints can only be added with a positive coefficient
				if (nominalConstrs[i].getSense() == GRB.LESS_EQUAL) {
					coeffConstr[i] = sepModel.addVar(0, Double.POSITIVE_INFINITY, 0, GRB.CONTINUOUS, "");
				}
				//ax>=b constraints can only be added with a negative coefficient (equivalent to adding -ax<=-b with positive coefficient)
				else if (nominalConstrs[i].getSense() == GRB.GREATER_EQUAL) {
					coeffConstr[i] = sepModel.addVar(Double.NEGATIVE_INFINITY, 0, 0, GRB.CONTINUOUS, "");
				}
				//ax=b constraints can be added with positive or negative coefficients
				else {
					coeffConstr[i] = sepModel.addVar(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 0, GRB.CONTINUOUS, "");
				}
			}
			
			//If the lower bound exists (LB > -10^100):
			//Initializes coeff variables of x>=LB, which can be added with negative coeffs
			coeffLB = new GRBVar[robustProblem.getNominalVariables().length];
			for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
				double LB = robustProblem.getNominalModelVariables()[i].get(GRB.DoubleAttr.LB);
				if (Math.abs(LB) < Math.pow(10, 100)) {
					coeffLB[i] = sepModel.addVar(Double.NEGATIVE_INFINITY, 0, 0, GRB.CONTINUOUS, "");
				}
			}
			
			//If the upper bound exists (UB < 10^100):
			//Initializes coeff variables of x<=UB, which can be added with positive coeffs
			coeffUB = new GRBVar[robustProblem.getNominalVariables().length];
			for (int i = 0; i < robustProblem.getNominalModelVariables().length; i++) {
				double UB = robustProblem.getNominalModelVariables()[i].get(GRB.DoubleAttr.UB);
				if (Math.abs(UB) < Math.pow(10, 100)) {
					coeffUB[i] = sepModel.addVar(0, Double.POSITIVE_INFINITY, 0, GRB.CONTINUOUS, "");
				}
			}
			
			//The weighted sum of the right hand sides has to be equal to pizero
			GRBLinExpr piZeroConstrExpr = new GRBLinExpr();
			for (int i = 0; i < nominalConstrs.length; i++) {
				piZeroConstrExpr.addTerm(nominalConstrs[i].getRHS(), coeffConstr[i]);
			}
			for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
				double UB = robustProblem.getNominalModelVariables()[i].get(GRB.DoubleAttr.UB);
				if (coeffUB[i] != null) {
					piZeroConstrExpr.addTerm(UB, coeffUB[i]);
				}
				double LB = robustProblem.getNominalModelVariables()[i].get(GRB.DoubleAttr.LB);
				if (coeffLB[i] != null) {
					piZeroConstrExpr.addTerm(LB, coeffLB[i]);
				}
			}
			sepModel.addConstr(piZeroConstrExpr, GRB.EQUAL, piZero, "");
			
			//Constraints on the weighted sum of the left hand sides
			Map<GRBVar, GRBLinExpr> piConstrExprMap = new HashMap<GRBVar, GRBLinExpr>();
			for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
				GRBLinExpr grbLinExpr = new GRBLinExpr();
				if (coeffUB[i] != null) {
					grbLinExpr.addTerm(1, coeffUB[i]);
				}
				if (coeffLB[i] != null) {
					grbLinExpr.addTerm(1, coeffLB[i]);
				}
				piConstrExprMap.put(robustProblem.getNominalModelVariables()[i], grbLinExpr);
			}
			for (int j = 0; j < nominalConstrs.length; j++) {
				for (int i = 0; i < nominalConstrs[j].size(); i++) {
					piConstrExprMap.get(nominalConstrs[j].getNominalVariable(i).getModelVariable()).addTerm(nominalConstrs[j].getCoefficient(i), coeffConstr[j]);
				}
			}
			//The weighted sum of the left hand sides has to match the pi variables for all uncertain variables
			for (int i = 0; i < robustProblem.getUncertainVariables().length; i++) {
				sepModel.addConstr(piConstrExprMap.get(robustProblem.getUncertainModelVariables()[i]), GRB.EQUAL, pi[i], "");
			}
			//For all certain variables, the weighted sum has to be >=0 in order to obtain a valid inequality that can be recycled
			for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
				if (robustProblem.getNominalVariables()[i].getDeviation() == 0) {
					sepModel.addConstr(piConstrExprMap.get(robustProblem.getNominalModelVariables()[i]), GRB.GREATER_EQUAL, 0, "");
				}
			}
		}
		
		/**
		 * First tries to separate recycled inequalities from already recyclable ones.
		 * If this is not successful, solves the separation LP.
		 */
		@Override
		public Collection<RobustInequality> computeRobustInequalities(double[] xValues, double[] pValues, double zValue) throws GRBException, IOException {
			//Tries to separate using recyclable inequalities
			Collection<RobustInequality> recycledInequalities = recyclableSeparator.getViolatedRobustInequalities(xValues, pValues, zValue, algorithmParameters.getFeasibilityTolerance());
			//Only proceeds if none are violated
			if (!recycledInequalities.isEmpty()) {
				return recycledInequalities;
			}
			
			//Sets the objective function to be the violation of the recycled constraint
			GRBLinExpr objExpr = new GRBLinExpr();
			for (Variable variable : robustProblem.getUncertainVariables()) {
				objExpr.addTerm(xValues[variable.getNominalIndex()]*variable.getDeviation()/scalingQuotient-pValues[variable.getUncertainIndex()], pi[variable.getUncertainIndex()]);
			}
			objExpr.addTerm(-zValue, piZero);
			sepModel.setObjective(objExpr, GRB.MAXIMIZE);
			
			//Solves the separation problem for each set in the partition of combinable constraints
			for (Set<NominalConstr> combinableConstraints : combinableConstraintsSets) {
				//Sets variables of all constraints outside the set to zero
				List<GRBConstr> tempCons = new ArrayList<GRBConstr>();
				if (combinableConstraints.size() != nominalConstrs.length) {
					for (int i = 0; i < nominalConstrs.length; i++) {
						if (!combinableConstraints.contains(nominalConstrs[i])) {
							tempCons.add(sepModel.addConstr(coeffConstr[i], GRB.EQUAL, 0, ""));
						}
					}
				}
				
				//Solves the separation LP in the remaining time limit 
				if (getRemainingTime().isPresent()) {
					if (getRemainingTime().get() <= 0) {
						break;
					}
					sepModel.set(GRB.DoubleParam.TimeLimit, getRemainingTime().get());
				}
				sepModel.update();
				sepModel.optimize();
				
				//Constructs and stores the cut if the violation is positive
				if (sepModel.get(GRB.IntAttr.Status) == GRB.OPTIMAL && sepModel.get(GRB.DoubleAttr.ObjVal) > algorithmParameters.getFeasibilityTolerance()) {
					double[] piValues = sepModel.get(GRB.DoubleAttr.X, pi);
					double piZeroValue = piZero.get(GRB.DoubleAttr.X);
					List<Variable> variables = new ArrayList<Variable>();
					List<Double> coeffs = new ArrayList<Double>();
					for (int i = 0; i < piValues.length; i++) {
						if (piValues[i] > algorithmParameters.getFeasibilityTolerance()) {
							coeffs.add(piValues[i]);
							variables.add(robustProblem.getUncertainVariables()[i]);
						}
					}
					recycledInequalities.add(new RecyclableInequality(variables, coeffs, piZeroValue).getRecycledInequality());
				}
				
				//Removes temporary constraints fixing variables
				for (GRBConstr tempCon : tempCons) {
					sepModel.remove(tempCon);
				}
				sepModel.update();
			}
			return recycledInequalities;
		}
		
		/**
		 * Partitions the set of constraints into subsets that are reasonable to be combined into recyclable inequalities.
		 * The partitioning corresponds to the connected components in a graph, in which there exists an edge between two constraint-nodes
		 * if the constraints share at least one variable with coefficients of opposite signs.
		 */
		@SuppressWarnings("unchecked")
		private void computeCombinableConstraints() throws GRBException {
			combinableConstraintsSets = new ArrayList<Set<NominalConstr>>();
						
			//We compute a directed forest whose roots are representatives of the connected constraints
			//At the start, each node is its own component and thus representative
			//We also store the rank of each tree to append shorter trees to the root of longer ones
			int[] componentTreePredecessors = new int[nominalConstrs.length];
			int[] rank = new int[nominalConstrs.length];
			for (int j = 0; j < componentTreePredecessors.length; j++) {
				componentTreePredecessors[j] = j;
				rank[j] = 1;
			}
			
			//Function computing the representative and compressing the path
			Function<Integer, Integer> findRepresentativeAndCompress = constraintIndex -> {
				//Finds representative
				int representative = componentTreePredecessors[constraintIndex];
				while (componentTreePredecessors[representative] != representative) {
					representative = componentTreePredecessors[representative];
				}
				//Compresses path
				int currentIndex = constraintIndex;
				int oldPredecessor;
				while ((oldPredecessor = componentTreePredecessors[currentIndex]) != representative) {
					componentTreePredecessors[currentIndex] = representative;
					currentIndex = oldPredecessor;
				}
				return representative;
			};
			
			//Stores for each variable the indices of the constraints in which they occur positively or negatively
			List<Integer>[] constraintsSupportingVarPositive = new ArrayList[robustProblem.getNominalVariables().length];
			List<Integer>[] constraintsSupportingVarNegative = new ArrayList[robustProblem.getNominalVariables().length];
			for (int i = 0; i < constraintsSupportingVarPositive.length; i++) {
				constraintsSupportingVarPositive[i] = new ArrayList<Integer>();
				constraintsSupportingVarNegative[i] = new ArrayList<Integer>();
			}
			//Iterates over all constraints and adds them to the above lists
			for (int j = 0; j < nominalConstrs.length; j++) {
				NominalConstr nominalConstr = nominalConstrs[j];
				List<Integer> signs = new ArrayList<Integer>(2);
				if (nominalConstr.getSense() == GRB.LESS_EQUAL || nominalConstr.getSense() == GRB.EQUAL) {
					signs.add(1);
				}
				if (nominalConstr.getSense() == GRB.GREATER_EQUAL || nominalConstr.getSense() == GRB.EQUAL) {
					signs.add(-1);
				}
				for (int i = 0; i < nominalConstr.size(); i++) {
					int varIndex = nominalConstr.getNominalVariable(i).getNominalIndex();
					for (Integer sign : signs) {
						if (sign*nominalConstr.getCoefficient(i) > 0) {
							constraintsSupportingVarPositive[varIndex].add(j);
						}
						else if (sign*nominalConstr.getCoefficient(i) < 0) {
							constraintsSupportingVarNegative[varIndex].add(j);
						}
					}
				}
			}
			
			//Counts the number of components such that we can terminate if the number is 1
			int numberComponents = nominalConstrs.length;
			
			//Iterates over all variables to connect constraints.
			//Two constraints are connected if there exists a variable that is contained in both
			//and the variable occurs positively and negatively in some constraints.
			for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
				if (numberComponents == 1) {
					break;
				}
				//If the variable is contained positively and negatively in some constraints
				if (!constraintsSupportingVarPositive[i].isEmpty() && !constraintsSupportingVarNegative[i].isEmpty()) {
					//Connects all supporting constraints into one component
					//Set of all representatives of the previous components
					Set<Integer> oldRepresentatives = new HashSet<>();
					//Computes the new representative of this component such that the rank is maximum
					int representative = -1;
					int maxRank = -1;
					for (int j : constraintsSupportingVarPositive[i]) {
						int oldRepresentative = findRepresentativeAndCompress.apply(j);
						int oldRank = rank[oldRepresentative];
						if (maxRank < oldRank) {
							representative = oldRepresentative;
							maxRank = oldRank;
						}
						oldRepresentatives.add(oldRepresentative);
					}
					for (int j : constraintsSupportingVarNegative[i]) {
						int oldRepresentative = findRepresentativeAndCompress.apply(j);
						int oldRank = rank[oldRepresentative];
						if (maxRank < oldRank) {
							representative = oldRepresentative;
							maxRank = oldRank;
						}
						oldRepresentatives.add(oldRepresentative);
					}

					//Sets the representative for all constraints that support the current variable
					for (int oldrep : oldRepresentatives) {
						//Update representative if the old one was different
						if (oldrep != representative) {
							componentTreePredecessors[oldrep] = representative;
							//Updates rank if necessary
							if (rank[oldrep] == rank[representative]) {
								rank[representative]++;
							}
						}
					}
					numberComponents -= oldRepresentatives.size()-1;
				}
			}
			
			//Builds components consisting of nominal constraints
			Map<Integer, Set<NominalConstr>> representativeToComponent = new HashMap<Integer, Set<NominalConstr>>();
			for (int j = 0; j < componentTreePredecessors.length; j++) {
				int representative = findRepresentativeAndCompress.apply(j);
				if (!representativeToComponent.containsKey(representative)) {
					representativeToComponent.put(representative, new HashSet<NominalConstr>());
				}
				representativeToComponent.get(representative).add(nominalConstrs[j]);
			}
			//Only adds those to combinableConstraintsSets that have at least two constraints.
			//Otherwise, there is nothing to combine.
			combinableConstraintsSets = new ArrayList<Set<NominalConstr>>();
			for (Set<NominalConstr> component : representativeToComponent.values()) {
				if (component.size() > 1) {
					combinableConstraintsSets.add(component);
				}
			}
		}
	}
	
	/**
	 * Considers non-recyclable constraints and fixes the problematic variables to zero.
	 * Then recycles the resulting inequality and lifts the fixed variables into the recycled inequality.
	 */
	private class LiftingSeparator implements Separator{
		/**
		 * Constraints that are either already recyclable or are suitable for lifting.
		 */
		private List<LiftableConstr> liftableConstrs;
		
		/**
		 * Constructor computing liftable constraints.
		 * Constraints ax<=b are considered suitable for lifting if they support
		 * at least two uncertain variables with positive coefficients.
		 */
		public LiftingSeparator() throws GRBException, IOException {
			liftableConstrs = new ArrayList<LiftableConstr>();
			for (NominalConstr nominalConstr : nominalConstrs) {
				List<Integer> signs = new ArrayList<Integer>();
				if (nominalConstr.getSense() == GRB.LESS_EQUAL || nominalConstr.getSense() == GRB.EQUAL) {
					signs.add(1);
				}
				if (nominalConstr.getSense() == GRB.GREATER_EQUAL || nominalConstr.getSense() == GRB.EQUAL) {
					signs.add(-1);
				}
				for (Integer sign : signs) {
					int numPosVars = 0;
					for (int i = 0; i < nominalConstr.size(); i++) {
						if (nominalConstr.getCoefficient(i)*sign > 0) {
							numPosVars++;
						}
					}
					if (numPosVars >= 2) {
						liftableConstrs.add(new LiftableConstr(nominalConstr, sign));
					}
				}
			}
		}
		
		@Override
		public Collection<RobustInequality> computeRobustInequalities(double[] xValues, double[] pValues, double zValue) throws GRBException, IOException {
			List<RobustInequality> robustInequalities = new ArrayList<RobustInequality>();
			for (LiftableConstr liftableConstr : liftableConstrs) {
				//Only lifts if the violation before lifting is positive
				if (liftableConstr.computePreLiftingViolation(xValues, pValues, zValue) > algorithmParameters.getFeasibilityTolerance()) {
					robustInequalities.add(liftableConstr.computeRobustInequality(xValues, pValues, zValue));
				}
			}
			return robustInequalities;
		}
		
		/**
		 * Constraints that are either already recyclable or are suitable for lifting
		 * Constraints ax<=b are considered suitable for lifting if they support
		 * at least two uncertain variables with positive coefficients.
		 * Otherwise, the resulting recycled inequality will be dominated by the robustness constraints p+z>=dx.
		 */
		class LiftableConstr extends NominalConstr{
			/**
			 * Uncertain variables with positive coeff in the inequality.
			 * Paired with their index in the row for faster lookup of coefficients.
			 * Sorted non-increasing with respect to the deviation.
			 */
			private List<VariableIndexPair> positiveUncertainVarPairs;
			
			/**
			 * Variables with negative coeff in the inequality.
			 * Paired with their index in the row for faster lookup of coefficients.
			 */
			private List<VariableIndexPair> negativeVarPairs;
			
			/**
			 * Sign used to bring the constraint into the form ax <= b.
			 */
			private int sign;
			
			/**
			 * Constructor reading positive/negative variables and sorting positive variables.
			 */
			LiftableConstr(NominalConstr nominalConstr, int sign) throws GRBException {
				super(nominalConstr);
				this.sign = sign;
				positiveUncertainVarPairs = new ArrayList<VariableIndexPair>();
				negativeVarPairs = new ArrayList<VariableIndexPair>();
				for (int i = 0; i < this.size(); i++) {
					Variable variable = this.getNominalVariable(i);
					if (this.getCoefficient(i) >= 0) {
						if (variable.getDeviation() > 0) {
							positiveUncertainVarPairs.add(new VariableIndexPair(variable, i));
						}
					}
					else {
						negativeVarPairs.add(new VariableIndexPair(variable, i));
					}
				}
				//Sorts non-increasing
				positiveUncertainVarPairs.sort((var1, var2) -> Double.compare(var2.getVariable().getDeviation(), var1.getVariable().getDeviation()));
			}

			@Override
			public Double getCoefficient(int index) throws GRBException {
				return super.getCoefficient(index)*sign;
			}
			@Override
			public double getRHS() throws GRBException {
				return sign*super.getRHS();
			}
			@Override
			public char getSense() throws GRBException {
				return GRB.LESS_EQUAL;
			}
			
			/**
			 * Computes the violation of the (potentially invalid) recycled inequality before lifting.
			 */
			double computePreLiftingViolation(double[] xValues, double[] pValues, double zValue) throws GRBException {
				double violation = -this.getRHS() * zValue;
				for (VariableIndexPair positiveVarPair : positiveUncertainVarPairs) {
					Variable variable = positiveVarPair.getVariable();
					double coeff = this.getCoefficient(positiveVarPair.getIndex());
					double impact = coeff*(variable.getDeviation()/scalingQuotient*xValues[variable.getNominalIndex()] - coeff*pValues[variable.getUncertainIndex()]);
					if (impact > 0) {
						violation += impact;
					}
				}
				return violation;
			}
			
			/**
			 * Decides which variables to lift and computes their lifting coefficients.
			 */
			RobustInequality computeRobustInequality(double[] xValues, double[] pValues, double zValue) throws GRBException, IOException {
				//List of variables with positive coeff and positive impact on the violation that will be part of the recycled inequality 
				List<VariableIndexPair> chosenPositiveVarPairs = new ArrayList<VariableIndexPair>(positiveUncertainVarPairs.size());
				for (VariableIndexPair positiveVarPair : positiveUncertainVarPairs) {
					Variable variable = positiveVarPair.getVariable();
					double coeff = this.getCoefficient(positiveVarPair.getIndex());
					double impact = coeff*(variable.getDeviation()/scalingQuotient*xValues[variable.getNominalIndex()] - coeff*pValues[variable.getUncertainIndex()]);
					if (impact >= 0) {
						chosenPositiveVarPairs.add(positiveVarPair);
					}
				}
				RobustInequality robustInequality = new RobustInequality();
				for (VariableIndexPair positiveVarPair : chosenPositiveVarPairs) {
					Variable variable = positiveVarPair.getVariable();
					double coeff = this.getCoefficient(positiveVarPair.getIndex());
					robustInequality.addVariable(variable, coeff*variable.getDeviation()/scalingQuotient, coeff);
				}
				
				//Right-hand side that also contains terms of possible estimated variables with negative coeff
				double rhs = this.getRHS();
				
				boolean lifted = false;
				if (!negativeVarPairs.isEmpty()) {
					//List of variables with negative coeff that will be lifted
					//Initialized to be all such variables
					List<VariableIndexPair> liftingPairs = new ArrayList<VariableIndexPair>(negativeVarPairs);
					
					//Sorts the negative vars non-decreasing w.r.t. their solution value
					liftingPairs.sort((pair1, pair2) -> Double.compare(xValues[pair1.getVariable().getNominalIndex()], xValues[pair2.getVariable().getNominalIndex()]));
					
					//While the rhs is negative, we estimate the variables with the highest solution values
					while (rhs < 0) {
						rhs = rhs - this.getCoefficient(liftingPairs.remove(0).getIndex());
					}
					
					//Solves the fractional knapsack subproblem for lifting with the current rhs as capacity
					double[] knapsackSolution = solveFractionalKnapsack(chosenPositiveVarPairs, rhs, 0, 0, -1);
					//Value for the current rhs-capacity
					double knapsackValueRHS = knapsackSolution[0];
					//Weight of the items completely fitting into the knapsack
					double weightSumLastFitRHS = knapsackSolution[1];
					//Last index of variable that completely fits into the knapsack
					int indexLastFitRHS = (int) knapsackSolution[2];
										
					//Sum of (lifting coefficients * solution values) for all lifted variables 
					double liftedTermsSum = 0;
					boolean updateLiftedTermsSum = true;
					
					//Checks for all lifting variables whether they should better not be lifted
					for (int l = liftingPairs.size()-1; l >= 0; l--) {
						VariableIndexPair liftingPair = liftingPairs.get(l);
						if (xValues[liftingPair.getVariable().getNominalIndex()] < algorithmParameters.getIntegerFeasibilityTolerance()) {
							break;
						}
						
						//Updates liftedTermsSum if necessary
						if (updateLiftedTermsSum) {
							liftedTermsSum = 0;
							for (VariableIndexPair liftingPair2 : liftingPairs) {
								double increasedCapacity = rhs - this.getCoefficient(liftingPair2.getIndex());
								double knapsackValueIncreasedCapacity = solveFractionalKnapsack(chosenPositiveVarPairs, increasedCapacity, knapsackValueRHS, weightSumLastFitRHS, indexLastFitRHS)[0];
								
								double liftingCoeff = knapsackValueRHS - knapsackValueIncreasedCapacity;
								liftedTermsSum += liftingCoeff*xValues[liftingPair2.getVariable().getNominalIndex()];
							}
							updateLiftedTermsSum = false;
						}
						
						//New rhs if the variable is estimated
						double estimatedRHS = rhs - this.getCoefficient(liftingPair.getIndex());
						knapsackSolution = solveFractionalKnapsack(chosenPositiveVarPairs, estimatedRHS, knapsackValueRHS, weightSumLastFitRHS, indexLastFitRHS);
						double knapsackValueEstimatedRHS = knapsackSolution[0];
						double weightSumLastFitEstimatedRHS = knapsackSolution[1];
						int indexLastFitEstimatedRHS = (int) knapsackSolution[2];
						
						//Value if the variable is estimated instead of lifted
						double valueEstimation = zValue*this.getCoefficient(liftingPair.getIndex());
						for (VariableIndexPair liftingPair2 : liftingPairs) {
							if (liftingPair2 == liftingPair) {
								continue;
							}
							double increasedCapacity = estimatedRHS - this.getCoefficient(liftingPair2.getIndex());
							double knapsackValueIncreasedCapacity = solveFractionalKnapsack(chosenPositiveVarPairs, increasedCapacity, knapsackValueEstimatedRHS, weightSumLastFitEstimatedRHS, indexLastFitEstimatedRHS)[0];
							
							double liftingCoeff = knapsackValueEstimatedRHS - knapsackValueIncreasedCapacity;
							valueEstimation += liftingCoeff*xValues[liftingPair2.getVariable().getNominalIndex()];
						}
						
						//We estimate and remove from the list of lifting vars if the estimation value is greater 
						if (valueEstimation >= liftedTermsSum) {
							rhs = estimatedRHS;
							knapsackValueRHS = knapsackValueEstimatedRHS;
							weightSumLastFitRHS = weightSumLastFitEstimatedRHS;
							indexLastFitRHS = indexLastFitEstimatedRHS;
							updateLiftedTermsSum = true;
							liftingPairs.remove(l);
						}
					}
					
					//Adds the lifted variables to the inequality
					for (VariableIndexPair liftingPair : liftingPairs) {
						double increasedCapacity = rhs - this.getCoefficient(liftingPair.getIndex());
						double knapsackValueIncreasedCapacity = solveFractionalKnapsack(chosenPositiveVarPairs, increasedCapacity, knapsackValueRHS, weightSumLastFitRHS, indexLastFitRHS)[0];
						double liftingCoeff = knapsackValueRHS - knapsackValueIncreasedCapacity;
						robustInequality.addVariable(liftingPair.getVariable(), liftingCoeff, 0);
					}
					if (!liftingPairs.isEmpty()) {
						lifted = true;
					}
				}
				robustInequality.zCoeff = rhs;
				if (lifted) {
					double violation =  robustInequality.computeViolation(xValues, pValues, zValue);
					if (violation > algorithmParameters.getFeasibilityTolerance()) {
						writeOutput("##### Lifted constraint with violation "+violation);						
					}
				}
				return robustInequality;
			}
			
			/**
			 * Solves the fractional knapsack problem on the given (already sorted) vars
			 * currentValue, currentWeightSum, and currentIndexLastFit are given by an already computed sub-solution
			 * and are used for faster solving.
			 * Returns an array with solution value, weight sum of all completely fitting vars, and the index of the last fitting var.
			 */
			double[] solveFractionalKnapsack(List<VariableIndexPair> knapsackVars, double capacity, double currentValue, double currentWeightSum, int currentIndexLastFit) throws GRBException {
				while (currentIndexLastFit < knapsackVars.size()-1) {
					Variable variable = knapsackVars.get(currentIndexLastFit+1).getVariable();
					double coeff = this.getCoefficient(knapsackVars.get(currentIndexLastFit+1).getIndex());
					if (currentWeightSum + coeff <= capacity) {
						currentWeightSum += coeff;
						currentValue += coeff*variable.getDeviation()/scalingQuotient;
						currentIndexLastFit++;
					}
					else {
						currentValue += (capacity - currentWeightSum)*variable.getDeviation()/scalingQuotient;
						break;
					}
				}
				double[] solution = {currentValue, currentWeightSum, currentIndexLastFit};
				return solution;
			}
			

			/**
			 * Pair of variable and index combined for sorting.
			 */
			private class VariableIndexPair{
				Variable variable;
				int index;
				
				public VariableIndexPair(Variable variable, int index) {
					this.variable = variable;
					this.index = index;
				}
				Variable getVariable() {
					return variable;
				}
				int getIndex() {
					return index;
				}
			}
		}
	}
	
	
	/**
	 * Transforms each constraint of the constraint matrix into a recyclable inequality and stores them in a list if they are relevant.
	 * That is, the sum of the left-hand side coefficients of the original inequality exceeds the right-hand side value.
	 * Otherwise, the corresponding recycled inequality is dominated by the robustness constraints x*deviation <= p+z. 
	 */
	private List<RecyclableInequality> computeRecyclableInequalities() throws GRBException {
		List<RecyclableInequality> preparedRecycledInequalities = new ArrayList<RecyclableInequality>();
		
		//Iterates over all nominal constraints
		for (NominalConstr constraint : nominalConstrs) {
			//We transform the row into the standard form ax<=b.
			//If we have ax>=b then we multiply everything with -1 to obtain -ax<=-b
			//If we have ax=b then we consider ax<=b and -ax<=-b.
			List<Integer> signs = new ArrayList<Integer>(2);
			if (constraint.getSense() == GRB.LESS_EQUAL || constraint.getSense() == GRB.EQUAL) {
				signs.add(1);
			}
			if (constraint.getSense() == GRB.GREATER_EQUAL || constraint.getSense() == GRB.EQUAL) {
				signs.add(-1);
			}
			
			//Obtains the lhs of the current constraint
			for (int sign : signs) {
				//Adapts the constraint such that we have ax<=b or -ax<=-b 
				double rhs = sign*constraint.getRHS();
				
				//Indices of the variables in the recycled inequality.
				List<Variable> recVariables = new ArrayList<Variable>();
				//Coefficients of the variables in the recycled inequality.
				List<Double> recCoeffs = new ArrayList<Double>();
				
				//The recycled inequality is only useful if the sum over all xCoeffs is greater than rhs.
				double coeffSum = 0;
				
				//We transform the constraint into a recyclable inequality that only contains
				//uncertain variables with positive coeffs at the lhs. All other variables
				//cannot be used and have thus to be moved to the rhs. Setting these at their bounds
				//gives a valid inequality with the desired property.
				for (int i = 0; i < constraint.size(); i++) {
					Variable variable = constraint.getNominalVariable(i);
					double coeff = sign*constraint.getCoefficient(i);
					if (coeff > 0){
						//If the variable is uncertain and the coefficient is positive, it can be added to the recycled inequality.
						if (variable.getDeviation() > 0) {
							recVariables.add(variable);
							recCoeffs.add(coeff);
							coeffSum += coeff;
						}
						//If the variable is certain and has a positive coefficient, we set it to its lower bound
						//in order to maximize satisfaction of the constraint.
						else {
							rhs -= coeff*variable.getModelVariable().get(GRB.DoubleAttr.LB);
						}
					}
					//If the variable has a negative coefficient, we set it to its upper bound
					//in order to maximize satisfaction of the constraint.
					else {
						rhs -= coeff*variable.getModelVariable().get(GRB.DoubleAttr.UB);
					}
				}
				//The recycled inequality is only useful if the sum over all xCoeffs is greater than rhs.
				if (coeffSum > rhs) {
					preparedRecycledInequalities.add(new RecyclableInequality(recVariables, recCoeffs, rhs));
				}
			}
		}
		return preparedRecycledInequalities;
	}
	
	/**
	 * Models robust inequalities of the form
	 * zCoeff*z + sum(pCoeff*p) - sum(xCoeff*x) >= constant 
	 */
	private class RobustInequality {
		List<Variable> variables = new ArrayList<Variable>();
		List<Double> xCoeffs = new ArrayList<Double>();
		List<Double> pCoeffs = new ArrayList<Double>();
		double zCoeff;
		double constant;

		public List<Variable> getVariables() {
			return variables;
		}

		public List<Double> getXCoeffs() {
			return xCoeffs;
		}

		public List<Double> getPCoeffs() {
			return pCoeffs;
		}

		public double getZCoeff() {
			return zCoeff;
		}

		public double getConstant() {
			return constant;
		}
		
		void addVariable(Variable variable, double xCoeff, double pCoeff) {
			variables.add(variable);
			xCoeffs.add(xCoeff);
			pCoeffs.add(pCoeff);
		}
		
		/**
		 * Returns lhs of the robust inequality: lhsExpr>=constant.
		 */
		GRBLinExpr getLHSexpr() {
			GRBLinExpr lhsExpr = new GRBLinExpr();
			lhsExpr.addTerm(getZCoeff(), z);
			for (int i = 0; i < getVariables().size(); i++) {
				Variable variable = getVariables().get(i);
				if (variable.getDeviation() > 0) {
					lhsExpr.addTerm(getPCoeffs().get(i), p[variable.getUncertainIndex()]);
				}
				lhsExpr.addTerm(-getXCoeffs().get(i), variable.getModelVariable());
			}
			return lhsExpr;
		}
		
		/**
		 * Computes the violation of the recycled inequality.
		 */
		double computeViolation(double[] xValues, double[] pValues, double zValue) {
			double violation = constant - zCoeff * zValue;
			for (int i = 0; i < variables.size(); i++) {
				Variable variable = variables.get(i);
				violation += xCoeffs.get(i)*xValues[variable.getNominalIndex()];
				if (variable.getDeviation() > 0) {
					violation -= pCoeffs.get(i)*pValues[variable.getUncertainIndex()];
				}
			}
			return violation;
		}
		
		void writeSeparationMessage(double violation) throws IOException {
			writeOutput("##### Found robust inequality with violation "+violation);
		}
	}
	
	
	/**
	 * Nominal recyclable inequalities
 	 * Offers method for computing a corresponding recycled (sub-)inequality,
 	 * obtained by dropping some x, that maximizes the violation for given solutions x,p.
	 */
	private class RecyclableInequality{
		List<Variable> variables;
		List<Double> coeffs;
		double rhs;

		/**
		 * Constructor directly obtaining variables, coefficients and rhs
		 */
		RecyclableInequality(List<Variable> variables, List<Double> coeffs, double rhs) throws GRBException {
			this.variables = variables;
			this.coeffs = coeffs;
			this.rhs = rhs;
		}

		/**
		 * Computes the corresponding recycled (sub-)inequality, obtained by dropping some x,
		 * that maximizes the violation for given solutions x,p.
		 */
		RobustInequality getRecycledInequality(double[] xValues, double[] pValues) {
			RobustInequality robustInequality = new RobustInequality();
			robustInequality.zCoeff = rhs;
			for (int i = 0; i < variables.size(); i++) {
				Variable variable = variables.get(i);
				if (coeffs.get(i)*pValues[variable.getUncertainIndex()] <= coeffs.get(i)*variable.getDeviation()/scalingQuotient*xValues[variable.getNominalIndex()]) {
					robustInequality.addVariable(variable, coeffs.get(i)*variable.getDeviation()/scalingQuotient, coeffs.get(i));
				}
			}
			return robustInequality;
		}
		
		/**
		 * Computes the corresponding recycled inequality
		 */
		RobustInequality getRecycledInequality() {
			RobustInequality robustInequality = new RobustInequality();
			robustInequality.zCoeff = rhs;
			for (int i = 0; i < variables.size(); i++) {
				Variable variable = variables.get(i);
				robustInequality.addVariable(variable, coeffs.get(i)*variable.getDeviation()/scalingQuotient, coeffs.get(i));
			}
			return robustInequality;
		}
	}
	
	/**
	 * Constraints of the nominal problem with getters for simplicity
	 */
	private class NominalConstr {
		private GRBConstr grbConstr;
		private GRBLinExpr row;
		
		public NominalConstr(NominalConstr nominalConstr) throws GRBException {
			this(nominalConstr.grbConstr);
		}
		
		public NominalConstr(GRBConstr grbConstr) throws GRBException {
			this.grbConstr = grbConstr;
			this.row = robustProblem.getModel().getRow(grbConstr);
		}
		
		public Variable getNominalVariable(int index) throws GRBException {
			return modelVarToVariable.get(row.getVar(index));
		}
		public Double getCoefficient(int index) throws GRBException {
			return row.getCoeff(index);
		}
		public double getRHS() throws GRBException {
			return grbConstr.get(GRB.DoubleAttr.RHS);
		}
		public char getSense() throws GRBException {
			return grbConstr.get(GRB.CharAttr.Sense);
		}
		public int size() throws GRBException {
			return row.size();
		}
	}
	
	
	/**
	 * Specifies separation and Gurobi cuts strategies;
	 */
	public static class RecyclingStrategies extends AlgStrategies {
		/**
		 * Enum type specifying whether we solve the integer program or the relaxation.
		 */
		public enum Relaxation{
			SOLVE_INTEGER,
			SOLVE_RELAXED;
		}
		
		/**
		 * Enum type specifying whether and how to add recycled inequalities directly to the formulation.
		 * DIRRECYCLE_NONE adds no recycled inequalities.
		 * DIRRECYCLE_CONSTRAINTS adds recycled nominal constraints.
		 */
		public enum DirectRecyclingStrategy{
			DIRRECYCLE_NONE,
			DIRRECYCLE_CONSTRAINTS;
		}
		
		/**
		 * Enum type specifying how to separate recycled inequalities.
		 * SEPRECYCLE_NONE separates no recycled inequalities.
		 * SEPARATE_CONSTRAINTS separates violated recycled nominal (sub-)constraints within the root node.
		 * SEPARATE_CLIQUES works on the conflict graph and separates violated recycled clique inequalities.
		 * SEPARATE_LP combines nominal inequalities into a recyclable inequality by solving an LP.
		 * SEPARATE_LIFT fixes variables with negative coefficients before recycling and lifts them into the recycled inequality.
		 */
		public enum SeparationRecyclingStrategy{
			SEPRECYCLE_NONE,
			SEPRECYCLE_CONSTRAINTS,
			SEPRECYCLE_CLIQUES,
			SEPRECYCLE_LP,
			SEPRECYCLE_LIFT;
		}

		/**
		 * Enum type specifying whether we allow GUROBI to add its own cuts.
		 */
		public enum GurobiCutStrategy{
			GCUTS_ENABLE,
			GCUTS_DISABLE;
		}
		
		Relaxation relaxation;
		DirectRecyclingStrategy directRecyclingStrategy;
		SeparationRecyclingStrategy separationStrategy;
		GurobiCutStrategy gurobiCutStrategy;
		
		/**
		 * Constructor obtaining arguments which are matched to the enums defining strategies.
		 */
		public RecyclingStrategies(List<String> argList, AlgorithmParameters algorithmParameters) throws IOException {
			super(argList, algorithmParameters);
		}
		
		@Override
		void setDefaultStrategies() {
			relaxation = Relaxation.SOLVE_INTEGER;
			directRecyclingStrategy = DirectRecyclingStrategy.DIRRECYCLE_NONE;
			separationStrategy = SeparationRecyclingStrategy.SEPRECYCLE_CONSTRAINTS;
			gurobiCutStrategy = GurobiCutStrategy.GCUTS_ENABLE;			
		}
		
		public Relaxation getRelaxation() {
			return relaxation;
		}
		public void setRelaxation(Relaxation relaxation) {
			this.relaxation = relaxation;
		}
		
		public DirectRecyclingStrategy getDirectRecyclingStrategy() {
			return directRecyclingStrategy;
		}
		public void setDirectRecyclingStrategy(DirectRecyclingStrategy directRecyclingStrategy) {
			this.directRecyclingStrategy = directRecyclingStrategy;
		}

		public SeparationRecyclingStrategy getSeparationStrategy() {
			return separationStrategy;
		}
		public void setSeparationStrategy(SeparationRecyclingStrategy separationStrategy) {
			this.separationStrategy = separationStrategy;
		}

		public GurobiCutStrategy getGurobiCutStrategy() {
			return gurobiCutStrategy;
		}
		public void setGurobiCutStrategy(GurobiCutStrategy gurobiCutStrategy) {
			this.gurobiCutStrategy = gurobiCutStrategy;
		}
	}
}