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
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

/**
 * This class implements the standard reformulation by Bertsimas and Sim together with the recycling of valid inequalities
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
	RecyclingStrategies recyclingStrategies;
	
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
	private GRBConstr[] nominalConstrs;
	
	/**
	 * A list of recyclable inequalities, prepared for faster separation.
	 */
	private List<RecyclableInequality> recyclableInequalities;
	
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
		
		//If we solve the integer problem and separate cuts then we add the callback for this
		//In this case, we set the precrush value to 1 such that cuts can always be used by gurobi
		if (recyclingStrategies.relaxation == Relaxation.SOLVE_INTEGER && separator != null){
			model.setCallback(new Callback());
			model.set(GRB.IntParam.PreCrush, 1);
		}
		
		//Solves the problem in the remaining time
		model.update();
		robustProblem.solve(getRemainingTime());
		primalBound = robustProblem.getPrimalBound();
		dualBound = robustProblem.getDualBound();
		
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
			double[] xValues = model.get(GRB.DoubleAttr.X, robustProblem.getUncertainModelVariables());
			double[] pValues = model.get(GRB.DoubleAttr.X, p);
			Double zValue = z.get(GRB.DoubleAttr.X);
			
			Collection<RecyclableInequality> separatedRecyclableInequalities = separator.getRecyclableInequalities(xValues, pValues, zValue);

			boolean foundViolated = false;
			for (RecyclableInequality recyclableInequality : separatedRecyclableInequalities) {
				double violation = recyclableInequality.computeViolation(scalingQuotient, xValues, pValues, zValue);
				if (violation > algorithmParameters.getFeasibilityTolerance()) {
					foundViolated = true;
					writeOutput("##### Found cut with violation = "+violation+"\n");
					model.addConstr(recyclableInequality.getRecycledExpr(scalingQuotient, xValues, pValues), GRB.GREATER_EQUAL, 0, "violatedRecycledIneq");
				}
			}
			if (!foundViolated) {
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
		
		//Stores best solution found
		solution = new LinkedHashMap<Variable, Double>(robustProblem.getNominalVariables().length);
		for (int i = 0; i < robustProblem.getNominalVariables().length; i++) {
			solution.put(robustProblem.getNominalVariables()[i], robustProblem.getNominalVariablesSolutionValues()[i]);
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
		nominalConstrs = model.getConstrs();
		
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
		//The size of p may vary depending on whether we merge variables into cliques.
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
		
		//Prepares the list of recyclable inequalities.
		if (recyclingStrategies.directRecyclingStrategy == RecyclingStrategies.DirectRecyclingStrategy.DIRRECYCLE_CONSTRAINTS
				|| recyclingStrategies.separationStrategy == RecyclingStrategies.SeparationRecyclingStrategy.SEPRECYCLE_CONSTRAINTS) {
			recyclableInequalities = prepareRecyclableInequalities();
		}
		
		//Adds recycled inequalities to the model if the corresponding strategy is chosen.
		if (recyclingStrategies.directRecyclingStrategy == RecyclingStrategies.DirectRecyclingStrategy.DIRRECYCLE_CONSTRAINTS) {
			for (RecyclableInequality recyclableInequality: recyclableInequalities) {
				model.addConstr(recyclableInequality.getRecycledExpr(scalingQuotient), GRB.GREATER_EQUAL, 0, "robustness_cut");
			}
		}
		
		//Sets separators if the corresponding options are chosen 
		if (recyclingStrategies.separationStrategy == RecyclingStrategies.SeparationRecyclingStrategy.SEPRECYCLE_CONSTRAINTS) {
			separator = new PreparedSeparator();
		}
		else if (recyclingStrategies.separationStrategy == RecyclingStrategies.SeparationRecyclingStrategy.SEPRECYCLE_CLIQUES) {
			separator = new CliqueSeparator();
		}
		
		//Disables gurobis cuts if the option is chosen
		if (recyclingStrategies.gurobiCutStrategy == RecyclingStrategies.GurobiCutStrategy.GCUTS_DISABLE) {
			model.set(GRB.IntParam.Cuts, 0);
		}
		
		//Relaxes variables if the option is chosen
		if (recyclingStrategies.relaxation == Relaxation.SOLVE_RELAXED) {
			for (GRBVar grbVar : robustProblem.getNominalModelVariables()) {
				grbVar.set(GRB.CharAttr.VType, GRB.CONTINUOUS);
			}
		}
		model.update();
	}
	
	/**
	 * Callback that is called at the root node of the branch and bound.
	 * Passes the current solution x,p,z to a separator and adds the recycled inequality if it is violated.
	 */
	private class Callback extends GRBCallback {
		
		@Override
		protected void callback() {
			try {
				//Checks whether we are in the root node and are allowed to add user cuts
				if (where == GRB.CB_MIPNODE && getDoubleInfo(GRB.CB_MIPNODE_NODCNT) == 0 && getIntInfo(GRB.CB_MIPNODE_STATUS) == GRB.Status.OPTIMAL) {
					//Queries the values of the current relaxed solution.
					double[] xValues = getNodeRel(robustProblem.getUncertainModelVariables());
					double[] pValues = getNodeRel(p);
					Double zValue = getNodeRel(z);
					
					//Passes the current solution to the separator 
					Collection<RecyclableInequality> recyclableInequalities = separator.getRecyclableInequalities(xValues, pValues, zValue);
					
					//Adds recycled inequalities if violated
					for (RecyclableInequality recyclableInequality : recyclableInequalities) {
						double violation = recyclableInequality.computeViolation(scalingQuotient, xValues, pValues, zValue);
						if (violation > algorithmParameters.getFeasibilityTolerance()) {
							writeOutput("##### Found cut with violation = "+violation);
							addCut(recyclableInequality.getRecycledExpr(scalingQuotient, xValues, pValues), GRB.GREATER_EQUAL, 0);
						}
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
	 * Interface for separators which obtain a solution x,p,z and return recyclable inequalities.
	 */
	private interface Separator {
		Collection<RecyclableInequality> getRecyclableInequalities(double[] xValues, double[] pValues, double zValue) throws GRBException;
	}
	
	/**
	 * Provides the already preprocessed recyclable inequalities.
	 */
	private class PreparedSeparator implements Separator {
		@Override
		public List<RecyclableInequality> getRecyclableInequalities(double[] xValues, double[] pValues, double zValue) {
			return recyclableInequalities;
		}
	}
	
	/**
	 * Computes a conflict graph for the nominal problem and uses it to separate recycled clique inequalities.
	 */
	private class CliqueSeparator implements Separator{
		/**
		 * Maps variables to the corresponding neighborhood in the conflict graph.
		 */
		Map<Variable, Set<Variable>> uncertainVariableToConflictingNeighborhood;
		
		/**
		 * Constructor computing the conflict graph and neighborhoods of variables.
		 */
		CliqueSeparator() throws GRBException, IOException {
			ConflictGraph conflictGraph = new ConflictGraph(robustProblem.getModel(), robustProblem.getUncertainModelVariables(), algorithmParameters);
			Map<GRBVar, Variable> modelVarToVariableMap = new HashMap<GRBVar, Variable>();
			for (Variable uncertainVariable : robustProblem.getUncertainVariables()) {
				modelVarToVariableMap.put(uncertainVariable.getModelVariable(), uncertainVariable);
			}
			
			uncertainVariableToConflictingNeighborhood = new HashMap<Variable, Set<Variable>>();
			for (GRBVar uncertainModelVariable : robustProblem.getUncertainModelVariables()) {
				Set<Variable> neighborhood = new HashSet<Variable>();
				for (GRBVar neighborVar : conflictGraph.getNeighbors(uncertainModelVariable)) {
					neighborhood.add(modelVarToVariableMap.get(neighborVar));
				}
				uncertainVariableToConflictingNeighborhood.put(modelVarToVariableMap.get(uncertainModelVariable), neighborhood);
			}
			
		}
		
		@Override
		public Set<RecyclableInequality> getRecyclableInequalities(double[] xValues, double[] pValues, double zValue) {
			//Assigns each variable its impact on the violation in the recycled inequality.
			Function<Variable, Double> violationImpactFunc = variable -> xValues[variable.getUncertainIndex()]*variable.getDeviation()/scalingQuotient - pValues[variable.getUncertainIndex()];
			//The clique inequalities to be recycled
			Set<RecyclableInequality> cliqueInequalities = new HashSet<RecyclableInequality>();
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
					cliqueInequalities.add(new RecyclableInequality(new ArrayList<Variable>(clique), 1, 1));
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
	 * Transforms each constraint of the constraint matrix into a recyclable inequality
	 * and stores them in a list if they are interesting for recycling.
	 * That is, the sum of the left-hand side coefficients exceeds the right-hand side value.
	 * Otherwise, the corresponding recycled inequality is dominated by the robustness constraints x*deviation <= p+z. 
	 */
	private List<RecyclableInequality> prepareRecyclableInequalities() throws GRBException {
		//For construction of recyclable inequalities, we need to map the uncertain modelvars to their indices.
		Variable[] uncertainVariables = robustProblem.getUncertainVariables();
		Map<GRBVar, Integer> uncertainModelVariableToIndex = new HashMap<GRBVar, Integer>();
		for (int i = 0; i < uncertainVariables.length; i++) {
			uncertainModelVariableToIndex.put(uncertainVariables[i].getModelVariable(), i);
		}
		
		List<RecyclableInequality> preparedCuts = new ArrayList<RecyclableInequality>();		
		
		//Iterates over all nominal inequalities
		for (GRBConstr constraint : nominalConstrs) {
			//We transform the row into the standard form ax<=b.
			//If we have ax>=b then we multiply everything with -1 to obtain -ax<=-b
			//If we have ax=b then we consider ax<=b and -ax<=-b.
			List<Integer> signs = new ArrayList<Integer>(2);
			if (constraint.get(GRB.CharAttr.Sense) == GRB.LESS_EQUAL || constraint.get(GRB.CharAttr.Sense) == GRB.EQUAL) {
				signs.add(1);
			}
			if (constraint.get(GRB.CharAttr.Sense) == GRB.GREATER_EQUAL || constraint.get(GRB.CharAttr.Sense) == GRB.EQUAL) {
				signs.add(-1);
			}
			
			//Obtains the lhs of the current constraint
			GRBLinExpr nominalExpr = robustProblem.getModel().getRow(constraint);
			for (int sign : signs) {
				//Adapts the constraint such that we have ax<=b or -ax<=-b 
				double rhs = sign*constraint.get(GRB.DoubleAttr.RHS);
				Map<GRBVar, Double> varToCoeff = new HashMap<GRBVar,Double>();
				for (int i = 0; i< nominalExpr.size(); i++) {
					varToCoeff.put(nominalExpr.getVar(i), sign*nominalExpr.getCoeff(i));
				}
				
				//Indices of the variables in the recyclable inequality.
				List<Variable> variables = new ArrayList<Variable>();
				//Coefficients of the variables in the recyclable inequality.
				List<Double> xCoeffs = new ArrayList<Double>();
				
				//The recyclable inequality is only useful if the sum over all xCoeffs is greater than rhs.
				double coeffSum = 0;
				
				//We transform the constraint into a recyclable inequality that only contains
				//uncertain variables with positive coeffs at the lhs. All other variables
				//cannot be used and have thus be moved to the rhs. Setting these at their bounds
				//gives a valid inequality with the desired property.
				for (GRBVar grbVar : varToCoeff.keySet()) {
					double coeff = varToCoeff.get(grbVar);
					if (coeff > 0){
						//If the variable is uncertain and the coefficient is positive, it can be added to the recyclable inequality.
						if (uncertainModelVariableToIndex.containsKey(grbVar)) {
							int varIndex = uncertainModelVariableToIndex.get(grbVar);
							variables.add(uncertainVariables[varIndex]);
							xCoeffs.add(coeff);
							coeffSum += coeff;
						}
						//If the variable is certain and has a positive coefficient, we set it to its lower bound
						//in order to maximize satisfaction of the constraint.
						else {
							rhs -= coeff*grbVar.get(GRB.DoubleAttr.LB);
						}
					}
					//If the variable has a negative coefficient, we set it to its upper bound
					//in order to maximize satisfaction of the constraint.
					else {
						rhs -= coeff*grbVar.get(GRB.DoubleAttr.UB);
					}
				}
				//The recyclable inequality is only useful if the sum over all xCoeffs is greater than rhs.
				if (coeffSum > rhs) {
					//If all xCoeffs are same then only pass this coeff to the constructor.
					if (xCoeffs.stream().allMatch(val -> val == xCoeffs.get(0))) {
						preparedCuts.add(new RecyclableInequality(variables, xCoeffs.get(0), rhs));
					}
					else {
						preparedCuts.add(new RecyclableInequality(variables, xCoeffs, rhs));
					}
				}
			}
		}
		return preparedCuts;
	}
	
	/**
	 * Preprocessed recyclable inequalities used for fast construction of recycled inequalities.
 	 * Contains the right-hand side of the inequality as well as the indices of the supported variables x and their coefficients in the recyclable inequality.
 	 * Offers method for computing a corresponding recycled (sub-)inequality that maximizes the violation for given solutions x,p.
 	 * Also offers a method for computing the violation of this inequality.
	 */
	private class RecyclableInequality {
		private List<Variable> variables;
		private List<Double> xCoeffs;
		private Double commonXCoeff;
		private double RHS;
		
		/**
		 * Constructor explicitly stating the coefficients of variables.
		 */
		RecyclableInequality(List<Variable> variables, List<Double> xCoeffs, double RHS) {
			this.variables = variables;
			this.xCoeffs = xCoeffs;
			this.RHS = RHS;
		}
		
		/**
		 * Constructor implicitly stating the coefficients of variables by a common coefficient for all variables. 
		 */
		RecyclableInequality(List<Variable> variables, double commonXCoeff, double RHS) {
			this.variables = variables;
			this.commonXCoeff = commonXCoeff;
			this.RHS = RHS;
		}
		
		/**
		 * Returns lhs of the recycled inequality: recycledExpr>=0.
		 * Only includes variables that do not negatively impact the violation for a given solution of x,p. That is, p <= deviation/scalingQuotient*x.
		 */
		GRBLinExpr getRecycledExpr(double scalingQuotient, double[] xValues, double[] pValues) {
			GRBLinExpr recycledExpr = new GRBLinExpr();
			recycledExpr.addTerm(RHS, z);
			for (int i = 0; i < variables.size(); i++) {
				Variable variable = variables.get(i);
				if (xValues == null || pValues[variable.getUncertainIndex()] <= variable.getDeviation()/scalingQuotient * xValues[variable.getUncertainIndex()]) {
					Double xCoeff = commonXCoeff;
					if (xCoeffs != null) {
						xCoeff = xCoeffs.get(i);
					}
					recycledExpr.addTerm(xCoeff, p[variable.getUncertainIndex()]);
					recycledExpr.addTerm(-xCoeff*variable.getDeviation()/scalingQuotient, variable.getModelVariable());
				}
			}
			return recycledExpr;
		}
		
		/**
		 * Returns lhs of the recycled inequality: recycledExpr>=0.
		 */
		GRBLinExpr getRecycledExpr(double scalingQuotient) {
			return getRecycledExpr(scalingQuotient, null, null);
		}
		
		/**
		 * Computes the violation of the optimal recycled sub-inequality, only consisting of variables that positively contribute to the violation.
		 */
		double computeViolation(double scalingQuotient, double[] xValues, double[] pValues, double zValue) {
			double violation = - RHS * zValue;
			for (int i = 0; i < variables.size(); i++) {
				Variable variable = variables.get(i);
				if (pValues[variable.getUncertainIndex()] < variable.getDeviation()/scalingQuotient * xValues[variable.getUncertainIndex()]) {
					Double xCoeff = commonXCoeff;
					if (xCoeffs != null) {
						xCoeff = xCoeffs.get(i);
					}
					violation += xCoeff * (variable.getDeviation()/scalingQuotient * xValues[variable.getUncertainIndex()] - pValues[variable.getUncertainIndex()]);
				}
			}
			return violation;
		}
	}
	
	
	/**
	 * Specifies separation and gurobi cuts strategies;
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
		 * DIRRECYCLE_CLIQUES works on the conflict graph and adds recycled maximal clique inequalities.
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
		 */
		public enum SeparationRecyclingStrategy{
			SEPRECYCLE_NONE,
			SEPRECYCLE_CONSTRAINTS,
			SEPRECYCLE_CLIQUES;
			
//			@Override
			static String getSimpleName(){
				return "";
			}
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
			directRecyclingStrategy = DirectRecyclingStrategy.DIRRECYCLE_CONSTRAINTS;
			separationStrategy = SeparationRecyclingStrategy.SEPRECYCLE_CONSTRAINTS;
			gurobiCutStrategy = GurobiCutStrategy.GCUTS_ENABLE;			
		}
		
		public Relaxation getRelaxation() {
			return relaxation;
		}
		public void setRelaxation(Relaxation relaxation) {
			this.relaxation = relaxation;
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