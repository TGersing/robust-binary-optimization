package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import alg.AbstractAlgorithm.AlgorithmParameters;
import alg.AlgBranchAndBound.BnBCallbackIntegerSubproblems;
import alg.SubproblemBoundedGurobi.BoundingSubproblemsStrategies.BoundingTerminationStrategy;
import alg.SubproblemBoundedGurobi.BoundingSubproblemsStrategies.BoundingLPWarmStartStrategy;
import alg.SubproblemGurobi.SubproblemsStrategies.ImprovingZStrategy;
import gurobi.GRB;
import gurobi.GRBCallback;
import gurobi.GRBConstr;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBVar;
import util.BnBNode;
import util.CliquePartitioning;
import util.PossibleZ;
import util.Variable;
import util.Variable.Type;

/**
 * This class implements the integer and relaxed robust subproblems solved during @AlgBranchAndBound.
 * 
 * @author Timo Gersing
 */
public class SubproblemBoundedGurobi extends SubproblemGurobi{
	/**
	 * The current lower bound on z.
	 */
	private PossibleZ lowerBoundZ;
	
	/**
	 * The current upper bound on z.
	 */
	private PossibleZ upperBoundZ;
	
	/**
	 * Contains the robustness constraints using the lower bound for each uncertain variable or clique.
	 */
	private GRBConstr[] robustnessConstraintsLowerBound;
	
	/**
	 * Contains the robustness constraints using the upper bound for each uncertain variable or clique.
	 */
	private GRBConstr[] robustnessConstraintsUpperBound;
	
	/**
	 * Objects representing constraints whose slack were non-basic in the parent model.
	 * Used to add the same constraints to the current model.
	 */
	private List<NonBasicSlackConstraint> parentNonBasicSlackConstraints = new ArrayList<NonBasicSlackConstraint>();
	
	/**
	 * The List of GRBConstr added to the model to include non-basic slack constraints from the parent model.
	 */
	private List<GRBConstr> parentNonBasicSlackGRBConstraints = new ArrayList<GRBConstr>();
	
	/**
	 * Indicates whether an uncertain variable cannot be used for substituting p using the current bound.
	 * Not possible for all, as they need to remain unchanged for warm starting.
	 */
	private boolean[] isIndexNonSubstitutable;
	
	/**
	* Upper bounds with which the p variables have been substituted.
	* Can vary for different variables, as some need to be unchanged for warm starting.
	*/
	private double[] upperBoundsUsedForSubstituting;

	/**
	 * Constraints of the nominal problem that need to be stored such that we can report whether they are basic or not.
	 */
	private GRBConstr[] nominalConstrs;

	/**
	 * The variables p from the (clique) reformulation.
	 */
	private GRBVar[] p;
	
	/**
	 * The variable z from the (clique) reformulation.
	 */
	private GRBVar z;
	
	/**
	 * The value of the variable z.
	 * Only available if z is part of the model and was not omitted by Lagrange relaxation
	 */
	private Optional<Double> zSolutionValue;
	
	/**
	 * The values of the p variables in the computed solution (after resubstituting).
	 * Only available if p is part of the model and was not omitted by Lagrange relaxation
	 */
	private double[] pSolutionValues;
	
	/**
	 * The value of the variable z for the improved solution.
	 */
	private Optional<Double> improvedZ;
	
	/**
	 * The primal bound provided by a potentially improved solution.
	 */
	private double improvedPrimalBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
	
	/**
	 * Values of variables in a new incumbent. Needs to be stored because an improved sub-optimal solution
	 * might be better than the optimal solution of the subproblem. 
	 */
	private double[] improvedNominalVariablesSolutionValues;
	
	/**
	 * Decides whether we solve the LP relaxation with robustness constraints added or not.
	 * If the constraints are omitted, solving the LP is much faster but yields a worse dual bound.
	 */
	private boolean applyLagrangeRelaxation;
	
	/**
	 * The Lagrange multiplier used for the relaxation.
	 */
	double lagrangeMultipier;
	
	/**
	 * Specifies strategies.
	 */
	private BoundingSubproblemsStrategies boundingStrategies;
	
	/**
	 * Callback passed to the integer subproblem.
	 * Reports bounds from the integer subproblem to the master and asks for termination.
	 */
	private BnBCallbackIntegerSubproblems bnbCallback;
	
	/**
	 * States whether variables are currently continuous or integer.
	 */
	private boolean continuous = false;
	
	/**
	 * Constructor receiving paths to problem files as well as parameters and strategies.
	 */
	protected SubproblemBoundedGurobi(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, BoundingSubproblemsStrategies boundingStrategies) throws IOException, GRBException {
		super(problemPath, robustPath, algorithmParameters);
		this.boundingStrategies = boundingStrategies;
	}
	
	/**
	 * Sets the clique partitioning and adds merged clique constraints to the nominal constraints.
	 */
	@Override
	void setCliquePartitioning(CliquePartitioning cliquePartitioning) {
		this.cliquePartitioning = cliquePartitioning;
		try {
			for (int index : cliquePartitioning.getMergedCliqueIndices()) {
				List<Integer> clique = cliquePartitioning.getCliques().get(index);
				GRBLinExpr cliqueExpr = new GRBLinExpr();
				for (Integer varIndex : clique) {
					cliqueExpr.addTerm(1, uncertainModelVariables[varIndex]);
				}
					model.addConstr(cliqueExpr, GRB.LESS_EQUAL, 1, "MergedCliqueConstraint"+index);
			}
			model.update();
		} catch (GRBException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Resets the problem to an unsolved state.
	 */
	@Override
	protected void resetProblem() throws GRBException {
		super.resetProblem();
		zSolutionValue = Optional.empty();
		improvedZ = Optional.empty();
		pSolutionValues = null;
	}
	
	/**
	 * Solves the linear relaxation of the subproblem.
	 */
	void solveRelaxed(Optional<Double> timeLimit, BnBNode node) throws IOException, GRBException {
		//Sets all variables to continuous
		if (!continuous) {
			for (GRBVar grbVar : nominalModelVariables) {
				grbVar.set(GRB.CharAttr.VType, GRB.CONTINUOUS);
			}
			//Removes callback.
			model.setCallback(null);
			continuous = true;
		}
		
		//Warm starts the LP if the option is chosen
		if (boundingStrategies.warmStartStrategy == BoundingLPWarmStartStrategy.WARMSTART_ENABLE) {
			//Checks whether there exists a basis from which we can start
			if (node.getParentVBasisNominal() != null) {
				//Sets basis for nominal variables and slack variables
				model.set(GRB.IntAttr.VBasis, nominalModelVariables, node.getParentVBasisNominal());
				model.set(GRB.IntAttr.CBasis, nominalConstrs, node.getParentCBasisNominal());
				

				//All new constraints that were not in the parent model are set to basic.
				//This ensures that a basis remains a basis.
				for (GRBConstr tempConstr : temporaryConstraints) {
					tempConstr.set(GRB.IntAttr.CBasis, 0);
				}

				//If slacks of constraints were non-basic in the parent model, then we set them here also to non-basic.
				//This ensures that a basis remains a basis.
				for (GRBConstr parentConstr : parentNonBasicSlackGRBConstraints) {
					parentConstr.set(GRB.IntAttr.CBasis, -1);
				}
				
				if (applyLagrangeRelaxation) {
					//If they were previously added to the model (now not contributing to constraints or objective), 
					//set p and z to non-basic such that we obtain a primal basis.
					if (z != null) {
						z.set(GRB.IntAttr.VBasis, -1);
						for (int i = 0; i < p.length; i++) {
							p[i].set(GRB.IntAttr.VBasis, -1);
						}
					}
					
					//As this is a primal feasible basis if we apply lagrange (in case we don't use optimality-cuts), we choose the primal simplex 
					if (algorithmParameters.getNumberThreads().isPresent()) {
						if (algorithmParameters.getNumberThreads().get() <= 1) {
							//Primal simplex
							model.set(GRB.IntParam.Method, 0);
						}
						else if (algorithmParameters.getNumberThreads().get() <= 2) {
							//Concurrent Simplex
							model.set(GRB.IntParam.Method, 5);
						}
					}
				}
				//If lagrange relaxation is not applied, we need to set values for p, z, and the robustness constraints
				else {
					//If robustness constraints were already part of the previous model, then we construct a dual feasible basis
					if (node.getParentVBasisZ() != null) {
						//p and z are set as before
						z.set(GRB.IntAttr.VBasis, node.getParentVBasisZ());
						for (int i = 0; i < p.length; i++) {
							p[i].set(GRB.IntAttr.VBasis, node.getParentVBasisP()[i]);
						}
					}
					//Otherwise, we construct a basis
					else {
						//We set z to non-basic
						z.set(GRB.IntAttr.VBasis, -1);
						for (int i = 0; i < p.length; i++) {
							//There will at most be one robustness constraint in the model
							//We set the slack variable to non-basic and p basic
							if (robustnessConstraintsLowerBound[i] != null) {
								p[i].set(GRB.IntAttr.VBasis, 0);
								robustnessConstraintsLowerBound[i].set(GRB.IntAttr.CBasis, -1);
							}
							//If no constraint is in the model, we set p to non-basic
							else {
								p[i].set(GRB.IntAttr.VBasis, -1);
							}
						}

					}
				}
			}
		}
		
		this.solve(timeLimit);
		
		//Sets LP solve method to default
		model.set(GRB.IntParam.Method, -1);
		
		//Writes information to the node log.
		String output = "\n##### Finished LP"
				+ "\nDual bound = "+dualBound;
		AbstractAlgorithm.writeOutput(output, algorithmParameters);
	}
	
	/**
	 * Returns an integer array showing which nominal variable are basic in the computed solution.
	 */
	protected int[] getVBasisNominal () {
		try {
			return model.get(GRB.IntAttr.VBasis, nominalModelVariables);
		} catch (GRBException e) {
			return null;
		}

	}
	
	/**
	 * Returns an integer array showing which p variable are basic in the computed solution.
	 */
	protected int[] getVBasisP () {
		try {
			if (p == null) {
				return null;
			}
			else {
				return model.get(GRB.IntAttr.VBasis, p);
			}
		} catch (GRBException e) {
			return null;
		}

	}
	
	/**
	 * Returns an integer showing whether z is basic in the computed solution.
	 */
	protected Integer getVBasisZ () {
		try {
			if (z == null) {
				return null;
			}
			else {
				return z.get(GRB.IntAttr.VBasis);
			}
		} catch (GRBException e) {
			return null;
		}

	}
	
	/**
	 * Returns an integer array showing which slack variables of nominal constraints
	 * are basic in the computed solution.
	 */
	protected int[] getCBasisNominal () {
		try {
			return model.get(GRB.IntAttr.CBasis, nominalConstrs);
		} catch (GRBException e) {
			return null;
		}

	}
	
	/**
	* Returns the upper bounds with which the p variables have been substituted.
	* Can vary for different variables, as some need to be unchanged for warm starting.
	*/
	public double[] getUpperBoundsUsedForSubstituting() {
		return upperBoundsUsedForSubstituting;
	}
	
	/**
	 * Returns objects representing constraints whose slack variables are non-basic in the current solution.
	 */
	protected List<NonBasicSlackConstraint> getNonBasicSlackConstraints () {
		try {
			List<NonBasicSlackConstraint> nonBasicSlackConstraints = new ArrayList<>();
			if (robustnessConstraintsLowerBound != null) {
				for (int i = 0; i < robustnessConstraintsLowerBound.length; i++) {
					GRBConstr lbConstr = robustnessConstraintsLowerBound[i];
					if (lbConstr != null && lbConstr.get(GRB.IntAttr.CBasis) == -1) {
						nonBasicSlackConstraints.add(new NonBasicSlackConstraint(lbConstr, "lb", Optional.of(lowerBoundZ.getValue()), Optional.of(i)));
					}
					GRBConstr ubConstr = robustnessConstraintsUpperBound[i];
					if (ubConstr != null && ubConstr.get(GRB.IntAttr.CBasis) == -1) {
						nonBasicSlackConstraints.add(new NonBasicSlackConstraint(ubConstr, "ub", Optional.of(upperBoundZ.getValue()), Optional.of(i)));
					}
				}
			}

			if (parentNonBasicSlackConstraints != null) {
				for (int i = 0; i < parentNonBasicSlackConstraints.size(); i++) {
					if (parentNonBasicSlackGRBConstraints.get(i).get(GRB.IntAttr.CBasis) == -1) {
						nonBasicSlackConstraints.add(parentNonBasicSlackConstraints.get(i));
					}
				}
			}
			
			if (optimalityCuts != null) {
				for (GRBConstr optimalityCut : optimalityCuts) {
					if (optimalityCut.get(GRB.IntAttr.CBasis) == -1) {
						nonBasicSlackConstraints.add(new NonBasicSlackConstraint(optimalityCut));
					}
				}
			}
			
			return nonBasicSlackConstraints;
		} catch (GRBException e) {
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * Objects representing constraints whose slack are non-basic in the current solution.
	 * Can be used to add the constraint again into a child model.
	 * If this constraint contains a p variable, then this object holds this information.
	 * Important for verifying whether p can be substituted.
	 * 
	 * Also stores the typ of the constraint and the used bound if it is a robustness constraint.
	 * Important for not adding the same constraint twice in the child problem.
	 */
	public class NonBasicSlackConstraint {
		
		private GRBLinExpr lhsExpr;
		private char sense;
		private double rhs;
		private boolean isLowerBoundConstraint = false;
		private boolean isUpperBoundConstraint = false;
		private Optional<Double> bound = Optional.empty();
		private Optional<Integer> correspondingPIndex = Optional.empty();
		
		/**
		 * Constructor for optimality cuts
		 */
		NonBasicSlackConstraint(GRBConstr constr) throws GRBException {
			this.lhsExpr = model.getRow(constr);
			this.sense = constr.get(GRB.CharAttr.Sense);
			this.rhs = constr.get(GRB.DoubleAttr.RHS);
		}
		
		/**
		 * Constructor for robustness constraints
		 */
		NonBasicSlackConstraint(GRBConstr constr, String lbOrUB, Optional<Double> bound, Optional<Integer> correspondingPIndex) throws GRBException {
			this.lhsExpr = model.getRow(constr);
			this.sense = constr.get(GRB.CharAttr.Sense);
			this.rhs = constr.get(GRB.DoubleAttr.RHS);
			if (lbOrUB.toLowerCase().equals("lb")) {
				isLowerBoundConstraint = true;
			}
			if (lbOrUB.toLowerCase().equals("ub")) {
				isUpperBoundConstraint = true;
			}
			this.bound = bound;
			this.correspondingPIndex = correspondingPIndex;
		}
		
		public GRBLinExpr getLhsExpr() {
			return lhsExpr;
		}
		public char getSense() {
			return sense;
		}
		public double getRhs() {
			return rhs;
		}
		public boolean isLowerBoundConstraint() {
			return isLowerBoundConstraint;
		}
		public boolean isUpperBoundConstraint() {
			return isUpperBoundConstraint;
		}
		public Optional<Double> getBound() {
			return bound;
		}
		public Optional<Integer> getCorrespondingPIndex() {
			return correspondingPIndex;
		}
	}
	
	/**
	 * Solves the integer subproblem.
	 */
	void solveInteger(Optional<Double> timeLimit, BnBCallbackIntegerSubproblems bnbCallback) throws IOException, GRBException  {
		this.bnbCallback = bnbCallback;
		
		//Sets variables to their original type
		if (continuous) {
			for (Variable var : nominalVariables) {
				if (var.getType() == Type.BINARY) {
					var.getModelVariable().set(GRB.CharAttr.VType, GRB.BINARY);
				}
				else if (var.getType() == Type.INTEGER) {
					var.getModelVariable().set(GRB.CharAttr.VType, GRB.INTEGER);
				}
				else if (var.getType() == Type.CONTINUOUS) {
					var.getModelVariable().set(GRB.CharAttr.VType, GRB.CONTINUOUS);
				}
			}
			//Adds callback
			model.setCallback(new Callback());
			continuous = false;
		}
		
		this.solve(timeLimit);
		
		//Writes information to the node log
		String output = "\n#####Finished MILP\n"
				+ "Primal Bound: "+primalBound+"\n"
				+ "Dual Bound: "+dualBound;
		AbstractAlgorithm.writeOutput(output, algorithmParameters);
		
		//Updates the improved primal bound if necessary
		if (primalBound < improvedPrimalBound) {
			improvedPrimalBound = primalBound;
			try {
				improvedNominalVariablesSolutionValues = model.get(GRB.DoubleAttr.X, nominalModelVariables);
			} catch (GRBException e) { }
		}
		bnbCallback.updatePrimalDualBounds(improvedPrimalBound, dualBound);
	}
	
	/**
	 * Updates bounds on z and the reformulation for a given node of the branching tree.
	 */
	void updateBoundsAndFormulation(BnBNode node) throws GRBException, IOException {
		if (nominalConstrs == null) {
			nominalConstrs = model.getConstrs();
		}
		
		//Removes all temporary constraints added for the previous subproblem
		for (GRBConstr grbConstr : temporaryConstraints) {
			model.remove(grbConstr);
		}
		//Clears the list of temporary constraints and non-basic slack constraints from the parent model
		temporaryConstraints = new ArrayList<GRBConstr>();
		optimalityCuts = new ArrayList<GRBConstr>();
		parentNonBasicSlackGRBConstraints = new ArrayList<GRBConstr>();
		
		//Gets the non-basic slack constraints from the parent model
		parentNonBasicSlackConstraints = node.getParentNonBasicSlackConstraints();
		
		//Stores whether we apply Lagrange relaxation.
		applyLagrangeRelaxation = node.isApplyLagrangeRelaxation();
		
		//In case we do not apply lagrange relaxation and the robustness vars are not yet added, we do add them to the model
		if (!applyLagrangeRelaxation && z == null) {
			z = model.addVar(0, Double.MAX_VALUE, Gamma, GRB.CONTINUOUS, "z");
			
			if (cliquePartitioning != null) {
				p = new GRBVar[cliquePartitioning.getCliques().size()];
			}
			else {
				p = new GRBVar[uncertainVariables.length];
			}
			for (int i = 0; i < p.length; i++) {
				p[i] = model.addVar(0, Double.MAX_VALUE, 1, GRB.CONTINUOUS, "p"+i);
			}
		}
		
		//Updates bounds
		this.lowerBoundZ=node.getLowerBoundZ();
		this.upperBoundZ=node.getUpperBoundZ();
		
		//Sets bounds for z
		if (!applyLagrangeRelaxation) {
			z.set(GRB.DoubleAttr.LB, lowerBoundZ.getValue());
			z.set(GRB.DoubleAttr.UB, upperBoundZ.getValue());
		}
		
		//Indicates whether a lower/upper bound constraint with the same bound is already inherited from the parent problem.
		//In this case, we don't have to add it again.
		boolean[] isLowerBoundConstraintAlreadyAdded = null;
		boolean[] isUpperBoundConstraintAlreadyAdded = null;
		if (p != null) {
			isLowerBoundConstraintAlreadyAdded = new boolean[p.length];
			isUpperBoundConstraintAlreadyAdded = new boolean[p.length];
		}
		//Adds constraints whose slack were non-basic in the parent model
		//Only relevant if we don't solve integer and do warm start.
		if (!node.isSolveInteger()
				&& boundingStrategies.warmStartStrategy == BoundingLPWarmStartStrategy.WARMSTART_ENABLE
				&& parentNonBasicSlackConstraints != null) {
			for (int j = 0; j < parentNonBasicSlackConstraints.size(); j++) {
				NonBasicSlackConstraint parentNonBasicSlackConstraint = node.getParentNonBasicSlackConstraints().get(j);
				GRBConstr constr = model.addConstr(parentNonBasicSlackConstraint.getLhsExpr(), 
						parentNonBasicSlackConstraint.getSense(),
						parentNonBasicSlackConstraint.getRhs(), "");
				parentNonBasicSlackGRBConstraints.add(constr);
				temporaryConstraints.add(constr);
				
				//A robustness constraint is already added, if it is with respect to the same bound
				if (parentNonBasicSlackConstraint.isLowerBoundConstraint()
						&& parentNonBasicSlackConstraint.getBound().get() == lowerBoundZ.getValue()) {
					isLowerBoundConstraintAlreadyAdded[parentNonBasicSlackConstraint.getCorrespondingPIndex().get()] = true;
				}
				else if (parentNonBasicSlackConstraint.isUpperBoundConstraint()
						&& parentNonBasicSlackConstraint.getBound().get() == upperBoundZ.getValue()) {
					isUpperBoundConstraintAlreadyAdded[parentNonBasicSlackConstraint.getCorrespondingPIndex().get()] = true;
				}
			}
		}
		
		//Indices of uncertain variables that cannot be used for substituting p with the current upper bound,
		//as they are basic themselves or part of a non-basic slack constraint in the previous model.
		//Only relevant if we don't solve integer, do warm start, did not consider lagrange for the parent model.
		//If the latter case does not apply, then we are not given upper bounds for substituting in the parent model.
		isIndexNonSubstitutable = new boolean[uncertainVariables.length];
		if (!node.isSolveInteger()
				&& boundingStrategies.warmStartStrategy == BoundingLPWarmStartStrategy.WARMSTART_ENABLE
				&& node.getParentVBasisNominal() != null
				&& node.getParentUpperBoundsUsedForSubstituting() != null) {
			//Sets the uncertain variables within non-basic slack constraints to non-substitutable
			for (int j = 0; j < parentNonBasicSlackConstraints.size(); j++) {
				NonBasicSlackConstraint parentNonBasicSlackConstraints = node.getParentNonBasicSlackConstraints().get(j);
				if (parentNonBasicSlackConstraints.correspondingPIndex.isPresent()) {
					if (cliquePartitioning != null) {
						for (Integer varIndex : cliquePartitioning.getCliques().get(parentNonBasicSlackConstraints.correspondingPIndex.get())) {
							isIndexNonSubstitutable[varIndex] = true;
						}
					}
					else {
						isIndexNonSubstitutable[parentNonBasicSlackConstraints.correspondingPIndex.get()] = true;
					}
				}
			}
			
			//Adds all basic x variables.
			for (int i = 0; i < uncertainVariables.length; i++) {
				if (node.getParentVBasisNominal()[uncertainVariables[i].getNominalIndex()] != -1) {
					isIndexNonSubstitutable[i] = true;
				}
			}
		}
		
		//Alters the objective function.
		GRBLinExpr objLinExpr = new GRBLinExpr();
		//Adds nominal objective coeffs
		for (Variable nominalVariable : nominalVariables) {
			objLinExpr.addTerm(nominalVariable.getObjectiveCoefficient(), nominalVariable.getModelVariable());
		}
		
		//Adds the constant obtained from the lower bound if we apply Lagrange relaxation
		if (applyLagrangeRelaxation) {
			objLinExpr.addConstant(Gamma*lowerBoundZ.getValue());
		}
		//Adds the terms obtained from the upper bound on z for substitution
		upperBoundsUsedForSubstituting = new double[uncertainVariables.length];
		for (int i = 0; i < uncertainVariables.length; i++) {
			Variable uncertainVariable = uncertainVariables[i];
			//If the variable can be substituted, we use the current upper bound 
			if (!isIndexNonSubstitutable[i]) {
				upperBoundsUsedForSubstituting[i] = upperBoundZ.getValue();
			}
			//Otherwise, we take the last upper bound that has been used for substituting this variable
			else {
				upperBoundsUsedForSubstituting[i] = node.getParentUpperBoundsUsedForSubstituting()[i];
			}
			if (uncertainVariable.getDeviation() > upperBoundsUsedForSubstituting[i]) {
				objLinExpr.addTerm(uncertainVariable.getDeviation() - upperBoundsUsedForSubstituting[i], uncertainVariable.getModelVariable());
			}
		}
		
		//Adds p and z to the objective in case we do not apply Lagrange relaxation
		if (!applyLagrangeRelaxation) {
			for (int i = 0; i < p.length; i++) {
				objLinExpr.addTerm(1, p[i]);
			}
			objLinExpr.addTerm(Gamma, z);
		}
		//Otherwise applies Lagrange relaxation to the robustness constraints
		else {
			//Counts the number of robustness constraints that would be in the model
			int constraintCount = 0;
			if (cliquePartitioning != null) {
				for (List<Integer> clique : cliquePartitioning.getCliques()) {
					for (int varIndex : clique) {
						if (uncertainVariables[varIndex].getDeviation() > lowerBoundZ.getValue()) {
							constraintCount++;
							break;
						}
					}
				}
			}
			else {
				for (Variable uncertainVariable : uncertainVariables) {
					if (uncertainVariable.getDeviation() > lowerBoundZ.getValue()){
						constraintCount++;
					}
				}
			}
			//Relaxes each constraint with Lagrange multiplier = Gamma/constraintCount
			lagrangeMultipier = Math.min(1, 1.0*Gamma/constraintCount);
			for (Variable uncertainVariable : uncertainVariables) {
				if (uncertainVariable.getDeviation() > lowerBoundZ.getValue()) {
					objLinExpr.addTerm(lagrangeMultipier*(Math.min(uncertainVariable.getDeviation(), upperBoundZ.getValue())-lowerBoundZ.getValue()), uncertainVariable.getModelVariable());
				}
			}
		}
		model.setObjective(objLinExpr);
		
		//Adds the constraints from the (clique) reformulation
		if (!applyLagrangeRelaxation) {
			if (cliquePartitioning != null) {
				robustnessConstraintsLowerBound = new GRBConstr[cliquePartitioning.getCliques().size()];
				robustnessConstraintsUpperBound = new GRBConstr[cliquePartitioning.getCliques().size()];
				for (int i = 0; i < cliquePartitioning.getCliques().size(); i++) {
					//Only adds the constraint if it is not already inherited from the parent problem
					if (!isLowerBoundConstraintAlreadyAdded[i]) {
						//Tests whether the lower bound constraint should be added
						//We only add constraints if they contain a deviation that is above the lower bound
						boolean addConstraint = false;
						for (int varIndex : cliquePartitioning.getCliques().get(i)) {
							if (Math.min(uncertainVariables[varIndex].getDeviation(), upperBoundsUsedForSubstituting[varIndex]) - lowerBoundZ.getValue() > algorithmParameters.getFeasibilityTolerance()) {
								addConstraint = true;
								break;
							}
						}
						//Adds the lower bound constraint
						if (addConstraint) {
							GRBLinExpr robustnessExprLowerBound = new GRBLinExpr();
							robustnessExprLowerBound.addTerm(1, p[i]);
							robustnessExprLowerBound.addTerm(1, z);
							for (int varIndex : cliquePartitioning.getCliques().get(i)) {
								Variable var = uncertainVariables[varIndex];
								if (Math.min(var.getDeviation(), upperBoundsUsedForSubstituting[varIndex]) > lowerBoundZ.getValue()) {
									robustnessExprLowerBound.addTerm(-(Math.min(var.getDeviation(), upperBoundsUsedForSubstituting[varIndex])-lowerBoundZ.getValue()), var.getModelVariable());
								}
							}
							GRBConstr constrLowerBound = model.addConstr(robustnessExprLowerBound, GRB.GREATER_EQUAL, lowerBoundZ.getValue(), "");
							robustnessConstraintsLowerBound[i] = constrLowerBound;
							temporaryConstraints.add(constrLowerBound);
						}
					}
					//Only adds the constraint if it is not already inherited from the parent problem
					if (!isUpperBoundConstraintAlreadyAdded[i]) {
						//Tests whether the upper bound constraint should be added
						boolean addConstraint = false;
						for (int varIndex : cliquePartitioning.getCliques().get(i)) {
							if (Math.min(uncertainVariables[varIndex].getDeviation(), upperBoundsUsedForSubstituting[varIndex]) - upperBoundZ.getValue() > algorithmParameters.getFeasibilityTolerance()) {
								addConstraint = true;
								break;
							}
						}
						if (addConstraint) {
							GRBLinExpr robustnessExprUpperBound = new GRBLinExpr();
							robustnessExprUpperBound.addTerm(1, p[i]);
							for (int varIndex : cliquePartitioning.getCliques().get(i)) {
								Variable var = uncertainVariables[varIndex];
								if (Math.min(var.getDeviation(), upperBoundsUsedForSubstituting[varIndex]) > upperBoundZ.getValue()) {
									robustnessExprUpperBound.addTerm(-(Math.min(var.getDeviation(), upperBoundsUsedForSubstituting[varIndex])-upperBoundZ.getValue()), var.getModelVariable());
								}
							}
							GRBConstr constrUpperBound = model.addConstr(robustnessExprUpperBound, GRB.GREATER_EQUAL, 0, "");
							robustnessConstraintsUpperBound[i] = constrUpperBound;
							temporaryConstraints.add(constrUpperBound);
						}
					}
				}
			}
			else {
				robustnessConstraintsLowerBound = new GRBConstr[uncertainVariables.length];
				robustnessConstraintsUpperBound = new GRBConstr[uncertainVariables.length];
				for (int i = uncertainVariables.length-1; i >= 0; i--) {
					Variable var = uncertainVariables[i];
					//Tests whether the lower bound constraint should be added
					if (Math.min(var.getDeviation(), upperBoundsUsedForSubstituting[i]) - lowerBoundZ.getValue() > algorithmParameters.getFeasibilityTolerance()) {
						GRBLinExpr robustnessExprLowerBound = new GRBLinExpr();
						robustnessExprLowerBound.addTerm(1, p[i]);
						robustnessExprLowerBound.addTerm(1, z);
						robustnessExprLowerBound.addTerm(-(Math.min(var.getDeviation(), upperBoundsUsedForSubstituting[i])-lowerBoundZ.getValue()), var.getModelVariable());
						GRBConstr constrLowerBound;
						constrLowerBound = model.addConstr(robustnessExprLowerBound, GRB.GREATER_EQUAL, lowerBoundZ.getValue(), "");
						robustnessConstraintsLowerBound[i] = constrLowerBound;
						temporaryConstraints.add(constrLowerBound);
					}
					
					//If we don't apply substitution to a variable and their deviation is high, then we add an upper bound constraint
					if (isIndexNonSubstitutable[i]) {
						//Tests whether the upper bound constraint should be added
						if (Math.min(var.getDeviation(), upperBoundsUsedForSubstituting[i]) - upperBoundZ.getValue() > algorithmParameters.getFeasibilityTolerance()) {
							GRBLinExpr robustnessExprUpperBound = new GRBLinExpr();
							robustnessExprUpperBound.addTerm(1, p[i]);
							robustnessExprUpperBound.addTerm(-(Math.min(var.getDeviation(), upperBoundsUsedForSubstituting[i])-upperBoundZ.getValue()), var.getModelVariable());
							GRBConstr constrUpperBound = model.addConstr(robustnessExprUpperBound, GRB.GREATER_EQUAL, 0, "");
							robustnessConstraintsUpperBound[i] = constrUpperBound;
							temporaryConstraints.add(constrUpperBound);
						}
					}
				}
			}
		}
		//No robustness constraints if we apply lagrange
		else {
			robustnessConstraintsLowerBound = null;
			robustnessConstraintsUpperBound = null;
		}
		model.update();
	}
	
	/**
	 * Obtains an incumbent solution, computes the optimal z and improves the primal bound if
	 * the new solution is the new best.
	 */
	private void improveZ (double[] incumbentNominalVariablesValues, double[] incumbentUncertainVariablesValues) {
		//Computes optimal value for z
		double optimalZ = computeOptimalBilinearZ(incumbentUncertainVariablesValues);
		//Computes the optimal objective value
		double objectiveValue = computeBilinearSolutionValue(incumbentNominalVariablesValues, incumbentUncertainVariablesValues, optimalZ);
		//Updates the improved primal bound if the new incumbent is better
		if (objectiveValue < improvedPrimalBound) {
			improvedPrimalBound = objectiveValue;
			improvedZ = Optional.of(optimalZ);
			improvedNominalVariablesSolutionValues = incumbentNominalVariablesValues;
		}
	}
	
	/**
	 * Implements the callbacks for termination and improving of incumbent solutions.
	 */
	private class Callback extends GRBCallback {
				
		@Override
		protected void callback() {
			//Stores incumbent solutions and potentially tries to improve an incumbent solution computing an optimal value for z.
			try {
				if (where == GRB.CB_MIPSOL) {
					if (boundingStrategies.getImprovingZStrategy() == ImprovingZStrategy.IMPROVINGZ_ENABLE) {
						double[] incumbentNominalVariablesValues = getSolution(nominalModelVariables);
						double[] incumbentUncertainVariablesValues = getSolution(uncertainModelVariables);
						improveZ(incumbentNominalVariablesValues, incumbentUncertainVariablesValues);
					}
					else {
						if (getDoubleInfo(GRB.CB_MIPSOL_OBJ) < improvedPrimalBound) {
							improvedPrimalBound = getDoubleInfo(GRB.CB_MIPSOL_OBJ);
							improvedZ = Optional.of(getSolution(z));
							improvedNominalVariablesSolutionValues = getSolution(nominalModelVariables);
						}
					}
				}
			} catch (GRBException e) {
				e.printStackTrace();
			}

			
			//Reports current primal and dual bounds to the master and possibly asks for termination
			try {
				if (where == GRB.CB_MIP) {
					dualBound = getDoubleInfo(GRB.CB_MIP_OBJBND);
					primalBound = getDoubleInfo(GRB.CB_MIP_OBJBST);
					
					if (boundingStrategies.getTerminationStrategy() == BoundingTerminationStrategy.TERMINATION_ENABLE) {
						if (bnbCallback.updateBoundsAndDecideTermination(improvedPrimalBound, primalBound, dualBound)) {
							abort();
						}
					}
					else {
						bnbCallback.updatePrimalDualBounds(improvedPrimalBound, dualBound);
					}
				}
			} catch (GRBException e1) {
				e1.printStackTrace();
			}

		}
	}
	
	/**
	 * Solves the model using Gurobi.
	 */
	@Override
	void solve(Optional<Double> timeLimit) throws IOException, GRBException  {
		super.solve(timeLimit);
		//Tries to query the solution values for z.
		//Value is null if this fails.
		try {
			if (!applyLagrangeRelaxation) {
				zSolutionValue = Optional.of(z.get(GRB.DoubleAttr.X));
				pSolutionValues = model.get(GRB.DoubleAttr.X, p);
				
				//Potentially resubstitutes the variables p.
				if (cliquePartitioning != null) {
					for (int i = 0; i < p.length; i++) {
						for (int varIndex : cliquePartitioning.getCliques().get(i)) {
							pSolutionValues[i] += Math.max(0,(uncertainVariables[varIndex].getDeviation() - upperBoundsUsedForSubstituting[varIndex]) * uncertainVariablesSolutionValues[varIndex]);
						}
					}
				}
				else {
					for (int i = 0; i < p.length; i++) {
						pSolutionValues[i] += Math.max(0,(uncertainVariables[i].getDeviation() - upperBoundsUsedForSubstituting[i]) * uncertainVariablesSolutionValues[i]);
					}
				}
			}
		} catch (GRBException e) { }
	}
	

	/**
	 * Computes the lowest possible value of z in the current node that would be an efficient branching point.
	 */
	public double getLowerBoundEffectiveBranching(BnBNode chosenNode) {
		if (applyLagrangeRelaxation) {
			//If we apply Lagrangean relaxation and theta is the branching point,
			//then the sum of x*Lagrange for variables whose deviation is greater or equal than theta must be less than Gamma.
			
			//The sum for all x*Lagrange whose deviation is greater or equal than the minimum value for z
			double sum = 0;
			//The largest index of a variable that has not been counted in the sum
			int largestUncountedIndex= -1;
			for (int varIndex = uncertainVariables.length-1; varIndex >= 0; varIndex--) {
				if (uncertainVariables[varIndex].getDeviation() < lowerBoundZ.getValue()) {
					largestUncountedIndex = varIndex;
					break;
				}
				sum += uncertainVariablesSolutionValues[varIndex]*lagrangeMultipier;
			}
			
			//Checks for all possible z
			for (int zIndex = 0; zIndex < chosenNode.getPossibleZs().size()-1 ; zIndex++) {
				double possibleLowerBound = chosenNode.getPossibleZs().get(zIndex).getValue();
				//Reduces the sum and sets the new uncounted index
				for (int varIndex = largestUncountedIndex+1; varIndex < uncertainVariables.length; varIndex++) {
					if (uncertainVariables[varIndex].getDeviation() >= possibleLowerBound) {
						largestUncountedIndex = varIndex-1;
						break;
					}
					sum -= uncertainVariablesSolutionValues[varIndex]*lagrangeMultipier;
				}
				//Checks whether sum is smaller than Gamma
				if (sum < Gamma) {
					return possibleLowerBound;
				}
			}
			//The second largest possible z is always effective
			return chosenNode.getPossibleZs().get(chosenNode.getPossibleZs().size()-2).getValue();

		}
		else {
			//If we don't apply Lagrangean relaxation, then we compute the minimum branching point
			//such that one of the robustness constraints is violated.
			
			//The second largest possible z is always effective
			double lowerBoundEffectiveBranching = chosenNode.getPossibleZs().get(chosenNode.getPossibleZs().size()-2).getValue();
			//Checks all robustness constraints
			for (int pIndex = 0; pIndex < p.length; pIndex++) {
				//Breaks if the bound is minimum
				if (lowerBoundEffectiveBranching == lowerBoundZ.getValue()) {
					break;
				}
				//Computes the minimum effective bound for this constraint
				double lhs = zSolutionValue.get() + pSolutionValues[pIndex];
				for (int k = 0; k < chosenNode.getPossibleZs().size()-1; k++) {
					double possibleLowerBound = chosenNode.getPossibleZs().get(k).getValue();
					//Breaks if the currently tested value is not smaller than an already found bound.
					if (possibleLowerBound >= lowerBoundEffectiveBranching) {
						break;
					}
					double rhs = possibleLowerBound;
					if (cliquePartitioning != null) {
						for (int varIndex : cliquePartitioning.getCliques().get(pIndex)) {
							Variable var = uncertainVariables[varIndex];
							if (uncertainVariablesSolutionValues[varIndex] > 0 && var.getDeviation() > possibleLowerBound) {
								rhs += (var.getDeviation() - possibleLowerBound) * uncertainVariablesSolutionValues[varIndex];
							}
						}
					}
					else {
						Variable var = uncertainVariables[pIndex];
						if (uncertainVariablesSolutionValues[pIndex] > 0 && var.getDeviation() > possibleLowerBound) {
							rhs += (var.getDeviation() - possibleLowerBound) * uncertainVariablesSolutionValues[pIndex];
						}
					}
					//Updates the bound if we found a violation
					if (lhs < rhs - algorithmParameters.getFeasibilityTolerance()) {
						lowerBoundEffectiveBranching = possibleLowerBound;
					}
				}
			}
			return lowerBoundEffectiveBranching;
		}
	}
	
	
	/**
	 * Computes the highest possible value of z in the current node that would be an efficient branching point.
	 */
	public double getUpperBoundEffectiveBranching(BnBNode chosenNode) {
		if (applyLagrangeRelaxation) {
			//If we apply Lagrangean relaxation and theta is the branching point, then there must be a
			//variable with deviation > theta, positive solution value, and Lagrange multiplier < 1.
			//Since our choice of Lagrange multipliers is never equal to one (otherwise we would stop branching),
			//we search for the variable with largest deviation and positive solution value.
			double largestPosSolutionDeviation = 0;
			for (int varIndex = uncertainVariables.length-1; varIndex >=0 ; varIndex--) {
				if (uncertainVariablesSolutionValues[varIndex] > 0) {
					largestPosSolutionDeviation = uncertainVariables[varIndex].getDeviation();
					break;
				}
			}
			//Finds the largest possible branching point that is smaller than the deviation.
			for (int zIndex = chosenNode.getPossibleZs().size()-2; zIndex >=0 ; zIndex--) {
				if (chosenNode.getPossibleZs().get(zIndex).getValue() < largestPosSolutionDeviation) {
					return chosenNode.getPossibleZs().get(zIndex).getValue();
				}
			}
			//The smallest possible z is always effective.
			return lowerBoundZ.getValue();
		}
		else {
			//If we don't apply Lagrangean relaxation, then we compute the maximum branching point
			//such that one of the robustness constraints is violated.
			
			//The smallest possible z is always effective.
			double upperBoundEffectiveBranching = lowerBoundZ.getValue();
			//Checks all robustness constraints
			for (int pIndex = 0; pIndex < p.length; pIndex++) {
				//Breaks if the bound is maximum
				if (upperBoundEffectiveBranching == chosenNode.getPossibleZs().get(chosenNode.getPossibleZs().size()-2).getValue()) {
					break;
				}
				//Computes the maximum effective bound for this constraint
				double lhs = pSolutionValues[pIndex];
				for (int zIndex = chosenNode.getPossibleZs().size()-2; zIndex >= 0; zIndex--) {
					double possibleUpperBound = chosenNode.getPossibleZs().get(zIndex).getValue();
					//Breaks if we the currently tested value is not greater than an already found bound.
					if (possibleUpperBound <= upperBoundEffectiveBranching) {
						break;
					}
					double rhs = 0;
					if (cliquePartitioning != null) {
						for (int varIndex : cliquePartitioning.getCliques().get(pIndex)) {
							Variable var = uncertainVariables[varIndex];
							if (uncertainVariablesSolutionValues[varIndex] > 0 && var.getDeviation() > possibleUpperBound) {
								rhs += (var.getDeviation() - possibleUpperBound) * uncertainVariablesSolutionValues[varIndex];
							}
						}
					}
					else {
						Variable var = uncertainVariables[pIndex];
						if (uncertainVariablesSolutionValues[pIndex] > 0 && var.getDeviation() > possibleUpperBound) {
							rhs += (var.getDeviation() - possibleUpperBound) * uncertainVariablesSolutionValues[pIndex];
						}
					}
					//Updates the bound if we found a violation
					if (lhs < rhs - algorithmParameters.getFeasibilityTolerance()) {
						upperBoundEffectiveBranching = possibleUpperBound;
					}
				}
			}
			return upperBoundEffectiveBranching;
		}
	}
	
	/**
	 * Checks whether the current solution is integer feasible.
	 */
	boolean isSolutionIntegerFeasible() {
		boolean feasible = true;
		for (int i = 0; i < nominalVariables.length; i++) {
			Variable variable = nominalVariables[i];
			if (variable.getType() == Type.BINARY || variable.getType() == Type.INTEGER) {
				if (Math.abs(nominalVariablesSolutionValues[i] - Math.round(nominalVariablesSolutionValues[i])) > algorithmParameters.getIntegerFeasibilityTolerance()) {
					feasible = false;
					break;
				}
			}
		}
		return feasible;
	}
	
	/**
	 * Returns the value for z computed by Gurobi.
	 */
	Optional<Double> getComputedZ() {
		return zSolutionValue;
	}
	
	/**
	 * Returns the improved primal bound.
	 */
	double getImprovedPrimalBound() {
		return this.improvedPrimalBound;
	}
	
	/**
	 * Returns the improved value for z.
	 */
	Optional<Double> getImprovedZ() {
		return this.improvedZ;
	}
	
	/**
	 * Returns solution values of the incumbent.
	 */
	double[] getIncumbentValues() {
		return this.improvedNominalVariablesSolutionValues;
	}
	
	/**
	 * Specifies strategies for algorithms solving bounded subproblems.
	 */
	abstract static class BoundingSubproblemsStrategies extends SubproblemsStrategies{
		/**
		 * Enum type specifying whether we terminate robust subproblems prematurely.
		 */
		public enum BoundingTerminationStrategy {
	 		TERMINATION_ENABLE,
			TERMINATION_DISABLE;
		}
		
		/**
		 * Enum type specifying whether we use optimality-cuts.
		 */
		public enum BoundingLPOptimalityCutsStrategy {
			LPOPTCUTS_ENABLE,
			LPOPTCUTS_DISABLE;
		}
		
		/**
		 * Enum type specifying whether we use optimality-cuts.
		 */
		public enum BoundingIPOptimalityCutsStrategy {
			IPOPTCUTS_ENABLE,
			IPOPTCUTS_DISABLE;
		}
		
		/**
		 * Enum type specifying whether we use optimality-cuts.
		 */
		public enum BoundingLPLagrangeRelaxStrategy {
			LAGRANGERELAX_NEVER,
			LAGRANGERELAX_HYBRID,
			LAGRANGERELAX_ALWAYS;
		}
		
		/**
		 * Enum type specifying whether we use optimality-cuts.
		 */
		public enum BoundingLPWarmStartStrategy {
			WARMSTART_ENABLE,
			WARMSTART_DISABLE;
		}

		protected BoundingTerminationStrategy terminationStrategy;
		protected BoundingLPOptimalityCutsStrategy lpOptimalityCutsStrategy;
		protected BoundingIPOptimalityCutsStrategy ipOptimalityCutsStrategy;
		protected BoundingLPLagrangeRelaxStrategy lagrangeRelaxStrategy;
		protected BoundingLPWarmStartStrategy warmStartStrategy;

		/**
		 * Constructor obtaining arguments which are matched to the enums defining strategies.
		 */
		public BoundingSubproblemsStrategies(List<String> argList, AlgorithmParameters algorithmParameters) throws IOException {
			super(argList, algorithmParameters);
		}

		public BoundingTerminationStrategy getTerminationStrategy() {
			return terminationStrategy;
		}
		public void setTerminationStrategy(BoundingTerminationStrategy terminationStrategy) {
			this.terminationStrategy = terminationStrategy;
		}
		
		public BoundingLPOptimalityCutsStrategy getLpOptimalityCutsStrategy() {
			return lpOptimalityCutsStrategy;
		}
		public void setLpOptimalityCutsStrategy(BoundingLPOptimalityCutsStrategy lpOptimalityCutsStrategy) {
			this.lpOptimalityCutsStrategy = lpOptimalityCutsStrategy;
		}
		
		public BoundingIPOptimalityCutsStrategy getIPOptimalityCutsStrategy() {
			return ipOptimalityCutsStrategy;
		}
		public void setIPOptimalityCutsStrategy(BoundingIPOptimalityCutsStrategy ipOptimalityCutsStrategy) {
			this.ipOptimalityCutsStrategy = ipOptimalityCutsStrategy;
		}

		public BoundingLPLagrangeRelaxStrategy getLagrangeRelaxLPStrategy() {
			return lagrangeRelaxStrategy;
		}
		public void setLagrangeRelaxLPStrategy(BoundingLPLagrangeRelaxStrategy lagrangeRelaxStrategy) {
			this.lagrangeRelaxStrategy = lagrangeRelaxStrategy;
		}
		
		public BoundingLPWarmStartStrategy getWarmStartStrategy() {
			return warmStartStrategy;
		}
		public void setWarmStartStrategy(BoundingLPWarmStartStrategy warmStartStrategy) {
			this.warmStartStrategy = warmStartStrategy;
		}
	}
}
