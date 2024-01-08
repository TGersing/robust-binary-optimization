package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import alg.SubproblemBoundedGurobi.BoundingSubproblemsStrategies;
import alg.SubproblemBoundedGurobi.BoundingSubproblemsStrategies.BoundingIPOptimalityCutsStrategy;
import alg.SubproblemBoundedGurobi.BoundingSubproblemsStrategies.BoundingLPLagrangeRelaxStrategy;
import alg.SubproblemBoundedGurobi.BoundingSubproblemsStrategies.BoundingLPOptimalityCutsStrategy;
import alg.SubproblemBoundedGurobi.BoundingSubproblemsStrategies.BoundingLPWarmStartStrategy;
import alg.SubproblemBoundedGurobi.NonBasicSlackConstraint;
import gurobi.GRBException;
import util.BnBNode;
import util.CliquePartitioning;
import util.ConflictGraph;
import util.PossibleZ;
import util.Variable;

/**
 * This class implements the branch and bound algorithm as described
 * in the paper "A Branch and Bound Algorithm for Robust Binary Optimization with Budget Uncertainty".
 * 
 * @author Timo Gersing
 */
public class AlgBranchAndBound extends AbstractAlgorithm implements RobustAlgorithm {
	/**
	 * The robust subproblem that we adapt for different sets Z containing possible values for z.
	 */
	private SubproblemBoundedGurobi subproblemBounded;
	
	/**
	 * Specifies strategies.
	 */
	private BnBStrategies bnbStrategies;

	/**
	 * A partitioning of the uncertain variables into cliques.
	 */
	private CliquePartitioning cliquePartitioning;
	
	/**
	 * A list containing possible optimal values for z.
	 */
	private List<PossibleZ> possibleZs;
	
	/**
	 * A set of tree nodes containing sets of yet to be considered values of z.
	 */
	private Set<BnBNode> remainingNodes = new HashSet<BnBNode>();
	
	/**
	 * A set containing nodes whose possible values for z were already considered for an integer robust subproblem.
	 * Sorted with respect to the contained possible values for z.
	 * Important for deciding which bounds to use for the optimality-cuts.
	 */
	private TreeSet<BnBNode> alreadyConsideredNodesforRobustSubproblems = new TreeSet<BnBNode>((node1, node2) -> Double.compare(node1.getLowerBoundZ().getValue(), node2.getLowerBoundZ().getValue()));
	
	/**
	 * Tracks the number of yet to be considered values of z.
	 */
	private int numberRemainingPossibleZ;
	
	/**
	 * Counts the number of integer robust subproblems already considered.
	 */
	private int numberIntegerSubproblemsStarted = 0;
	
	/**
	 * Counts the cumulative time (in seconds) elapsed while solving integer robust subproblems.
	 */
	private double timeInIntegerSubproblems = 0;

	/**
	 * Counts the number of linear relaxations already considered.
	 */
	private int numberLinearSubproblemsStarted = 0;
	
	/**
	 * Counts the cumulative time (in seconds) elapsed while solving linear relaxations.
	 */
	private double timeInLinearizedSubproblems = 0;
	
	/**
	 * Stores the value of the primal bound after the first integer robust subproblem.
	 * Used for evaluating the quality of the node choice. 
	 */
	private Double primalBoundAfterFirstInteger = null;
	
	/**
	 * 	In order to report a general dual bound, we consider the minimum dual bound of all leaves in the branching tree.
	 *  Since we will prune nodes over the course of the algorithm, we store their minimum dual bound.
	 */
	private double minimumDualBoundPrunedNodes = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
	
	/**
	 * Stores values of the current best solution.
	 */
	double[] solutionValues;
	
	/**
	 * Constructor receiving problem paths as well as parameters and strategies.
	 */
	public AlgBranchAndBound(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, BnBStrategies bnbStrategies) throws IOException, GRBException {
		super(algorithmParameters);
		this.subproblemBounded = new SubproblemBoundedGurobi(problemPath, robustPath, algorithmParameters, bnbStrategies);
		this.bnbStrategies = bnbStrategies;
	}
	
	/**
	 * Executes the branch and bound algorithm.
	 */
	@Override
	protected void executeAlgorithm() throws IOException, GRBException {
		String output = "######################################################\n"
				+ "##### Solving Problem via Branch and Bound"
				+ "\n######################################################";
		writeOutput(output);
		
		output = "\n###########################\n"
				+ "##### Start Preprocessing"
				+ "\n###########################";
		writeOutput(output);
		
		//Computes a conflict graph and clique partitioning if the corresponding strategy is chosen.
		ConflictGraph conflictGraph = null;
		if (bnbStrategies.getCliqueStrategy() == RobustAlgorithmStrategies.CliqueStrategy.CLIQUES_ENABLE) {
			conflictGraph = new ConflictGraph(subproblemBounded.getModel(), subproblemBounded.getUncertainModelVariables(), algorithmParameters);
			cliquePartitioning = new CliquePartitioning(subproblemBounded.getUncertainVariables(), conflictGraph, algorithmParameters);
			subproblemBounded.setCliquePartitioning(cliquePartitioning);
		}
		
		//Computes a list of possible optimal choices for z respecting the chosen filtering and clique strategies.
		possibleZs = computePossibleZs(subproblemBounded.getUncertainVariables(), subproblemBounded.getGamma(), bnbStrategies, conflictGraph, algorithmParameters);
		
		//Decides whether the first LP should be solved with a Lagrange relaxation or not
		boolean applyLagrangeRelaxation = (bnbStrategies.getLagrangeRelaxLPStrategy() == BoundingLPLagrangeRelaxStrategy.LAGRANGERELAX_ALWAYS 
				|| bnbStrategies.getLagrangeRelaxLPStrategy() == BoundingLPLagrangeRelaxStrategy.LAGRANGERELAX_HYBRID);
		
		//Initialize the root node containing all possible values for z.
		remainingNodes.add(new BnBNode(possibleZs, applyLagrangeRelaxation));
		numberRemainingPossibleZ = possibleZs.size();
		
		output = "\n###########################\n"
				+ "##### Start Branching"
				+ "\n###########################";
		writeOutput(output);
		
		//While there exist active nodes in the branching tree.
		while (!remainingNodes.isEmpty()) {
			//Terminates the algorithm if the time limit is reached.
			if (deadline.isPresent() && System.nanoTime() > deadline.get()) {
				output = "\n###########################\n"
						+ "##### Reached Time Limit"
						+ "\n###########################";
				writeOutput(output);
				break;
			}
			
			//Chooses the node with the lowest dual bound.
			BnBNode currentNode = Collections.min(remainingNodes);
			setDualBound(Math.min(minimumDualBoundPrunedNodes, currentNode.getDualBound()));
			primalDualIntegral.update(true);
			
			//Writes information of the current state of the BnB
			output = "\n##### Branching Tree Information\n"
					+ "Elapsed Time = "+(System.nanoTime()-startTime)/Math.pow(10, 9)+"\n"
					+ "Number Remaining Nodes = "+remainingNodes.size()
					+ "\nNumber Remaining Possible Z = "+numberRemainingPossibleZ
					+ "\nCurrent Primal Bound = "+getPrimalBound()
					+ "\nCurrent Dual Bound = "+getDualBound()
					+ "\nRelative Optimality Gap = "+(getRelativeGap()*100)+"%"
					+ "\nCurrent Primal-Dual Integral = "+getPrimalDualIntegral();
			writeOutput(output);
			
			//Removes the chosen node from the set of remaining nodes to consider.
			remainingNodes.remove(currentNode);
			
			//Solves the integer robust subproblem if the current node is marked correspondingly.
			if (currentNode.isSolveInteger()) {
				numberIntegerSubproblemsStarted++;
				solveInteger(currentNode);
			}
			//Solves the linear relaxation otherwise. Either prunes or branches the node or marks it to be solved as an integer subproblem.
			else {
				numberLinearSubproblemsStarted++;
				solveLinear(currentNode);
			}

			//Updates the performance indicator if we did not yet considered the second integer robust subproblem.
			if (numberIntegerSubproblemsStarted <= 1) {
				primalBoundAfterFirstInteger = getPrimalBound();
			}
			
			//Prunes nodes if possible
			pruneNodes();
			//Prunes single possible values for z within nodes if possible
			prunePossibleZinNodes();
		}
		
		//Updates the dual bound to be the minimum of all pruned nodes and
		//the minimum dual bound of the not yet considered nodes.
		updateGlobalDualBound();
		
		//Stores best solution found
		if (solutionValues != null) {
			solution = new LinkedHashMap<Variable, Double>(solutionValues.length);
			for (int i = 0; i < solutionValues.length; i++) {
				solution.put(subproblemBounded.getNominalVariables()[i], solutionValues[i]);
			}
		}
	}
	
	/**
	 * Solves the integer robust subproblem for a given node.
	 */
	private void solveInteger(BnBNode chosenNode) throws IOException, GRBException {
		String output = "\n##### Solving MILP for "+chosenNode;
		writeOutput(output);

		//Updates the robust subproblem for the chosen node.
		subproblemBounded.updateBoundsAndFormulation(chosenNode);
		
		//Set of nodes, sorted ascending with respect to the contained possible values of z, that is [0,5] < [10,15].
		//Will contain all nodes for which we might apply estimators.
		TreeSet<BnBNode> surroundingNodes = new TreeSet<BnBNode>((node1, node2) -> Double.compare(node1.getLowerBoundZ().getValue(), node2.getLowerBoundZ().getValue()));
		
		//Determines estimators and optimality cuts if the corresponding options are chosen.
		if (bnbStrategies.getEstimatorStrategy() != BnBStrategies.BnBEstimatorStrategy.ESTIMATORS_NONE
				|| bnbStrategies.getIPOptimalityCutsStrategy() == BoundingIPOptimalityCutsStrategy.IPOPTCUTS_ENABLE) {
			
			//Initializes surrounding nodes with all remaining nodes.
			for (BnBNode bnBNode : remainingNodes) {
				surroundingNodes.add(bnBNode);
			}
			
			//The bounds for the optimality-cuts are given by the smallest and largest remaining values
			//around or within the chosen node such that there exists no already considered node in between.
			if (bnbStrategies.getIPOptimalityCutsStrategy() == BoundingIPOptimalityCutsStrategy.IPOPTCUTS_ENABLE) {
				//If we set optimality-cuts then we do not compute estimators for nodes outside of the corresponding bounds.				
				//Removes all nodes that are smaller than an already considered node that is smaller than the chosen node.
				BnBNode greatestSmallerConsideredNode = alreadyConsideredNodesforRobustSubproblems.floor(chosenNode);
				if (greatestSmallerConsideredNode != null) {
					surroundingNodes.headSet(greatestSmallerConsideredNode).clear();
				}
				//Removes all nodes that are greater than an already considered node that is larger than the chosen node.
				BnBNode leastGreaterConsideredNode = alreadyConsideredNodesforRobustSubproblems.ceiling(chosenNode);
				if (leastGreaterConsideredNode != null) {
					surroundingNodes.tailSet(leastGreaterConsideredNode).clear();
				}
				
				//Initializes the bounds for the optimality-cuts and adds these to the problem.
				TreeSet<BnBNode> optCutNodes = new TreeSet<BnBNode>(surroundingNodes);
				optCutNodes.add(chosenNode);
				PossibleZ lowerBoundOptCut = optCutNodes.first().getLowerBoundZ();
				PossibleZ upperBoundOptCut = optCutNodes.last().getUpperBoundZ();
				subproblemBounded.addOptimalityCuts(lowerBoundOptCut, upperBoundOptCut);
			}
			
			//Computes estimators according to the chosen option.
			if (bnbStrategies.getEstimatorStrategy() != BnBStrategies.BnBEstimatorStrategy.ESTIMATORS_NONE) {
				TreeSet<PossibleZ> comparedPossibleZ = new TreeSet<PossibleZ>();
				for (BnBNode node : surroundingNodes) {
					comparedPossibleZ.addAll(node.getPossibleZs());
				}
				if (bnbStrategies.getEstimatorStrategy() == BnBStrategies.BnBEstimatorStrategy.ESTIMATORS_IMPROVED) {
					//If the node consists of only one possible value for z then this value defines the estimators for all lower and higher values.
					if (chosenNode.getLowerBoundZ() == chosenNode.getUpperBoundZ()) {
						chosenNode.getUpperBoundZ().setImprovedEstimators(comparedPossibleZ, subproblemBounded.getGamma(), subproblemBounded.getUncertainVariables());
					}
					//Otherwise the estimators are defined by the lowest and highest possible z.
					else {
						chosenNode.getLowerBoundZ().setImprovedEstimatorsBelowZ(comparedPossibleZ, subproblemBounded.getGamma(), subproblemBounded.getUncertainVariables());
						chosenNode.getUpperBoundZ().setImprovedEstimatorsAboveZ(comparedPossibleZ, subproblemBounded.getUncertainVariables());
					}
				}
				//If we use the estimators by Hansknecht, Richter, Stiller then these are defined by the lowest value in the node.
				else if (bnbStrategies.getEstimatorStrategy() == BnBStrategies.BnBEstimatorStrategy.ESTIMATORS_HRS) {
					chosenNode.getLowerBoundZ().setHRSEstimators(comparedPossibleZ, subproblemBounded.getGamma());
				}
			}
		}
		writeOutput("");
		
		//Initializes the callback used for reporting primal and dual bounds from the subproblem to the master.
		BnBCallbackIntegerSubproblems bnbCallback = new BnBCallbackIntegerSubproblems(chosenNode, surroundingNodes);
		//Tries to solve the robust subproblem in the time limit.
		long starttimeInteger = System.nanoTime();
		subproblemBounded.solveInteger(getRemainingTime(), bnbCallback);
		timeInIntegerSubproblems+=(System.nanoTime()-starttimeInteger)/Math.pow(10, 9);
		
		//Adds the node to the set of nodes that were already considered for subproblems.
		alreadyConsideredNodesforRobustSubproblems.add(chosenNode);
		
		//Updates the number of remaining z as well as the minimum dual bound of pruned nodes.
		numberRemainingPossibleZ -= chosenNode.getPossibleZs().size();
		minimumDualBoundPrunedNodes = Math.min(minimumDualBoundPrunedNodes, chosenNode.getDualBound());
	}
	
	/**
	 * Callback passed to the integer subproblems.
	 * Reports bounds from the integer subproblem to the master and asks for termination.
	 */
	protected class BnBCallbackIntegerSubproblems {
		/**
		 * The node for which we solve the integer subproblem.
		 */
		private BnBNode currentNode;
		
		/**
		 * Set of nodes around the current node for which we compute estimators.
		 */
		private TreeSet<BnBNode> surroundingNodes;

		/**
		 * Saves the last time stamp at which we were able to prune a possible value of z via estimators. 
		 */
		private Long timeStampLastPrune = null;
		
		/**
		 * Maps the bounds of the nodes to not yet pruned values of z for which we have computed a estimators from the bounds.
		 */
		private Map<PossibleZ, List<PossibleZ>> boundToRemainingPossibleZwithEstimator;

		/**
		 * Constructor receiving the chosen node and its surrounding nodes.
		 */
		private BnBCallbackIntegerSubproblems(BnBNode currentNode, TreeSet<BnBNode> surroundingNodes) {
			boundToRemainingPossibleZwithEstimator = new HashMap<PossibleZ, List<PossibleZ>>();
			boundToRemainingPossibleZwithEstimator.put(currentNode.getUpperBoundZ(), new ArrayList<PossibleZ>(currentNode.getUpperBoundZ().getEstimators().keySet()));
			if (currentNode.getUpperBoundZ() != currentNode.getLowerBoundZ()) {
				boundToRemainingPossibleZwithEstimator.put(currentNode.getLowerBoundZ(), new ArrayList<PossibleZ>(currentNode.getLowerBoundZ().getEstimators().keySet()));
			}
			this.currentNode = currentNode;
			this.surroundingNodes = surroundingNodes;
		}
		
		/**
		 * Receives primal and dual bounds from the subproblem to update the bounds in the master problem.
		 */
		protected boolean updatePrimalDualBounds(double primalBound, double subproblemDualBound) {
			boolean boundChanged = false;
			
			//Updates the primal bound
			if (primalBound < getPrimalBound()) {
				setPrimalBound(primalBound);
				solutionValues = subproblemBounded.getIncumbentValues();
				boundChanged = true;
				
				try {
					writeOutput("###Found new global primal bound = "+getPrimalBound());
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			//Updates dual bounds of possible z contained in this node and potentially surrounding nodes via estimators.
			if (currentNode.getDualBound() < subproblemDualBound) {
				currentNode.updateDualBound(subproblemDualBound);
				boundChanged = true;
				
				if (bnbStrategies.getEstimatorStrategy() == BnBStrategies.BnBEstimatorStrategy.ESTIMATORS_IMPROVED) {
					for (BnBNode node : surroundingNodes.tailSet(currentNode, false)) {
						boolean individualBoundChanged = false;
						for (PossibleZ comparedPossibleZ : node.getPossibleZs()) {
							double estimator = currentNode.getUpperBoundZ().getEstimators().get(comparedPossibleZ);
							if (comparedPossibleZ.updateDualBound(subproblemDualBound - estimator)) {
								individualBoundChanged = true;
							}
						}
						if (individualBoundChanged) {
							node.updateDualBound();
						}
					}
				}
				if (bnbStrategies.getEstimatorStrategy() == BnBStrategies.BnBEstimatorStrategy.ESTIMATORS_IMPROVED
						|| bnbStrategies.getEstimatorStrategy() == BnBStrategies.BnBEstimatorStrategy.ESTIMATORS_HRS) {
					for (BnBNode node : surroundingNodes.headSet(currentNode, false)) {
						boolean individualBoundChanged = false;
						for (PossibleZ comparedPossibleZ : node.getPossibleZs()) {
							double estimator = currentNode.getLowerBoundZ().getEstimators().get(comparedPossibleZ);
							if (comparedPossibleZ.updateDualBound(subproblemDualBound - estimator)) {
								individualBoundChanged = true;
							}
						}
						if (individualBoundChanged) {
							node.updateDualBound();
						}
					}
				}
				updateGlobalDualBound(currentNode);
			}
			return boundChanged;
		}
		
		/**
		 * Receives primal and dual bounds for updating and returns whether the subproblem may be terminated.
		 */
		protected boolean updateBoundsAndDecideTermination(double improvedPrimalBound, double subproblemPrimalBound, double subproblemDualBound) {
			boolean boundChanged = this.updatePrimalDualBounds(improvedPrimalBound, subproblemDualBound);
			
			//Checks whether the subproblem can be terminated with respect to the global primal bound and its current dual bound.
			if (isOptimal(getPrimalBound(), subproblemDualBound)) {
				//Initializes the time stamp and remaining possible z for pruning.
				if (timeStampLastPrune == null) {
					timeStampLastPrune = System.nanoTime();
				}
				//Checks whether we can hope for pruning additional possible z.
				boolean pruningPossible = false;
				if (boundChanged) {
					for (PossibleZ boundingZ : boundToRemainingPossibleZwithEstimator.keySet()) {
						List<PossibleZ> remainingPossibleZwithEstimator = boundToRemainingPossibleZwithEstimator.get(boundingZ);
						for (int i = remainingPossibleZwithEstimator.size()-1; i >= 0; i--) {
							PossibleZ possibleZ = remainingPossibleZwithEstimator.get(i);
							//Checks whether the dual bound is good enough for pruning.
							if (isOptimal(getPrimalBound(), possibleZ.getDualBound())) {
								//Removes the possible z and sets the time stamp.
								remainingPossibleZwithEstimator.remove(i);
								timeStampLastPrune = System.nanoTime();
							}
							//If raising the dual bound up to the primal bound would be sufficient
							//to prune the possible z then we can hope for pruning it later and
							//don't terminate the subproblem if the last time stamp is not too far in the past.
							else if (isOptimal(getPrimalBound(), subproblemPrimalBound - boundingZ.getEstimators().get(possibleZ))) {
								pruningPossible = true;
							}
						}
					}
				}
				//If the bounds didn't change, then pruning was possible before and is still now
				else {
					pruningPossible= true;
				}
				//If there is hope for pruning further possible z and the last pruning
				//has happened in the last 10 seconds then we continue solving the subproblem.
				if (!pruningPossible || (System.nanoTime()-timeStampLastPrune)/Math.pow(10, 9) >= 10) {
					String output = "###Terminated due to global primal bound";
					try {
						writeOutput(output);
					} catch (IOException e) {}
					
					return true;
				}
			}
			return false;
		}
	}

	/**
	 * Solves the linear relaxation for a given node.
	 */
	private void solveLinear(BnBNode chosenNode) throws IOException, GRBException {
		String output;
		if (chosenNode.isApplyLagrangeRelaxation()) {
			output = "\n##### Solving Lagrangean relaxed LP for "+chosenNode;
		}
		else {
			output = "\n##### Solving robust LP for "+chosenNode;
		}
		writeOutput(output);
		
		//Updates bounds for the linear subproblem.
		subproblemBounded.updateBoundsAndFormulation(chosenNode);
		
		//Adds optimality-cuts if the corresponding option is chosen.
		if (bnbStrategies.getLpOptimalityCutsStrategy() == BoundingLPOptimalityCutsStrategy.LPOPTCUTS_ENABLE) {
			subproblemBounded.addOptimalityCuts(chosenNode.getLowerBoundZ(), chosenNode.getUpperBoundZ());
		}
		
		writeOutput("");
		
		//Tries to solve the problem within the time limit.
		long startTimeLinearizedSubproblem = System.nanoTime();
		subproblemBounded.solveRelaxed(getRemainingTime(), chosenNode);
		timeInLinearizedSubproblems += (System.nanoTime() - startTimeLinearizedSubproblem)/Math.pow(10, 9);
		
		//If the linear relaxation is infeasible then the node's dual bound is infinity and it can be pruned.
		if (subproblemBounded.isStatusInfeasible()) {
			output = "LP is infeasible\n"
					+ "Node is pruned";
			writeOutput(output);
			numberRemainingPossibleZ -= chosenNode.getPossibleZs().size();
			return;
		}
		
		//If the linear relaxation could not be solved to optimality then we branch the node in the middle of the interval or solve it as an MILP if it already consists of only one possible value for z.
		//As we obtain no new bounds from the linear relaxation, we use the old bound of the node for the new nodes emerging from branching.
		else if (!subproblemBounded.isStatusOptimal()) {
			output = "LP could not be solved to optimality.";
			if (chosenNode.getPossibleZs().size() >1) {
				output += "\nNode is branched at the middle of the interval, using bounds from parent node.";
				branchNodeMiddle(chosenNode, chosenNode.getDualBound(), chosenNode.isApplyLagrangeRelaxation());
			}
			//If the node consists only of one possible value for z then we reinsert it and solve it as an integer robust subproblem.
			else {
				output += "\nNode is considered as an MILP.";
				chosenNode.setSolveInteger(true);
				remainingNodes.add(chosenNode);
			}
			writeOutput(output);
			return;
		}
		
		//We solved the linear relaxation to optimality and use the objective value to update the nodes dual bound.
		chosenNode.updateDualBound(subproblemBounded.getDualBound());
		
		//If the solver found an integer solution to the linear relaxation then we obtain a new incumbent.
		if (subproblemBounded.isSolutionIntegerFeasible()) {
			output = "\nSolution is integer feasible.";
			double incumbentValue;
			//If robustness constraints are applied then the value of the LP is the value of the incumbent
			if (!chosenNode.isApplyLagrangeRelaxation()) {
				incumbentValue = subproblemBounded.getPrimalBound();
			}
			//Otherwise, we have to compute the solution value using the optimal value for z
			else {
				output += "\nSolution is converted from Lagrangean relaxation into a robust solution with proper values for p and z.";
				double bilinearZ = subproblemBounded.computeOptimalBilinearZ();
				incumbentValue = subproblemBounded.computeBilinearSolutionValue(bilinearZ);
			}
			if (incumbentValue < getPrimalBound()) {
				output += "\nSolution defines new primal bound of value "+incumbentValue;
				setPrimalBound(incumbentValue);
				solutionValues = subproblemBounded.getNominalVariablesSolutionValues();
			}
			else {
				output += "\nThe incumbent value of "+incumbentValue+" does not define a new primal bound.";
			}
			writeOutput(output);
		}
		
		//Tries to use the new dual bound for pruning the node.
		if (isOptimal(getPrimalBound(), subproblemBounded.getDualBound())) {
			output = "\nPrune node, as dual bound is within optimality gap.";
			writeOutput(output);
			minimumDualBoundPrunedNodes = Math.min(minimumDualBoundPrunedNodes, subproblemBounded.getDualBound());
			numberRemainingPossibleZ -= chosenNode.getPossibleZs().size();
		}
		//If the node cannot be pruned then we check whether the node will be branched or solved as an integer robust subproblem.
		else {
			//There is no need to branch if there is only one remaining value of z
			if (chosenNode.getPossibleZs().size() == 1) {
				chosenNode.setSolveInteger(true);
				remainingNodes.add(chosenNode);
				output = "Only one value remaining, solve as integer.";
				writeOutput(output);
				return;
			}
			//In theory, we have linear=bilinear if z is equal to its bounds
			//Checks this in case of numerical issues
			if (subproblemBounded.getComputedZ().isPresent()) {
				double computedZ = subproblemBounded.getComputedZ().get();
				if (computedZ <= chosenNode.getLowerBoundZ().getValue()+algorithmParameters.getFeasibilityTolerance()
					|| computedZ >= chosenNode.getUpperBoundZ().getValue()-algorithmParameters.getFeasibilityTolerance()) {
					chosenNode.setSolveInteger(true);
					remainingNodes.add(chosenNode);
					output = "Variable z is at its bounds, and thus we have linear=bilinear. Solve as integer.";
					writeOutput(output);
					return;
				}
			}
			
			//Otherwise, computes the optimal bilinear z for the solution and projects it into the interval for z defined by the node.
			double bilinearZ = subproblemBounded.computeOptimalBilinearZ();
			bilinearZ = Math.min(chosenNode.getUpperBoundZ().getValue(), bilinearZ);
			bilinearZ = Math.max(chosenNode.getLowerBoundZ().getValue(), bilinearZ);
			//Computes a primal bound on the optimal objective value of the bilinear problem with restricted z.
			double bilinearPrimalBound = subproblemBounded.computeBilinearSolutionValue(bilinearZ);
			output = "Bilinear primal bound = "+bilinearPrimalBound;
			writeOutput(output);
			//If the dual bound linear relaxation and the primal bound on the bilinear problem are close together then we stop branching.
			if (isOptimal(bilinearPrimalBound, subproblemBounded.getDualBound())) {
				chosenNode.setSolveInteger(true);
				remainingNodes.add(chosenNode);
				output = "Stop branching due to low bilinear primal bound.";
				writeOutput(output);
				return;
			}
			
			//Otherwise we branch the node according to the chosen option.
			//Decides whether we consider lagrange relaxation or not for the next LP
			boolean applyLagrangeRelaxation;
			//If hybrid, we might change from lagrange to non-lagrange in case the gap between the value of the lagrange relaxation and the bilinear primal bound is small 
			if (bnbStrategies.getLagrangeRelaxLPStrategy() == BoundingLPLagrangeRelaxStrategy.LAGRANGERELAX_HYBRID) {
				applyLagrangeRelaxation = chosenNode.isApplyLagrangeRelaxation() && getRelativeGap(bilinearPrimalBound, subproblemBounded.getPrimalBound()) > Math.sqrt(algorithmParameters.getRelativeGapTolerance());
			}
			else {
				applyLagrangeRelaxation = chosenNode.isApplyLagrangeRelaxation();
			}
			
			if (bnbStrategies.getBranchingStrategy() == BnBStrategies.BnBBranchingStrategy.BRANCHAROUND_Z) {
				branchNodeAroundZ(chosenNode, applyLagrangeRelaxation);
			}
			else if (bnbStrategies.getBranchingStrategy() == BnBStrategies.BnBBranchingStrategy.BRANCHAROUND_CONVEXCOMB){
				branchNodeConvexcombination(chosenNode, applyLagrangeRelaxation);
			}
		}
	}
	
	/**
	 * Branches a node around a given branching point.
	 * Sets a given value as dual bound for the child nodes.
	 * Defines whether child nodes should be considered with lagrange relaxation.
	 * Stores information for warm start if the corresponding option is chosen.
	 */
	private void branchNode(BnBNode chosenNode, double branchingPoint, double dualBoundLP, boolean applyLagrangeRelaxation) throws IOException {
		//Splits the interval of possible z around the branching point.
		List<PossibleZ> lowerPossibleZs = new ArrayList<PossibleZ>(chosenNode.getPossibleZs().size());
		int index = 0;
		while (index < chosenNode.getPossibleZs().size()-1 && chosenNode.getPossibleZs().get(index).getValue() <= branchingPoint) {
			lowerPossibleZs.add(chosenNode.getPossibleZs().get(index));
			index++;
		}
		BnBNode lowerNode = new BnBNode(lowerPossibleZs, dualBoundLP, applyLagrangeRelaxation);
		

		List<PossibleZ> upperPossibleZs = new ArrayList<PossibleZ>(chosenNode.getPossibleZs().size()-lowerPossibleZs.size());
		while (index < chosenNode.getPossibleZs().size()) {
			upperPossibleZs.add(chosenNode.getPossibleZs().get(index));
			index++;
		}
		BnBNode upperNode = new BnBNode(upperPossibleZs, dualBoundLP, applyLagrangeRelaxation);
		
		//Stores information for warm start if the corresponding option is chosen.
		if (bnbStrategies.warmStartStrategy == BoundingLPWarmStartStrategy.WARMSTART_ENABLE) {
			int[] VBasisNominal = subproblemBounded.getVBasisNominal();
			int[] VBasisP = subproblemBounded.getVBasisP();
			Integer VBasisZ = subproblemBounded.getVBasisZ();
			int[] CBasisNominal = subproblemBounded.getCBasisNominal();
			
			List<NonBasicSlackConstraint> parentNonBasicSlackConstraints = subproblemBounded.getNonBasicSlackConstraints();
			double[] upperBoundsUsedForSubstituting = null;
			//We only have to take care of substitution if we didn't apply lagrange to the parent node
			if (!chosenNode.isApplyLagrangeRelaxation()) {
				upperBoundsUsedForSubstituting = subproblemBounded.getUpperBoundsUsedForSubstituting();
			}
			
			for (BnBNode childNode : new BnBNode[] {lowerNode, upperNode}) {
				childNode.setParentVBasisNominal(VBasisNominal);
				childNode.setParentVBasisP(VBasisP);
				childNode.setParentVBasisZ(VBasisZ);
				childNode.setParentCBasisNominal(CBasisNominal);
				childNode.setParentNonBasicSlackConstraints(parentNonBasicSlackConstraints);
				childNode.setUpperBoundsUsedForSubstituting(upperBoundsUsedForSubstituting);

			}
		}

		String output = "New Nodes: "+lowerNode+" "+upperNode;
		writeOutput(output);
		
		//Adds the corresponding new nodes to the set of remaining nodes.
		remainingNodes.add(upperNode);
		remainingNodes.add(lowerNode);
	}
	
	/**
	 * Branches a node around the mean value of the interval.
	 */
	private void branchNodeMiddle(BnBNode chosenNode, double dualBoundLP, boolean applyLagrangeRelaxation) throws IOException {
		double branchingPoint = (chosenNode.getLowerBoundZ().getValue()+chosenNode.getUpperBoundZ().getValue())/2;
		branchNode(chosenNode, branchingPoint, dualBoundLP, applyLagrangeRelaxation);
	}
	
	/**
	 * Branches a node around the solution value of z. If the value is not present, we compute the optimal bilinear value of z.
	 */
	private void branchNodeAroundZ(BnBNode chosenNode, boolean applyLagrangeRelaxation) throws IOException {
		double branchingPoint;
		if (subproblemBounded.getComputedZ().isPresent()) {
			branchingPoint = subproblemBounded.getComputedZ().get();
		}
		else {
			branchingPoint = subproblemBounded.computeOptimalBilinearZ();
			branchingPoint = Math.min(chosenNode.getUpperBoundZ().getValue(), branchingPoint);
			branchingPoint = Math.max(chosenNode.getLowerBoundZ().getValue(), branchingPoint);
		}
		branchNode(chosenNode, branchingPoint, subproblemBounded.getDualBound(), applyLagrangeRelaxation);
	}
	
	/**
	 * Branches a node around a convex combination the solution value of z and the mean of the interval.
	 */
	private void branchNodeConvexcombination(BnBNode chosenNode, boolean applyLagrangeRelaxation) throws IOException {
		//Value of optimal bilinear z if we considered the Lagrangean relaxation.
		//Otherwise, the value of the computed z.
		double zValue;
		if (subproblemBounded.getComputedZ().isEmpty()) {
			zValue = subproblemBounded.computeOptimalBilinearZ();
			zValue = Math.min(chosenNode.getUpperBoundZ().getValue(), zValue);
			zValue = Math.max(chosenNode.getLowerBoundZ().getValue(), zValue);
		}
		else {
			zValue = subproblemBounded.getComputedZ().get();
		}
		
		//Bounds for which a branching point is effective.
		double lowerBoundEffectiveBranching = subproblemBounded.getLowerBoundEffectiveBranching(chosenNode);
		double upperBoundEffectiveBranching = subproblemBounded.getUpperBoundEffectiveBranching(chosenNode);
		
		
		//Sets the branching point to be the mean of zValue and the mean of the effective interval.
		double branchingPoint = zValue/2+(lowerBoundEffectiveBranching+upperBoundEffectiveBranching)/4;
		
		//Projects the branching point into the middle 80% of the effective interval
		branchingPoint = Math.min(upperBoundEffectiveBranching-0.1*(upperBoundEffectiveBranching-lowerBoundEffectiveBranching), branchingPoint);
		branchingPoint = Math.max(lowerBoundEffectiveBranching+0.1*(upperBoundEffectiveBranching-lowerBoundEffectiveBranching), branchingPoint);
		
		String output = "Branch around "+branchingPoint;
		writeOutput(output);
		
		branchNode(chosenNode, branchingPoint, subproblemBounded.getDualBound(), applyLagrangeRelaxation);
	}
	
	/**
	 * Updates the dual bound to be the minimum of all pruned nodes, the current node, and the minimum dual bound of the not yet considered nodes.
	 */
	private void updateGlobalDualBound(BnBNode currentNode) {
		double newDualBoud = Math.min(minimumDualBoundPrunedNodes, currentNode.getDualBound());
		for (BnBNode node : remainingNodes) {
			newDualBoud = Math.min(newDualBoud, node.getDualBound());
		}
		setDualBound(newDualBoud);
	}
	
	/**
	 * Updates the dual bound to be the minimum of all pruned nodes and the minimum dual bound of the not yet considered nodes.
	 */
	private void updateGlobalDualBound() {
		double newDualBoud = minimumDualBoundPrunedNodes;
		for (BnBNode node : remainingNodes) {
			newDualBoud = Math.min(newDualBoud, node.getDualBound());
		}
		setDualBound(newDualBoud);
	}
	
	/**
	 * Prunes all nodes whose dual bound is close enough to the current primal bound.
	 * Returns a list of all pruned nodes.
	 */
	private List<BnBNode> pruneNodes() throws IOException {
		List<BnBNode> prunedNodes = new ArrayList<BnBNode>();
		String output = "Pruned Node";
		Iterator<BnBNode> nodeIterator = remainingNodes.iterator();
		while (nodeIterator.hasNext()) {
			BnBNode node = nodeIterator.next();
			if (isOptimal(getPrimalBound(), node.getDualBound())) {
				output += " "+node;
				nodeIterator.remove();
				minimumDualBoundPrunedNodes = Math.min(minimumDualBoundPrunedNodes, node.getDualBound());
				numberRemainingPossibleZ -= node.getPossibleZs().size();
				prunedNodes.add(node);
			}
		}
		if (!prunedNodes.isEmpty()) {
			writeOutput(output);
		}
		return prunedNodes;
	}
	
	/**
	 * Tries to prune the possible values for z within nodes using the primal bound and the individual dual bounds.
	 * Returns a List of pruned values. 
	 */
	private List<PossibleZ> prunePossibleZinNodes() throws IOException {
		List<PossibleZ> prunedPossibleZs = new ArrayList<PossibleZ>();
		String output = "Pruned Possible z:\n";
		for (BnBNode node : remainingNodes) {
			for (int zIndex = node.getPossibleZs().size()-1; zIndex >= 0; zIndex--) {
				PossibleZ possibleZ = node.getPossibleZs().get(zIndex);
				if (isOptimal(getPrimalBound(), possibleZ.getDualBound())) {
					output += " "+possibleZ;
					minimumDualBoundPrunedNodes = Math.min(minimumDualBoundPrunedNodes, possibleZ.getDualBound());
					node.getPossibleZs().remove(zIndex);
					prunedPossibleZs.add(possibleZ);
					numberRemainingPossibleZ--;
				}
			}
		}
		if (!prunedPossibleZs.isEmpty()) {
			writeOutput(output);
		}
		return prunedPossibleZs;
	}
	
	/**
	 * Returns the number of uncertain variables in the robust problem.
	 */
	public Integer getNumberUncertainVariables() {
		return  subproblemBounded.getUncertainVariables().length;
	}
	
	/**
	 * Returns the number of cliques computed for the robust problem.
	 * Returns the number of uncertain variables if we did not compute cliques.
	 */
	public int getNumberCliques() {
		if (cliquePartitioning != null) {
			return cliquePartitioning.getCliques().size();
		}
		else {
			return subproblemBounded.getUncertainVariables().length;
		}
	}
	
	/**
	 * Returns the number of possible values of z computed for the robust problem.
	 */
	public Integer getNumberPossibleZ() {
		return possibleZs.size();
	}
	
	/**
	 * Returns the primal bound obtained after the first integer robust subproblem has been solved.
	 */
	public Double getPrimalBoundAfterFirstInteger() {
		return primalBoundAfterFirstInteger;
	}

	/**
	 * Returns the number of integer subproblems considered.
	 */
	public Integer getNumberIntegerStarted() {
		return numberIntegerSubproblemsStarted;
	}

	/**
	 * Returns the cumulative time (in seconds) elapsed while solving integer robust subproblems. 
	 */
	public Double getTimeSpendInInteger() {
		return timeInIntegerSubproblems;
	}

	/**
	 * Returns the number of linear relaxations considered.
	 */
	public Integer getNumberLPStarted() {
		return numberLinearSubproblemsStarted;
	}

	/**
	 * Returns the cumulative time (in seconds) elapsed while solving linear relaxations. 
	 */
	public Double getTimeSpendInLP() {
		return timeInLinearizedSubproblems;
	}
	
	/**
	 * Specifies all strategies for the BnB algorithm;
	 */
	public static class BnBStrategies extends BoundingSubproblemsStrategies{
		/**
		 * Specifies whether we use the improved estimators from the paper or the estimators
		 * of Hansknecht, Richter, Stiller to carry over dual bounds from one subproblem to another
		 * or whether we don't use estimators at all.
		 */
		public enum BnBEstimatorStrategy{
			ESTIMATORS_IMPROVED,
			ESTIMATORS_HRS,
			ESTIMATORS_NONE;
		}
		
		/**
		 * Specifies whether we branch around the linear relaxation value of z or
		 * whether we branch around a convex combination of z and the middle of the interval.
		 */
		public enum BnBBranchingStrategy{
			BRANCHAROUND_CONVEXCOMB,
			BRANCHAROUND_Z;
		}

		BnBEstimatorStrategy estimatorStrategy;
		BnBBranchingStrategy branchingStrategy;
		
		/**
		 * Constructor obtaining arguments which are matched to the enums defining strategies.
		 */
		public BnBStrategies(List<String> argList, AlgorithmParameters algorithmParameters) throws IOException {
			super(argList, algorithmParameters);
		}

		@Override
		void setDefaultStrategies() {
			super.setDefaultStrategies();
			improvingZStrategy = ImprovingZStrategy.IMPROVINGZ_ENABLE;
			terminationStrategy = BoundingTerminationStrategy.TERMINATION_ENABLE;
			estimatorStrategy = BnBEstimatorStrategy.ESTIMATORS_IMPROVED;
			branchingStrategy = BnBBranchingStrategy.BRANCHAROUND_CONVEXCOMB;
			lpOptimalityCutsStrategy = BoundingLPOptimalityCutsStrategy.LPOPTCUTS_DISABLE;
			ipOptimalityCutsStrategy = BoundingIPOptimalityCutsStrategy.IPOPTCUTS_DISABLE;
			lagrangeRelaxStrategy = BoundingLPLagrangeRelaxStrategy.LAGRANGERELAX_HYBRID;
			warmStartStrategy = BoundingLPWarmStartStrategy.WARMSTART_ENABLE;
		}
		
		public BnBEstimatorStrategy getEstimatorStrategy() {
			return estimatorStrategy;
		}
		public void setEstimatorStrategy(BnBEstimatorStrategy estimatorStrategy) {
			this.estimatorStrategy = estimatorStrategy;
		}
		
		public BnBBranchingStrategy getBranchingStrategy() {
			return branchingStrategy;
		}
		public void setBranchingStrategy(BnBBranchingStrategy branchingStrategy) {
			this.branchingStrategy = branchingStrategy;
		}
	}
}
