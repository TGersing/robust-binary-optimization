package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.TreeSet;

import alg.AlgBranchAndBound.BnBStrategies.BnBIntegerOptimalityCutsStrategy;
import alg.AlgBranchAndBound.BnBStrategies.BnBLPOptimalityCutsStrategy;
import alg.SubproblemBoundedGurobi.BoundingSubproblemsStrategies;
import gurobi.GRBException;
import util.BnBNode;
import util.CliquePartitioning;
import util.ConflictGraph;
import util.PossibleZ;
import util.Variable;

/**
 * This class implements the branch and bound algorithm as described
 * in the paper "A Branch & Bound Algorithm for Robust Binary Optimization with Budget Uncertainty".
 * 
 * @author Timo Gersing
 */
public class AlgBranchAndBound extends AbstractAlgorithm implements RobustAlgorithm {
	/**
	 * The robust subproblem that we adapt for different sets Z containing possible values for z.
	 */
	private SubproblemBoundedGurobi subproblemBounded;

	/**
	 * A partitioning of the uncertain variables into cliques.
	 */
	private CliquePartitioning cliquePartitioning;
	
	/**
	 * A list containing possible optimal values for z.
	 */
	private List<PossibleZ> possibleZs;
	
	/**
	 * A sorted set of tree nodes containing Sets of yet to be considered values of z.
	 * The first node is always the one that is chosen for the next iteration.
	 * The nodes are sorted with respect to their dual bound.
	 */
	private TreeSet<BnBNode> remainingNodes = new TreeSet<BnBNode>();
	
	/**
	 * A set containing possible values for z that we already considered for an integer robust subproblem.
	 * Important for deciding which bounds to use for the optimality-cuts.
	 */
	private TreeSet<PossibleZ> consideredZForRobustSubproblems = new TreeSet<PossibleZ>();
	
	/**
	 * Tracks the number of yet to be considered values of z.
	 */
	private int numberRemainingPossibleZ;
	
	/**
	 * Specifies strategies.
	 */
	private BnBStrategies bnbStrategies;
	
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
	 * Executes the BnB.
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
		if (bnbStrategies.getCliqueStrategy() == RobustAlgorithmStrategies.CliqueStrategy.CLIQUES) {
			conflictGraph = new ConflictGraph(subproblemBounded.getModel(), subproblemBounded.getUncertainModelVariables(), algorithmParameters);
			cliquePartitioning = new CliquePartitioning(subproblemBounded.getUncertainVariables(), conflictGraph, algorithmParameters);
			subproblemBounded.setCliquePartitioning(cliquePartitioning);
		}
		
		//Computes a list of possible optimal choices for z respecting the chosen filtering and clique strategies.
		possibleZs = computePossibleZs(subproblemBounded.getUncertainVariables(), subproblemBounded.getGamma(), bnbStrategies, conflictGraph, algorithmParameters);
		
		//Initialize the root node containing all possible values for z.
		remainingNodes.add(new BnBNode(possibleZs));
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
			BnBNode chosenNode = remainingNodes.first();
			
			dualBound = Math.min(minimumDualBoundPrunedNodes, chosenNode.getDualBound());
			output = "\n##### Branching Tree Information\n"
					+ "Elapsed Time = "+(System.nanoTime()-startTime)/Math.pow(10, 9)+"\n"
					+ "Number Remaining Nodes = "+remainingNodes.size()
					+ "\nNumber Remaining Possible Z = "+numberRemainingPossibleZ
					+ "\nCurrent Primal Bound = "+primalBound
					+ "\nCurrent Dual Bound = "+dualBound
					+ "\nRelative Optimality Gap = "+(AbstractAlgorithm.getRelativeGap(primalBound, dualBound)*100)+"%";
			writeOutput(output);
			
			//Solves the integer robust subproblem if the current node is marked correspondingly.
			if (chosenNode.isSolveInteger()) {
				solveInteger(chosenNode);
			}
			//Solves the linear relaxation otherwise. Either branches the node afterwards or marks it to be solved as an integer subproblem.
			else {
				solveLinear(chosenNode);
			}

			//Updates the primal bound if we did not yet considered the second integer robust subproblem.
			if (numberIntegerSubproblemsStarted <= 1) {
				primalBoundAfterFirstInteger = primalBound;
			}
			
			//Checks whether there exist nodes that can be pruned.
			//It is sufficient to consider the last node in the set, as they are sorted w.r.t. their dual bounds. 
			while (!remainingNodes.isEmpty()) {
				BnBNode highestDualBoundNode = remainingNodes.last();
				if (isOptimal(primalBound, highestDualBoundNode.getDualBound())) {
					output = "Pruned Node "+highestDualBoundNode;
					writeOutput(output);
					remainingNodes.remove(highestDualBoundNode);
					minimumDualBoundPrunedNodes = Math.min(minimumDualBoundPrunedNodes, highestDualBoundNode.getDualBound());
					numberRemainingPossibleZ -= highestDualBoundNode.getPossibleZs().size();
				}
				else {
					break;
				}
			}
		}
		
		//Updates the dual bound to be the minimum of all pruned nodes and
		//the minimum dual bound of the not yet considered nodes.
		dualBound = minimumDualBoundPrunedNodes;
		for (BnBNode node : remainingNodes) {
			dualBound = Math.min(dualBound, node.getDualBound());
		}
		
		//Stores best solution found
		solution = new LinkedHashMap<Variable, Double>(solutionValues.length);
		for (int i = 0; i < solutionValues.length; i++) {
			solution.put(subproblemBounded.getNominalVariables()[i], solutionValues[i]);
		}
	}
	
	/**
	 * Solves the integer robust subproblem for a given node.
	 */
	private void solveInteger(BnBNode chosenNode) throws IOException, GRBException {
		String output = "\n##### Solving MILP for "+chosenNode;
		writeOutput(output);
		numberIntegerSubproblemsStarted++;
		
		//Updates the robust subproblem for the chosen node and sets the current primal bound for the termination.
		subproblemBounded.updateBounds(chosenNode);
		subproblemBounded.setGlobalPrimalBound(primalBound);
		
		//Determines estimators and optimality cuts if the corresponding options are chosen.
		if (bnbStrategies.getEstimatorStrategy() != BnBStrategies.BnBEstimatorStrategy.NO_ESTIMATORS || bnbStrategies.getIntOptimalityCutsStrategy() == BnBIntegerOptimalityCutsStrategy.CUTS) {
			//The set of remaining values.
			TreeSet<PossibleZ> remainingPossibleZs = new TreeSet<PossibleZ>();
			
			//Tries to prune remaining values. All not yet pruned values are then added to the set of remaining values.
			List<PossibleZ> prunedPossibleZs = new ArrayList<PossibleZ>();
			for (BnBNode node : remainingNodes) {
				prunedPossibleZs.addAll(performPruning(node));
				remainingPossibleZs.addAll(node.getPossibleZs());
			}
			//Prints to log if pruning was successful.
			if (!prunedPossibleZs.isEmpty()) {
				output = "Pruned values via estimated dual bound:\n";
				for (PossibleZ possibleZ : prunedPossibleZs) {
					output += possibleZ+"  ";
				}
				writeOutput(output);
			}
			
			//The bounds for the optimality-cuts are given by the smallest and largest remaining values
			//around the chosen node such that there exists no already considered value of z in between.
			if (bnbStrategies.getIntOptimalityCutsStrategy() == BnBIntegerOptimalityCutsStrategy.CUTS) {
				//If we set optimality-cuts then we do not compute estimators for values outside of the corresponding bounds.				
				//Removes all remaining values that are smaller than an already considered value for z that is smaller than the values in the node.
				PossibleZ greatestSmallerConsideredZ = consideredZForRobustSubproblems.floor(chosenNode.getLowerBoundZ());
				if (greatestSmallerConsideredZ != null) {
					remainingPossibleZs.headSet(greatestSmallerConsideredZ).clear();
				}
				//Removes all remaining values that are greater than a already considered value for z that is larger than the values in the node.
				PossibleZ leastGreaterConsideredZ = consideredZForRobustSubproblems.ceiling(chosenNode.getUpperBoundZ());
				if (leastGreaterConsideredZ != null) {
					remainingPossibleZs.tailSet(leastGreaterConsideredZ).clear();
				}				
				
				//Initializes the bounds for the optimality-cuts and adds these to the problem.
				PossibleZ lowerBoundOptCut = remainingPossibleZs.first();
				PossibleZ upperBoundOptCut = remainingPossibleZs.last();
				subproblemBounded.addOptimalityCuts(lowerBoundOptCut, upperBoundOptCut);
			}
			
			//Computes estimators according to the chosen option.
			if (bnbStrategies.getEstimatorStrategy() == BnBStrategies.BnBEstimatorStrategy.IMPROVED_ESTIMATORS) {
				//If the node consists of only one possible value for z then this value defines the estimators for all lower and higher values.
				if (chosenNode.getLowerBoundZ() == chosenNode.getUpperBoundZ()) {
					chosenNode.getUpperBoundZ().setImprovedEstimators(remainingPossibleZs, subproblemBounded.getGamma(), subproblemBounded.getUncertainVariables());
				}
				//Otherwise the estimators are defined by the lowest and highest possible z.
				else {
					chosenNode.getLowerBoundZ().setImprovedEstimatorsBelowZ(remainingPossibleZs, subproblemBounded.getGamma(), subproblemBounded.getUncertainVariables());
					chosenNode.getUpperBoundZ().setImprovedEstimatorsAboveZ(remainingPossibleZs, subproblemBounded.getUncertainVariables());
				}
			}
			//If we use the estimators by Hansknecht, Stichter, Stiller then these are defined by the lowest value in the node.
			else if (bnbStrategies.getEstimatorStrategy() == BnBStrategies.BnBEstimatorStrategy.HRS_ESTIMATORS) {
				chosenNode.getLowerBoundZ().setHRSEstimators(remainingPossibleZs, subproblemBounded.getGamma());
			}
		}
		writeOutput("");
		
		//Tries to solve the robust subproblem in the time limit.
		long starttimeInteger = System.nanoTime();
		subproblemBounded.solveInteger(getRemainingTime());
		timeInIntegerSubproblems+=(System.nanoTime()-starttimeInteger)/Math.pow(10, 9);
		
		//Updates the primal bound. Considers the global primal bound of the subproblem as incumbent, as this respects solutions with improved z.
		double incumbentValue = subproblemBounded.getGlobalPrimalBound();
		if (incumbentValue < primalBound) {
			primalBound = incumbentValue;
			solutionValues = subproblemBounded.getIncumbentValues();
		}
		
		//Updates the individual dual bound for each possible z in the node. 
		//Also adds them to the set of already considered possible z.
		for (PossibleZ possibleZ : chosenNode.getPossibleZs()) {
			possibleZ.updateDualBound(subproblemBounded.getDualBound());
			consideredZForRobustSubproblems.add(possibleZ);
		}
		
		//Removes the current node from the set of remaining nodes and updates the dual bound.
		remainingNodes.remove(chosenNode);
		numberRemainingPossibleZ -= chosenNode.getPossibleZs().size();
		chosenNode.updateDualBound();
		minimumDualBoundPrunedNodes = Math.min(minimumDualBoundPrunedNodes, chosenNode.getDualBound());		
		
		//Uses estimators to improve the individual dual bounds of other possible values for z.
		if (bnbStrategies.getEstimatorStrategy() != BnBStrategies.BnBEstimatorStrategy.NO_ESTIMATORS) {
			//If the node consists of only one possible value for z then this value defines the estimators for all lower and higher values.
			if (chosenNode.getLowerBoundZ() == chosenNode.getUpperBoundZ()) {
				chosenNode.getLowerBoundZ().estimateDualBounds(subproblemBounded.getDualBound());
			}
			//Otherwise the estimators are defined by the lowest and highest possible z.
			else {
				chosenNode.getLowerBoundZ().estimateDualBounds(subproblemBounded.getDualBound());
				chosenNode.getUpperBoundZ().estimateDualBounds(subproblemBounded.getDualBound());	
			}
			
			//Updates the dual bound of all remaining nodes.
			//Reinserts a node into the tree set if the dual bound has changed to maintain the ordering.
			Iterator<BnBNode> iterator = remainingNodes.iterator();
			List<BnBNode> reinsertNodes = new ArrayList<BnBNode>(remainingNodes.size());
			while (iterator.hasNext()) {
				BnBNode node = iterator.next();
				double prevDualBound = node.getDualBound();
				node.updateDualBound();
				if (prevDualBound != node.getDualBound()) {
					reinsertNodes.add(node);
					iterator.remove();
				}
			}
			remainingNodes.addAll(reinsertNodes);
		}
	}

	/**
	 * Solves the linear relaxation for a given node.
	 */
	private void solveLinear(BnBNode chosenNode) throws IOException, GRBException {
		String output = "\n##### Solving LP for "+chosenNode;
		writeOutput(output);
		numberLinearSubproblemsStarted++;
		
		//Removes the node from the set remaining nodes.
		remainingNodes.remove(chosenNode);
		
		//Tries to prune possible values for z in the node.
		List<PossibleZ> prunedPossibleZs = performPruning(chosenNode);
		if (!prunedPossibleZs.isEmpty()) {
			output = "Pruned values via estimated dual bound:\n";
			for (PossibleZ possibleZ : prunedPossibleZs) {
				output += possibleZ+"  ";
			}
			output += "\nChosen Node After Individual Pruning: "+chosenNode;
			writeOutput(output);
		}
		
		//Updates bounds for the linear subproblem.
		subproblemBounded.updateBounds(chosenNode);
		
		//Adds optimality-cuts if the corresponding option is chosen.
		if (bnbStrategies.getLpOptimalityCutsStrategy() == BnBLPOptimalityCutsStrategy.CUTS) {
			subproblemBounded.addOptimalityCuts(chosenNode.getLowerBoundZ(), chosenNode.getUpperBoundZ());
		}
		
		writeOutput("");
		
		//Tries to solve the problem within the time limit.
		long startTimeLinearizedSubproblem = System.nanoTime();
		subproblemBounded.solveRelaxed(getRemainingTime());
		timeInLinearizedSubproblems += (System.nanoTime() - startTimeLinearizedSubproblem)/Math.pow(10, 9);
		
		//If the linear relaxation is infeasible then the node's dual bound is infinity and it can be pruned.
		if (subproblemBounded.isStatusInfeasible()) {
			output = "LP is infeasible\n"
					+ "Node is pruned";
			writeOutput(output);
			numberRemainingPossibleZ -= chosenNode.getPossibleZs().size();
			return;
		}
		//If the linear relaxation could not be solved to optimality then we branch the node in the middle of the interval.
		//As we obtain no new bounds from the linear relaxation, we use the old bound of the node for the new nodes emerging from branching.
		else if (!subproblemBounded.isStatusOptimal()) {
			output = "LP could not be solved to optimality.";
			if (chosenNode.getPossibleZs().size() >1) {
				output += "\nNode is branched at the middle of the interval, using bounds from parent node.";
				branchNodeMiddle(chosenNode, chosenNode.getDualBound());
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
		chosenNode.updateDualBoundLP(subproblemBounded.getDualBound());

		//If the solver found an integer solution to the linear relaxation then we obtain a new incumbent.
		if (subproblemBounded.isSolutionIntegerFeasible()) {
			output = "\nSolution is integer feasible.";
			double incumbentValue = subproblemBounded.getPrimalBound();
			if (incumbentValue < primalBound) {
				output += "\nSolution defines new best primal bound.";
				primalBound = incumbentValue;
				solutionValues = subproblemBounded.getNominalVariablesSolutionValues();
			}
			writeOutput(output);
		}
		
		//Tries to use the new dual bound for pruning the node.
		if (isOptimal(primalBound, subproblemBounded.getDualBound())) {
			output = "\nPrune node, as dual bound is within optimality gap.";
			writeOutput(output);
			minimumDualBoundPrunedNodes = Math.min(minimumDualBoundPrunedNodes, subproblemBounded.getDualBound());
			numberRemainingPossibleZ -= chosenNode.getPossibleZs().size();
		}
		//If the node cannot be pruned then we check whether the node will be branched or solved as an integer robust subproblem.
		else {
			//Computes the optimal bilinear z for the solution and projects it into the interval for z defined by the node.
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
			}
			//Otherwise we branch the node according to the chosen option.
			else if (bnbStrategies.getBranchingStrategy() == BnBStrategies.BnBBranchingStrategy.BRANCHAROUND_Z) {
				branchNodeAroundZ(chosenNode);
			}
			else if (bnbStrategies.getBranchingStrategy() == BnBStrategies.BnBBranchingStrategy.BRANCHAROUND_CONVEXCOMB){
				branchNodeConvexcombination(chosenNode);
			}
		}
	}
	
	/**
	 * Branches a node around a given branching point. Sets a given value as dual bound for the child nodes.
	 */
	private void branchNode(BnBNode chosenNode, double branchingPoint, double dualBoundLP) throws IOException {		
		//Splits the interval of possible z around the branching point.
		List<PossibleZ> lowerPossibleZs = new ArrayList<PossibleZ>(chosenNode.getPossibleZs().size());
		int index = 0;
		while (chosenNode.getPossibleZs().get(index).getValue() <= branchingPoint) {
			lowerPossibleZs.add(chosenNode.getPossibleZs().get(index));
			index++;
		}
		BnBNode lowerNode = new BnBNode(lowerPossibleZs, dualBoundLP);

		List<PossibleZ> upperPossibleZs = new ArrayList<PossibleZ>(chosenNode.getPossibleZs().size()-lowerPossibleZs.size());
		while (index < chosenNode.getPossibleZs().size()) {
			upperPossibleZs.add(chosenNode.getPossibleZs().get(index));
			index++;
		}
		BnBNode upperNode = new BnBNode(upperPossibleZs, dualBoundLP);
		String output = "New Nodes: "+lowerNode+" "+upperNode;
		writeOutput(output);
		
		//Adds the corresponding new nodes to the set of remaining nodes.
		remainingNodes.add(upperNode);
		remainingNodes.add(lowerNode);
	}
	
	/**
	 * Branches a node around the mean value of the interval. Sets a given value as dual bound for the child nodes.
	 */
	private void branchNodeMiddle(BnBNode chosenNode, double dualBoundLP) throws IOException {
		double branchingPoint = (chosenNode.getLowerBoundZ().getValue()+chosenNode.getUpperBoundZ().getValue())/2;
		branchNode(chosenNode, branchingPoint, dualBoundLP);
	}
	
	/**
	 * Branches a node around the solution value of z.
	 */
	private void branchNodeAroundZ(BnBNode chosenNode) throws IOException {
		double branchingPoint = subproblemBounded.getComputedZ();
		branchNode(chosenNode, branchingPoint, subproblemBounded.getDualBound());
	}
	
	/**
	 * Branches a node around a convex combination the solution value of z and the mean of the interval.
	 */
	private void branchNodeConvexcombination(BnBNode chosenNode) throws IOException {
		//Sets the branching point to be the mean of the interval.
		double branchingPoint = (chosenNode.getLowerBoundZ().getValue()+chosenNode.getUpperBoundZ().getValue())/2;
		
		//While the branching does not cut-off the current optimal solution,
		//we adjust the branching point towards the solution value of z.
		//We obtain a convex combination of the mean and the solution value as we enter the loop at least once.
		boolean branchingHasEffect = false;
		//Adjusts the branching point until branching is effective.
		while (!branchingHasEffect) {
			//Chooses the middle between the current branching point and the solution value of z.
			branchingPoint = (subproblemBounded.getComputedZ()+branchingPoint)/2;
			
			//Determines the new bounds on z around the branching point.
			int LBindex = 0;
			while (chosenNode.getPossibleZs().get(LBindex+1).getValue() <= branchingPoint) {
				LBindex++;
			}
			double newUpperBound = chosenNode.getPossibleZs().get(LBindex).getValue();
			double newLowerBound = chosenNode.getPossibleZs().get(LBindex+1).getValue();
			
			//If the new upper bound is not greater than the current solution value then the branching is effective in theory.
			//Necessary to prevent infinity loops occurring because of numerical issues. 
			if (newUpperBound <= subproblemBounded.getComputedZ()) {
				branchingHasEffect = true;
			}
			//Tests if the branching is effective, i.e., the current solution will be cut-off.
			else {
				branchingHasEffect = subproblemBounded.isBranchingEffective(subproblemBounded.getGamma(), newUpperBound, newLowerBound, bnbStrategies.getLpOptimalityCutsStrategy());
			}
		}
		String output = "Branch around "+branchingPoint;
		writeOutput(output);
		
		branchNode(chosenNode, branchingPoint, subproblemBounded.getDualBound());
	}
	
	/**
	 * Tries to prune the possible values for z in a node using the primal bound and the individual dual bounds.
	 * Returns a set of pruned values. 
	 */
	private List<PossibleZ> performPruning(BnBNode chosenNode) throws IOException {
		List<PossibleZ> prunedPossibleZs = new ArrayList<PossibleZ>();
		for (int zIndex = chosenNode.getPossibleZs().size()-1; zIndex >= 0; zIndex--) {
			PossibleZ possibleZ = chosenNode.getPossibleZs().get(zIndex);
			if (isOptimal(primalBound, possibleZ.getDualBound())) {
				minimumDualBoundPrunedNodes = Math.min(minimumDualBoundPrunedNodes, possibleZ.getDualBound());
				chosenNode.getPossibleZs().remove(zIndex);
				prunedPossibleZs.add(possibleZ);
				numberRemainingPossibleZ--;
			}
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
			IMPROVED_ESTIMATORS,
			HRS_ESTIMATORS,
			NO_ESTIMATORS;
		}
		
		/**
		 * Specifies whether we branch around the linear relaxation value of z or
		 * whether we branch around a convex combination of z and the middle of the interval.
		 */
		public enum BnBBranchingStrategy{
			BRANCHAROUND_CONVEXCOMB,
			BRANCHAROUND_Z;
		}
		
		/**
		 * Enum type specifying whether we use optimality-cuts.
		 */
		public enum BnBLPOptimalityCutsStrategy {
			NO_CUTS,
			CUTS;
		}
		
		/**
		 * Enum type specifying whether we use optimality-cuts.
		 */
		public enum BnBIntegerOptimalityCutsStrategy {
			NO_CUTS,
			CUTS,
		}

		private BnBEstimatorStrategy estimatorStrategy;
		private BnBBranchingStrategy branchingStrategy;
		private BnBLPOptimalityCutsStrategy lpOptimalityCutsStrategy;
		private BnBIntegerOptimalityCutsStrategy intOptimalityCutsStrategy;
		
		public BnBStrategies() {
			super();
			improvingZStrategy = BoundingImprovingZStrategy.IMPROVE_Z;
			terminationStrategy = BoundingTerminationStrategy.TERMINATE;
			estimatorStrategy = BnBEstimatorStrategy.IMPROVED_ESTIMATORS;
			branchingStrategy = BnBBranchingStrategy.BRANCHAROUND_CONVEXCOMB;
			lpOptimalityCutsStrategy = BnBLPOptimalityCutsStrategy.NO_CUTS;
			intOptimalityCutsStrategy = BnBIntegerOptimalityCutsStrategy.CUTS;
		}
		
		public BnBEstimatorStrategy getEstimatorStrategy() {
			return estimatorStrategy;
		}
		public void setEstimatorStrategy(BnBEstimatorStrategy estimatorStrategy) {
			this.estimatorStrategy = estimatorStrategy;
		}
		public BnBLPOptimalityCutsStrategy getLpOptimalityCutsStrategy() {
			return lpOptimalityCutsStrategy;
		}
		public void setLpOptimalityCutsStrategy(BnBLPOptimalityCutsStrategy lpOptimalityCutsStrategy) {
			this.lpOptimalityCutsStrategy = lpOptimalityCutsStrategy;
		}
		public BnBIntegerOptimalityCutsStrategy getIntOptimalityCutsStrategy() {
			return intOptimalityCutsStrategy;
		}
		public void setIntOptimalityCutsStrategy(BnBIntegerOptimalityCutsStrategy intOptimalityCutsStrategy) {
			this.intOptimalityCutsStrategy = intOptimalityCutsStrategy;
		}
		public BnBBranchingStrategy getBranchingStrategy() {
			return branchingStrategy;
		}
		public void setBranchingStrategy(BnBBranchingStrategy branchingStrategy) {
			this.branchingStrategy = branchingStrategy;
		}
	}
}
