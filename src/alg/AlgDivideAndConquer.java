package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import alg.AlgDivideAndConquer.DnCStrategies.DnCOptimalityCutsStrategy;
import alg.SubproblemNominalGurobi.NOSStrategies;
import gurobi.GRBException;
import util.CliquePartitioning;
import util.ConflictGraph;
import util.DnCNode;
import util.PossibleZ;
import util.Variable;

/**
 * This class implements an enhanced version of the Divide and Conquer of Hansknecht, Richter, Stiller 
 * from "Fast robust shortest path computations".  The enhancements are as described
 * in "A Branch & Bound Algorithm for Robust Binary Optimization with Budget Uncertainty".
 * 
 * @author Timo Gersing
 */
public class AlgDivideAndConquer extends AbstractAlgorithm implements RobustAlgorithm{
	/**
	 * The nominal subproblem that we adapt for different choices of z.
	 */
	private SubproblemNominalGurobi subproblemNominal;
	
	/**
	 * Specifies strategies.
	 */
	private DnCStrategies dncStrategies;
	
	/**
	 * A sorted set of tree nodes containing yet to be considered values of z.
	 * The first node is always the one that is chosen for the next iteration.
	 * The nodes are sorted with respect to the primal bound provided by the values of z defining the boundaries of the node (interval),
	 * as described in the paper by Hansknecht, Richter, Stiller.
	 */
	private TreeSet<DnCNode> remainingNodes = new TreeSet<DnCNode>();
	
	/**
	 * Tracks the number of yet to be considered values of z.
	 */
	private int numberRemainingPossibleZ;
	
	/**
	 * The final dual bound is the minimum over the computed dual bound for all nominal subproblems.
	 */
	private double minimumDualBoundConsideredSubproblems = DEFAULT_PRIMAL_BOUND;
	
	/**
	 * Stores values of the current best solution.
	 */
	double[] solutionValues;
	
	/**
	 * Constructor receiving problem paths as well as parameters and strategies.
	 */
	public AlgDivideAndConquer(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, DnCStrategies dncStrategies) throws IOException, GRBException {
		super(algorithmParameters);
		this.subproblemNominal = new SubproblemNominalGurobi(problemPath, robustPath, algorithmParameters, dncStrategies);
		this.dncStrategies = dncStrategies;
	}
	
	/**
	 * Executes the algorithm.
	 */
	@Override
	protected void executeAlgorithm() throws GRBException, IOException {
		String output = "##########################################################\n"
				+ "##### Solving Problem via Divide and Conquer"
				+ "\n##########################################################";
		writeOutput(output);
		
		output = "\n###############################\n"
				+ "##### Start Preprocessing"
				+ "\n###############################";
		writeOutput(output);

		//Computes a conflict graph if the corresponding strategy is chosen.
		ConflictGraph conflictGraph = null;
		if (dncStrategies.getCliqueStrategy() == RobustAlgorithmStrategies.CliqueStrategy.CLIQUES) {
			conflictGraph = new ConflictGraph(subproblemNominal.getModel(), subproblemNominal.getUncertainModelVariables(), algorithmParameters);
			subproblemNominal.setCliquePartitioning(new CliquePartitioning(subproblemNominal.getUncertainVariables(), conflictGraph, algorithmParameters));
		}
		//Computes a list of possible optimal choices for z respecting the chosen filtering and clique strategies.
		List<PossibleZ> possibleZs = computePossibleZs(subproblemNominal.getUncertainVariables(), subproblemNominal.getGamma(), dncStrategies, conflictGraph, algorithmParameters);

		output = "\n###############################\n"
				+ "##### Start Divide and Conquer"
				+ "\n###############################";
		writeOutput(output);
		
		//We first solve the nominal subproblems in which z is equal to the lowest and highest value of all possible z.
		//For this, we initialize a node containing all possible z.
		//The node is needed, as the solveNOS method uses the chosen node and all remaining nodes to compute estimators for the remaining values of z and determine bounds for the optimality cuts.
		DnCNode initialNode = new DnCNode(possibleZs);
		remainingNodes.add(initialNode);
		numberRemainingPossibleZ = possibleZs.size();
		
		Set<PossibleZ> initialZs = new LinkedHashSet<PossibleZ>();
		initialZs.add(initialNode.getPossibleZs().get(0));
		initialZs.add(initialNode.getPossibleZs().get(possibleZs.size()-1));
		
		for (PossibleZ chosenZ : initialZs) {
			//Terminates the algorithm if the timelimit is reached
			if (deadline.isPresent() && System.nanoTime() > deadline.get()) {
				output = "\n###########################\n"
						+ "##### Reached Time Limit"
						+ "\n###########################";
				writeOutput(output);
				break;
			}
			//Solves the nominal subproblem for the chosen z.
			solveNOS(initialNode, chosenZ);
			//Removes the chosen z from the set of remaining possible values of z.
			initialNode.getPossibleZs().remove(chosenZ);
			numberRemainingPossibleZ--;
			//Updates the minimum dual bound computed for all subproblems.
			minimumDualBoundConsideredSubproblems = Math.min(chosenZ.getDualBound(), minimumDualBoundConsideredSubproblems);
		}
		
		//Removes the initial node if the set of remaining possible z is empty. This will terminate the algorithm.
		if (initialNode.getPossibleZs().isEmpty()) {
			remainingNodes.clear();
		}
						
		//Processes nodes while the set of nodes is not empty.
		while (!remainingNodes.isEmpty()) {
			//Terminates the algorithm if the time limit is reached.
			if (deadline.isPresent() && System.nanoTime() > deadline.get()) {
				output = "\n###########################\n"
						+ "##### Reached Time Limit"
						+ "\n###########################";
				writeOutput(output);
				break;
			}
			
			output = "\n##### DnC Tree Information\n"
					+ "Elapsed Time = "+(System.nanoTime()-startTime)/Math.pow(10, 9)+"\n"
					+ "Number Remaining Nodes = "+remainingNodes.size()
					+ "\nNumber Remaining Possible Z = "+numberRemainingPossibleZ
					+ "\nCurrent Primal Bound = "+primalBound;
			writeOutput(output);
			
			//We select the first node in our tree.
			DnCNode chosenNode = remainingNodes.first();
			remainingNodes.remove(chosenNode);
			
			//Tries to prune possible values for z in the node using their individual dual bounds.
			String chosenNodeString = chosenNode.toString();
			List<PossibleZ> prunedPossibleZs = performPruning(chosenNode);
			if (!prunedPossibleZs.isEmpty()) {
				//Prints the original node and the pruned values for z to the log.
				output = "\nChosen Node Before Individual Pruning: "+chosenNodeString
						+ "\nPruned values via estimated dual bound:\n";
				for (PossibleZ possibleZ : prunedPossibleZs) {
					output += possibleZ+"  ";
				}
				writeOutput(output);
				//The node can be pruned if there are no possible values left.
				if (chosenNode.getPossibleZs().isEmpty()) {
					output = "Node is empty after pruning";
					writeOutput(output);
					continue;
				}
				else {
					//Prints the node after pruning to the log.
					output = "Chosen Node After Individual Pruning: "+chosenNode;
					writeOutput(output);
				}
			}
			else {
				//Prints the original node to the log if we could not prune individual values.
				output = "\nChosen Node: "+chosenNode;
				writeOutput(output);
			}
			
			//If the node could not be pruned, we choose the middle value of the possible values of z.
			int indexZ = chosenNode.getPossibleZs().size()/2;
			PossibleZ chosenPossibleZ = chosenNode.getPossibleZs().get(indexZ);
			
			//Solves the nominal subproblem for the chosen z.
			solveNOS(chosenNode, chosenPossibleZ);
			numberRemainingPossibleZ--;
			
			//Splits the interval of possible z into two subsets, excluding the already considered middle value of z.
			//If the subsets are not empty, we add new nodes.
			output = "\n#####Divide Node";
			if (indexZ > 0) {
				List<PossibleZ> possibleZList = new ArrayList<PossibleZ>();
				for (int index = 0; index < indexZ; index++) {
					possibleZList.add(chosenNode.getPossibleZs().get(index));
				}
				DnCNode newNode = new DnCNode(possibleZList);
				//Sets the primal bounds of the values of z defining the boundaries of the node that are used for the sorting of nodes.
				newNode.setLowerPrimalBound(chosenNode.getLowerPrimalBound());
				newNode.setUpperPrimalBound(subproblemNominal.getPrimalBound());
				remainingNodes.add(newNode);
				output += "\nAdded new node "+newNode;
			}
			if (indexZ+1 < chosenNode.getPossibleZs().size()) {
				List<PossibleZ> possibleZList = new ArrayList<PossibleZ>();
				for (int index = indexZ+1; index < chosenNode.getPossibleZs().size(); index++) {
					possibleZList.add(chosenNode.getPossibleZs().get(index));
				}
				DnCNode newNode = new DnCNode(possibleZList);
				//Sets the primal bounds of the values of z defining the boundaries of the node that are used for the sorting of nodes.
				newNode.setLowerPrimalBound(subproblemNominal.getPrimalBound());
				newNode.setUpperPrimalBound(chosenNode.getUpperPrimalBound());
				remainingNodes.add(newNode);
				output += "\nAdded new node "+newNode;
			}
			writeOutput(output);
		}
		
		//Sets the dual bound as the minimum dual bound of all pruned and remaining nodes.
		dualBound = minimumDualBoundConsideredSubproblems;
		for (DnCNode node : remainingNodes) {
			dualBound = Math.min(dualBound, node.getDualBound());
		}
		
		//Stores best solution found
		solution = new LinkedHashMap<Variable, Double>(solutionValues.length);
		for (int i = 0; i < solutionValues.length; i++) {
			solution.put(subproblemNominal.getNominalVariables()[i], solutionValues[i]);
		}
	}

	/**
	 * Solves the nominal subproblem for the chosen z, contained in the chosen node.
	 */
	private void solveNOS(DnCNode chosenNode, PossibleZ chosenPossibleZ) throws IOException, GRBException {
		String output = "\n##### Solving NOS(z) for z="+chosenPossibleZ;
		writeOutput(output);
		
		//Alters the nominal subproblem with respect to z.
		subproblemNominal.updateZ(chosenPossibleZ);
		subproblemNominal.setGlobalPrimalBound(primalBound);
		
		//Determines the set of remaining values for z for which we compute estimators.
		//If we use optimality cuts, then these are defined by the lowest and highest value in the node.
		//Then we only compute estimators for the values of z in the node.
		//If we do not use optimality cuts then we compute estimators for all remaining values for z.
		TreeSet<PossibleZ> restrictedRemainingPossibleZs = new TreeSet<PossibleZ>();
		if (dncStrategies.getOptimalityCutsStrategy() == DnCOptimalityCutsStrategy.CUTS) {
			restrictedRemainingPossibleZs.addAll(chosenNode.getPossibleZs());
			PossibleZ lowerBoundOptCut = restrictedRemainingPossibleZs.first();
			PossibleZ upperBoundOptCut = restrictedRemainingPossibleZs.last();
			subproblemNominal.addOptimalityCuts(lowerBoundOptCut, upperBoundOptCut);
		}
		else {
			for (DnCNode node : remainingNodes) {
				restrictedRemainingPossibleZs.addAll(node.getPossibleZs());
			}
		}
		
		//Computes estimators for obtaining dual bounds on other nominal subproblems.
		if (dncStrategies.getEstimatorStrategy() == DnCStrategies.DnCEstimatorStrategy.IMPROVED_ESTIMATORS) {
			chosenPossibleZ.setImprovedEstimators(restrictedRemainingPossibleZs, subproblemNominal.getGamma(), subproblemNominal.getUncertainVariables());
		}
		else if (dncStrategies.getEstimatorStrategy() == DnCStrategies.DnCEstimatorStrategy.HRS_ESTIMATORS) {
			chosenPossibleZ.setHRSEstimators(restrictedRemainingPossibleZs, subproblemNominal.getGamma());
		}
		
		writeOutput("");
		//Tries to solve the nominal subproblem within the remaining time.
		subproblemNominal.solve(getRemainingTime());
		
		//Updates the dual bound of the chosen value for z.
		chosenPossibleZ.updateDualBound(subproblemNominal.getDualBound());
		
		//Updates the primal bound.
		double incumbentValue = subproblemNominal.getGlobalPrimalBound();
		if (incumbentValue < primalBound) {
			primalBound = incumbentValue;
			solutionValues = subproblemNominal.getIncumbentValues();
		}
		
		//Updates the individual dual bounds and the dual bounds for all nodes.
		chosenPossibleZ.estimateDualBounds(subproblemNominal.getDualBound());
		for (DnCNode node : remainingNodes) {
			node.updateDualBound();
		}
	}
	
	/**
	 * Tries to prune possible values for z in a node using the primal bound and the individual dual bounds.
	 * Returns a set of pruned values. 
	 */
	private List<PossibleZ> performPruning(DnCNode chosenNode) throws IOException {
		List<PossibleZ> prunedPossibleZs = new ArrayList<PossibleZ>();
		for (int zIndex = chosenNode.getPossibleZs().size()-1; zIndex >= 0; zIndex--) {
			PossibleZ possibleZ = chosenNode.getPossibleZs().get(zIndex);
			if (isOptimal(primalBound, possibleZ.getDualBound())) {
				//Updates the dual bound of all considered nominal subproblems.
				minimumDualBoundConsideredSubproblems = Math.min(minimumDualBoundConsideredSubproblems, possibleZ.getDualBound());
				chosenNode.getPossibleZs().remove(zIndex);
				prunedPossibleZs.add(possibleZ);
				numberRemainingPossibleZ--;
			}
		}
		return prunedPossibleZs;
	}
	

	
	/**
	 * Specifies all strategies for the DnC algorithm;
	 */
	public static class DnCStrategies extends NOSStrategies{
		/**
		 * Enum type specifying whether we use the estimators of Hansknecht, Richter, Stiller to obtain
		 * dual bounds from one nominal subproblem to another or whether we use the improved estimators.
		 */
		public enum DnCEstimatorStrategy{
			HRS_ESTIMATORS,
			IMPROVED_ESTIMATORS; 
		}
		
		/**
		 * Enum type specifying whether we use optimality-cuts.
		 */
		public enum DnCOptimalityCutsStrategy {
			NO_CUTS,
			CUTS;
		}
		
		private DnCEstimatorStrategy estimatorStrategy;
		private DnCOptimalityCutsStrategy optimalityCutsStrategy;
		
		public DnCStrategies() {
			super();
			improvingZStrategy = NOSImprovingZStrategy.IMPROVE_Z;
			terminationStrategy = NOSTerminationStrategy.TERMINATE_ESTIMATORS;
			estimatorStrategy = DnCEstimatorStrategy.IMPROVED_ESTIMATORS;
			optimalityCutsStrategy = DnCOptimalityCutsStrategy.CUTS;
		}

		public DnCEstimatorStrategy getEstimatorStrategy() {
			return estimatorStrategy;
		}
		public void setEstimatorStrategy(DnCEstimatorStrategy estimatorStrategy) {
			this.estimatorStrategy = estimatorStrategy;
		}

		public DnCOptimalityCutsStrategy getOptimalityCutsStrategy() {
			return optimalityCutsStrategy;
		}
		public void setOptimalityCutsStrategy(DnCOptimalityCutsStrategy optimalityCutsStrategy) {
			this.optimalityCutsStrategy = optimalityCutsStrategy;
		}		
	}
}