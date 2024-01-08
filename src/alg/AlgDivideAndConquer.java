package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import alg.AlgDivideAndConquer.DnCStrategies.DnCOptimalityCutsStrategy;
import alg.SubproblemNominalGurobi.MasterCallback;
import alg.SubproblemNominalGurobi.NOSStrategies;
import alg.SubproblemNominalGurobi.NOSStrategies.NOSTerminationStrategy;
import gurobi.GRBException;
import util.CliquePartitioning;
import util.ConflictGraph;
import util.DnCNode;
import util.PossibleZ;
import util.Variable;

/**
 * This class implements an enhanced version of the Divide and Conquer of Hansknecht, Richter, Stiller 
 * from "Fast robust shortest path computations".  The enhancements are as described
 * in "A Branch and Bound Algorithm for Robust Binary Optimization with Budget Uncertainty".
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
		if (dncStrategies.getCliqueStrategy() == RobustAlgorithmStrategies.CliqueStrategy.CLIQUES_ENABLE) {
			conflictGraph = new ConflictGraph(subproblemNominal.getModel(), subproblemNominal.getUncertainModelVariables(), algorithmParameters);
			//Computes cliques, which are stored with the variables for computing better estimators
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
			//Terminates the algorithm if the time limit is reached
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
		}
		
		//Prunes values of z from the node 
		prunePossibleZinNodes();
		
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
			
			//Sets the dual bound as the minimum dual bound of all pruned and remaining nodes.
			updateGlobalDualBound();
			primalDualIntegral.update(true);
			
			output = "\n##### DnC Tree Information\n"
					+ "Elapsed Time = "+(System.nanoTime()-startTime)/Math.pow(10, 9)+"\n"
					+ "Number Remaining Nodes = "+remainingNodes.size()
					+ "\nNumber Remaining Possible Z = "+numberRemainingPossibleZ
					+ "\nCurrent Primal Bound = "+getPrimalBound()
					+ "\nCurrent Dual Bound = "+getDualBound()
					+ "\nRelative Optimality Gap = "+(getRelativeGap()*100)+"%"
					+ "\nCurrent Primal-Dual Integral = "+getPrimalDualIntegral();
			writeOutput(output);
			
			//We select and remove the first node in our tree.
			DnCNode currentNode = remainingNodes.pollFirst();
			
			//We choose the middle value of the possible values of z.
			int indexZ = currentNode.getPossibleZs().size()/2;
			PossibleZ chosenPossibleZ = currentNode.getPossibleZs().get(indexZ);
			
			//Solve the nominal subproblem for the chosen z.
			solveNOS(currentNode, chosenPossibleZ);
			
			//The chosen z doesn't need to be considered any longer.
			numberRemainingPossibleZ--;
			
			//Splits the interval of possible z into two subsets, excluding the already considered middle value of z.
			//If the subsets are not empty, we add new nodes.
			output = "\n#####Divide Node";
			if (indexZ > 0) {
				List<PossibleZ> possibleZList = new ArrayList<PossibleZ>();
				for (int index = 0; index < indexZ; index++) {
					possibleZList.add(currentNode.getPossibleZs().get(index));
				}
				DnCNode newNode = new DnCNode(possibleZList);
				//Sets the primal bounds of the values of z defining the boundaries of the node that are used for the sorting of nodes.
				newNode.setLowerPrimalBound(currentNode.getLowerPrimalBound());
				newNode.setUpperPrimalBound(subproblemNominal.getPrimalBound());
				remainingNodes.add(newNode);
				output += "\nAdded new node "+newNode;
			}
			if (indexZ+1 < currentNode.getPossibleZs().size()) {
				List<PossibleZ> possibleZList = new ArrayList<PossibleZ>();
				for (int index = indexZ+1; index < currentNode.getPossibleZs().size(); index++) {
					possibleZList.add(currentNode.getPossibleZs().get(index));
				}
				DnCNode newNode = new DnCNode(possibleZList);
				//Sets the primal bounds of the values of z defining the boundaries of the node that are used for the sorting of nodes.
				newNode.setLowerPrimalBound(subproblemNominal.getPrimalBound());
				newNode.setUpperPrimalBound(currentNode.getUpperPrimalBound());
				remainingNodes.add(newNode);
				output += "\nAdded new node "+newNode;
			}
			writeOutput(output);
			
			//Prunes whole nodes w.r.t. to their dual bound if possible
			pruneNodes();
			//Prunes single possible values for z within nodes w.r.t. their individual dual bound if possible
			prunePossibleZinNodes();
		}
		
		//Updates the dual bound to be the minimum of all pruned nodes and
		//the minimum dual bound of the not yet considered nodes.
		updateGlobalDualBound();
		
		//Stores best solution found, if available
		if (solution != null) {
			solution = new LinkedHashMap<Variable, Double>(solutionValues.length);
			for (int i = 0; i < solutionValues.length; i++) {
				solution.put(subproblemNominal.getNominalVariables()[i], solutionValues[i]);
			}
		}
	}

	/**
	 * Solves the nominal subproblem for the chosen z, contained in the chosen node.
	 */
	private void solveNOS(DnCNode currentNode, PossibleZ chosenPossibleZ) throws IOException, GRBException {
		String output = "\n##### Solving NOS(z) for z="+chosenPossibleZ;
		writeOutput(output);
		
		//Alters the nominal subproblem with respect to z.
		subproblemNominal.updateZ(chosenPossibleZ);
		
		//Determines the set of remaining values for z for which we compute estimators.
		//If we use optimality cuts, then these are defined by the lowest and highest value in the node.
		//Then we only compute estimators for the values of z in the node.
		//If we do not use optimality cuts then we compute estimators for all remaining values for z.
		TreeSet<PossibleZ> restrictedRemainingPossibleZs = new TreeSet<PossibleZ>();
		if (dncStrategies.getOptimalityCutsStrategy() == DnCOptimalityCutsStrategy.IPOPTCUTS_ENABLE) {
			restrictedRemainingPossibleZs.addAll(currentNode.getPossibleZs());
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
		if (dncStrategies.getEstimatorStrategy() == DnCStrategies.DnCEstimatorStrategy.ESTIMATORS_IMPROVED) {
			chosenPossibleZ.setImprovedEstimators(restrictedRemainingPossibleZs, subproblemNominal.getGamma(), subproblemNominal.getUncertainVariables());
		}
		else if (dncStrategies.getEstimatorStrategy() == DnCStrategies.DnCEstimatorStrategy.ESTIMATORS_HRS) {
			chosenPossibleZ.setHRSEstimators(restrictedRemainingPossibleZs, subproblemNominal.getGamma());
		}
		writeOutput("");
		
		//Initializes the callback used for reporting primal and dual bounds from the subproblem to the master.
		subproblemNominal.setMasterCallback(new DnCCallbackIntegerSubproblems(currentNode, chosenPossibleZ));
		//Tries to solve the nominal subproblem within the remaining time.
		subproblemNominal.solve(getRemainingTime());
		
		//Updates the dual bound of the chosen value for z.
		chosenPossibleZ.updateDualBound(subproblemNominal.getDualBound());
		minimumDualBoundConsideredSubproblems = Math.min(chosenPossibleZ.getDualBound(), minimumDualBoundConsideredSubproblems);
	}
	

	
	/**
	 * Callback passed to the integer subproblems.
	 * Reports bounds from the integer subproblem to the master and asks for termination.
	 */
	protected class DnCCallbackIntegerSubproblems implements MasterCallback{
		/**
		 * The node for which we solve the integer subproblem.
		 */
		private DnCNode currentNode;
		
		/**
		 * Value of z for which we compute the nominal subproblem.
		 */
		private PossibleZ chosenPossibleZ;
		
		/**
		 * Saves the last time stamp at which we were able to prune a possible value of z via estimators. 
		 */
		private Long timeStampLastPrune = null;
		
		/**
		 * Remaining values of z that are not yet pruned for which we have computed an estimator.
		 */
		private ArrayList<PossibleZ> remainingComparedPossibleZs;

		/**
		 * Constructor receiving the chosen node and its surrounding nodes.
		 */
		private DnCCallbackIntegerSubproblems(DnCNode currentNode, PossibleZ chosenPossibleZ) {
			remainingComparedPossibleZs = new ArrayList<PossibleZ>(chosenPossibleZ.getEstimators().keySet());
			this.currentNode = currentNode;
			this.chosenPossibleZ = chosenPossibleZ;
		}
		
		/**
		 * Receives primal and dual bounds from the subproblem to update the bounds in the master problem.
		 */
		public boolean updatePrimalDualBounds(double primalBound, double subproblemDualBound) {
			boolean boundChanged = false;
			
			//Updates the primal bound
			if (primalBound < getPrimalBound()) {
				setPrimalBound(primalBound);
				solutionValues = subproblemNominal.getImprovedNominalSolutionValues();
				boundChanged = true;
				
				try {
					writeOutput("###Found new global primal bound = "+getPrimalBound());
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			//Updates dual bounds of possible z contained in this node and potentially surrounding nodes via estimators.
			if (chosenPossibleZ.getDualBound() < subproblemDualBound) {
				chosenPossibleZ.updateDualBound(subproblemDualBound);
				boundChanged = true;
				
				Set<DnCNode> estimatedNodes = new HashSet<>();
				estimatedNodes.add(currentNode);
				if (dncStrategies.optimalityCutsStrategy == DnCOptimalityCutsStrategy.IPOPTCUTS_DISABLE) {
					estimatedNodes.addAll(remainingNodes);
				}

				for (DnCNode node : estimatedNodes) {
					boolean individualBoundChanged = false;
					for (PossibleZ comparedPossibleZ : node.getPossibleZs()) {
						if (chosenPossibleZ.getEstimators().containsKey(comparedPossibleZ)) {
							double estimator = chosenPossibleZ.getEstimators().get(comparedPossibleZ);
							if (comparedPossibleZ.updateDualBound(subproblemDualBound - estimator)) {
								individualBoundChanged = true;
							}
						}
					}
					if (individualBoundChanged) {
						node.updateDualBound();
					}
				}
				
				updateGlobalDualBound(currentNode);
			}
			return boundChanged;
		}
		
		/**
		 * Receives primal and dual bounds for updating and returns whether the subproblem may be terminated.
		 */
		public boolean updateBoundsAndDecideTermination(double improvedPrimalBound, double subproblemPrimalBound, double subproblemDualBound) {
			boolean boundChanged = this.updatePrimalDualBounds(improvedPrimalBound, subproblemDualBound);
			
			//Checks whether the subproblem can be terminated with respect to the global primal bound and its current dual bound.
			if (isOptimal(getPrimalBound(), subproblemDualBound)) {
				
				if (dncStrategies.terminationStrategy == NOSTerminationStrategy.TERMINATION_DIRECT) {
					return true;
				}
				else if (dncStrategies.terminationStrategy == NOSTerminationStrategy.TERMINATION_ESTIMATORS) {
					//Initializes the time stamp and remaining possible z for pruning.
					if (timeStampLastPrune == null) {
						timeStampLastPrune = System.nanoTime();
					}
					//Checks whether we can hope for pruning additional possible z.
					boolean pruningPossible = false;
					if (boundChanged) {
						for (int i = remainingComparedPossibleZs.size()-1; i >= 0; i--) {
							PossibleZ possibleZ = remainingComparedPossibleZs.get(i);
							//Checks whether the dual bound is good enough for pruning.
							if (isOptimal(getPrimalBound(), possibleZ.getDualBound())) {
								//Removes the possible z and sets the time stamp.
								remainingComparedPossibleZs.remove(i);
								timeStampLastPrune = System.nanoTime();
							}
							//If raising the dual bound up to the primal bound would be sufficient
							//to prune the possible z then we can hope for pruning it later and
							//don't terminate the subproblem if the last time stamp is not too far in the past.
							else if (isOptimal(getPrimalBound(), subproblemPrimalBound - chosenPossibleZ.getEstimators().get(possibleZ))) {
								pruningPossible = true;
							}
						}
					}
					//If the bounds didn't change, then pruning was possible before and is still now
					else {
						pruningPossible= true;
					}
					//If there is hope for pruning further possible z and the last pruning
					//has happened in the last 10 seconds, then we continue solving the subproblem.
					//Otherwise we return true and terminate.
					if (!pruningPossible || (System.nanoTime()-timeStampLastPrune)/Math.pow(10, 9) >= 10) {
						return true;
					}
				}
			}
			return false;
		}
	}
	
	/**
	 * Updates the dual bound to be the minimum of all pruned nodes, the current node, and the minimum dual bound of the not yet considered nodes.
	 */
	private void updateGlobalDualBound(DnCNode currentNode) {
		double newDualBoud = Math.min(minimumDualBoundConsideredSubproblems, currentNode.getDualBound());
		for (DnCNode node : remainingNodes) {
			newDualBoud = Math.min(newDualBoud, node.getDualBound());
		}
		setDualBound(newDualBoud);
	}
	
	/**
	 * Updates the dual bound to be the minimum of all pruned nodes and the minimum dual bound of the not yet considered nodes.
	 */
	private void updateGlobalDualBound() {
		double newDualBoud = minimumDualBoundConsideredSubproblems;
		for (DnCNode node : remainingNodes) {
			newDualBoud = Math.min(newDualBoud, node.getDualBound());
		}
		setDualBound(newDualBoud);
	}
	
	/**
	 * Prunes all nodes whose dual bound is close enough to the current primal bound.
	 * Returns a list of all pruned nodes.
	 */
	private List<DnCNode> pruneNodes() throws IOException {
		List<DnCNode> prunedNodes = new ArrayList<DnCNode>();
		String output = "Pruned Node";
		Iterator<DnCNode> nodeIterator = remainingNodes.iterator();
		while (nodeIterator.hasNext()) {
			DnCNode node = nodeIterator.next();
			if (isOptimal(getPrimalBound(), node.getDualBound())) {
				output += " "+node;
				nodeIterator.remove();
				minimumDualBoundConsideredSubproblems = Math.min(minimumDualBoundConsideredSubproblems, node.getDualBound());
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
		for (DnCNode node : remainingNodes) {
			for (int zIndex = node.getPossibleZs().size()-1; zIndex >= 0; zIndex--) {
				PossibleZ possibleZ = node.getPossibleZs().get(zIndex);
				if (isOptimal(getPrimalBound(), possibleZ.getDualBound())) {
					output += " "+possibleZ;
					minimumDualBoundConsideredSubproblems = Math.min(minimumDualBoundConsideredSubproblems, possibleZ.getDualBound());
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
	 * Specifies all strategies for the DnC algorithm;
	 */
	public static class DnCStrategies extends NOSStrategies{
		/**
		 * Enum type specifying whether we use the estimators of Hansknecht, Richter, Stiller to obtain
		 * dual bounds from one nominal subproblem to another or whether we use the improved estimators.
		 */
		public enum DnCEstimatorStrategy{
			ESTIMATORS_IMPROVED,
			ESTIMATORS_HRS;
		}
		
		/**
		 * Enum type specifying whether we use optimality-cuts.
		 */
		public enum DnCOptimalityCutsStrategy {
			IPOPTCUTS_ENABLE,
			IPOPTCUTS_DISABLE;
		}
		
		DnCEstimatorStrategy estimatorStrategy;
		DnCOptimalityCutsStrategy optimalityCutsStrategy;
		
		/**
		 * Constructor obtaining arguments which are matched to the enums defining strategies.
		 */
		public DnCStrategies(List<String> argList, AlgorithmParameters algorithmParameters) throws IOException {
			super(argList, algorithmParameters);
		}

		@Override
		void setDefaultStrategies() {
			super.setDefaultStrategies();
			improvingZStrategy = ImprovingZStrategy.IMPROVINGZ_ENABLE;
			terminationStrategy = NOSTerminationStrategy.TERMINATION_ESTIMATORS;
			estimatorStrategy = DnCEstimatorStrategy.ESTIMATORS_IMPROVED;
			optimalityCutsStrategy = DnCOptimalityCutsStrategy.IPOPTCUTS_ENABLE;
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