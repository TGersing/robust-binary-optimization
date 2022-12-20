package util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import alg.AbstractAlgorithm;
import alg.AbstractAlgorithm.AlgorithmParameters;
import gurobi.GRBVar;
import util.ConflictGraph.HyperEdge;

/**
 * Models a partitioning of uncertain variables into cliques within the conflict graph.
 * 
 * @author Timo Gersing
 */
public class CliquePartitioning {
	/**
	 * The uncertain variables for which we define the cliques.
	 */
	private Variable[] uncertainVariables;
	
	/**
	 * A Clique is modeled to be a list of indices corresponding to the array of uncertain variables.
	 * The partitioning is modeled as a list of cliques.
	 */
	private List<List<Integer>> cliques;
		
	/**
	 * Initializes and constructs a partitioning of the array of uncertain variables using a given conflict graph.
	 */
	public CliquePartitioning(Variable[] uncertainVariables, ConflictGraph conflictGraph, AlgorithmParameters algorithmParameters) throws IOException {
		this.uncertainVariables = uncertainVariables;
		this.cliques = new ArrayList<List<Integer>>(uncertainVariables.length);
		this.partitionCliques(uncertainVariables, conflictGraph, algorithmParameters);
	}
	
	/**
	 * Returns the list of cliques.
	 */
	public List<List<Integer>> getCliques() {
		return cliques;
	}
	
	/**
	 * Adds a clique to the partitioning and sets the clique for the included variables.
	 */
	public void addClique(List<Integer> clique) {
		cliques.add(clique);
		for (Integer index : clique) {
			uncertainVariables[index].setClique(clique);
		}
	}
	
	/**
	 * Partitions an array of variables into cliques for a given conflict graph.
	 */
	private void partitionCliques(Variable[] uncertainVariables, ConflictGraph conflictGraph, AlgorithmParameters algorithmParameters) throws IOException {
		long startTimePartitioning = System.nanoTime();
		String output = "\nStart Partitioning Variables into Cliques";
		AbstractAlgorithm.writeOutput(output, algorithmParameters);
		
		//The nodes in the conflict graph are defined on the variables of the model and have to be mapped to their index.
		Map<GRBVar, Integer> modelVariableToIndex = new HashMap<GRBVar, Integer>();
		for (int i = 0; i < uncertainVariables.length; i++) {
			modelVariableToIndex.put(uncertainVariables[i].getModelVariable(), i);
		}
		
		//Initialize the set of not yet merged variables.
		Set<GRBVar> remainingVariables = new LinkedHashSet<GRBVar>(uncertainVariables.length);
		for (Variable variable : uncertainVariables) {
			remainingVariables.add(variable.getModelVariable());
		}
		
		//Merge variables while there are remaining variables.
		while (!remainingVariables.isEmpty()) {
			//The first remaining variable is chosen to be the pivot variable around which we build the next clique.
			GRBVar pivotCliqueMember = remainingVariables.iterator().next();
			List<GRBVar> modelVarClique = new ArrayList<GRBVar>();
			
			//If the pivot variable is part of a hyper edge in the conflict graph then this already forms a clique.
			//We choose the largest hyper edge, only considering remaining variables, in which the pivot variable is contained.
			//This saves time since we do not have to test whether the added variables are pairwise neighbors.
			List<GRBVar> largestHyperEdge = null;
			for (HyperEdge hyperEdge : conflictGraph.getIncidentHyperEdges(pivotCliqueMember)) {
				List<GRBVar> remainingHyperEdge = new ArrayList<GRBVar>(hyperEdge.getVariables());
				remainingHyperEdge.retainAll(remainingVariables);
				if (largestHyperEdge == null || remainingHyperEdge.size() > largestHyperEdge.size()) {
					largestHyperEdge = remainingHyperEdge;
				}
			}
			if (largestHyperEdge == null) {
				modelVarClique.add(pivotCliqueMember);
			}
			else {
				modelVarClique.addAll(largestHyperEdge);
			}
			
			//All potential clique members are in the neighborhood of the pivot variable, intersected with the set of remaining variables.
			Set<GRBVar> potentialCliqueMembers = conflictGraph.getNeighbors(pivotCliqueMember);
			potentialCliqueMembers.removeAll(modelVarClique);
			potentialCliqueMembers.retainAll(remainingVariables);
			
			//All potential clique members must be neighbors of the already added clique members.
			for (GRBVar cliqueMember : modelVarClique) {
				if (potentialCliqueMembers.isEmpty()) {
					break;
				}
				if (cliqueMember == pivotCliqueMember) {
					continue;
				}
				potentialCliqueMembers.retainAll(conflictGraph.getNeighbors(cliqueMember));
			}
			
			//While there exists a potential member, which is guaranteed to be a neighbor of all clique members,
			//we add the new member and update the set of potential members.
			while (!potentialCliqueMembers.isEmpty()) {
				GRBVar cliqueMember = potentialCliqueMembers.iterator().next();
				modelVarClique.add(cliqueMember);
				potentialCliqueMembers.remove(cliqueMember);
				potentialCliqueMembers.retainAll(conflictGraph.getNeighbors(cliqueMember));
			}
			
			//The clique is removed from the set of remaining variables and transformed into a list of indices.
			remainingVariables.removeAll(modelVarClique);
			List<Integer> clique = new ArrayList<Integer>(modelVarClique.size());
			for (GRBVar modelVar : modelVarClique) {
				clique.add(modelVariableToIndex.get(modelVar));
			}
			this.addClique(clique);
		}
		output = "Finished Partitioning within "+(Math.pow(10, -9)*(System.nanoTime()-startTimePartitioning))+" sec"
				+"\nMerged "+uncertainVariables.length+" Uncertain Variables in "+this.cliques.size()+" Cliques.";
		AbstractAlgorithm.writeOutput(output, algorithmParameters);
	}
}
