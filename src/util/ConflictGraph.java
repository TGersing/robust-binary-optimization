package util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import alg.AbstractAlgorithm;
import alg.AbstractAlgorithm.AlgorithmParameters;
import gurobi.GRB;
import gurobi.GRBConstr;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

/**
 * Models a conflict graph consisting of variables with an edge between two variables if
 * they cannot both be equal to one.
 * The graph is modeled using a combination of hyper edges, containing multiple variables that are
 * all adjacent to each other, and edges connecting two variables.
 * This structure is used, since hyper edges are more memory efficient if they contain many variables.
 * A mapping maps a variable to a list of hyper edges.
 * We also store a mapping that maps a variable to all adjacent variables connected by a simple edge.
 * 
 * @author Timo Gersing
 */
public class ConflictGraph {
	/**
	 * Maps a variable to a list of hyper edges, all containing the variable.
	 * The list may contain duplicates of hyper edges.
	 */
	private Map<GRBVar, List<HyperEdge>> incidentHyperEdges;
	
	/**
	 * Maps a variable to a list of adjacent variables, connected via a "simple edge". The list may contain duplicate variables.
	 * It may, but does not necessarily, contain variables that are within an incident hyper edge. 
	 */
	private Map<GRBVar, List<GRBVar>> adjacentVariables;
	
	/**
	 * Constructor obtaining a model with its uncertain variables and computing the conflict graph.
	 */
	public ConflictGraph(GRBModel model, GRBVar[] uncertainModelVariables, AlgorithmParameters algorithmParameters) throws GRBException, IOException {
		incidentHyperEdges = new HashMap<GRBVar, List<HyperEdge>>();
		adjacentVariables = new HashMap<GRBVar, List<GRBVar>>();
		computeConflictGraph(model, uncertainModelVariables, algorithmParameters);
	}
	
	/**
	 * Adds a hyper edge to the graph by appending it to the list of incident hyper edges
	 * for all variables contained within the edge.
	 */
	private void addHyperEdge(HyperEdge hyperEdge) {
		for (GRBVar variable : hyperEdge.getVariables()) {
			if (!incidentHyperEdges.containsKey(variable)) {
				incidentHyperEdges.put(variable, new ArrayList<HyperEdge>(1));
			}
			incidentHyperEdges.get(variable).add(hyperEdge);
		}
	}
	
	/**
	 * Adds a simple edge to the graph by adding the end nodes to the lists of adjacent variables.
	 */
	private void addEdge(GRBVar variable1, GRBVar variable2) {
		if (!adjacentVariables.containsKey(variable1)) {
			adjacentVariables.put(variable1, new ArrayList<GRBVar>(1));
		}
		adjacentVariables.get(variable1).add(variable2);
		
		if (!adjacentVariables.containsKey(variable2)) {
			adjacentVariables.put(variable2, new ArrayList<GRBVar>(1));
		}
		adjacentVariables.get(variable2).add(variable1);
	}
	
	/**
	 * Constructs and returns a set containing all neighbors of a variable,
	 * either contained in the list of adjacent variables or connected via a hyper edge.
	 * The set is a LinkedHashSet to guarantee the ordering.
	 */
	public Set<GRBVar> getNeighbors(GRBVar variable) {
		Set<GRBVar> neighbors = new LinkedHashSet<GRBVar>();
		if (adjacentVariables.containsKey(variable)) {
			neighbors.addAll(adjacentVariables.get(variable));
		}
		if (incidentHyperEdges.containsKey(variable)) {
			for (HyperEdge hyperEdge : incidentHyperEdges.get(variable)) {
				neighbors.addAll(hyperEdge.getVariables());
			}
		}
		neighbors.remove(variable);
		return neighbors;
	}
	
	/**
	 * Returns the list of incident hyper edges of a variable.
	 */
	public List<HyperEdge> getIncidentHyperEdges(GRBVar variable) {
		if (incidentHyperEdges.containsKey(variable)) {
			return incidentHyperEdges.get(variable);
		}
		else {
			return new ArrayList<HyperEdge>();
		}
	}
	
	/**
	 * Iterates over the constraints of the model to find conflicts between variables.
	 */
	private void computeConflictGraph(GRBModel model, GRBVar[] uncertainModelVariables, AlgorithmParameters algorithmParameters) throws GRBException, IOException {
		long startTimeConflictGraph = System.nanoTime();
		String output = "\nStart Computing Conflict Graph";
		AbstractAlgorithm.writeOutput(output, algorithmParameters);
		
		//Writes all uncertain variables of the model in a set such that we can quickly check whether
		//a variable is uncertain.
		Set<GRBVar> uncertainModelVariablesSet = new HashSet<GRBVar>(Arrays.asList(uncertainModelVariables));
		
		//Iterates over all constraints of the model.
		GRBConstr[] constraints = model.getConstrs();
		for (GRBConstr constraint : constraints) {
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
			GRBLinExpr lhsExpr = model.getRow(constraint);
			for (int sign : signs) {
				//The largest max{coeff,0} of an uncertain variable.
				//Useful for checking whether there exist any conflicts.
				double maxCoeff = 0;

				//We compute the maximum satisfaction of the constraint
				//by setting variables to their appropriate bound.
				double maxSatisfaction = sign*constraint.get(GRB.DoubleAttr.RHS);
				for (int i = 0; i< lhsExpr.size(); i++) {
					double coeff = sign*lhsExpr.getCoeff(i);
					//If the coefficient is positive then we maximize the satisfaction
					//by setting the variable to its lower bound.
					if (coeff > 0){
						maxSatisfaction -= coeff*lhsExpr.getVar(i).get(GRB.DoubleAttr.LB);
						if (coeff > maxCoeff && uncertainModelVariablesSet.contains(lhsExpr.getVar(i))) {
							maxCoeff = coeff;
						}
					}
					//If the coefficient is negative then we maximize the satisfaction
					//by setting the variable to its upper bound.
					else if (coeff < 0) {
						maxSatisfaction -= coeff*lhsExpr.getVar(i).get(GRB.DoubleAttr.UB);
					}
				}
				
				//maxSatisfaction should be >=0, as the problem is otherwise infeasible.
				//If an uncertain variables coefficient is smaller than zero then we already considered
				//it at its upper bound 1. Hence, setting this variable to 1 does not change the satisfaction.
				//If a variable's coefficient is >=0 then it is at its lower bound 0. Setting this variable
				//to 1 reduces the satisfaction by the coefficient.
				//Accordingly, two variables are in conflict if max{coeff1,0}+max{coeff2,0} is greater than
				//the maximum satisfaction.
				//This can only be the case if 2*maxCoeff is greater than maxSatisfaction.
				if (2*maxCoeff > maxSatisfaction+algorithmParameters.getFeasibilityTolerance()) {
					//Constructs a list of candidate variables together with their max{coeff,0}.
					//A variable is a candidate iff it has a conflict with the variable defining maxCoeff.
					//That is iff maxCoeff + max{coeff,0} is greater than maxSatisfaction.
					List<Tuple> candidates = new ArrayList<Tuple>(lhsExpr.size());
					for (int i = 0; i< lhsExpr.size(); i++) {
						double coeff = Math.max(0, sign*lhsExpr.getCoeff(i));
						if (coeff + maxCoeff > maxSatisfaction+algorithmParameters.getFeasibilityTolerance()) {
							GRBVar variable = lhsExpr.getVar(i);
							if (uncertainModelVariablesSet.contains(variable)) {
								candidates.add(new Tuple(variable, coeff));
							}
						}
					}
					//If there are only two candidates then we add the edge to the conflict graph.
					if (candidates.size() == 2) {
						addEdge(candidates.get(0).getVar(), candidates.get(1).getVar());
					}
					else if (candidates.size() > 2) {
						//Otherwise we sort the candidates non-decreasing w.r.t. coefficients.
						Collections.sort(candidates);
						//If coeff[i]+coeff[i+1] is greater than maxSatisfaction then {i,i+1,...} form a hyper edge.
						//We first search for the smallest i with that property and add the variables to a list.
						List<GRBVar> variablesOriginalHyperEdge = new ArrayList<GRBVar>();
						//The smallest index i is at most candidates.size()-2.
						variablesOriginalHyperEdge.add(candidates.get(candidates.size()-1).getVar());
						variablesOriginalHyperEdge.add(candidates.get(candidates.size()-2).getVar());
						int candidateIndex = candidates.size()-3;
						while (candidateIndex >= 0) {
							Tuple candidate = candidates.get(candidateIndex);
							if (candidate.getCoeff() + candidates.get(candidateIndex+1).getCoeff() > maxSatisfaction+algorithmParameters.getFeasibilityTolerance()) {
								variablesOriginalHyperEdge.add(candidates.get(candidateIndex).getVar());
								candidateIndex--;
							}
							else {
								break;
							}
						}
						//If the constructed hyper edge contains more than two variables then we add it to the graph.
						if (variablesOriginalHyperEdge.size() > 2) {
							addHyperEdge(new HyperEdge(variablesOriginalHyperEdge, variablesOriginalHyperEdge.size(), null));
						}
						//Otherwise we add it as a simple edge.
						else {
							addEdge(variablesOriginalHyperEdge.get(0), variablesOriginalHyperEdge.get(1));
						}
						
						//If there are still candidates that are not in the above hyper edge than they are neighbors to the
						//largest variable (index 0), but not the smallest (last index) in the hyper edge.
						//We find the last index in the hyper edge such that all variables up to this index are neighbors to the candidate.
						//Together with the original edge and this index, the variable forms a new hyper edge.
						int maxIndexSubEdge = variablesOriginalHyperEdge.size()-2;
						while (candidateIndex >= 0) {
							Tuple candidate = candidates.get(candidateIndex);
							while (candidate.getCoeff() + candidates.get(candidates.size()-1-maxIndexSubEdge).getCoeff() <= maxSatisfaction+algorithmParameters.getFeasibilityTolerance()) {
								maxIndexSubEdge--;
							}
							if (maxIndexSubEdge ==  0) {
								addEdge(candidate.getVar(), candidates.get(0).getVar());
							}
							else {
								addHyperEdge(new HyperEdge(variablesOriginalHyperEdge, maxIndexSubEdge, candidate.getVar()));
							}
							candidateIndex--;
						}
					}
				}
			}
		}
		output = "Finished Computing Conflict Graph within "+((System.nanoTime()-startTimeConflictGraph)/Math.pow(10, 9))+" sec";
		AbstractAlgorithm.writeOutput(output, algorithmParameters);
	}
	
	/**
	 * This class models hyper edges in the conflict graph.
	 * A hyper edge is usually simply a list of variables that are all neighbors in the conflict graph.
	 * For some constraints, we might add an extra variable to the hyper edge when removing other variables
	 * that are not neighbors. Since the variables in the hyper edge are sorted non-decreasing
	 * with respect to their coefficient, the new hyper edge consists of the new variable and all variables
	 * of the original hyper edge up to a certain index. Storing hyper edges this way is more memory efficient
	 * if the list of original variables is shared as a reference by all these hyper edges.
	 * 
	 */
	protected class HyperEdge {
		/**
		 * Variables in the original hyper edge.
		 */
		private List<GRBVar> variablesOriginalHyperEdge;
		
		/**
		 * The extra variable that is added by dropping other variables.
		 * This field is null if the current HyperEdge models the original one. 
		 */
		private GRBVar extraVar;

		/**
		 * Maximum index up to which the variables in the original hyper edge are
		 * neighbors to the extra variable.
		 */
		private int maxIndexSubEdge;
		
		/**
		 * Constructor setting the original variable, the maximum index of the su hyper edge and the extra variable.
		 */
		HyperEdge(List<GRBVar> variablesOriginalHyperEdge, int maxIndexSubEdge, GRBVar extraVar) {
			this.variablesOriginalHyperEdge = variablesOriginalHyperEdge;
			this.maxIndexSubEdge = maxIndexSubEdge;
			this.extraVar = extraVar;
		}
		
		/**
		 * Returns the variables of the hyper edge.
		 * If we have an extra variable, then this is added to the appropriate sub hyper edge.
		 * Otherwise, we return the variables of the original hyper edge.
		 */
		List<GRBVar> getVariables() {
			if (extraVar == null) {
				return variablesOriginalHyperEdge;
			}
			else {
				List<GRBVar> variables = new ArrayList<GRBVar>(1+variablesOriginalHyperEdge.size()-maxIndexSubEdge);
				variables.add(extraVar);
				variables.addAll(variablesOriginalHyperEdge.subList(0, maxIndexSubEdge+1));
				return variables;
			}
		}
	}

	
	/**
	 * Auxiliary class for storing variables together with their coefficient in a tuple.
	 * To allow for sorting, the tuple is comparable with respect to the coefficient.  
	 */
	private class Tuple implements Comparable<Tuple>{
		private  GRBVar var;
		private  double coeff;
		
		public Tuple(GRBVar var, double coeff) {
			this.var = var;
			this.coeff = coeff;
		}

		public GRBVar getVar() {
			return var;
		}
		public double getCoeff() {
			return coeff;
		}

		@Override
		public int compareTo(Tuple tuple) {
			return Double.compare(this.coeff, tuple.coeff);
		}
	}
}
