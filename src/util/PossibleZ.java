package util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import alg.AbstractAlgorithm;

/**
 * This class models a possibly optimal value for z. An instance contains the value of z,
 * a dual bound on the objective value of all solutions fulfilling the optimality criterion from the paper,
 * and a map of estimators which we use to obtain dual bounds for other possible values for z.
 * Possible z can be compared with respect to their value.
 * 
 * @author Timo Gersing
 */
public class PossibleZ implements Comparable<PossibleZ>{
	/**
	 * The possible value of z. Equals the deviation of a variable or 0.
	 */
	private double value;
	
	/**
	 * A dual bound on the objective value of all solutions fulfilling the optimality criterion from the paper.
	 */
	private double dualBound = AbstractAlgorithm.DEFAULT_DUAL_BOUND;
	
	/**
	 * A map of estimators which we use to obtain dual bounds for other possible values for z.
	 */
	private Map<PossibleZ, Double> estimators = new HashMap<PossibleZ, Double>();
	
	/**
	 * Constructor setting the value of z.
	 */
	public PossibleZ(double z) {
		this.value = z;
	}
	
	/**
	 * Returns the possible value of z.
	 */
	public double getValue() {
		return value;
	}
	
	/**
	 * Returns the dual bound.
	 */
	public double getDualBound() {
		return dualBound;
	}
	
	/**
	 * Updates the dual bound to be the maximum of the current and a new dual bound.
	 */
	public void updateDualBound(double dualBound) {
		this.dualBound = Math.max(this.dualBound, dualBound);
	}
	
	/**
	 * Returns the map of estimators.
	 */
	public Map<PossibleZ, Double> getEstimators() {
		return estimators;
	}	
	
	/**
	 * Clears the map of estimators and then sets the estimators used by Hansknecht, Richter, Stiller
	 * for possible values for z in a given set.
	 * The estimators are only defined for possible z that are smaller than the current z.
	 */
	public void setHRSEstimators(TreeSet<PossibleZ> comparedZs, double Gamma) {
		estimators = new HashMap<PossibleZ, Double>();
		for (PossibleZ possibleZ : comparedZs.headSet(this, false)) {
			this.estimators.put(possibleZ, Gamma*(this.getValue() - possibleZ.getValue()));
		}
	}
	
	/**
	 * Clears the map of estimators and then sets improved estimators, as described in the paper,
	 * for all values for z in a given set.
	 */
	public <V> void setImprovedEstimators(TreeSet<PossibleZ> comparedZs, double Gamma, Variable[] uncertainVariables) {
		estimators = new HashMap<PossibleZ, Double>();
		computeImprovedEstimatorsAboveZ(comparedZs, uncertainVariables);
		computeImprovedEstimatorsBelowZ(comparedZs, Gamma, uncertainVariables);
	}
	
	/**
	 * Clears the map of estimators and then sets improved estimators, as described in the paper,
	 * for all values for z in a given set that are greater than the current z.
	 */
	public <V> void setImprovedEstimatorsAboveZ(TreeSet<PossibleZ> comparedZs, Variable[] uncertainVariables) {
		estimators = new HashMap<PossibleZ, Double>();
		computeImprovedEstimatorsAboveZ(comparedZs, uncertainVariables);
	}

	/**
	 * Clears the map of estimators and then sets improved estimators, as described in the paper,
	 * for all values for z in a given set that are smaller than the current z.
	 */
	public <V> void setImprovedEstimatorsBelowZ(TreeSet<PossibleZ> comparedZs, double Gamma, Variable[] uncertainVariables) {
		estimators = new HashMap<PossibleZ, Double>();
		computeImprovedEstimatorsBelowZ(comparedZs, Gamma, uncertainVariables);
	}
	
	/**
	 * Computes improved estimators, as described in the paper, for all values for z in a given set
	 * that are greater than the current z and adds them to the map of estimators.
	 */
	private <V> void computeImprovedEstimatorsAboveZ(TreeSet<PossibleZ> comparedZs, Variable[] uncertainVariables) {
		//Value of the estimator.
		double estimator = 0;
		
		//Searches for the index of the smallest variable whose deviation is greater then the value of the current z.
		//If this index does not exist then there exist no greater possible value for z.
		int variableIndex = 0;
		while (variableIndex < uncertainVariables.length && uncertainVariables[variableIndex].getDeviation() <= this.getValue()) {
			variableIndex++;
		}
		
		//Identifies a clique with the contained variable that is currently considered for the estimator. 
		Map<List<Integer>, Variable> variablesUsedForCliques = new HashMap<List<Integer>, Variable>();
		
		//Iterates in ascending order over the set of possible z that are greater than the current z.
		for (PossibleZ comparedZ : comparedZs.tailSet(this, false)) {
			//Iterates over all uncertain variables whose deviation is smaller than or equal to the value of the compared z.
			while (variableIndex < uncertainVariables.length && uncertainVariables[variableIndex].getDeviation() <= comparedZ.getValue()) {
				Variable variable = uncertainVariables[variableIndex];
				//If the variable is contained in a clique then we identify its clique with this variable
				//and potentially replace the variable with which the clique was identified before. 
				if (variable.getClique() != null && variable.getClique().size() > 1) {
					if (variablesUsedForCliques.containsKey(variable.getClique())) {
						estimator -= variablesUsedForCliques.get(variable.getClique()).getDeviation() - this.getValue();
					}
					variablesUsedForCliques.put(variable.getClique(), variable);
				}
				//Adds the difference in the variables deviation and the current possible z to the estimator.
				estimator += uncertainVariables[variableIndex].getDeviation() - this.getValue();
				variableIndex++;
			}
			//Adds the estimator to the map.
			this.estimators.put(comparedZ, estimator);
		}
	}
	
	/**
	 * Computes improved estimators, as described in the paper, for all values for z in a given set
	 * that are smaller than the current z and adds them to the map of estimators.
	 */
	private <V> void computeImprovedEstimatorsBelowZ(TreeSet<PossibleZ> comparedZs, double Gamma, Variable[] uncertainVariables) {
		//Value of the estimator.
		double dualBoundEstimator = 0;
		
		//Searches for the index of the greatest variable whose deviation is smaller then the value of the current z.
		//If this index does not exist then there exist no greater possible value for z.
		int variableIndex = uncertainVariables.length-1;
		while (variableIndex >= 0 && uncertainVariables[variableIndex].getDeviation() >= this.getValue()) {
			variableIndex--;
		}
		
		//Identifies a clique with the contained variable that is currently considered for the estimator. 
		Map<List<Integer>, Variable> variablesUsedForCliques = new HashMap<List<Integer>, Variable>();
		
		//Keeps track of the variables that are currently considered in non-ascending order w.r.t. their deviation.
		List<Variable> variablesUsed = new ArrayList<Variable>();

		//Iterates in descending order over the set of possible z that are smaller than the current z.
		Iterator<PossibleZ> descendingIterator = comparedZs.headSet(this, false).descendingIterator();
		while (descendingIterator.hasNext()) {
			PossibleZ comparedZ = descendingIterator.next();
			//If the compared value of z is zero then the estimation is Gamma*(this.value - compared.value)=Gamma*this.value
			if (comparedZ.getValue() == 0) {
				this.estimators.put(comparedZ, Gamma*this.getValue());
				break;
			}
			//Iterates over all uncertain variables whose deviation is greater than or equal to the value of the compared z.
			while (variableIndex >= 0 && uncertainVariables[variableIndex].getDeviation() >= comparedZ.getValue()) {
				Variable variable = uncertainVariables[variableIndex];
				//If the variable is contained in a clique then we identify its clique with this variable
				//and potentially replace the variable with which the clique was identified before. 
				if (variable.getClique() != null && variable.getClique().size() > 1) {
					if (variablesUsedForCliques.containsKey(variable.getClique())) {
						Variable removableVariable = variablesUsedForCliques.get(variable.getClique());
						dualBoundEstimator -= this.getValue() - removableVariable.getDeviation();
						variablesUsed.remove(removableVariable);
					}
					variablesUsedForCliques.put(variable.getClique(), variable);
				}
				//Marks the variable to be used.
				variablesUsed.add(variable);
				//Adds the difference in the variables deviation and the current possible z to the estimator.
				dualBoundEstimator += this.getValue() - variable.getDeviation();
				variableIndex--;
			}
			//If we consider more than Gamma variables then we remove the ones with the highest deviation,
			//which are at the start of the list.
			while (variablesUsed.size() > Gamma) {
				Variable removableVariable = variablesUsed.remove(0);
				dualBoundEstimator -= this.getValue() - removableVariable.getDeviation();
				variablesUsedForCliques.remove(removableVariable.getClique());
			}
			//Adds the estimator to the map.
			this.estimators.put(comparedZ, dualBoundEstimator);
		}
	}

	/**
	 * Obtains a dual bound that is used to update the dual bounds of all possible z for which we
	 * have computed an estimator. The new potential dual bound is dualBound-estimator.
	 */
 	public void estimateDualBounds(double dualBound) {
		for (PossibleZ possibleZ : this.getEstimators().keySet()) {
			possibleZ.updateDualBound(dualBound - this.getEstimators().get(possibleZ));
		}
	}
	
 	/**
 	 * Possible z are compared with respect to their value in the natural ordering.
 	 */
	@Override
	public int compareTo(PossibleZ possibleZ) {
		return Double.compare(this.value, possibleZ.getValue());
	}
	
	/**
	 * Possible z are represented as a sting using their value.
	 */
	@Override
	public String toString() {
		return Double.toString(value);
	}
}
