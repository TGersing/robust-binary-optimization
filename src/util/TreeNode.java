package util;

import java.util.ArrayList;
import java.util.List;

import alg.AbstractAlgorithm;

/**
 * This class is parent to the classes modeling nodes in the tree of the branch and bound and the divide and conquer.
 * 
 * @author Timo Gersing
 */
public abstract class TreeNode {
	/**
	 * A list of possibly optimal values for z.
	 */
	protected List<PossibleZ> possibleZs;

	/**
	 * Dual bound on the optimal objective value corresponding to the possible values for z.
	 */
	protected double dualBound = AbstractAlgorithm.DEFAULT_DUAL_BOUND;
	
	/**
	 * Constructor obtaining a list (should be sorted) of possibly optimal values for z.
	 * The dual bound is initialized to be the minimum over the individual dual bounds of all z.
	 */
	public TreeNode(List<PossibleZ> possibleZs) {
		this.possibleZs = new ArrayList<PossibleZ>(possibleZs);
		updateDualBound();
	}
	
	/**
	 * Returns the list of possible values for z.
	 */
	public List<PossibleZ> getPossibleZs() {
		return this.possibleZs;
	}
	
	/**
	 * Returns the lowest possible z from the (sorted) list. 
	 */
	public PossibleZ getLowerBoundZ() {
		return this.possibleZs.get(0);
	}
	
	/**
	 * Returns the highest possible z from the (sorted) list. 
	 */
	public PossibleZ getUpperBoundZ() {
		return this.possibleZs.get(possibleZs.size()-1);
	}
	
	/**
	 * Returns the node's dual bound.
	 */
	public double getDualBound() {
		return dualBound;
	}
	
	/**
	 * Updates the node's dual bound to be the maximum of the worst individual dual bound of possible z
	 * and the current dual bound, which may not only be based on individual dual bounds, but also LP values.
	 */
	public void updateDualBound() {
		double minimumIndividualDualBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
		for (PossibleZ possibleZ : possibleZs) {
			minimumIndividualDualBound = Math.min(minimumIndividualDualBound, possibleZ.getDualBound());
		}
		dualBound = Math.max(dualBound, minimumIndividualDualBound);
	}
	
	/**
	 * Writes a node to a string representing the interval defined by the possible values for z.
	 */
	@Override
	public String toString() {
		if (getLowerBoundZ() == getUpperBoundZ()) {
			return "["+getLowerBoundZ()+"]";
		}
		else {
			return "["+getLowerBoundZ()+","+getUpperBoundZ()+"]";
		}
	}
}