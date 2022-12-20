package util;

import java.util.List;

/**
 * This class models nodes in the branching tree used for the branch and bound algorithm.
 * Nodes are compared with respect to their dual bounds.
 * 
 * @author Timo Gersing
 */
public class BnBNode extends TreeNode implements Comparable<BnBNode>{
	/**
	 * A node can be marked to be considered for an MILP the next time it is chosen.
	 * If this field is false then we solve the corresponding LP relaxation and either branch
	 * the node or set the field to true for the next time the node is called.
	 */
	private boolean solveInteger = false;
		
	/**
	 * Constructor obtaining a list (should be sorted) of possibly optimal values for z.
	 * The dual bound is initialized to be the minimum over the individual dual bounds of all z.
	 */
	public BnBNode(List<PossibleZ> possibleZs) {
		super(possibleZs);
	}
	
	/**
	 * Constructor obtaining a list (should be sorted) of possibly optimal values for z.
	 * The dual bound is initialized to be the maximum of a given dual bound and
	 * the minimum over the individual dual bounds of all z.
	 */
	public BnBNode(List<PossibleZ> possibleZs, double dualBoundLP) {
		super(possibleZs);
		updateDualBoundLP(dualBoundLP);
	}
		
	/**
	 * Updates the node's dual bound using a bound obtained from an LP relaxation.
	 * Also updates the bound for all possible z within the node.  
	 */
	public void updateDualBoundLP(double dualBoundLP) {
		if (dualBoundLP > dualBound) {
			dualBound = dualBoundLP;
			for (PossibleZ possibleZ : possibleZs) {
				possibleZ.updateDualBound(dualBoundLP);
			}
		}
	}

	/**
	 * Returns the decision whether to solve an MILP corresponding to the node or not. 
	 */
	public boolean isSolveInteger() {
		return solveInteger;
	}
	
	/**
	 * Sets the decision whether to solve an MILP corresponding to the node the next time the node is chosen. 
	 */
	public void setSolveInteger(boolean solveInteger) {
		this.solveInteger = solveInteger;
	}
	
	/**
	 * A node is smaller than another node if its dual bounds is smaller.
	 * If the dual bounds are equal, we break ties comparing the lowest possible values in the nodes.
	 */
	@Override
	public int compareTo(BnBNode node) {
		if (this.getDualBound() > node.getDualBound()) {
			return 1;
		}
		else if (this.getDualBound() == node.getDualBound()) {
			return this.getLowerBoundZ().compareTo(node.getLowerBoundZ());
		}
		else {
			return -1;
		}
	}
}