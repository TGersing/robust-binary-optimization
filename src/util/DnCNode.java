package util;

import java.util.List;

import alg.AbstractAlgorithm;

/**
 * This class models nodes in the tree used for the divide and conquer algorithm.
 * 
 * The set of remaining possible z is partitioned into intervals
 * defined by the already considered values for z as splitting points.
 * 
 * A node consists of a list containing the possible z in such an interval.
 * 
 * Nodes are compared with respect to the primal bound computed for the already considered possible z
 * defining the lower and upper bound of the node.
 * 
 * @author Timo Gersing
 */
public class DnCNode extends TreeNode implements Comparable<DnCNode>{
	/**
	 * Primal bound computed for the already considered possible z defining the lower bound of the node.
	 */
	private double lowerPrimalBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
	
	/**
	 * Primal bound computed for the already considered possible z defining the upper bound of the node.
	 */
	private double upperPrimalBound = AbstractAlgorithm.DEFAULT_PRIMAL_BOUND;
	
	/**
	 * Constructor obtaining a list (should be sorted) of possibly optimal values for z.
	 * The dual bound is initialized to be the minimum over the individual dual bounds of all z.
	 */
	public DnCNode(List<PossibleZ> possibleZs) {
		super(possibleZs);
	}
	
	/**
	 * Sets the primal bound computed for the already considered possible z
	 * defining the lower bound of the node.
	 */
	public void setLowerPrimalBound(double lowerPrimalBound) {
		this.lowerPrimalBound = lowerPrimalBound;
	}
	
	/**
	 * Returns the primal bound computed for the already considered possible z
	 * defining the lower bound of the node.
	 */
	public double getLowerPrimalBound() {
		return lowerPrimalBound;
	}

	/**
	 * Sets the primal bound computed for the already considered possible z
	 * defining the upper bound of the node.
	 */
	public void setUpperPrimalBound(double upperPrimalBound) {
		this.upperPrimalBound = upperPrimalBound;
	}

	/**
	 * Returns the primal bound computed for the already considered possible z
	 * defining the upper bound of the node.
	 */
	public double getUpperPrimalBound() {
		return upperPrimalBound;
	}
	
	/**
	 * Returns the minimum of the primal bounds computed for the already considered possible z
	 * defining the node.
	 */
	public double getPrimalBound() {
		return Math.min(lowerPrimalBound, upperPrimalBound);
	}
	
	/**
	 * A node is smaller than another node if the minimum of the primal bounds is smaller.
	 * If the primal bounds are equal, we break ties comparing the lowest possible values in the nodes.
	 */
	@Override
	public int compareTo(DnCNode node) {
		if (this.getPrimalBound() < node.getPrimalBound()) {
			return -1;
		}
		else if (this.getPrimalBound() == node.getPrimalBound()) {
			return this.getLowerBoundZ().compareTo(node.getLowerBoundZ());
		}
		else {
			return 1;
		}
	}
}