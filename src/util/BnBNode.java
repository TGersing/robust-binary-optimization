package util;

import java.util.List;

import alg.SubproblemBoundedGurobi.NonBasicSlackConstraint;

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
	 * Decides whether we solve the LP relaxation with robustness constraints added or not.
	 * If the constraints are omitted, we add their violation to the objective function.
	 * In this case, solving the LP is much faster but yields a worse dual bound.
	 */
	private boolean applyLagrangeRelaxation;
	
	/**
	 * Objects representing constraints whose slack were non-basic in the parent model.
	 * Used to add the same constraints to the current model.
	 */
	private List<NonBasicSlackConstraint> parentNonBasicSlackConstraints;
	
	/**
	* Upper bounds with which the p variables have been substituted in the parent model.
	* Can vary for different variables, as some need to be unchanged for warm starting.
	*/
	private double[] parentUpperBoundsUsedForSubstituting;
	
	/**
	 * Integer array showing which nominal variables were basic in the computed solution of the parent model.
	 */
	private int[] parentVBasisNominal;
	
	/**
	 * Integer array showing which p variables were basic in the computed solution of the parent model.
	 */
	private int[] parentVBasisP;
	
	/**
	 * Integer showing whether z was basic in the computed solution of the parent model.
	 */
	private Integer parentVBasisZ;
	
	/**
	 * Integer array showing which nominal constraints were basic in the computed solution of the parent model.
	 */
	private int[] parentCBasisNominal;
	
		
	/**
	 * Constructor obtaining a list (should be sorted) of possibly optimal values for z.
	 * The dual bound is initialized to be the minimum over the individual dual bounds of all z.
	 */
	public BnBNode(List<PossibleZ> possibleZs, boolean applyLagrangeRelaxation) {
		super(possibleZs);
		this.applyLagrangeRelaxation = applyLagrangeRelaxation;
	}
	
	/**
	 * Constructor obtaining a list (should be sorted) of possibly optimal values for z.
	 * The dual bound is initialized to be the maximum of a given dual bound and
	 * the minimum over the individual dual bounds of all z.
	 */
	public BnBNode(List<PossibleZ> possibleZs, double dualBound, boolean applyLagrangeRelaxation) {
		super(possibleZs);
		updateDualBound(dualBound);
		this.applyLagrangeRelaxation = applyLagrangeRelaxation;
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
		if (solveInteger) {
			applyLagrangeRelaxation = false;
			parentNonBasicSlackConstraints = null;
		}
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

	public boolean isApplyLagrangeRelaxation() {
		return applyLagrangeRelaxation;
	}
	
	public int[] getParentVBasisNominal() {
		return parentVBasisNominal;
	}
	public void setParentVBasisNominal(int[] parentVBasisNominal) {
		this.parentVBasisNominal = parentVBasisNominal;
	}

	public int[] getParentVBasisP() {
		return parentVBasisP;
	}
	public void setParentVBasisP(int[] parentVBasisP) {
		this.parentVBasisP = parentVBasisP;
	}

	public Integer getParentVBasisZ() {
		return parentVBasisZ;
	}
	public void setParentVBasisZ(Integer parentVBasisZ) {
		this.parentVBasisZ = parentVBasisZ;
	}

	public int[] getParentCBasisNominal() {
		return parentCBasisNominal;
	}
	public void setParentCBasisNominal(int[] parentCBasisNominal) {
		this.parentCBasisNominal = parentCBasisNominal;
	}

	public List<NonBasicSlackConstraint> getParentNonBasicSlackConstraints() {
		return parentNonBasicSlackConstraints;
	}
	public void setParentNonBasicSlackConstraints(List<NonBasicSlackConstraint> parentNonBasicSlackConstraints) {
		this.parentNonBasicSlackConstraints = parentNonBasicSlackConstraints;
	}

	public double[] getParentUpperBoundsUsedForSubstituting() {
		return parentUpperBoundsUsedForSubstituting;
	}
	public void setUpperBoundsUsedForSubstituting(double[] parentUpperBoundsUsedForSubstituting) {
		this.parentUpperBoundsUsedForSubstituting = parentUpperBoundsUsedForSubstituting;
	}
	
}