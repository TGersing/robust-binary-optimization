package util;

import java.util.List;

import gurobi.GRB;
import gurobi.GRB.CharAttr;
import gurobi.GRB.DoubleAttr;
import gurobi.GRBException;
import gurobi.GRBVar;

/**
 * This class stores variables of the solver's instance together with other information.
 * Variables can be compared, and thus sorted, with respect to their deviation.
 * A variable is smaller than a different variable if it's deviation is smaller.
 * 
 * @author Timo Gersing
 */
public class Variable implements Comparable<Variable>{
	/**
	 * Enum specifying the type of a variable.
	 */
	public enum Type {
		BINARY,
		INTEGER,
		CONTINUOUS;
	}

	/**
	 * Index in the list of nominal variables.
	 */
	private int nominalIndex;

	/**
	 * Index in the list of uncertain variables. Is null if the variable is not uncertain.
	 */
	private Integer uncertainIndex;

	/**
	 * Variable in the solver's instance.
	 */
	private GRBVar modelVariable;
	
	/**
	 * The variable's objective coefficient.
	 */
	private double objectiveCoefficient;
	
	/**
	 * A potential deviation of the objective coefficient
	 * (equals 0 for certain variables).
	 */
	private double deviation;
	
	/**
	 * Type of the variable.
	 */
	private Type type;
	
	/**
	 * A clique represented by the indices in the array of uncertain variables
	 * (only used for uncertain variables).
	 */
	private List<Integer> clique;
		
	/**
	 * Constructor obtaining a variable of the solver's instance and an objective coefficient.
	 * Also sets the variable's type (BINARY, INTEGER, CONTINUOUS). 
	 */
	public Variable(int index, GRBVar modelVariable, double objectiveCoefficient) throws GRBException {
		this.nominalIndex = index;
		this.modelVariable = modelVariable;
		this.objectiveCoefficient = objectiveCoefficient;
		setType();
	}

	/**
	 * Returns the index of the variable in the list of nominal variables.
	 */
	public int getNominalIndex() {
		return nominalIndex;
	}
	
	/**
	 * Returns the index of the variable in the list of uncertain variables.
	 */
	public Integer getUncertainIndex() {
		return uncertainIndex;
	}

	/**
	 * Sets the index of the variable in the list of uncertain variables.
	 */
	public void setUncertainIndex(Integer uncertainIndex) {
		this.uncertainIndex = uncertainIndex;
	}
	
	/**
	 * Returns the variable of the solver's instance.
	 */
	public GRBVar getModelVariable() {
		return modelVariable;
	}
	
	/**
	 * Returns the variable's objective coefficient.
	 */
	public double getObjectiveCoefficient() {
		return objectiveCoefficient;
	}
	
	/**
	 * Returns the variable's deviation (0 for certain variables).
	 */
	public double getDeviation() {
		return deviation;
	}
	
	/**
	 * Sets the variable's deviation.
	 */
	public void setDeviation(double deviation) {
		this.deviation = deviation;
	}
	
	/**
	 * Queries and sets the variable's type from the model.
	 */
	private void setType() throws GRBException {
		if (modelVariable.get(CharAttr.VType) == GRB.BINARY ||
				(modelVariable.get(CharAttr.VType) == GRB.INTEGER && modelVariable.get(DoubleAttr.LB) == 0 && modelVariable.get(DoubleAttr.UB) == 1)) {
			type = Type.BINARY;
		}
		else if (modelVariable.get(CharAttr.VType) == GRB.INTEGER){
			type = Type.INTEGER;
		}
		else {
			type = Type.CONTINUOUS;
		}
	}
	
	/**
	 * Returns the variable's type.
	 */
	public Type getType() {
		return this.type;
	}
	
	/**
	 * Returns the clique of the conflict graph in which the variable was merged.
	 * The clique is null if the variable is certain or if no cliques have been merged. 
	 */
	public List<Integer> getClique() {
		return this.clique;
	}
	
	/**
	 * Sets a clique of the conflict graph in which the variable is contained.
	 */
	public void setClique(List<Integer> clique) {
		this.clique = clique;
	}

	/**
	 * A variable is smaller than another variable if its deviation is smaller.
	 */
	@Override
	public int compareTo(Variable var) {
		return Double.compare(this.deviation, var.getDeviation());
	}
}