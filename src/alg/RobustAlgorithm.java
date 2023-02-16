package alg;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import alg.AbstractAlgorithm.AlgStrategies;
import alg.AbstractAlgorithm.AlgorithmParameters;
import gurobi.GRBVar;
import util.ConflictGraph;
import util.PossibleZ;
import util.Variable;

/**
 * The RobustAlgorithm interface is implemented by all algorithms solving robust problems.
 * 
 * @author Timo Gersing
 */
public interface RobustAlgorithm {
	/**
	 * Computes for an array of sorted variables (non-decreasing w.r.t deviation) 
	 * the index of the variable defining highest possible optimal value for z.
	 * The index can be -1, indicating that the value of z is 0.
	 */
	private int computeIndexHighestPossibleZ(Variable[] uncertainVariables, double Gamma, RobustAlgorithmStrategies robustAlgorithmStrategies) {
		int indexHighestPossibleZ = uncertainVariables.length-1;

		if (robustAlgorithmStrategies.getFilterStrategy() == RobustAlgorithmStrategies.FilterStrategy.FILTERINGZ_ENABLE) {
			if (robustAlgorithmStrategies.getCliqueStrategy() == RobustAlgorithmStrategies.CliqueStrategy.CLIQUES_ENABLE) {
				indexHighestPossibleZ = uncertainVariables.length;
				int numberCliques = 0;
				Set<List<Integer>> usedCliques = new HashSet<List<Integer>>();
				while (indexHighestPossibleZ > -1 && numberCliques <= Gamma) {
					indexHighestPossibleZ--;
					if (indexHighestPossibleZ > -1) {
						List<Integer> clique = uncertainVariables[indexHighestPossibleZ].getClique();
						if (clique == null) {
							numberCliques++;
						}
						else if (!usedCliques.contains(clique)) {
							usedCliques.add(clique);
							numberCliques++;
						}
					}
				}
			}
			else {
				indexHighestPossibleZ = uncertainVariables.length - (int)Gamma - 1;
			}
		}
		return indexHighestPossibleZ;
	}
	
	/**
	 * Computes for an array of sorted variables (non-decreasing w.r.t deviation) the highest possible optimal value for z.
	 * If we do not filter then the value is the largest of all deviations.
	 * If we filter then the value is the Gamma+1 largest value, potentially improved using a partitioning into cliques.
	 */
	default double computeHighestPossibleZ(Variable[] uncertainVariables, double Gamma, RobustAlgorithmStrategies robustAlgorithmStrategies, AlgorithmParameters algorithmParameters) {
		try {
			String output = "\nStart Computing Highest Possible Value for z";
			AbstractAlgorithm.writeOutput(output, algorithmParameters);
		} catch (IOException e) {
			e.printStackTrace();
		}

		int indexHighestPossibleZ = computeIndexHighestPossibleZ(uncertainVariables, Gamma, robustAlgorithmStrategies);
		if (indexHighestPossibleZ < 0) {
			return 0;
		}
		else {
			return uncertainVariables[indexHighestPossibleZ].getDeviation();
		}
	}
	
	/**
	 * Computes for an array of sorted variables (non-decreasing w.r.t deviation) a sorted (increasing) list of possible optimal values for z.
	 * If we do not filter then the set contains all deviations and 0.
	 */
	default List<PossibleZ> computePossibleZs(Variable[] uncertainVariables, double Gamma, RobustAlgorithmStrategies robustAlgorithmStrategies, ConflictGraph conflictGraph, AlgorithmParameters algorithmParameters) {
		long startTimePossibleZ = System.nanoTime();
		
		try {
			String output = "\nStart Computing Possible Optimal Values for z";
			AbstractAlgorithm.writeOutput(output, algorithmParameters);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//Adds zero as a possible optimal value.
		List<PossibleZ> possibleZs = new ArrayList<PossibleZ>();
		possibleZs.add(new PossibleZ(0));

		
		//Computes the largest index that has to be considered. 
		int indexHighestPossibleZ = computeIndexHighestPossibleZ(uncertainVariables, Gamma, robustAlgorithmStrategies);

		if (indexHighestPossibleZ > -1) {
			//If we do not filter or if Gamma is not integer then we add all deviations up to the highest index once.
			if (robustAlgorithmStrategies.getFilterStrategy() == RobustAlgorithmStrategies.FilterStrategy.FILTERINGZ_DISABLE || (int)Gamma != Gamma) {
				for (int i = 0; i <= indexHighestPossibleZ; i++) {
					if (possibleZs.get(possibleZs.size()-1).getValue() != uncertainVariables[i].getDeviation()) {
						possibleZs.add(new PossibleZ(uncertainVariables[i].getDeviation()));
					}
				}
			}
			//Otherwise we check for all indices 0,...,highestPossibleLowerZIndex-1 if the corresponding deviation must be added
			else {
				for (int i = 0; i < indexHighestPossibleZ; i++) {
					if(uncertainVariables[i].getDeviation() != possibleZs.get(possibleZs.size()-1).getValue()) {
						//If the current deviation is not already considered, we store all variables with a lower index whose deviation is higher than the last added z in a list.
						List<GRBVar> intermediateVariables = new ArrayList<GRBVar>();
						for (int j = i-1; j >= 0; j--) {
							if (possibleZs.get(possibleZs.size()-1).getValue() == uncertainVariables[j].getDeviation()) {
								break;
							}
							else {
								intermediateVariables.add(uncertainVariables[j].getModelVariable());
							}
						}
						//The deviation has to be added if there exists a variable in the list that is not in the neighborhood of the current variable in the conflict graph.
						if (!intermediateVariables.isEmpty()) {
							Set<GRBVar> neighborhood = new HashSet<GRBVar>();
							if (conflictGraph != null) {
								neighborhood = conflictGraph.getNeighbors(uncertainVariables[i].getModelVariable());
							}
							if (!neighborhood.containsAll(intermediateVariables)) {
								possibleZs.add(new PossibleZ(uncertainVariables[i].getDeviation()));
							}	
						}
					}
				}
				
				//Add the mandatory deviation of the variable with index highestPossibleLowerZIndex, if not already done.
				if(uncertainVariables[indexHighestPossibleZ].getDeviation() != possibleZs.get(possibleZs.size()-1).getValue()) {
					possibleZs.add(new PossibleZ(uncertainVariables[indexHighestPossibleZ].getDeviation()));
				}
			}
		}
		
		try {
			String output = "Finished Computing Possible z within "+((System.nanoTime()-startTimePossibleZ)/Math.pow(10, 9))+" sec\n"
					+"Number Uncertain Variables = "+uncertainVariables.length+"\n"
					+"Number Possible z = "+possibleZs.size();
			AbstractAlgorithm.writeOutput(output, algorithmParameters);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return possibleZs;
	}
	
	/**
	 * Specifies clique and filtering strategies;
	 */
	public class RobustAlgorithmStrategies extends AlgStrategies{
		/**
		 * Specifies whether the set of possibly optimal values for z should be filtered in advance.
		 */
		public enum FilterStrategy {
			FILTERINGZ_ENABLE,
			FILTERINGZ_DISABLE;
		}
		
		/**
		 * Specifies whether variables should be merged into cliques of a conflict graph.
		 */
		public enum CliqueStrategy {
			CLIQUES_ENABLE,
			CLIQUES_DISABLE;
		}
		
		CliqueStrategy cliqueStrategy;
		FilterStrategy filterStrategy;
		
		/**
		 * Constructor obtaining arguments which are matched to the enums defining strategies.
		 */
		public RobustAlgorithmStrategies(List<String> argList, AlgorithmParameters algorithmParameters) throws IOException {
			super(argList, algorithmParameters);
		}

		@Override
		void setDefaultStrategies() {
			cliqueStrategy = CliqueStrategy.CLIQUES_ENABLE;
			filterStrategy = FilterStrategy.FILTERINGZ_ENABLE;
		}
		
		public CliqueStrategy getCliqueStrategy() {
			return cliqueStrategy;
		}
		public void setCliqueStrategy(CliqueStrategy cliqueStrategy) {
			this.cliqueStrategy = cliqueStrategy;
		}
		
		public FilterStrategy getFilterStrategy() {
			return filterStrategy;
		}
		public void setFilterStrategy(FilterStrategy filterStrategy) {
			this.filterStrategy = filterStrategy;
		}
	}
}
