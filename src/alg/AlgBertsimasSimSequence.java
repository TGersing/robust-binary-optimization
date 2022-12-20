package alg;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;

import alg.RobustAlgorithm.RobustAlgorithmStrategies.CliqueStrategy;
import alg.SubproblemNominalGurobi.NOSStrategies;
import gurobi.GRBException;
import util.ConflictGraph;
import util.PossibleZ;
import util.Variable;

/**
 * This class implements the algorithm of Bertsimas and Sim solving a sequence of nominal subproblems.
 * 
 * @author Timo Gersing
 */
public class AlgBertsimasSimSequence extends AbstractAlgorithm implements RobustAlgorithm{
	/**
	 * The nominal subproblem that we adapt for different choices of z.
	 */
	private SubproblemNominalGurobi subproblemNominal;
	
	/**
	 * Specifies strategies.
	 */
	private BSStrategies bsStrategies;
	
	/**
	 * Stores values of the current best solution.
	 */
	double[] solutionValues;

		
	/**
	 * Constructor storing the given robust problem and strategies.
	 */
	public AlgBertsimasSimSequence(String problemPath, String robustPath, AlgorithmParameters algorithmParameters, BSStrategies bsStrategies) throws IOException, GRBException {
		super(algorithmParameters);
		this.subproblemNominal = new SubproblemNominalGurobi(problemPath, robustPath, algorithmParameters, bsStrategies);
		this.bsStrategies = bsStrategies;
	}
	
	
	/**
	 * Executes the algorithm.
	 */
	@Override
	protected void executeAlgorithm() throws IOException, GRBException {
		String output = "##########################################################\n"
				+ "##### Solving Problem via Bertsimas Sim Sequence"
				+ "\n##########################################################";
		writeOutput(output);

		output = "\n###############################\n"
				+ "##### Start Preprocessing"
				+ "\n###############################";
		writeOutput(output);
		//Computes a conflict graph if the corresponding strategy is chosen.
		ConflictGraph conflictGraph = null;
		if (bsStrategies.getCliqueStrategy() == CliqueStrategy.CLIQUES) {
			conflictGraph = new ConflictGraph(subproblemNominal.getModel(), subproblemNominal.getUncertainModelVariables(), algorithmParameters);
		}
		//Computes a list of possible optimal choices for z respecting the chosen filtering and clique strategies.
		List<PossibleZ> possibleZs = computePossibleZs(subproblemNominal.getUncertainVariables(), subproblemNominal.getGamma(), bsStrategies, conflictGraph, algorithmParameters);
		//The final dual bound is the minimum over the computed dual bound for all nominal subproblems.
		double minimumDualBoundSubproblems = DEFAULT_PRIMAL_BOUND;
		
		//The index of the current possibleZ in the list possibleZs.
		//This is needed to determine whether the algorithm stopped after considering the last nominal subproblem or before, due to the time limit.
		int indexCurrentZ = 0;
		
		output = "\n###############################\n"
				+ "##### Start Sequence"
				+ "\n###############################";
		writeOutput(output);
		
		//Iterate over all possible values of z.
		while (indexCurrentZ < possibleZs.size()) {
			//Terminates the algorithm if the time limit is reached.
			if (deadline.isPresent() && System.nanoTime() > deadline.get()) {
				output = "\n###########################\n"
						+ "##### Reached Time Limit"
						+ "\n###########################";
				writeOutput(output);
				break;
			}
			//Output to the log file.
			output = "\n##### BS Sequence Information\n"
					+ "\nElapsed Time = "+(System.nanoTime()-startTime)/Math.pow(10, 9)
					+ "\nNumber Remaining Possible Z = "+(possibleZs.size()-indexCurrentZ)+" of "+possibleZs.size()
					+ "\nCurrent Primal Bound = "+primalBound;
			writeOutput(output);


			//Initializes the next z to consider.
			PossibleZ currentZ = possibleZs.get(indexCurrentZ);
			
			//Alters the nominal subproblem with respect to z and then tries to solve it within the remaining time.
			output = "\n##### Solving NOS(z) for z="+currentZ;
			writeOutput(output);
			
			subproblemNominal.updateZ(currentZ);
			subproblemNominal.setGlobalPrimalBound(primalBound);
			subproblemNominal.solve(getRemainingTime());
			
			//Updates the primal bound.
			double incumbentValue = subproblemNominal.getGlobalPrimalBound();
			if (incumbentValue < primalBound) {
				primalBound = incumbentValue;
				solutionValues = subproblemNominal.getIncumbentValues();
			}
			
			//Updates the minimum dual bound over all considered nominal subproblems.
			minimumDualBoundSubproblems = Math.min(minimumDualBoundSubproblems, subproblemNominal.getDualBound());
			indexCurrentZ++;
		}
		
		//If the algorithm terminated after considering all nominal subproblems then the dual bound is the minimum over the computed dual bound for all nominal subproblems.
		//Otherwise the dual bound the dual bound is infinity.
		if (indexCurrentZ == possibleZs.size()) {
			dualBound = minimumDualBoundSubproblems;
		}
		else {
			dualBound = DEFAULT_DUAL_BOUND;
		}
		
		//Stores best solution found
		solution = new LinkedHashMap<Variable, Double>(solutionValues.length);
		for (int i = 0; i < solutionValues.length; i++) {
			solution.put(subproblemNominal.getNominalVariables()[i], solutionValues[i]);
		}
	}
	
	/**
	 * Specifies all strategies for the BS algorithm;
	 */
	public static class BSStrategies extends NOSStrategies{
		
		public BSStrategies() {
			improvingZStrategy = NOSImprovingZStrategy.IMPROVE_Z;
			terminationStrategy = NOSTerminationStrategy.TERMINATE_DIRECT;
		}
		
		public NOSImprovingZStrategy getImprovingZStrategy() {
			return improvingZStrategy;
		}
		public void setImprovingZStrategy(NOSImprovingZStrategy improvingZStrategy) {
			this.improvingZStrategy = improvingZStrategy;
		}

		public NOSTerminationStrategy getTerminationStrategy() {
			return terminationStrategy;
		}
		public void setTerminationStrategy(NOSTerminationStrategy terminationStrategy) {
			this.terminationStrategy = terminationStrategy;
		}
	}

}