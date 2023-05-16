package main;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import alg.AlgBranchAndBound;
import alg.AbstractAlgorithm;
import alg.AbstractAlgorithm.AlgorithmParameters;
import alg.AlgBertsimasSimSequence;
import alg.AlgBertsimasSimSequence.BSStrategies;
import alg.AlgBranchAndBound.BnBStrategies;
import alg.AlgCuttingPlanesGurobi;
import alg.AlgDivideAndConquer;
import alg.AlgDivideAndConquer.DnCStrategies;
import alg.AlgRP1Gurobi;
import alg.AlgRP2Gurobi;
import alg.AlgRP3Gurobi;
import alg.AlgRP4Gurobi;
import alg.AlgRecycleInequalitiesGurobi;
import alg.AlgRecycleInequalitiesGurobi.RecyclingStrategies;
import alg.AlgReformulationGurobi;
import alg.RobustAlgorithm.RobustAlgorithmStrategies;
import gurobi.GRB;
import gurobi.GRBException;
import util.Variable;

/**
 * Receives paths to problem and robustness components, the name of the algorithm,
 * and optionally a time limit, log path, and a path to the solution file.
 * 
 * @author Timo Gersing
 */
public class AlgorithmExecuter {

	public static void main(String[] args) throws GRBException, IOException {
		//Path to the problem and the robustness components.
		String problemPath;
		String robustPath;
		//Algorithm to be used. Can be one of the following:
		//bnb, dnc, rec, ref, cut, bss, rp1, rp2, rp3, rp4
		String algorithmName;
		
		//Optionally also receives "timelimit=*", "logpath=*.txt", "solutionpath=*.txt"
		Optional<Integer> timeLimit = Optional.empty();
		AlgorithmParameters algorithmParameters = new AlgorithmParameters();
		Optional<String> solutionPath = Optional.empty();
		Optional<String> resultsPath = Optional.empty();
		
		List<String> argList = new ArrayList<String>(Arrays.asList(args));
		//If arguments are given
		if (args.length > 0) {
			//Obtains the path to the problem and the robustness components.
			problemPath = argList.get(0);
			argList.remove(0);
			robustPath = argList.get(0);
			argList.remove(0);
			//Obtains the algorithm to be used.
			algorithmName = argList.get(0);
			argList.remove(0);
			
			//Obtains optional arguments.
			Iterator<String> argIterator = argList.iterator();
			while (argIterator.hasNext()) {
				String arg = argIterator.next();
				if (arg.contains("=")) {
					String key = arg.substring(0, arg.indexOf('='));
					String value = arg.substring(arg.indexOf('=')+1);
					if (key.toLowerCase().equals("timelimit")) {
						timeLimit = Optional.of(Integer.parseInt(value));
						argIterator.remove();
					}
					else if (key.toLowerCase().equals("logpath")) {
						algorithmParameters.setLogPath(Optional.of(value));
						argIterator.remove();
					}
					else if (key.toLowerCase().equals("resultspath")) {
						resultsPath = Optional.of(value);
						argIterator.remove();
					}
					else if (key.toLowerCase().equals("solutionpath")) {
						solutionPath = Optional.of(value);
						argIterator.remove();
					}
				}
			}
		}
		//Reads from console if no arguments are given
		else {
			BufferedReader consoleReader = new BufferedReader(new InputStreamReader(System.in));
			System.out.println("State the path to the nominal problem:");
			problemPath = consoleReader.readLine();
			System.out.println("State the path to the robustness components:");
			robustPath = consoleReader.readLine();
			System.out.println("State the algorithm to be used (bnb, dnc, rec, ref, cut, bss, rp1, rp2, rp3, rp4):");
			algorithmName = consoleReader.readLine();
			
			System.out.println("State the time limit (leave blank for no time limit):");
			String input = consoleReader.readLine();
			if (!input.equals("")) {
				timeLimit = Optional.of(Integer.parseInt(input));
			}
			
			System.out.println("State the destination for the log file (leave blank for no log file):");
			input = consoleReader.readLine();
			if (!input.equals("")) {
				algorithmParameters.setLogPath(Optional.of(input));
			}
			
			System.out.println("State the destination for the results file (leave blank for no results file):");
			input = consoleReader.readLine();
			if (!input.equals("")) {
				resultsPath = Optional.of(input);
			}
			
			System.out.println("State the destination for the solution file (leave blank for no solution file):");
			input = consoleReader.readLine();
			if (!input.equals("")) {
				solutionPath = Optional.of(input);
			}
		}
		
		String resultsOutput = "";
		if (resultsPath.isPresent()) {
			//Clears results file
			FileWriter resultsWriter = new FileWriter(resultsPath.get(), true);
			//Writes results file
			resultsOutput = "Instance"
					+ ";Robustness Components"
					+ ";Algorithm"
					+ ";Strategies"
					+ ";Primal Bound"
					+ ";Dual Bound"
					+ ";Primal Dual Integral"
					+ ";Relative Gap"
					+ ";Optimal"
					+ ";Elapsed Time\n";
			
			resultsOutput += problemPath.split("/")[problemPath.split("/").length-1];
			resultsOutput += ";"+robustPath.split("/")[robustPath.split("/").length-1];
			resultsOutput += ";"+algorithmName;
			resultsOutput += ";";
			for (int i = 0; i < argList.size(); i++) {
				if (i==0) {
					resultsOutput += argList.get(i);
				}
				else {
					resultsOutput += "_"+argList.get(i);
				}
			}
			resultsWriter.write(resultsOutput);
			resultsWriter.write(";No Results Available\n");
			resultsWriter.close();
		}

		
		//Initializes the chosen algorithm
		AbstractAlgorithm algo = null;
		if (args.length == 0) {
			argList = null;
		}
		switch (algorithmName.toLowerCase()) {
			case "bnb":
				algo = new AlgBranchAndBound(problemPath, robustPath, algorithmParameters, new BnBStrategies(argList, algorithmParameters));
				break;
			case "dnc":
				algo = new AlgDivideAndConquer(problemPath, robustPath, algorithmParameters, new DnCStrategies(argList, algorithmParameters));
				break;
			case "rec":
				algo = new AlgRecycleInequalitiesGurobi(problemPath, robustPath, algorithmParameters, new RecyclingStrategies(argList, algorithmParameters));
				break;
			case "ref":
				algo = new AlgReformulationGurobi(problemPath, robustPath, algorithmParameters, new RobustAlgorithmStrategies(argList, algorithmParameters));
				break;
			case "cut":
				algo = new AlgCuttingPlanesGurobi(problemPath, robustPath, algorithmParameters);
				break;
			case "bss":
				algo = new AlgBertsimasSimSequence(problemPath, robustPath, algorithmParameters, new BSStrategies(argList, algorithmParameters));
				break;
			case "rp1":
				algo = new AlgRP1Gurobi(problemPath, robustPath, algorithmParameters, new RobustAlgorithmStrategies(argList, algorithmParameters));
				break;
			case "rp2":
				algo = new AlgRP2Gurobi(problemPath, robustPath, algorithmParameters);
				break;
			case "rp3":
				algo = new AlgRP3Gurobi(problemPath, robustPath, algorithmParameters);
				break;
			case "rp4":
				algo = new AlgRP4Gurobi(problemPath, robustPath, algorithmParameters, new RobustAlgorithmStrategies(argList, algorithmParameters));
				break;
			default:
				throw new IllegalArgumentException("Unexpected value: " + algorithmName.toLowerCase());
		}
		
		//Solves the problem using the chosen algorithm
		algo.solve(timeLimit);
		
		//Writes solution file
		if (solutionPath.isPresent()) {
			//Obtains solution
			Map<Variable, Double> solution = null;
			solution = algo.getSolution();

			//Clears solution file
			new FileWriter(solutionPath.get(), false).close();
			
			//Writes solution
			String output = "";
			for (Variable var : solution.keySet()) {
				output += var.getModelVariable().get(GRB.StringAttr.VarName) + "=" + solution.get(var)+"\n";
			}
			AbstractAlgorithm.writeOutput(output, solutionPath);
		}
		if (resultsPath.isPresent()) {
			//Clears results file
			FileWriter resultsWriter = new FileWriter(resultsPath.get(), false);
			resultsWriter.write(resultsOutput);
			
			resultsWriter.write(";"+algo.getPrimalBound()
					+ ";"+algo.getDualBound()
					+ ";"+algo.getPrimalDualIntegral()
					+ ";"+algo.getRelativeGap()
					+ ";"+algo.isOptimal()
					+ ";"+algo.getElapsedTime()+"\n");
			resultsWriter.close();
		}
	}
}