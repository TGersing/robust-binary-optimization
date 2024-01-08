package main;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

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
import alg.AlgSubmodularCuts;
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
		//bnb, dnc, rec, ref, sub, cut, bss, rp1, rp2, rp3, rp4
		String algorithmName = null;
		
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
					else if (key.toLowerCase().equals("numberthreads")) {
						algorithmParameters.setNumberThreads(Optional.of(Integer.parseInt(value)));
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
			JFileChooser fileChooser = new JFileChooser();
			File workingDirectory = new File(System.getProperty("user.dir"));
			fileChooser.setCurrentDirectory(workingDirectory);
			
			int result = fileChooser.showDialog(null, "Select the file stating the nominal problem");
			if (result == JFileChooser.APPROVE_OPTION) {
			    File selectedFile = fileChooser.getSelectedFile();
			    System.out.println("Selected file: " + selectedFile.getAbsolutePath());
				problemPath = selectedFile.getAbsolutePath();
			}
			else {
				return;
			}

			result = fileChooser.showDialog(null, "Select the file stating the robustness components");
			if (result == JFileChooser.APPROVE_OPTION) {
			    File selectedFile = fileChooser.getSelectedFile();
			    System.out.println("Selected file: " + selectedFile.getAbsolutePath());
			    robustPath = selectedFile.getAbsolutePath();
			}
			else {
				return;
			}
			
	        String response = (String) JOptionPane.showInputDialog(null,
		            "Choose One", "Select Algorithm",
		            JOptionPane.QUESTION_MESSAGE, null,
		            new String[] {"Branch and Bound", "Divide and Conquer", "Recycling", "Standard Reformulation", "Submodular Cuts", "Scenario Separation", "Bertsimas Sim Sequence", "Atamturk RP1", "Atamturk RP2", "Atamturk RP3", "Atamturk RP4"},
		            "Branch and Bound");
			switch (response) {
				case "Branch and Bound":
					algorithmName = "bnb";
					break;
				case "Divide and Conquer":
					algorithmName = "dnc";
					break;
				case "Recycling":
					algorithmName = "rec";
					break;
				case "Standard Reformulation":
					algorithmName = "ref";
					break;
				case "Submodular Cuts":
					algorithmName = "sub";
					break;
				case "Scenario Separation":
					algorithmName = "cut";
					break;
				case "Bertsimas Sim Sequence":
					algorithmName = "bss";
					break;
				case "Atamturk RP1":
					algorithmName = "rp1";
					break;
				case "Atamturk RP2":
					algorithmName = "rp2";
					break;
				case "Atamturk RP3":
					algorithmName = "rp3";
					break;
				case "Atamturk RP4":
					algorithmName = "rp4";
					break;
			}
			
			JPanel parameterPanel = new JPanel();
			parameterPanel.setLayout(new BoxLayout(parameterPanel, BoxLayout.Y_AXIS));
			
			JTextField timeLimitField = new JTextField();
			parameterPanel.add(new JLabel("State the time limit in seconds (leave blank for no limit)"));
			parameterPanel.add(timeLimitField);
			parameterPanel.add(Box.createVerticalStrut(20));
			
			JTextField threadsField = new JTextField();
			parameterPanel.add(new JLabel("State the thread limit (leave blank for no limit)"));
			parameterPanel.add(threadsField);
			parameterPanel.add(Box.createVerticalStrut(20));

			result = JOptionPane.showConfirmDialog(null, parameterPanel, "Enter Algorithm Parameters", JOptionPane.DEFAULT_OPTION);
			if (result == JOptionPane.OK_OPTION) {
				if (!timeLimitField.getText().equals("")) {
					timeLimit = Optional.of(Integer.parseInt(timeLimitField.getText()));
				}
				if (!threadsField.getText().equals("")) {
					algorithmParameters.setNumberThreads(Optional.of(Integer.parseInt(threadsField.getText())));
				}
			}
			
			result = fileChooser.showDialog(null, "Select the destination for the log file (cancel for no log file)");
			if (result == JFileChooser.APPROVE_OPTION) {
			    File selectedFile = fileChooser.getSelectedFile();
			    System.out.println("Selected file: " + selectedFile.getAbsolutePath());
			    algorithmParameters.setLogPath(Optional.of(selectedFile.getAbsolutePath()));
			}
			
			result = fileChooser.showDialog(null, "Select the destination for the results file (cancel for no results file)");
			if (result == JFileChooser.APPROVE_OPTION) {
			    File selectedFile = fileChooser.getSelectedFile();
			    System.out.println("Selected file: " + selectedFile.getAbsolutePath());
			    resultsPath = Optional.of(selectedFile.getAbsolutePath());
			}
			
			result = fileChooser.showDialog(null, "Select the destination for the solution file (cancel for no solution file)");
			if (result == JFileChooser.APPROVE_OPTION) {
			    File selectedFile = fileChooser.getSelectedFile();
			    System.out.println("Selected file: " + selectedFile.getAbsolutePath());
			    solutionPath = Optional.of(selectedFile.getAbsolutePath());
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
					+ ";Elapsed Time";
			if (algorithmName.toLowerCase().equals("bnb")) {
				resultsOutput += ";Number Uncertain Variables"
						+ ";Number Cliques"
						+ ";Number Possible Z"
						+ ";Number Integer Subproblems Started"
						+ ";Time in Integer Subproblems"
						+ ";Primal Bound after First Integer Subproblem"
						+ ";Gap Primal Bound and First Solution"
						+ ";Number LP Started"
						+ ";Time in LP";
			}
			resultsOutput+="\n";
			
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

		
		try {
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
				case "sub":
					algo = new AlgSubmodularCuts(problemPath, robustPath, algorithmParameters, new RobustAlgorithmStrategies(argList, algorithmParameters));
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
			if (solutionPath.isPresent() && algo.getSolution() != null) {
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
						+ ";"+algo.getElapsedTime());
				if (algo instanceof AlgBranchAndBound) {
					AlgBranchAndBound bnbAlgorithm = (AlgBranchAndBound) algo;
					resultsWriter.write(";"+bnbAlgorithm.getNumberUncertainVariables()
							+ ";"+bnbAlgorithm.getNumberCliques()
							+ ";"+bnbAlgorithm.getNumberPossibleZ()
							+ ";"+bnbAlgorithm.getNumberIntegerStarted()
							+ ";"+bnbAlgorithm.getTimeSpendInInteger()
							+ ";"+bnbAlgorithm.getPrimalBoundAfterFirstInteger()
							+ ";"+AbstractAlgorithm.getRelativeGap(bnbAlgorithm.getPrimalBoundAfterFirstInteger(), bnbAlgorithm.getPrimalBound())
							+ ";"+bnbAlgorithm.getNumberLPStarted()
							+ ";"+bnbAlgorithm.getTimeSpendInLP());
				}
				resultsWriter.write("\n");
				resultsWriter.close();
			}
		} catch (Throwable e) {
			AbstractAlgorithm.writeOutput(e.toString(), algorithmParameters);
			
			e.printStackTrace();
			
			if (resultsPath.isPresent()) {
				//Clears results file
				FileWriter resultsWriter = new FileWriter(resultsPath.get(), false);
				resultsWriter.write(resultsOutput);
				
				resultsWriter.write(";"+e.toString());
				resultsWriter.write("\n");
				resultsWriter.close();
			}
		}
	}
}