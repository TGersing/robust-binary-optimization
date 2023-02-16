package alg;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import gurobi.GRB;
import gurobi.GRBException;
import gurobi.GRBModel;
import util.Variable;


/**
 * The AbstractAlgorithm class is an abstract parent class to all algorithms solving optimization problems.
 * We only support minimization problems. If the problem is a maximization problem then we convert it.
 * 
 * @author Timo Gersing
 */
public abstract class AbstractAlgorithm {
	
	/**
	 * Primal bound before a feasible solution is found. As we consider minimization problems, the primal bound is initialized to be positive infinity.
	 */
	public static final double DEFAULT_PRIMAL_BOUND = Double.POSITIVE_INFINITY;
	/**
	 * Dual bound before an actual bound is found. As we consider minimization problems, the dual bound is initialized to be negative infinity.
	 */
	public static final double DEFAULT_DUAL_BOUND = Double.NEGATIVE_INFINITY;
	/**
	 * The best known solution value. 
	 */
	protected double primalBound = DEFAULT_PRIMAL_BOUND;
	/**
	 * The best known dual bound on the solution value.
	 */
	protected double dualBound = DEFAULT_DUAL_BOUND;
	/**
	 * Stores the best known solution by mapping each variable to its value.  
	 */
	protected LinkedHashMap<Variable, Double> solution;
	
	/**
	 * System time in nanoseconds at the start of the algorithm.
	 */
	protected long startTime;
	/**
	 * System time in nanoseconds when the algorithm should be finished.
	 */
	protected Optional<Long> deadline = Optional.empty();
	/**
	 * System time in nanoseconds at the end of the algorithm.
	 */
	protected long endTime;
	
	/**
	 * Parameters for the algorithm.
	 */
	protected AlgorithmParameters algorithmParameters = new AlgorithmParameters();
	
	/**
	 * Constructor for the algorithm receiving parameters and potentially clearing the log file.
	 */
	public AbstractAlgorithm(AlgorithmParameters algorithmParameters) throws IOException, GRBException {
		this.algorithmParameters = algorithmParameters;
		if (algorithmParameters.getLogPath().isPresent()) {
			FileWriter fileWriter = new FileWriter(algorithmParameters.getLogPath().get(), false);
			fileWriter.close();
		}
	}
	
	/**
	 * Executes the specific algorithm.
	 */
	protected abstract void executeAlgorithm() throws IOException, GRBException;
	
	/**
	 * Initializes the start time and potential deadline, executes the algorithm,
	 * and writes the computed bounds as well as the gap and the elapsed time to the log. 
	 */
	public void solve(Optional<Integer> timeLimit) throws IOException, GRBException {
		this.startTime = System.nanoTime();
		if (timeLimit.isPresent()) {
			this.deadline = Optional.of(System.nanoTime() + (long)Math.pow(10, 9)*timeLimit.get());
		}
		executeAlgorithm();
		endTime = System.nanoTime();
		
		String output = "\n##########################################################\n"
				+ "##### Finished Solving"
				+ "\n##########################################################\n"
				+ "Primal Bound: "+primalBound+"\n"
				+ "Dual Bound: "+dualBound+"\n"
				+ "Relative Gap: "+(100*getRelativeGap())+"%\n"
				+ "Elapsed Time: "+getElapsedTime()+" sec";
		writeOutput(output);
	}

	/**
	 * Initializes the start time and deadline, executes the algorithm,
	 * and writes the computed bounds as well as the gap and the elapsed time to the log. 
	 */
	public void solve(int timeLimit) throws IOException, GRBException {
		this.solve(Optional.of(timeLimit));
	}

	/**
	 * Initializes the start time, executes the algorithm,
	 * and writes the computed bounds as well as the gap and the elapsed time to the log. 
	 */
	public void solve() throws IOException, GRBException {
		this.solve(Optional.empty());
	}
	
	/**
	 * Returns the solution values of all variables.
	 */
	public Map<Variable, Double> getSolution() {
		return solution;
	}
	
	/**
	* Returns the primal bound computed by the algorithm.
	*/
	public double getPrimalBound() {
		return primalBound;
	}
	/**
	* Returns the dual bound computed by the algorithm.
	*/
	public double getDualBound() {
		return dualBound;
	}

	/**
	 * Returns the relative gap between a given primal and dual bound.
	 */
	public static double getRelativeGap(double primalBound, double dualBound) {
		if (primalBound <= dualBound) {
			return 0.;
		}
		else if (primalBound != 0 && primalBound != DEFAULT_PRIMAL_BOUND) {
			return (primalBound-dualBound)/Math.abs(primalBound);
		}
		else {
			return Double.POSITIVE_INFINITY;
		}
	}
	
	/**
	 * Returns the relative gap between the primal and dual bound computed by the algorithm.
	 */
	public double getRelativeGap() {
		return getRelativeGap(this.getPrimalBound(), this.getDualBound());
	}

	/**
	 * Decides whether a problem is solved to optimality for a given primal and dual bound as well as optimality tolerances given in the parameters.
	 */
	public static boolean isOptimal(double primalBound, double dualBound, AlgorithmParameters algorithmParameters) {
		return (getRelativeGap(primalBound, dualBound) <= algorithmParameters.getRelativeGapTolerance() 
				|| primalBound - dualBound <= algorithmParameters.getAbsoluteGapTolerance());
	}

	/**
	 * Decides whether a pair of primal/dual bounds is within tolerance with respect to this.algorithmParameters.
	 */
	public boolean isOptimal(double primalBound, double dualBound) {
		return isOptimal(primalBound, dualBound, this.algorithmParameters);
	}
	
	/**
	 * Decides whether the primal/dual bounds computed by the algorithm are within tolerance with respect to this.algorithmParameters.
	 */
	public boolean isOptimal() {
		return isOptimal(this.getPrimalBound(), this.getDualBound());
	}
	
	/**
	 * Returns the time in seconds until the algorithm should be finished, if a deadline is present.
	 */
	protected Optional<Double> getRemainingTime() {
		if (deadline.isPresent()) {
			return Optional.of((deadline.get()-System.nanoTime())/Math.pow(10, 9));
		}
		else {
			return Optional.empty();
		}
	}
	
	/**
	 * Returns the time elapsed during the execution of the algorithm in seconds.
	 */
	public double getElapsedTime() {
		return (endTime-startTime)/Math.pow(10, 9);
	}
	
	/**
	 * Writes output to the console and log file as specified in this.algorithmParameters, if present.
	 */
	protected void writeOutput(String output) throws IOException {
		writeOutput(output, this.algorithmParameters);
	}
	
	/**
	 * Writes output to the console and log file as specified by the parameters, if present.
	 */
	public static void writeOutput(String output, AlgorithmParameters algorithmParameters) throws IOException {
		writeOutput(output, algorithmParameters.getLogPath());
	}
	
	/**
	 * Writes output to the console and log file as specified by logPath, if present.
	 */
	public static void writeOutput(String output, Optional<String> logPath) throws IOException {
		System.out.println(output);
		if (logPath.isPresent()) {
			FileWriter fileWriter = new FileWriter(logPath.get(), true);
			fileWriter.write(output+"\n");
			fileWriter.close();
		}
	}
	
	
	/**
	 * Parameters for the algorithm.
	 * Specifies logging, number threads, and tolerances.  
	 */
	public static class AlgorithmParameters {
		/**
		 * Path to the file to which the log is written.
		 */
		private Optional<String> logPath = Optional.empty();
		/**
		 * Decides whether Gurobi output is printed to console.
		 */
		private boolean gurobiConsoleOutput = true;
		/**
		 * Decides whether Gurobi output is printed to log file.
		 */
		private boolean gurobiLogOutput = true;
		/**
		 * Number of threads used by MILP solver
		 */
		private Optional<Integer> numberThreads = Optional.empty();
		/**
		 * Numerical tolerance for checking whether a value is integer. (default equals Gurobi default)
		 */
		private double integerFeasibilityTolerance = Math.pow(10, -5);
		/**
		 * Numerical tolerance for checking whether a constraint is satisfied. (default equals Gurobi default)
		 */
		private double feasibilityTolerance = Math.pow(10, -6);
		/**
		 * Tolerance indicating for which relative difference between primal and dual bound a problem is considered to be solved (default equals Gurobi default)
		 */
		private double relativeGapTolerance = Math.pow(10, -4);
		/**
		 * Tolerance indicating for which absolute difference between primal and dual bound a problem is considered to be solved (default equals Gurobi default)
		 */
		private double absoluteGapTolerance = Math.pow(10, -10);
		
		/**
		 * Constructor specifying the log path and using default for other parameters.
		 */
		public AlgorithmParameters(String logPath) {
			this.logPath = Optional.of(logPath);
		}
		/**
		 * Constructor initializing default parameters.
		 */
		public AlgorithmParameters() {}
		
		public Optional<String> getLogPath() {
			return logPath;
		}
		public void setLogPath(Optional<String> logPath) {
			this.logPath = logPath;
		}
		public boolean isGurobiConsoleOutput() {
			return gurobiConsoleOutput;
		}
		public void setGurobiConsoleOutput(boolean gurobiConsoleOutput) {
			this.gurobiConsoleOutput = gurobiConsoleOutput;
		}
		public boolean isGurobiLogOutput() {
			return gurobiLogOutput;
		}
		public void setGurobiLogOutput(boolean gurobiLogOutput) {
			this.gurobiLogOutput = gurobiLogOutput;
		}
		public Optional<Integer> getNumberThreads() {
			return numberThreads;
		}
		public void setNumberThreads(Optional<Integer> numberThreads) {
			this.numberThreads = numberThreads;
		}
		public double getIntegerFeasibilityTolerance() {
			return integerFeasibilityTolerance;
		}
		public void setIntegerFeasibilityTolerance(double integerFeasibilityTolerance) {
			this.integerFeasibilityTolerance = integerFeasibilityTolerance;
		}
		public double getFeasibilityTolerance() {
			return feasibilityTolerance;
		}
		public void setFeasibilityTolerance(double feasibilityTolerance) {
			this.feasibilityTolerance = feasibilityTolerance;
		}
		public double getRelativeGapTolerance() {
			return relativeGapTolerance;
		}
		public void setRelativeGapTolerance(double relativeGapTolerance) {
			this.relativeGapTolerance = relativeGapTolerance;
		}
		public double getAbsoluteGapTolerance() {
			return absoluteGapTolerance;
		}
		public void setAbsoluteGapTolerance(double absoluteGapTolerance) {
			this.absoluteGapTolerance = absoluteGapTolerance;
		}
		
		/**
		 * Applies parameters to the MILP model of a problem instance.
		 */
		public void applyParameters(ProblemGurobi problem) throws GRBException {
			GRBModel model = problem.getModel();
			if (gurobiLogOutput && logPath.isPresent()) {
				model.set(GRB.StringParam.LogFile, logPath.get());
			}
			if (gurobiConsoleOutput) {
				model.set(GRB.IntParam.LogToConsole, 1);
			}
			else {
				model.set(GRB.IntParam.LogToConsole, 0);
			}
			if (numberThreads.isPresent()) {
				model.set(GRB.IntParam.Threads, numberThreads.get());
			}
			model.set(GRB.DoubleParam.IntFeasTol, integerFeasibilityTolerance);
			model.set(GRB.DoubleParam.FeasibilityTol, feasibilityTolerance);
			model.set(GRB.DoubleParam.MIPGap, relativeGapTolerance);
			model.set(GRB.DoubleParam.MIPGapAbs, absoluteGapTolerance);
		}
	}
	
	static abstract class AlgStrategies {
		/**
		 * Constructor obtaining arguments which are matched to the enums defining strategies.
		 */
		public AlgStrategies(List<String> argList, AlgorithmParameters algorithmParameters) throws IOException {
			setDefaultStrategies();
			
			List<Field> fields = new ArrayList<Field>();
			Class<?> superClass = this.getClass();
			while (superClass != null) {
				fields.addAll(Arrays.asList(superClass.getDeclaredFields()));
				superClass = superClass.getSuperclass();
			}
			
			if (argList != null) {
				//Obtains optional arguments.
				for (String arg : argList) {
					boolean successful = false;
					for (Field field : fields) {
						try {
							field.set(this, Enum.valueOf(((Enum<?>) field.get(this)).getDeclaringClass(), arg.toUpperCase()));
							successful = true;
						} catch (ClassCastException | IllegalArgumentException | IllegalAccessException e) {}
					}
					if (!successful) {
						writeOutput("Argument "+arg+" could not be matched to an enum of "+this.getClass().getName()+" and will thus be ignored.", algorithmParameters);
					}
				}
				writeOutput("", algorithmParameters);
			}
			else {
				int counter = 0;
				while (true) {
					BufferedReader consoleReader = new BufferedReader(new InputStreamReader(System.in));
					String request = "\nDo you want to specify ";
					if (counter > 0) {
						request += "further ";
					}
					request += "strategies for the algorithm?\n"
							+ "Type in one of the following strategies to be defined (leave blank to not change strategies):\n";
					for (Field field : fields) {
						try {
							request += field.get(this).getClass().getSimpleName()+" ";
						} catch (IllegalArgumentException | IllegalAccessException e) {	}
					}
					System.out.println(request);
					String fieldClassName = consoleReader.readLine();
					if (fieldClassName.isEmpty()) {
						break;
					}
					Field chosenField = null;
					for (Field field : fields) {
						try {
							if (field.get(this).getClass().getSimpleName().toLowerCase().equals(fieldClassName.toLowerCase())) {
								chosenField = field;
								fieldClassName = field.get(this).getClass().getSimpleName();
							}
						} catch (IllegalArgumentException | IllegalAccessException e) {	}
					}
					if (chosenField != null) {
						try {
							request = "Type in one of the following options:\n";
							for (Enum<?> strategy : EnumSet.allOf(((Enum<?>) chosenField.get(this)).getDeclaringClass())) {
								request += strategy.name()+" ";
							}
							System.out.println(request);
							String input = consoleReader.readLine().toUpperCase();
							chosenField.set(this, Enum.valueOf(((Enum<?>) chosenField.get(this)).getDeclaringClass(), input));
							System.out.println("Set "+fieldClassName+" to "+input);
							counter++;
						} catch (IllegalArgumentException | IllegalAccessException e) {
							System.out.println("The input could not be matched to a strategy.");
						}
					}
					else {
						System.out.println("The input could not be matched to a strategy.");
					}
				}
			}
		}

		abstract void setDefaultStrategies();
	}
}