package dro.stat;

import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

/**
 * This is the abstract superclass for the bi- and multivariate GrangerCausality following C.J.W.
 * Grangers idea for causality. It basically describes, how useful one variable x is to predict
 * the future development of variable y (see wikipedia for more information). This is done by
 * creating 2 hypotheses H0 and H1 and applying an F-Test on both to analyse the statistical
 * evidence, provided in the data models, whether one can safely reject H0 or not.
 * <br>
 *
 * Please also have a look at Sergey Edunov's GNU implementation for the Granger-Causality.
 * https://code.google.com/p/jquant/source/browse/trunk/src/ru/algorithmist/jquant/math/GrangerTest.java?r=8
 *
 * @author MD
 * @version 1.0
 * @since 26/04/15
 */
public abstract class GrangerCausality {

	// The critical value used to evaluate, whether "X granger-causes Y" or not
	// if the pValue > criticalValue => "X granger-causes Y"
	public static double CRITICAL_VALUE = 0.95;
	// The lag size used to compute the Granger-causality
	protected int lagSize	= 1;
	private NumberFormat FORMATTER = new DecimalFormat("#.###");
	
	/**
	 * Super-constructor for Granger bivariate and multivariate causality computation.
	 * @param aLagSize The lag size that defines the amount of historical information incorporated
	 * into the Granger-Causality computation. Has to be > 0. <b>HINT: Check the method
     *                 #getMaximumLagSize to be on the safe side.</b>
	 * @param aCriticalValue The critical value is used to check, if the computed p-Value indicates
	 * that "X granger-causes Y". Has to be a value in the interval [0, 1.0], with 1.0 being the
	 * highest possible value. 
	 */
	public GrangerCausality(int aLagSize, double aCriticalValue) {
		this.lagSize = aLagSize;
		CRITICAL_VALUE = aCriticalValue;
	}

    /**
     * This method applies the Granger-Causality to the given time series parameters. The following
     * parameter configuration is expected:<br>
     * <ul>
     *      <li>1st vector: The time series for the y variable</li>
     *      <li>2nd vector: The time series for the x variable</li>
     *      <li>All following vectors: The time series for the remaining universe</li>
     * </ul>
     * Additionally, all given time series must have the same length, which has to be greater than
     * the provided lag window size at construction time.
     * @param tss The time series for the entire universe (containing at least y and x, in this
     *            order).
     * @return The computed GrangerCausalIndicator to answer the question: "Does x granger-cause y
     * in the given universe u?"
     */
    public abstract GrangerCausalIndicator apply(double[] ... tss);

	/**
	 * Suppose y is the target to be checked, x is the impact to be tested, u is the universe
	 * of all variables (including y, x and all other time series) and L is the given lag. The
	 * question is then:<br>
	 * <b>Does x Granger-cause y in the given universe u with the given lag L?</b><br>
	 * <br>
	 * The H0-test is then the following equation, which'll be solved via OLS:
	 * <pre>
	 *  y_t 	= (y_t-1 ... y_y-L u0_t-1 ... u0_t-L ... un_t-1 ... un_t-L) * (beta_1 ... beta_(n*L))^T + (epsilon_1 ... epsilon_(n*L))^T
	 *  ...
	 *  y_t-d	= (y_t-d-1 ... y_y-d-L u0_t-d-1 ... u0_t-d-L ... un_t-d-1 ... un_t-d-L) * (beta_1 ... beta_(n*L))^T + (epsilon_1 ... epsilon_(n*L))^T
	 * 	which translates to	
	 *  Y = U(without x) * BETA + EPSILON
	 * </pre>
	 * The H1-test is then the exact same as H0 but U contains also x.
	 * 
	 * This translates to the following variable names:
	 * <ul>
	 *   <li>Y => strippedY</li>
	 *   <li>X(without x) => laggedH0</li> 
	 *   <li>X => laggedH1</li>
	 * </ul>
	 * 
	 * @param laggedH0Rows The number of rows in U without x (the laggedH0 matrix). This is the
	 * d index in the above defined equation (which can also be called "number of observations").
	 * @param strippedY The Y vector that is cut off with the lag size.
	 * @param laggedH0 The universe matrix without the x variable whose influence is to be tested.
	 * The Matrix is also lagged with the given lag L.
     * @param variablesH0 The number of variables in the universe without the variable x.
     * @param laggedH1 The universe matrix with the x variable. The Matrix is also lagged with the
     * given lag L.
     * @param variablesH1 The number of variables in the universe with the variable x (naturally,
     *                    this should be variablesH0+1).
     * @return The computed GrangerCausalIndicator object, in case of success.
     * Null - If anything went wrong (e.g., the OLS can't compute the parameters due to a
     * SingularMatrixException)
	 */
	protected GrangerCausalIndicator performGranger(int laggedH0Rows, double[] strippedY,
			double[][] laggedH0, int variablesH0, double[][] laggedH1, int variablesH1) {
		try {
	        OLSMultipleLinearRegression h0 = new OLSMultipleLinearRegression();
	        OLSMultipleLinearRegression h1 = new OLSMultipleLinearRegression();
	        h0.newSampleData(strippedY, laggedH0);
	//        print(laggedH0);
	        h1.newSampleData(strippedY, laggedH1);
	//		print(laggedH1);
	        double rs0[] = h0.estimateResiduals();
	        double rs1[] = h1.estimateResiduals();
	//        System.out.print("b = {");
	//        for(int i=0; i<b0.length;i++){System.out.print(FORMATTER.format(b0[i])+", ");};
	//        System.out.println();
	//        System.out.print("residuals = {");
	//        for(int i=0; i<rs0.length;i++){System.out.print(FORMATTER.format(rs0[i])+", ");};
	//        System.out.println();
	        double RSS0 = sqrSum(rs0);
	        double RSS1 = sqrSum(rs1);
	        int nbrParametersModel0 = this.lagSize * variablesH0;
	        int nbrParametersModel1 = this.lagSize * variablesH1;
	        // 			(RSS1-RSS2) / (p2-p1)
	        // f	= -------------------------
	        //			    RSS2 / (n-p2)
//	        double ftest = ((RSS0 - RSS1)/this.lagSize) / (RSS1 / (laggedH0Rows - 2*this.lagSize - 1));
	        double ftest = ((RSS0 - RSS1) / (nbrParametersModel1-nbrParametersModel0))
	        		/ (RSS1 / (laggedH0Rows-nbrParametersModel1));
	
	//        System.out.println(RSS0 + " " + RSS1);
	//        System.out.println("F-test " + ftest);
	
	        FDistribution fDist = new FDistribution(this.lagSize, laggedH0Rows-2*this.lagSize-1);
	        double pValue = 1.0 - fDist.cumulativeProbability(ftest);
	        //        System.out.println("P-value " + pValue);
	        
	        return new GrangerCausalIndicator(pValue, CRITICAL_VALUE, this.lagSize);
		} catch(SingularMatrixException smex) {
			smex.printStackTrace();
			return null;
		}
	}
	
    /**
     * This method creates a lagged matrix from all given vectors by applying the lag window size
     * to each vector.
     *
     * <b>PRE-CONDITION:<br></b>
     * <ul>
     *     <li>Each vector has to be bigger than the lag window size, provided at creation time of
     *     this object.</li>
     * </ul>
     * @param a A number of input vectors.
     * @return The lagged matrix for all given vectors.
     */
    protected double[][] createLaggedSide(double[]... a) {
        int n = a[0].length - this.lagSize;
        double[][] res = new double[n][this.lagSize*a.length];
        for(int i=0; i<a.length; i++)
        {
            double[] ai = a[i];
            for(int l=0; l<this.lagSize; l++)
            {
                for(int j=0; j<n; j++)
                {
                    res[j][i*this.lagSize+l] = ai[l+j];
                }
            }
        }
        return res;
    }

    /**
     * @param a The original array for which a sum of squares shall be computed.
     * @return The sum of all squares for each element in the given array: result = a1*a1 + a2*a2 +
     *  ... + an*an 
     *  TODO: Use stream -> map and sum
     */
    protected double sqrSum(double[] a){
        double res = 0;
        for(double v : a) { res+=v*v; }
        return res;
    }
    
    /**
     * Checks, whether the GrangerTest can be applied for the given time series data size and
     * the given lag size.
     * @param dataSize Is the number of observations.
     * @param nbrOfVariables The number of unique variables in the given universe (e.g., 2 in case
     *                       of only y and x). ATTENTION: This is not the number of parameters in
     *                       the Granger-model.
     * @param lagSize The lag window size (the amount of historical information).
     * @return True - The Granger computation can be applied for the given parameters
     * False - The Granger-computation (and all underlying computations, e.g., OLS, F-Test) will
     * most likely fail due to specific constraints.
     */
    protected boolean checkDataSizeConstraints(int dataSize, int nbrOfVariables, int lagSize) {
    	// 1) Constraint 1: There must be more data available than predictors to be computed
    	//  Available data resolves to the overall data length - the lag size
    	int observedDataSize = dataSize-lagSize;
    	// Nbr of predictors is computed from the lagged XY matrix length (this is the number
    	// of parameters in the VAR model
    	int nbrOfPredictors = nbrOfVariables*lagSize;
    	if (observedDataSize <= nbrOfPredictors) {
    		// Keep in mind that the referenced 'data' in the following error message DOES NOT
    		// refer to the given time series data size, but rather to the lagged matrices size!!
    		throw new RuntimeException(String.format("There is "
    			+ "not enough data availabe (%d) to satisfy all predictors (%d) computed for the "
    			+ "given lag size (%d)",
    			observedDataSize, nbrOfPredictors, lagSize));
    	}
    	// 2) Constraint 2: The CDF for the GrangerTest requires the lag size to be smaller than:
//    	int maxLagSize = (int)Math.ceil((dataSize - 1) / 3.0);
    	int maxLagSize = getMaximumLagSize(observedDataSize, nbrOfVariables);
    	if (lagSize > maxLagSize) throw new RuntimeException(String.format("The lag size (%d) is "
    			+ "to big for the given data size (%d) and the GrangerTest contraints "
    			+ "(max lag size = %d).",
    			lagSize, dataSize, maxLagSize));
    	return true;
    }

    /**
     * Computes the maximum possible lag window size (maximum amount of historical information) that
     * can be used for the given size of your universe.
     * @param nbrOfObservations The number of observations for the given data (data size).
     * @param nbrOfVariables The number of unique variables in the current universe. 
     * @return The maximum lag that can be used to compute the Granger-Causality for the given
     * data size and number of unique variables. 
     */
    public static int getMaximumLagSize(int nbrOfObservations, int nbrOfVariables) {
    	// 1) Check if there's enough data for at least a lag of 1.
    	//    If we can't even support a lag of 1 throw an error.
    	if (nbrOfObservations-1 < nbrOfVariables)
    		throw new RuntimeException(String.format("ERROR: There's not even enough data (%s) "
    				+ "available for a lag size of 1 for the given number of variables (%s)",
    				nbrOfObservations, nbrOfVariables));
    	
    	// 2) Find the maximum lag, that satisfies the following condition:
    	int maxLag = (int)Math.floor(nbrOfObservations/(nbrOfVariables+1));
//    	int maxLag = 1;
//    	for (int curLag=1; curLag<Integer.MAX_VALUE; curLag++) {
//    		if (nbrOfObservations-curLag > nbrOfVariables*curLag) {
//    			maxLag = curLag;
//    		} else {
//    			break;
//    		}
//    	}
    	return maxLag;
    }
    
//    private Pair<Integer, Integer> computeGrangerTestContraints(int dataSize, int lagSize) {
//    	// 1) Constraint 1: There must be more data available than predictors to be computed
//    	//  Available data resolves to the overall data length - the lag size
//    	int availableData = dataSize - lagSize;
//    	//  Nbr of predictors is computed from the lagged XY matrix length, which resolves to:
//    	//  2*lagSize (cause X and Y variable with lagged columns)
//    	int nbrOfPredictors = 2*lagSize;
//    	// 2) Constraint 2: The CDF for thactuale GrangerTest requires the lag size to be smaller than:
//    	int maxLagSize = (int)Math.ceil((dataSize - 1) / 3.0);
//    	return new Pair<Integer, Integer>(value0, value1);
//    }
    
    /**
     * Cuts of the first 'l' leading elements in the given array a and returns the remaining array
     * with a length of a.length-l.
     * @param l The number of leading values to cut of the given array.
     * @param a The original array that should be shortened.
     * @return An array that contains the remaining elements of a, after the first 'l' elements
     * have been cut off.
     */
    protected double[] strip(double[] a){
        double[] res = new double[a.length-this.lagSize];
        System.arraycopy(a, this.lagSize, res, 0, res.length);
        return res;
    }
	
	public void print(double[][] matrix) {
		String str = "|\t";
        for (int i=0; i<matrix.length; i++)
        {
        	double[] row = matrix[i];
            for (int j=0; j<row.length; j++)
            {
                str += FORMATTER.format(matrix[i][j]);
                str += "\t";
            }
            System.out.println(str + "|");
            str = "|\t";
        }
	}
}
