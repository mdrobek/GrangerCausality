package dro.stat;

/**
 * This class represents the Bivariate Granger causality test for two given stationary time series,
 * as described by C. W. J. Granger.
 * 
 * @author MD
 * @version 1.0
 * @since 26/04/15
 */
public class GrangerCausalityStrategy_Bivariate extends GrangerCausality {

    /**
     * Creates a bivariate Granger computation strategy, with a default critical value (0.95).
     * @param aLagSize The lag window size (= the historical information that is used from each
     *                 given time series to compute the Granger-Causality.
     */
	public GrangerCausalityStrategy_Bivariate(int aLagSize) {
		this(aLagSize, CRITICAL_VALUE);
	}
	
	/**
     * Creates a new bivariate Granger computation strategy.
     * @param aLagSize The lag window size (= the historical information that is used from each
     *                 given time series to compute the Granger-Causality.
	 * @param aCriticalValue The critical values used to check, if the computed p-Value indicates
	 * that "X granger-causes Y". Has to be a value in the interval [0, 1.0], with 1.0 being the
	 * highest possible value.
	 */ 
	public GrangerCausalityStrategy_Bivariate(int aLagSize, double aCriticalValue) {
		super(aLagSize, aCriticalValue);
	}
	
    @Override
    public GrangerCausalIndicator apply(double[] ... tss) {
        return apply(tss[0], tss[1]);
    }

    /**
     * Applies the Granger-Causality computation to the 2 given timeseries t1 and t2
     * PRE-CONDITION:
     * <ul>
     *     <li>t1 and t2 must be different to each other, cause it is commonly believed, that in any
     *     case, the historical information of a stationary time series is beneficial in predicting
     *     its future development. y therefore always granger-causes itself.</li>
     * </ul>
     * @param y The time series of the target variable.
     * @param x The time series of the impacting variable.
     * @return The GrangerCausalIndicator object that contains all information to answer the
     * question: "Does x granger-cause y?"
     * Null, in case y and x are equal (see the above statement for explanation).
     */
	public GrangerCausalIndicator apply(double[] y, double[] x) {
		// Granger-Causality can't be computed for the same time series (which have the same values)
		if (y.equals(x)) return null;
		// Check whether the GrangerTest can be applied or not
		if (checkDataSizeConstraints(y.length, 2, super.lagSize)) {
			// Apply the actual test with y = tsValues[0], x = tsValues[1]
			double[][] laggedY = createLaggedSide(y);
			double[][] laggedXY = createLaggedSide(x, y);
			double[] strippedY = strip(y);
			return super.performGranger(laggedY.length, strippedY, laggedY, 1, laggedXY, 2);
		}
		return null;
	}

	@Override
	public String toString() {
		return "Bivariate Granger-Causality";
	}

}