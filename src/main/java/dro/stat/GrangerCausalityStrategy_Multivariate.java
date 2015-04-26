package dro.stat;

/**
 * The created WeightMatrix and their elements Cell(rowIndex, columnIndex) have to be read
 * as follows:
 * <ul>
 * 	<li>y = Variable associated with the row index</li>  
 * 	<li>x = Variable associated with the column index</li>  
 * </ul>
 * Cell(rowIndex, columnIndex) means therefore:
 * <b>How strong columnIndex granger-causes rowIndex.</b><br>
 * Or in other words:
 * Each row shows how strong the associated variable is granger-caused by the respective column.
 * 
 * @author MD
 * @version 1.0
 * @since 26/04/15
 */
public class GrangerCausalityStrategy_Multivariate extends GrangerCausality {

    /**
     * Creates a new object for the multivariate GrangerCausality strategy with a default critical
     * value (0.95).
     * @param aLagSize The lag window size (amount of historical information to be incorporated in
     *                 the Granger computation).
     */
	public GrangerCausalityStrategy_Multivariate(int aLagSize) {
		this(aLagSize, CRITICAL_VALUE);
	}

	/**
     * Creates a new object for the multivariate GrangerCausality strategy.
     * @param aLagSize The lag window size (amount of historical information to be incorporated in
     *                 the Granger computation).
	 * @param aCriticalValue The critical values used to check, if the computed p-Value indicates
	 * that "X granger-causes Y". Has to be a value in the interval [0, 1.0], with 1.0 being the
	 * highest possible value.
	 */ 
	public GrangerCausalityStrategy_Multivariate(int aLagSize, double aCriticalValue) {
		super(aLagSize, aCriticalValue);
	}


    @Override
    public GrangerCausalIndicator apply(double[]... tss) {
        if (checkDataSizeConstraints(tss[0].length, tss.length, super.lagSize)) {
            apply(0, 1, tss, tss.length);
        }
        return null;
    }

	/**
	 * Applies the Granger-Causality computation to the given universe.
	 * @param yVectorRow The index of the row for the y variable in the given universe matrix.
	 * @param xVectorRow The index of the row for the x variable in the given universe matrix.
	 * @param universe The universe matrix. Each row reflects the time series of one variable.
	 * @param nbrVariablesInUniverse The number of unique time series in the universe (not the
	 * number of columns in the universe matrix).
	 * @return The GrangerCausalIndicator that contains all information to answer the question:
     * "Does x granger-cause y in the given universe U?"
     * Null, if some constraints for the underlying computations (OLS, F-Test) are not met.
	 */
	private GrangerCausalIndicator apply(int yVectorRow, int xVectorRow,
			double[][] universe, int nbrVariablesInUniverse) {
		// Order the universe: y variable will be in row 0 and x variable will be in row 1
		double[][] orderedUniverse = orderUniverse(universe, yVectorRow, xVectorRow);
		// The x variable is after the ordering process always in row 1 
		double[][] orderedUniverseWithoutX = removeRow(orderedUniverse, 1);
		double[][] laggedH0 = super.createLaggedSide(orderedUniverseWithoutX);
		double[][] laggedH1 = super.createLaggedSide(orderedUniverse);
		double[] strippedY = super.strip(universe[yVectorRow]);
		return super.performGranger(laggedH0.length, strippedY, laggedH0,
				nbrVariablesInUniverse-1, laggedH1, nbrVariablesInUniverse);
	}
	
	protected double[][] removeRow(double[][] orig, int rowIndex) {
		int observations = orig[0].length;
		double[][] smallerMat = new double[orig.length-1][observations];
		int smallerMatRowCounter = 0;
		for (int row=0; row<orig.length; row++) {
			if (row != rowIndex) {
				System.arraycopy(orig, row, smallerMat, smallerMatRowCounter++, 1);
			}
		}
		return smallerMat;
	}
	
	protected double[][] orderUniverse(double[][] universe, int yRow, int xRow) {
		double[][] ordered = new double[universe.length][universe[0].length];
		// Copy y vector to first row
		System.arraycopy(universe, yRow, ordered, 0, 1);
		// Copy x vector to second row
		System.arraycopy(universe, xRow, ordered, 1, 1);
		// Copy all remaining vectors in their appearing order
		int targetRow = 2;
		for (int i=0; i<universe.length; i++) {
			if (i != yRow && i != xRow)
				System.arraycopy(universe, i, ordered, targetRow++, 1);
		}
		return ordered;
	}

	@Override
	public String toString() { return "Multivariate Granger-Causality";	}

}
