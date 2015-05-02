import dro.stat.GrangerCausalityStrategy_Bivariate;
import org.junit.Assert;
import org.junit.Test;

/**
 * Class that provides all unit tests for the default Granger-Causality methods and the specific
 * bivariate Granger-Methods.
 *
 * @author MD
 * @version 1.0
 * @since 02/05/15
 */
public class BivariateGrangerTest {

    /**
     * Tests the square sum method.
     */
    @Test
    public void testSqrSumMethod() {
        // T1) Square sum of an empty array is .0
        BivariateGrangerImpl big = new BivariateGrangerImpl(22);
        double sqrSum = big.sqrSum(new double[]{});
        Assert.assertEquals(.0, sqrSum, .0);

        // T2) Square sum of an example has to be correct
        double[] testArr = new double[] { 1, 2, 3, 4, 5 };
        sqrSum = big.sqrSum(testArr);
        Assert.assertEquals(55, sqrSum, .0);
    }


    /**
     * Tests the strip method that cuts off l given elements from an array.
     */
    @Test
    public void testStripMethod() {
        // T1) Can't strip an array of l elements, if the array length is smaller than l
        //     => expect an exception that states this
        final int lagSize = 5;
        BivariateGrangerImpl big = new BivariateGrangerImpl(lagSize);
        double[] testArr = new double[] { .0, -2.1, Math.PI, 22.2 };
        try {
            big.strip(testArr);
            Assert.fail();
        } catch (Exception e) { /* Nothing to do here */ }

        // T2) Stripping of an array should always create a new array
        testArr = new double[] { 0, 0, 0, 0, 0 };
        double[] testArr1 = new double[] { 0, 0, 0, 0, 0 };
        double[] stripped = big.strip(testArr);
        Assert.assertFalse(stripped.equals(testArr));
        BivariateGrangerImpl big0 = new BivariateGrangerImpl(0);
        stripped = big0.strip(testArr);
        Assert.assertFalse(stripped.equals(testArr));

        // T3) Striping l elements from array with length n, leaves a new array of length n-l
        testArr = new double[] { .0, -2.1, Math.PI, 22.2, -343434 };
        BivariateGrangerImpl big3 = new BivariateGrangerImpl(3);
        stripped = big3.strip(testArr);
        Assert.assertArrayEquals(stripped, new double[]{ 22.2, -343434 }, .0);
    }


    /**
     * Private test class that extends the bivariate Granger-Causality class to gain access to all
     * protected methods.
     */
    private class BivariateGrangerImpl extends GrangerCausalityStrategy_Bivariate {
        public BivariateGrangerImpl(int aLagSize) { super(aLagSize); }
        @Override
        protected double[] strip(double[] a) { return super.strip(a); }
        @Override
        protected double sqrSum(double[] a) { return super.sqrSum(a); }
    }
}
