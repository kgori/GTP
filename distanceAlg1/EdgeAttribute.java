package distanceAlg1;

import java.text.DecimalFormat;
import java.util.Arrays;

public class EdgeAttribute {
    double[] vect;

    public static final double TOLERANCE = 0.000000000000001;  // compares double values up to TOLERANCE (currently 15 decimal places)

    public EdgeAttribute() {
        this.vect = null;
    }

    public EdgeAttribute(double[] vect) {
        this.vect = vect;
    }

    public EdgeAttribute(String s) {
        try {
            // if s is a vector
            if ((s.charAt(0) == '[') && (s.charAt(s.length() - 1) == ']')) {
                // remove [ and ]
                s = s.substring(1, s.length() - 1);
                // split the string at ','
                String[] elements = s.split(" ");
                this.vect = new double[elements.length];
                // convert vector element strings into doubles
                for (int i = 0; i < elements.length; i++) {
                    vect[i] = Double.parseDouble(elements[i]);
                }
            }
            // else s is just a double (vector of length 1)
            // (or s is neither, but then will be caught by the NumberFormatException)
            else {
                this.vect = new double[1];
                vect[0] = Double.parseDouble(s);
            }
        } catch (NumberFormatException e) {
            System.err.println("Error creating new edge attribute: input string does not have double where expected or bracket problem: " + e.getMessage());
            System.exit(1);
        }
    }


    public double getAttribute() {
        if (vect == null) {
            return 0;
        }
        return vect[0];
    }

    public void setEdgeAttribute(EdgeAttribute attrib) {
        this.vect = attrib.vect;
    }

    public EdgeAttribute clone() {
        return new EdgeAttribute(Arrays.copyOf(vect, vect.length));
    }

    public String toString() {
        if (vect == null) {
            return "";
        }

        // 10 decimals
        DecimalFormat df = new DecimalFormat("#0.##########");
        if (vect.length == 1) {
            return df.format(vect[0]);
        }

        String str = "[" + df.format(vect[0]);
        for (int i = 1; i < vect.length; i++) {
            str = str + " " + df.format(vect[i]);
        }
        return str + "]";
    }

    // TODO:  only set up to handle the attribute being a vector
    @Override
    public boolean equals(Object e) {
        if (e == null) {
            return false;
        }
        if (this == e) {
            return true;
        }

        if (!(e instanceof EdgeAttribute)) {
            return false;
        }

        // we cannot just use Arrays.equal, since we need to compare the double values with a tolerance.
        for (int i = 0; i < vect.length; i++) {
            if (Math.abs(vect[i] - ((EdgeAttribute) e).vect[i]) > TOLERANCE) {
                return false;
            }
        }

        return true;
    }

    /**
     * Compute the L2 norm or equivalent for this kind of attribute.
     *
     * @param attrib
     * @return
     */
    public double norm() {
        if (vect == null) {
            return 0.0;
        }

        double norm = 0;

        for (int i = 0; i < vect.length; i++) {
            norm = norm + Math.pow(vect[i], 2);
        }

        return Math.sqrt(norm);
    }

    /**
     * Compute the difference between two attributes,
     * by subtracting the second from the first.
     *
     * @param a1
     * @param a2
     * @return
     */
    public static EdgeAttribute difference(EdgeAttribute a1, EdgeAttribute a2) {
        if (a1 == null && a2 == null) {
            System.out.println("Calculating difference between two null edge attributes; returning null");
            return null;
        }

        if (a1 == null) {
            return a2;
        }

        if (a2 == null) {
            return a1;
        }

        // XXX: might be able to get rid of the next two if statements.
        if (a1.vect == null) {
            return a2;
        }
        if (a2.vect == null) {
            return a1;
        }

        if (a1.vect.length != a2.vect.length) {
            System.err.println("Error:  vectors different lengths when computing difference of EdgeAttributes " + a1 + " and " + a2);
            System.exit(1);
        }
        int diffLength = a1.vect.length;

        EdgeAttribute diff = new EdgeAttribute();
        diff.vect = new double[diffLength];

        for (int i = 0; i < diffLength; i++) {
            diff.vect[i] = a1.vect[i] - a2.vect[i];
        }

        return diff;
    }

    /**
     * Adds the two attributes.
     *
     * @param a1
     * @param a2
     */
    public static EdgeAttribute add(EdgeAttribute a1, EdgeAttribute a2) {
        if (a1.vect == null) {
            return a2;
        }
        if (a2.vect == null) {
            return a1;
        }

        if (a1.vect.length != a2.vect.length) {
            System.err.println("Error:  vectors different lengths when adding EdgeAttributes " + a1 + " and " + a2);
            System.exit(1);
        }


        int vectLength = a1.vect.length;

        EdgeAttribute sum = new EdgeAttribute();
        sum.vect = new double[vectLength];

        for (int i = 0; i < vectLength; i++) {
            sum.vect[i] = a1.vect[i] + a2.vect[i];
        }

        return sum;

    }

    /**
     * Find the specified point (given by position) on the line between start and target.
     * position is between 0 and 1, where position = 0 returns the start and
     * position = 1 returns the target.
     * Pass in null to compute distance to/from the 0 attribute.
     *
     * @param start
     * @param target
     * @param position
     * @return
     */
    public static EdgeAttribute weightedPairAverage(EdgeAttribute start, EdgeAttribute target, double position) {
        if (start == null && target == null) {
            System.out.println("Calculating point between two null edge attributes; returning null");
            return null;
        }

        if (start == null) {
            start = EdgeAttribute.zeroAttribute(target.size());
        }

        if (target == null) {
            target = EdgeAttribute.zeroAttribute(start.size());
        }

        if (start.vect.length != target.vect.length) {
            System.err.println("Error calculating point between edge attributes.  Edge attributes are different lenghts: " + start + " and " + target);
            System.exit(1);
        }

        if (position < 0 || position > 1) {
            System.err.println("Error calculating point between edge attributes:  position " + position + " must be between 0 and 1");
        }

        if (start.equals(target)) {
            return start;
        }

        EdgeAttribute point = new EdgeAttribute(new double[start.vect.length]);

        for (int i = 0; i < start.vect.length; i++) {
            point.vect[i] = (1 - position) * start.vect[i] + position * target.vect[i];
        }
        return point;
    }


    /**
     * Scales each of the elements of vect by a.
     *
     * @param s
     */
    public void scaleBy(double a) {
        for (int i = 0; i < vect.length; i++) {
            vect[i] = vect[i] * a;
        }
    }

    /**
     * Returns of the "size" of the EdgeAttribute.
     * In this case (EdgeAttribute is double[]), the size is the length of vect.
     *
     * @return
     */
    public int size() {
        if (vect == null) {
            return 0;
        }
        return vect.length;
    }

    /**
     * If vect has length 1, make sure it is positive.
     */
    public void ensurePositive() {
        if (vect.length == 1) {
            vect[0] = Math.abs(vect[0]);
        }
    }

    /**
     * Returns the EdgeAttribute with vector length length,
     * with all entries 0.
     *
     * @param length
     * @return
     */
    public static EdgeAttribute zeroAttribute(int size) {
        if (size < 1) {
            System.err.println("Error creating zero edge attribute of size " + size + "; invalid size");
            System.exit(1);
        }

        EdgeAttribute zero = new EdgeAttribute(new double[size]);
        Arrays.fill(zero.vect, 0.0);
        return zero;
    }
}