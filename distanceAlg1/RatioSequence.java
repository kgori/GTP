/** This file is part of GeodeMAPS and GTP, programs for computing the geodesic distance between phylogenetic trees.
 Copyright (C) 2008, 2009  Megan Owen, Scott Provan

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>. */


package distanceAlg1;

import java.text.DecimalFormat;
import java.util.*;

/**
 * An instance of this class represents a sequence of ratios.
 *
 * @author maowen
 */
public class RatioSequence extends Vector<Ratio> {
    private static final long serialVersionUID = 42L;  // to get rid of warning that serialVersionUID is not declared,
    // despite RatioSequence being a serializable class (inherited from Vector)

    private int combineCode = 0;    // represents a binary number where a 1 in the i-th position means
    // the i-th and (i+1)-th ratios have been combined.

    /**
     * Constructor
     */
    public RatioSequence() {
    }

    /**
     * Constructor
     *
     * @param ptA = string representing the coordinates of the first point ie. "1,2,3,"
     * @param ptB = string representing the negation of the coordinates of the target point "1,1,1,"
     */
    public RatioSequence(String ptA, String ptB) {
        // convert the two points into a ratio sequence
        while (!ptA.equals("")) {
            this.add(new Ratio(Double.valueOf(ptA.substring(0, ptA.indexOf(","))), Double.valueOf(ptB.substring(0, ptB.indexOf(",")))));
            ptA = ptA.substring(ptA.indexOf(",") + 1);
            ptB = ptB.substring(ptB.indexOf(",") + 1);
        }
    }

/*	public RatioSequence clone(RatioSequence v) {
        if (v == null) {
			return null;
		}

		Vector newV = new Vector();
		for (int i = 0; i < v.size(); i++) {
			if (v.get(i) == null) {
				newV.add(null);
			}
			else {
				newV.add(v.get(i).clone());
			}
		}
		return newV;
	}
*/

    public void setCombineCode(int c) {
        combineCode = c;
    }

    public int getCombineCode() {
        return combineCode;
    }

    /**
     * Returns the i-th ratio in the ratio sequence.
     *
     * @param i
     * @return
     */
    public Ratio getRatio(int i) {
        return this.get(i);
    }

    /**
     * Sets the i-th ratio in the ratio sequence to
     * public Ratio setRatio(int i, Ratio r) {
     * <p/>
     * }
     * <p/>
     * <p/>
     * /** Returns true if the ratio sequence rs is ascending.  Namely, the first ratio is <= to
     * the second ratio, which is <= the third ratio, etc.
     *
     * @param rs
     * @return
     */
    public boolean isAscending() {
        for (int i = 0; i < this.size() - 1; i++) {
            if (((Ratio) this.get(i)).getRatio() > ((Ratio) this.get(i + 1)).getRatio()) {
                return false;
            }
        }
        return true;
    }

    /**
     * Calculates the distance for the ratio sequence rs:
     * e1/f1, e2/f2, ... en/fn gives the distance sqrt( (e1 + f1)^2 + (e2 + f2)^2 + ... + (en + fn)^2 )
     *
     * @param rs
     * @return
     */
    public double getDistance() {
        double distSqrd = 0;

        for (int i = 0; i < size(); i++) {
            Ratio r = getRatio(i);
            distSqrd = distSqrd + Math.pow(r.getELength() + r.getFLength(), 2);
        }

/*		Iterator<Ratio> rsIter = this.iterator();

		while (rsIter.hasNext()) {
			Ratio r = rsIter.next();
			distSqrd = distSqrd + Math.pow(r.getELength() + r.getFLength(), 2);
		}*/
        return Math.sqrt(distSqrd);
    }

    /**
     * Gets the distance for this ratio sequence after it has been combined to be non-descending.
     *
     * @return
     */
    public double getMinNonDesRSDistance() {
        return getNonDesRSWithMinDist().getDistance();
    }


    /**
     * Interleaves the ratio sequences rs1 and rs2 after combining them to get
     * the ascending ratio sequence with the min distance,
     * to make a new ratio sequence.
     *
     * @param rs1
     * @param rs2
     * @return
     */
    public static RatioSequence interleave(RatioSequence rs1, RatioSequence rs2) {
        RatioSequence combined1 = rs1.getNonDesRSWithMinDist();
        RatioSequence combined2 = rs2.getNonDesRSWithMinDist();

        RatioSequence interleavedRS = new RatioSequence();
        int index1 = 0;  // index for stepping through combined1
        int index2 = 0;  // index for stepping through combined2
        while (index1 < combined1.size() && index2 < combined2.size()) {
            if (combined1.get(index1).getRatio() <= combined2.get(index2).getRatio()) {
                interleavedRS.add(combined1.get(index1));
                index1++;
            } else {
                interleavedRS.add(combined2.get(index2));
                index2++;
            }
        }
        // if we have finished adding rs2 but not rs1
        while (index1 < combined1.size()) {
            interleavedRS.add(combined1.get(index1));
            index1++;
        }
        // if we have finished adding rs1 but not rs2
        while (index2 < combined2.size()) {
            interleavedRS.add(combined2.get(index2));
            index2++;
        }

        return interleavedRS;
    }

    public RatioSequence clone() {
        if (this == null) {
            return null;
        }

        RatioSequence newRS = new RatioSequence();
        for (int i = 0; i < this.size(); i++) {
            if (this.get(i) == null) {
                newRS.add(null);
            } else {
                newRS.add(this.get(i).clone());
            }
        }
        newRS.setCombineCode(this.getCombineCode());
        return newRS;
    }

    /**
     * Returns a randomly generated RatioSequence with d ratios,
     * each with a numerator and denominator randomly generated between 0 and 1.
     *
     * @param d
     * @return
     */
    public static RatioSequence getRandomRS(int dim) {
        String ptA = "";
        String ptB = "";

        //generate points randomly.  # of coordinates = dimension
        for (int i = 0; i < dim; i++) {
            ptA = ptA + (Math.random() + 0.1) / 1.1 + ",";
            ptB = ptB + (Math.random() + 0.1) / 1.1 + ",";
        }

        return new RatioSequence(ptA, ptB);
    }

    /**
     * Returns the ratio sequence derived from rs by combining the ratios in rs according to combineCode.
     * In particular, consider combineCode in binary.  A 1 in the i-th position from the right means the
     * (i-1)th and i-th ratios should be combined.
     * TODO: add error checking
     *
     * @param rs
     * @param combineCode
     * @return
     */
    public RatioSequence getCombinedRS(int combineCode) {
        RatioSequence combinedRS = new RatioSequence(); // new ratio sequence formed by combining.
        Ratio ratioToCombineWith = this.get(0); // the ratio that we will combine the next ratio with if it is indicated we should combine the ratios

        // loop though the rs.size() - 1 ways to combine the ratios
        // each is represented as 2^j
        for (int j = 0; j < this.size() - 1; j++) {
            // if i contains j
            if ((combineCode & (int) Math.pow(2, j)) == (int) Math.pow(2, j)) {
                //combine ratio j and ratio j+1
                ratioToCombineWith = Ratio.combine(ratioToCombineWith, this.get(j + 1));
            } else {
                // move ratio j (and whatever earlier ratios it is combined with) into the vector combinedRS
                combinedRS.add(ratioToCombineWith);
                // ratio j+1 becomes the ratioToCombineWith
                ratioToCombineWith = this.get(j + 1);
            }
        }
        // add last ratio to the ratio sequence
        combinedRS.add(ratioToCombineWith);

        return combinedRS;
    }


    /**
     * Combines the ratios in rs so they are non-descending ie. so e1/f1 <= e2/f3 <= ... <= en/fn
     * Use the following algorithm: At the first pair of ratios which are descending, combine.
     * Compare the combined ratio with the previous one, and combine if necessary, etc.
     * <p/>
     * If the ratio sequence we are running this on has size < 2, return that ratio sequence.
     * XXX: combine code calculated wrong when have to compare with previously processed ratios (I think I fixed this???)
     *
     * @param rs
     * @return
     */
    public RatioSequence getNonDesRSWithMinDist() {
        if (this.size() < 2) {
            return this;
        }
        RatioSequence combinedRS = this.clone();
        int i = 0;   // index for stepping through rs
        int combineCode = 0;
        int ccArray[] = new int[this.size() - 1];   // array storing which ratios have been combined
        Arrays.fill(ccArray, 2);  // initialize array to 0.  Will mark a 1 if ratios combined.

        int a = 0; // array index

        while (i < combinedRS.size() - 1) {
            if (combinedRS.get(i).getRatio() > combinedRS.get(i + 1).getRatio()) {
                Ratio combinedRatio = Ratio.combine(combinedRS.get(i), combinedRS.get(i + 1));
                combinedRS.remove(i);
                combinedRS.remove(i);
                combinedRS.add(i, combinedRatio);
                ccArray[a] = 1;
                if (i > 0) {
                    i--;
                    // go back in the array to the last non-combined ratio
                    while (ccArray[a] == 1) {
                        a--;
                    }
                } else {
                    // we must advance a
                    while ((a < this.size() - 1) && (ccArray[a] != 2)) {
                        a++;
                    }
                }
            } else {
                ccArray[a] = 0;
                // the ratios are not-descending, so go on to the next pair
                i++;
                // we must jump ahead to the next ratio not considered for combination
                while ((a < this.size() - 1) && (ccArray[a] != 2)) {
                    a++;
                }
            }
        }

        for (int k = 0; k < this.size() - 1; k++) {
            if (ccArray[k] == 1) {
                combineCode = combineCode + (int) Math.pow(2, k);
            }
        }
        combinedRS.setCombineCode(combineCode);

        return combinedRS;
    }

    /**
     * Combines the ratios in rs so they are strictly ascending ie. so e1/f1 < e2/f3 < ... < en/fn
     * Use the following algorithm: At the first pair of ratios which are descending, combine.
     * Compare the combined ratio with the previous one, and combine if necessary, etc.
     * <p/>
     * If the ratio sequence we are running this on has size < 2, return that ratio sequence.
     * XXX: combine code calculated wrong when have to compare with previously processed ratios (I think I fixed this???)
     *
     * @param rs
     * @return
     */
    public RatioSequence getAscRSWithMinDist() {
        if (this.size() < 2) {
            return this;
        }
        RatioSequence combinedRS = this.clone();
        int i = 0;   // index for stepping through rs
        int combineCode = 0;
        int ccArray[] = new int[this.size() - 1];   // array storing which ratios have been combined
        Arrays.fill(ccArray, 2);  // initialize array to 0.  Will mark a 1 if ratios combined.

        int a = 0; // array index

        while (i < combinedRS.size() - 1) {
            if (combinedRS.get(i).getRatio() >= combinedRS.get(i + 1).getRatio()) {
                Ratio combinedRatio = Ratio.combine(combinedRS.get(i), combinedRS.get(i + 1));
                combinedRS.remove(i);
                combinedRS.remove(i);
                combinedRS.add(i, combinedRatio);
                ccArray[a] = 1;
                if (i > 0) {
                    i--;
                    // go back in the array to the last non-combined ratio
                    while (ccArray[a] == 1) {
                        a--;
                    }
                } else {
                    // we must advance a
                    while ((a < this.size() - 1) && (ccArray[a] != 2)) {
                        a++;
                    }
                }
            } else {
                ccArray[a] = 0;
                // the ratios are not-descending, so go on to the next pair
                i++;
                // we must jump ahead to the next ratio not considered for combination
                while ((a < this.size() - 1) && (ccArray[a] != 2)) {
                    a++;
                }
            }
        }

        for (int k = 0; k < this.size() - 1; k++) {
            if (ccArray[k] == 1) {
                combineCode = combineCode + (int) Math.pow(2, k);
            }
        }
        combinedRS.setCombineCode(combineCode);

        return combinedRS;
    }


    /**
     * Reverses the ratio sequence, including flipping the ratios.
     * So if it represents the  ratio sequence between T1 and T2,
     * we return the ratio sequence between T2 and T1.
     *
     * @return
     */
    public RatioSequence reverse() {
        RatioSequence revRS = new RatioSequence();
        for (int i = size() - 1; i >= 0; i--) {
            revRS.add(get(i).reverse());
        }
        return revRS;
    }

    /**
     * Returns the ratio sequence as just the ratio values.
     *
     * @return
     */
    public String toStringValueAndRatio() {
        DecimalFormat d4o = new DecimalFormat("#0.####");
        String s = "[ ";

        for (int i = 0; i < this.size(); i++) {
            s = s + d4o.format(this.get(i).getRatio()) + " = " + this.get(i) + " , ";
        }
        // we take a substring of s to not return the last " ,"
        return s.substring(0, s.length() - 2) + "]";
    }

    public String toStringValue() {
        DecimalFormat d4o = new DecimalFormat("#0.####");
        String s = "[ ";

        for (int i = 0; i < this.size(); i++) {
            s = s + d4o.format(this.get(i).getRatio()) + " , ";
        }
        // we take a substring of s to not return the last " ,"
        return s.substring(0, s.length() - 2) + "]";
    }

    public String toStringVerbose(Vector<String> leaf2NumMap) {
        String s = "";
        for (int i = 0; i < this.size(); i++) {
            s = s + "\nRatio " + i + ":  ";
            s = s + get(i).toStringVerbose(leaf2NumMap);
        }

        return s;
    }

    /*  Returns the combinatorial type of the geodesic, in a condensed format.
     *
     */
    public String toStringCombType() {
        String s = "";

        for (int i = 0; i < this.size(); i++) {
            s = s + get(i).toStringCombType() + ";";
        }
        return s;
    }

    public String toStringCombTypeAndValue() {
        String s = "";

        for (int i = 0; i < this.size(); i++) {
            s = s + get(i).toStringCombTypeAndValue() + ";";
        }
        return s;
    }
}

