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

// TODO:  create sameAs method and fix getCommonEdges

import java.util.*;

public class Bipartition implements Cloneable {

    protected BitSet partition;

    public Bipartition() {
        System.out.println("Bipartition constructor 1");
        partition = new BitSet();
    }

    public Bipartition(BitSet edge) {
        System.out.println("Bipartition constructor 2");
        this.partition = edge;
    }

    public Bipartition(String s) {
        System.out.println("Bipartition constructor 3");
        partition = new BitSet();

        for (int i = 0; i < s.length(); i++)
            if (s.charAt(i) == '1') {
                partition.set(i);
            } else if (s.charAt(i) != '0') {
                throw new RuntimeException("Error creating bipartition: input string " + s + " should only contain 0s and 1s");
            }
    }

    public BitSet getPartition() {
        System.out.println("getPartition");
        return partition;
    }

    public boolean isEmpty() {
        System.out.println("isEmpty");
/*		if (split.isEmpty()) {
            System.out.println("in isZero, split is " + this.edge + " and will return " + split.isEmpty());
		}*/
        return partition.isEmpty();
    }

    public void setPartition(BitSet edge) {
        System.out.println("setPartition");
        this.partition = edge;
    }

    /**
     * Adds a 1 to the split, or changes an element from one block to the other in the partition.
     * Adds the 1 at coordinate one, from the right.
     * Precondition:  this split contains a 0 at position one
     *
     * @param one
     */
    public void addOne(int one) {
        System.out.println("addOne");
//		System.out.println("split is " + split);
        partition.set(one);
//		System.out.println("split with flipped bit in " + one + " is " + split);
    }

    /**
     * Removes a 1 from the split.  In other words, transfers the element one to the other (zero)
     * block in the partition.
     * Precondition:  this split contains a 1 at position one
     *
     * @param leafNum
     */
    public void removeOne(int one) {
        System.out.println("removeOne");
        partition.clear(one);
    }

    /**
     * = e1 and e2 are disjoint edges
     * Precondition:  the leaf2NumMap for the trees containing e1 and e2 are the same.
     *
     * @param e1
     * @param e2
     * @return
     */
    public boolean disjointFrom(Bipartition e) {
        System.out.println("disjointFrom");
//		return (this.edge.and(e.edge).compareTo(BigInteger.ZERO) == 0 );
//		System.out.println("" + this.edge + " and " + e.edge + " are disjoint: " + !this.edge.intersects(e.edge));
        return !this.partition.intersects(e.partition);
//		return (this.edge & e.edge) == 0;
    }

    /**
     * = this split contains e
     * Precondition:  the leaf2NumMap for the tree containing this split and the tree containing e are the same.
     * XXX:  are we changing e outside of this method?
     *
     * @param e
     * @return
     */
    public boolean contains(Bipartition e) {
        System.out.println("contains v 1");
        BitSet edgeClone = (BitSet) e.partition.clone();
        edgeClone.and(this.partition);
//		System.out.println("edgeClone is " + edgeClone);
//		System.out.println("" + this.edge + " contains " + e.edge + "? " + edgeClone.equals(e.edge) );
        return edgeClone.equals(e.partition);
//		return (this.edge.and(e.edge).compareTo(e.edge) == 0);
//		return (split & e.edge) == e.edge;
    }

    /**
     * = this split contains a 1 in the bit corresponding to the i-th column
     *
     * @param i
     * @return
     */
    public boolean contains(int i) {
        System.out.println("contains v 2");
//		System.out.println("" + this.edge + " contains split " + i + ": " + this.edge.get(i));
        return this.partition.get(i);
        // return this.edge.testBit(i);
//		return this.contains( new Bipartition((long) Math.pow(2,i)) );
    }

    /**
     * = this split properly contains e  (ie. contains e, but is not equal to e)
     *
     * @param e
     * @return
     */
    public boolean properlyContains(Bipartition e) {
        System.out.println("properlyContains");
        return this.contains(e) && !e.contains(this);
    }

    /**
     * = this split crosses e
     *
     * @param e
     * @return
     */
    public boolean crosses(Bipartition e) {
        System.out.println("crosses");
        return !(disjointFrom(e) || this.contains(e) || e.contains(this));
    }

    /** Returns true if e contains ones exactly where this edges doesn't.
     * numLeaves gives the total number of ones each should contain.
     * @param e
     * @param numLeaves
     * @return
     */
/*	public boolean isComplementOf(Bipartition e, int numLeaves) {
//		System.out.println("Sum of edges " + e.edge + " and " + this.edge + " is " + (e.edge + this.edge) + " and want " + ((long)Math.pow(2,numLeaves) - 1));
		if ( (e.edge + this.edge) == ((long)Math.pow(2,numLeaves) - 1) )  {
			return true;
		}
		else {
			return false;
		}
	}*/

    /**
     * Changes e to be its complement.
     * XXX:  I think this is a slow method; make faster
     *
     * @param e
     * @param numLeaves
     * @return
     */
    public void complement(int numLeaves) {
        System.out.println("complement");
//		this.edge = (long)Math.pow(2,numLeaves) - 1 - this.edge;
        for (int i = 0; i < numLeaves; i++) {
            this.partition.flip(i);
        }
    }

    /**
     * Coverts split into a 0-1 vector.  XXX: leading zeros not shown.
     * XXX:  maybe fix???
     */
    public String toString() {
        System.out.println("toString");
//		return Long.toBinaryString(split);
//		return "" + this.edge.toByteArray();
        return this.partition.toString();
    }

    public static String toStringVerbose(BitSet edge, Vector<String> leaf2NumMap) {
        System.out.println("toStringVerbose");
/*		String toDisplay = "";

		String s = Long.toBinaryString(split);
		// we now have a 0-1 string, but it is reversed.  ie. 100 should become {2} and 1100 should become {3,2}
		for (int i = s.length()-1; i >=0; i--) {
			if (s.charAt(i) == '1') {
				toDisplay = toDisplay + leaf2NumMap.get(s.length() -1 - i) + "*";
			}
		}
		// remove the last *
		return toDisplay.substring(0, toDisplay.length() - 1); */
        String toDisplay = "";
        for (int i = 0; i < edge.length(); i++) {
            if (edge.get(i)) {
                toDisplay = toDisplay + leaf2NumMap.get(i) + ",";
            }
        }
        // remove the last ,
        return toDisplay.substring(0, toDisplay.length() - 1);
//		return toDisplay;
    }

    /**
     * XXX:  check that replacing the old equals with the new one that actually overrides the object equals
     * doesn't affect the code.
     *
     * @param e
     * @return
     */

    @Override
    public boolean equals(Object e) {
        System.out.println("equals");
        if (e == null) {
            return false;
        }
        if (this == e) {
            return true;
        }

        if (!(e instanceof Bipartition) || e instanceof PhyloTreeEdge) {
            return false;
        }

        return partition.equals(((Bipartition) e).partition);
    }

    //TODO:  currently this method does NOT override the object one, which has header "public Object clone()"
    // Do we care?  It seems like there are problems with making clones, and one could just use a copy constructor.
    public Bipartition clone() {
        System.out.println("clone");
        return new Bipartition((BitSet) partition.clone());
    }


    public boolean isCompatibleWith(Vector<Bipartition> splits) {
        boolean compatible = true;
        for (int i = 0; i < splits.size(); i++) {
            if (this.crosses(splits.get(i))) {
                compatible = false;
            }
        }
        return compatible;
    }

}
