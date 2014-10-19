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

import java.util.*;

public class PhyloTreeEdge extends Bipartition {
    private EdgeAttribute attribute;
    private Bipartition originalEdge;
    private int originalID;   // unique identifier based on order split is read in from tree string;
    // -1 means it hasn't been assigned and if the split is derived from an original split,
    // it will have the same ID as the original split

    // TODO: set length and originalID to be null if not given.  But need to make sure this doesn't break anything else.


    /**
     * Constructor: creates a new split with no child leaves and leave length as null
     */
    public PhyloTreeEdge() {
        System.out.println("PhyloTreeEdge::constructor 1");
        super();
        originalEdge = new Bipartition();
        originalID = -1;
    }

    /**
     * Constructor: creates a new split with the specified child leaves and leave length as null
     *
     * @param partition
     */
    public PhyloTreeEdge(BitSet edge) {
        System.out.println("PhyloTreeEdge::constructor 2");
        super(edge);
        originalEdge = new Bipartition();
        originalID = -1;
//		System.out.println("split " + split);
    }

    public PhyloTreeEdge(EdgeAttribute attrib) {
        System.out.println("PhyloTreeEdge::constructor 3");
        super();
        this.attribute = attrib;
        originalEdge = new Bipartition();
        originalID = -1;
//		System.out.println("double length ");
    }


    public PhyloTreeEdge(EdgeAttribute attrib, int originalID) {
        System.out.println("PhyloTreeEdge::constructor 4");
        super();
        this.attribute = attrib;
        originalEdge = new Bipartition();
        this.originalID = originalID;
    }

    public PhyloTreeEdge(EdgeAttribute attrib, Bipartition originalEdge, int originalID) {
        System.out.println("PhyloTreeEdge::constructor 5");
        super();
        this.attribute = attrib;
        this.originalEdge = new Bipartition(originalEdge.partition);
        this.originalID = originalID;
    }


    public PhyloTreeEdge(Bipartition edge, EdgeAttribute attrib, int originalID) {
        System.out.println("PhyloTreeEdge::constructor 6");
        super(edge.partition);
        this.attribute = attrib;
        originalEdge = new Bipartition(edge.partition);
        this.originalID = originalID;
    }

    public PhyloTreeEdge(BitSet edge, EdgeAttribute attrib, BitSet originalEdge, int originalID) {
        System.out.println("PhyloTreeEdge::constructor 7");
        super(edge);
        this.attribute = attrib;
        this.originalEdge = new Bipartition(originalEdge);
        this.originalID = originalID;
    }

    // Getters and Setters
    // Return
    public double getLength() {
        System.out.println("PhyloTreeEdge::getLength");
        return attribute.norm();
    }

    /**
     * Returns true if the length of the edge is 0.
     */
    public boolean isZero() {
        System.out.println("PhyloTreeEdge::isZero");
        return (this.getLength() == 0);        //XXX: Change to make dependent on tolerance.
    }

    /** String representation of an split.
     * Returns the length followed by a 0-1 vector, representing the children.
     */
/*	public String toString() {
        return "" + TreeDistance.truncate(length,6) + " " + super.toString();
	}
*/

    /**
     * String representation of an split.
     */
    public String toString() {
        System.out.println("PhyloTreeEdge::toString");
        return "" + attribute + " " + partition;
    }

    public String toStringVerbose(Vector<String> leaf2NumMap) {
        System.out.println("PhyloTreeEdge::toStringVerbose");
        return "" + originalID + "\t\t" + attribute + "\t\t" + Bipartition.toStringVerbose(this.partition, leaf2NumMap);
    }


    /**
     * Returns a string of what it printed out.
     *
     * @param edges
     * @param leaf2NumMap
     * @return
     */
    public static String printEdgesVerbose(Vector<PhyloTreeEdge> edges, Vector<String> leaf2NumMap, Boolean originalEdges) {
        System.out.println("PhyloTreeEdge::printEdgesVerbose");
        String output = "";

        System.out.println("Edge ID\t\tLength\t\tLeaves Below");
        output = output + "Edge ID\t\tLength\t\tLeaves Below\n";
        for (int i = 0; i < edges.size(); i++) {
            if (edges.get(i).getOriginalEdge().getPartition() != null && originalEdges) {
//				System.out.println(edges.get(i).getOriginalID() + "\t\t" + TreeDistance.truncate(edges.get(i).length, 8) + "\t\t" +  Bipartition.toStringVerbose(edges.get(i).getOriginalEdge().getEdge(), leaf2NumMap));
                System.out.println(edges.get(i).getOriginalID() + "\t\t" + edges.get(i).attribute + "\t\t" + Bipartition.toStringVerbose(edges.get(i).getOriginalEdge().getPartition(), leaf2NumMap));
                //			output = output + edges.get(i).getOriginalID() + "\t\t" + TreeDistance.truncate(edges.get(i).length, 8) + "\t\t" +  Bipartition.toStringVerbose(edges.get(i).getOriginalEdge().getEdge(), leaf2NumMap) + "\n";
                output = output + edges.get(i).getOriginalID() + "\t\t" + edges.get(i).attribute + "\t\t" + Bipartition.toStringVerbose(edges.get(i).getOriginalEdge().getPartition(), leaf2NumMap) + "\n";
            } else {
                System.out.println(edges.get(i).toStringVerbose(leaf2NumMap));
                output = output + edges.get(i).toStringVerbose(leaf2NumMap) + "\n";
            }
        }
        return output;
    }

    // TODO:  currently not overriding the object clone because doesn't return type Object.
    // Also, can not use constructor.
    public PhyloTreeEdge clone() {
        System.out.println("PhyloTreeEdge::clone");
        return new PhyloTreeEdge((BitSet) this.partition.clone(), (EdgeAttribute) this.attribute.clone(), (BitSet) this.originalEdge.getPartition().clone(), this.originalID);
    }

    @Override
    public boolean equals(Object e) {
        System.out.println("PhyloTreeEdge::equals");
        if (e == null) {
            return false;
        }
        if (this == e) {
            return true;
        }

        if (!(e instanceof PhyloTreeEdge)) {
            return false;
        }

        return partition.equals(((Bipartition) e).partition) && attribute.equals(((PhyloTreeEdge) e).attribute);
    }

    public boolean sameBipartition(PhyloTreeEdge e) {
        System.out.println("PhyloTreeEdge::sameBipartition v1");
        return this.partition.equals(e.partition);
    }

    public boolean sameBipartition(Bipartition e) {
        System.out.println("PhyloTreeEdge::sameBipartition v2");
        return this.partition.equals(e.partition);
    }

    /**
     * Already "cloned".
     *
     * @return
     */
    public Bipartition asSplit() {
        System.out.println("PhyloTreeEdge::asSplit");
        return new Bipartition((BitSet) this.partition.clone());
    }


    public Bipartition getOriginalEdge() {
        System.out.println("PhyloTreeEdge::getOriginalEdge");
        return originalEdge;
    }

    public void setOriginalEdge(Bipartition originalEdge) {
        System.out.println("PhyloTreeEdge::setOriginalEdge");
        this.originalEdge = originalEdge;
    }

    public int getOriginalID() {
        System.out.println("PhyloTreeEdge::getOriginalID");
        return originalID;
    }

    public void setOriginalID(int originalID) {
        System.out.println("PhyloTreeEdge::setOriginalID");
        this.originalID = originalID;
    }

    public void setAttribute(EdgeAttribute attrib) {
        System.out.println("PhyloTreeEdge::setAttribute");
        this.attribute = attrib;
    }

    public EdgeAttribute getAttribute() {
        System.out.println("PhyloTreeEdge::getAttribute");
        return attribute;
    }
}
