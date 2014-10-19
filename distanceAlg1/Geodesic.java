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

import java.io.IOException;
import java.util.*;

import polyAlg.PolyMain;

public class Geodesic {
    private RatioSequence rs;
    private Vector<PhyloTreeEdge> commonEdges;   // stores the common edges where their length is the difference in lengths of that edge between the two trees
    private double leafContributionSquared = 0;

    // constructors
    public Geodesic(RatioSequence rs) {
        System.out.println("Geodesic::constructor1");
        this.rs = rs;
        //		this.rs = rs.getAscRSWithMinDist();
        commonEdges = new Vector<PhyloTreeEdge>();
//		dist = this.rs.getDistance();
    }

    public Geodesic(RatioSequence rs, Vector<PhyloTreeEdge> cEdges) {
        System.out.println("Geodesic::constructor2");
        this.rs = rs;
        //		this.rs = rs.getAscRSWithMinDist();
        commonEdges = cEdges;
    }

    public Geodesic(RatioSequence rs, Vector<PhyloTreeEdge> cEdges, double leafContributionSquared) {
        System.out.println("Geodesic::constructor3");
        this.rs = rs;
        //		this.rs = rs.getAscRSWithMinDist();
        commonEdges = cEdges;
        this.leafContributionSquared = leafContributionSquared;
    }

    public RatioSequence getRS() {
        System.out.println("Geodesic::getRS");
        return rs;
    }

    public void setRS(RatioSequence rs) {
        System.out.println("Geodesic::setRS");
        this.rs = rs;
    }

    public double getDist() {
        System.out.println("Geodesic::getDist");
        double commonEdgeDistSquared = 0;
        for (int i = 0; i < commonEdges.size(); i++) {
            commonEdgeDistSquared = commonEdgeDistSquared + Math.pow(commonEdges.get(i).getLength(), 2);
        }
/*		if (TreeDistance.verbose > 0) {
            System.out.println("rs.getAscRSWithMinDist().getDistance() is " + rs.getAscRSWithMinDist().getDistance() + "; commonEdgeDistSquared is " +commonEdgeDistSquared + "; and leafContributionSquared is " + leafContributionSquared);
		} */
        return Math.sqrt(Math.pow(rs.getNonDesRSWithMinDist().getDistance(), 2) + commonEdgeDistSquared + leafContributionSquared);
    }


    /**
     * Adds split to list of common edges.  Assumes the length is the change in this common split
     * between the two trees.
     */
    public void addCommonEdge(PhyloTreeEdge e) {
        System.out.println("Geodesic::addCommonEdge");
        commonEdges.add(e);
    }

    /**
     * Displays the geodesic in user-friendly form.
     *
     * @return
     */
    public String toStringVerboseOld(PhyloTree t1, PhyloTree t2) {
        System.out.println("Geodesic::toStringVerboseOld");
        Vector<PhyloTreeEdge> commonEdges = this.getCommonEdges();
        Boolean cEdge = false;

        // display T1 only splits
        String toDisplay = "\nSplits only in T1:\n";
        for (int i = 0; i < t1.getEdges().size(); i++) {
            cEdge = false;
            for (int j = 0; j < commonEdges.size(); j++) {
                if (commonEdges.get(j).sameBipartition(t1.getEdge(i))) {
                    cEdge = true;
                    break;
                }
            }
            if (!cEdge) {
                toDisplay = toDisplay + t1.getEdge(i).toStringVerbose(t1.getLeaf2NumMap()) + "\n" + t1.getEdge(i) + "\n";
            }
        }

        //display T2 only splits
        toDisplay = toDisplay + "\n\nSplits only in T2:\n";
        for (int i = 0; i < t2.getEdges().size(); i++) {
            if (!commonEdges.contains(t2.getEdge(i))) {
                // if this is not a common split, display
//				toDisplay = toDisplay + t2.getEdge(i).toStringVerbose(t2.getLeaf2NumMap()) + "\n";
                toDisplay = toDisplay + t2.getEdge(i).toStringVerbose(t2.getLeaf2NumMap()) + "\n" + t2.getEdge(i) + "\n";
            }
        }

        // display common splits
        toDisplay = toDisplay + "\n\nCommon splits:\n";
        for (int i = 0; i < commonEdges.size(); i++) {
//			toDisplay = toDisplay + commonEdges.get(i).toStringVerbose(t1.getLeaf2NumMap()) + "\n";
            toDisplay = toDisplay + commonEdges.get(i).toStringVerbose(t1.getLeaf2NumMap()) + "\n" + commonEdges.get(i) + "\n";
        }

        return toDisplay;
    }

    public Geodesic clone() {
        System.out.println("Geodesic::clone");
        return new Geodesic(rs.clone(), TreeDistance.myVectorClonePhyloTreeEdge(commonEdges));
    }

    public String toString() {
        System.out.println("Geodesic::toString");
        return "" + getDist() + "; " + rs.getNonDesRSWithMinDist();
//		return "" + getDist() + "; " + rs;
    }

    public Vector<PhyloTreeEdge> getCommonEdges() {
        System.out.println("Geodesic::getCommonEdges");
        return commonEdges;
    }

    public void setCommonEdges(Vector<PhyloTreeEdge> commonEdges) {
        System.out.println("Geodesic::setCommonEdges");
        this.commonEdges = commonEdges;
    }

    public int numCommonEdges() {
        System.out.println("Geodesic::numCommonEdges");
        return commonEdges.size();
    }

    /**
     * Returns the number of orthants/topologies that the geodesic passes through (not including boundaries between orthants).
     * = # ratios in the strictly ascending ratio sequences + 1
     *
     * @return
     */
    public int numTopologies() {
        System.out.println("Geodesic::numTopologies");
        return rs.getAscRSWithMinDist().size() + 1;
    }

    /**
     * Returns the geodesic with the ratio sequence (and ratios) reversed.
     *
     * @return
     */
    public Geodesic reverse() {
        System.out.println("Geodesic::reverse");
        return new Geodesic(rs.reverse(), commonEdges, leafContributionSquared);
    }

    public double getLeafContributionSquared() {
        System.out.println("Geodesic::getLeafContributionSquared");
        return leafContributionSquared;
    }

    public void setLeafContributionSquared(double leafContributionSquared) {
        System.out.println("Geodesic::setLeafContributionSquared");
        this.leafContributionSquared = leafContributionSquared;
    }

    /**
     * Returns the tree at the given position (a number between 0 and 1) along the geodesic
     * between tree1 and tree2.
     *
     * @param tree1
     * @param tree2
     * @param position
     * @return
     */
    public static PhyloTree getTreeAt(PhyloTree t1, PhyloTree t2, double position) throws IOException {
        System.out.println("Geodesic::getTreeAt");
        int lowerRatioIndex = -2;    // the index of the ratio containing all f edges in the tree we want
        // i.e. the index of the ratio with time < position, but such that the next ratio has time >= position
        // if position is in the starting orthant, we don't want any f edges.
        int higherRatioIndex = -2;  // the index of the ratio containing all e edges in the tree we want
        // i.e. the index of the ratio with time >= position, but such that the next ratio has time >= position
        // if position is in the target orthant, we don't want any e edges
        Vector<PhyloTreeEdge> eEdges = new Vector<PhyloTreeEdge>();
        Vector<PhyloTreeEdge> fEdges = new Vector<PhyloTreeEdge>();

        // set the commonEdges
        Vector<PhyloTreeEdge> commonEdges = getCommonEdges(t1, t2, position);

        // remove all common edges from the two trees before finding the geodesic between them
        // (including ones that are not in the other tree, but compatible with all edges in that tree)
        // (but first clone t1 and t2, so we don't change the originals)
        PhyloTree t1NoCommonEdges = t1.clone();
        PhyloTree t2NoCommonEdges = t2.clone();

        for (int i = 0; i < commonEdges.size(); i++) {
            t1NoCommonEdges.removeSplit(commonEdges.get(i));
            t2NoCommonEdges.removeSplit(commonEdges.get(i));
        }

        //	System.out.println("t1NoCommonEdges is " + t1NoCommonEdges);
        //	System.out.println("t2NoCommonEdges is " + t2NoCommonEdges);

        Geodesic geo = PolyMain.getGeodesic(t1NoCommonEdges, t2NoCommonEdges, null);

        PhyloTree tree = new PhyloTree(commonEdges, t1.getLeaf2NumMap());

        // set the leaf lengths
        EdgeAttribute[] newLeafEdgeAttribs = new EdgeAttribute[t1.getLeafEdgeAttribs().length];
        for (int i = 0; i < newLeafEdgeAttribs.length; i++) {
            //			newLeafEdgeAttribs[i] = ((1-position)*t1.getLeafEdgeAttribs()[i] + position*t2.getLeafEdgeAttribs()[i]);
            newLeafEdgeAttribs[i] = EdgeAttribute.weightedPairAverage(t1.getLeafEdgeAttribs()[i], t2.getLeafEdgeAttribs()[i], position);

            //			System.out.println("t1 leaf length: " +  t1.getLeafEdgeLengths()[i] + "; t2 leaf length: " +  t2.getLeafEdgeLengths()[i] + "; new length: " + ((1-position)*t1.getLeafEdgeLengths()[i] + position*t2.getLeafEdgeLengths()[i]));

        }
        tree.setLeafEdgeAttribs(newLeafEdgeAttribs);

        if (geo.getRS().size() == 0) {
            // then we are done, because the two trees are in the same orthant
            return tree;
        }
        // figure out what orthant the new tree is in
        // first check if the new tree is in the starting orthant
        if (geo.getRS().getRatio(0).getTime() > position) {
            // new tree is in the interior of the starting orthant
            lowerRatioIndex = -1;
            //			System.out.println("in starting orthant: setting lower ratio index to -1.  First ratio time is " + geo.getRS().getRatio(0).getTime() + " and position is " + position);

            higherRatioIndex = 0;
            //			System.out.println("in starting orthant: setting higher ratio index to 0.  First ratio time is " + geo.getRS().getRatio(0).getTime() + " and position is " + position);

        }
        // if the new tree is in the last orthant
        else if (geo.getRS().getRatio(geo.getRS().size() - 1).getTime() < position) {
            lowerRatioIndex = geo.getRS().size() - 1;
            higherRatioIndex = geo.getRS().size();
            //			System.out.println("in target orthant: setting lower ratio index to " + (geo.getRS().size()-1) + ".  Final ratio time is " + geo.getRS().getRatio(geo.getRS().size()-1).getTime() + " and position is " + position);
            //			System.out.println("in target orthant: setting higehr ratio index to " + (geo.getRS().size()) + ".  Final ratio time is " + geo.getRS().getRatio(geo.getRS().size()-1).getTime() + " and position is " + position);

        }
        // the new tree is in an intermediate orthant
        else {
            for (int i = 0; i < geo.getRS().size(); i++) {
                // note:  want < instead of <= so we are in an orthant and not still on the boundary,
                // if we have a string of equalities
                double ratioTime = geo.getRS().getRatio(i).getTime();
                //				System.out.println("ratio: " + geo.getRS().getRatio(i));
                //				System.out.println("time: " + ratioTime);

                if ((lowerRatioIndex == -2) && (ratioTime >= position)) {
                    lowerRatioIndex = i - 1;
                    //					System.out.println("setting lower ratio index to " + (i-1) + ".  Time is " + ratioTime + " and position is " + position);
                }
                if ((higherRatioIndex == -2) && (lowerRatioIndex != -2) && (ratioTime > position)) {
                    higherRatioIndex = i;
                    //					System.out.println("setting higher ratio index to " + i + ".  Time is " + ratioTime + " and position is " + position);t hot
                }
            }
        }
        // if we didn't set the higherRatioIndex, then we are on the boundary with the target orthant.
        // we want all no e edges, so set higherRatioIndex to
        if (higherRatioIndex == -2) {
            higherRatioIndex = geo.getRS().size();
        }

        // add the edges for all f edges in ratios indexed <= lowerRatioIndex
        for (int i = 0; i <= lowerRatioIndex; i++) {
            fEdges = geo.getRS().getRatio(i).getFEdges();

            for (PhyloTreeEdge f : fEdges) {
                //				double newLength = ((position*geo.getRS().getRatio(i).getFLength() - (1-position)*geo.getRS().getRatio(i).getELength())/geo.getRS().getRatio(i).getFLength())*fEdges.get(j).getLength();
                EdgeAttribute newAttrib = f.getAttribute().clone();
                newAttrib.scaleBy((position * geo.getRS().getRatio(i).getFLength() - (1 - position) * geo.getRS().getRatio(i).getELength()) / geo.getRS().getRatio(i).getFLength());
                // don't have to clone newAttrib because a new object is created each iteration of this loop
                tree.addEdge(new PhyloTreeEdge(f.asSplit(), newAttrib, f.getOriginalID()));
                //				System.out.println("Added split with length " + newLength );
            }
        }

        // to the new tree, add the e edges in the ratios indexed >= higherRatioIndex
        for (int i = higherRatioIndex; i < geo.getRS().size(); i++) {
            eEdges = geo.getRS().getRatio(i).getEEdges();

            for (PhyloTreeEdge e : eEdges) {
                //				double newLength = ((1-position)*geo.getRS().getRatio(i).getELength() - position*geo.getRS().getRatio(i).getFLength())/geo.getRS().getRatio(i).getELength()*eEdges.get(j).getLength();
                EdgeAttribute newAttrib = e.getAttribute().clone();
                newAttrib.scaleBy(((1 - position) * geo.getRS().getRatio(i).getELength() - position * geo.getRS().getRatio(i).getFLength()) / geo.getRS().getRatio(i).getELength());
                tree.addEdge(new PhyloTreeEdge(e.asSplit(), newAttrib, e.getOriginalID()));
            }
        }
        return tree;
    }

    /**
     * Returns a vector containing edges with the same partition in t1 and t2.
     * If t1 and t2 do not have the same leaf2NumMap, then error.
     * If t1 and t2 have no common edges, then returns a vector of size 0.
     * We set the length of the returned edges to be (1-position)*length_in_t1 + position*length_in_t2.
     * Doesn't assume t1 and t2 have same # of edges.
     *
     * @param t1
     * @param t2
     * @return
     */
    public static Vector<PhyloTreeEdge> getCommonEdges(PhyloTree t1, PhyloTree t2, double position) {
        System.out.println("Geodesic::getCommonEdges");
        EdgeAttribute commonEdgeAttribute;
        Bipartition commonSplit;

        Vector<PhyloTreeEdge> commonEdges = new Vector<PhyloTreeEdge>();

        // if the two trees do not have the same leaf2NumMap
        if (!(t1.getLeaf2NumMap().equals(t2.getLeaf2NumMap()))) {
            System.err.println("Error: the two trees do not have the same leaves!");
            System.err.println("First tree's leaves are " + t1.getLeaf2NumMap());
            System.err.println("Second tree's leaves are " + t2.getLeaf2NumMap());
            throw new RuntimeException();
        }

        if (position < 0 || position > 1) {
            System.err.println("Error:  position " + position + " must be between 0 and 1");
            throw new RuntimeException();
        }

        // end error checking

        for (PhyloTreeEdge e : t1.getEdges()) {
            if (t2.getSplits().contains(e.asSplit())) {
                // then we have found the same split in each tree
                commonSplit = e.asSplit();
                commonEdgeAttribute = EdgeAttribute.weightedPairAverage(e.getAttribute(), t2.getAttribOfSplit(commonSplit), position);
                commonEdges.add(new PhyloTreeEdge(commonSplit.clone(), commonEdgeAttribute.clone(), -1));
            }
            // otherwise check if the split is compatible with all splits in t2
            else if (e.isCompatibleWith(t2.getSplits())) {
                commonEdgeAttribute = EdgeAttribute.weightedPairAverage(e.getAttribute(), null, position);
                commonEdges.add(new PhyloTreeEdge(e.asSplit(), commonEdgeAttribute.clone(), -1));
            }
        }
        // check for splits in t2 that are compatible with all splits in t1
        for (PhyloTreeEdge e : t2.getEdges()) {
            if (e.isCompatibleWith(t1.getSplits()) && !(t1.getSplits().contains(e.asSplit()))) {
                commonEdgeAttribute = EdgeAttribute.weightedPairAverage(null, e.getAttribute(), position);
                commonEdges.add(new PhyloTreeEdge(e.asSplit(), commonEdgeAttribute.clone(), -1));
            }
        }

        return commonEdges;
    }
}
