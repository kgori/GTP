/** This file is part of GTP, a program for computing the geodesic distance between phylogenetic trees.
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

package polyAlg;

import distanceAlg1.*;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

public class PolyMain {
    // stores pairs of trees with no common edges.  Should be reset at each new distance calculation
    public static Vector<PhyloTree> aTreesNoCommonEdges = new Vector<PhyloTree>();
    public static Vector<PhyloTree> bTreesNoCommonEdges = new Vector<PhyloTree>();

    //	public static boolean rooted = true;  //holds if the trees are rooted or not.
    public static boolean normalize = false;  // holds if we should normalize the tree split lengths

    public static int verbose = 0;

    public static String LEAF_CONTRIBUTION_SQUARED_DESCRIPTION = "(Leaf contribution squared = square of the length of the vector" +
            " whose i-th element is the absolute value of the difference between the length of the split ending in leaf i in the first tree" +
            " and the length of the split ending in leaf i in the second tree.)";

    /**
     * Stores subtrees with no common edges in the global variables aTreesNoCommonEdges (from t1)
     * and in bTreesNoCommonEdges (from t2).
     * <p/>
     * If one of the trees has 0 edges, then all edges in the other tree will be compatible with it.
     * Thus we will not get any more subtree pairs with disjoint leaves, and should return.
     */
    public static void splitOnCommonEdge(PhyloTree t1, PhyloTree t2) {
        int numEdges1 = t1.getEdges().size(); // number of edges in tree 1
        int numEdges2 = t2.getEdges().size(); /// number of edges in tree 2

        if (numEdges1 == 0 || numEdges2 == 0) {
            return;
        }
        // look for common edges
        Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1, t2);

        // if there are no common edges
        // XXX: need to check the following methods don't require the trees to have the same number of edges
        if (commonEdges.size() == 0) {
            aTreesNoCommonEdges.add(t1);
            bTreesNoCommonEdges.add(t2);
            return;
        }

        // else if there exists a common split: split the trees along the first split in commonEdges
        // and recursively call getDistance for the two new pairs of trees.
        PhyloTreeEdge commonEdge = commonEdges.get(0);

        // A will be the tree with leaves corresponding to 1's in commonEdge
        Vector<String> leaf2NumMapA = new Vector<String>();
        Vector<String> leaf2NumMapB = new Vector<String>();

        Vector<PhyloTreeEdge> edgesA1 = new Vector<PhyloTreeEdge>();
        Vector<PhyloTreeEdge> edgesA2 = new Vector<PhyloTreeEdge>();
        Vector<PhyloTreeEdge> edgesB1 = new Vector<PhyloTreeEdge>();
        Vector<PhyloTreeEdge> edgesB2 = new Vector<PhyloTreeEdge>();


        for (PhyloTreeEdge e : t1.getEdges()) {
            edgesA1.add(new PhyloTreeEdge(e.getAttribute().clone(), e.getOriginalEdge(), e.getOriginalID()));
            edgesB1.add(new PhyloTreeEdge(e.getAttribute().clone(), e.getOriginalEdge(), e.getOriginalID()));
        }

        for (PhyloTreeEdge e : t2.getEdges()) {
            edgesA2.add(new PhyloTreeEdge(e.getAttribute().clone(), e.getOriginalEdge(), e.getOriginalID()));
            edgesB2.add(new PhyloTreeEdge(e.getAttribute().clone(), e.getOriginalEdge(), e.getOriginalID()));

        }

        Boolean aLeavesAdded = false;  // if we have added a leaf in B representing the A tree
        int indexAleaves = 0;  // the index we are at in  the vectors holding the leaves in the A and B subtrees
        int indexBleaves = 0;
        // step through the leafs represented in commonEdge
        // (there should be two more leaves than edges)
        for (int i = 0; i < t1.getLeaf2NumMap().size(); i++) {
            if (commonEdge.contains(i)) {
                // commonEdge contains leaf i

                leaf2NumMapA.add((String) t1.getLeaf2NumMap().get(i));
                // these leaves must be added as a group to the B trees
                if (!aLeavesAdded) {
                    leaf2NumMapB.add(t1.getLeaf2NumMap().get(i) + "*");    // add a one of the leaves of the A tree to represent all the A trees leaves
//					 add the column corresponding to this leaf to the B edges vector (for the corresponding trees)
                    for (int j = 0; j < numEdges1; j++) {
                        if (t1.getEdge(j).properlyContains(commonEdge)) {
                            edgesB1.get(j).addOne(indexBleaves);
                        }
                    }
                    for (int j = 0; j < numEdges2; j++) {
                        if (t2.getEdge(j).properlyContains(commonEdge)) {
                            edgesB2.get(j).addOne(indexBleaves);
                        }
                    }
                    indexBleaves++;
                    aLeavesAdded = true;
                }
                // add the column corresponding to this leaf to the A edges vector (for the corresponding trees)
                // XXX: problem: might be adding edges which contain leaves in A but also
                for (int j = 0; j < numEdges1; j++) {
                    if (commonEdge.properlyContains(t1.getEdge(j)) && t1.getEdge(j).contains(i)) {
                        edgesA1.get(j).addOne(indexAleaves);
                    }
                }
                for (int j = 0; j < numEdges2; j++) {
                    if (commonEdge.properlyContains(t2.getEdge(j)) && t2.getEdge(j).contains(i)) {
                        edgesA2.get(j).addOne(indexAleaves);
                    }
                }
                indexAleaves++;
            } else {
                // commonEdge does not contain leaf i
                leaf2NumMapB.add((String) t1.getLeaf2NumMap().get(i));
//				 add the column corresponding to this leaf to the B edges vector (for the corresponding trees)
                for (int j = 0; j < numEdges1; j++) {
                    if (t1.getEdge(j).contains(i)) {
                        edgesB1.get(j).addOne(indexBleaves);
                    }
                }
                for (int j = 0; j < numEdges2; j++) {
                    if (t2.getEdges().get(j).contains(i)) {
                        edgesB2.get(j).addOne(indexBleaves);
                    }
                }
                indexBleaves++;
            }
        }
        edgesA1 = Tools.deleteEmptyEdges(edgesA1);
        edgesA2 = Tools.deleteEmptyEdges(edgesA2);
        edgesB1 = Tools.deleteEmptyEdges(edgesB1);
        edgesB2 = Tools.deleteEmptyEdges(edgesB2);

        // make the 4 trees
        PhyloTree tA1 = new PhyloTree(Tools.myVectorClonePhyloTreeEdge(edgesA1), Tools.myVectorCloneString(leaf2NumMapA));
        PhyloTree tB1 = new PhyloTree(Tools.myVectorClonePhyloTreeEdge(edgesB1), Tools.myVectorCloneString(leaf2NumMapB));
        PhyloTree tA2 = new PhyloTree(Tools.myVectorClonePhyloTreeEdge(edgesA2), Tools.myVectorCloneString(leaf2NumMapA));
        PhyloTree tB2 = new PhyloTree(Tools.myVectorClonePhyloTreeEdge(edgesB2), Tools.myVectorCloneString(leaf2NumMapB));
        splitOnCommonEdge(tA1, tA2);
        splitOnCommonEdge(tB1, tB2);
    }


    public static double getRobinsonFouldsDistance(PhyloTree tree1, PhyloTree tree2, boolean normalize) {
        Vector<PhyloTreeEdge> enic = tree1.getEdgesNotInCommonWith(tree2);
        enic.addAll(tree2.getEdgesNotInCommonWith(tree1));
        double rf_value = enic.size();
        if (normalize) rf_value /= tree1.numEdges() + tree2.numEdges();
        return rf_value;
    }

    public static double getWeightedRobinsonFouldsDistance(PhyloTree tree1, PhyloTree tree2, boolean normalize) {
        double wrf_value = 0;

        // Collect edges-in-common and edges-not-in-common...
        Vector<PhyloTreeEdge> eic = PhyloTree.getCommonEdges(tree1, tree2);
        Vector<PhyloTreeEdge> enic = tree1.getEdgesNotInCommonWith(tree2);
        enic.addAll(tree2.getEdgesNotInCommonWith(tree1));
        // ... and leaves
        EdgeAttribute[] leaves1 = tree1.getLeafEdgeAttribs();
        EdgeAttribute[] leaves2 = tree2.getLeafEdgeAttribs(); // Assuming these are the same

        // Collect length differences for internal edges...
        for (PhyloTreeEdge pte : eic) {
            wrf_value += pte.getLength();
        }

        for (PhyloTreeEdge pte : enic) {
            wrf_value += pte.getLength();
        }

        // ... and leaves
        for (int i = 0; i < leaves1.length; i++) {
            wrf_value += Math.abs(leaves1[i].norm() - leaves2[i].norm());
        }

        if (normalize) return wrf_value / (tree1.getBranchLengthSum() + tree2.getBranchLengthSum());
        return wrf_value;
    }

    public static double getEuclideanDistance(PhyloTree tree1, PhyloTree tree2, boolean normalize) {
        double euc_value = 0;

        // Collect edges-in-common and edges-not-in-common...
        Vector<PhyloTreeEdge> eic = PhyloTree.getCommonEdges(tree1, tree2);
        Vector<PhyloTreeEdge> enic = tree1.getEdgesNotInCommonWith(tree2);
        enic.addAll(tree2.getEdgesNotInCommonWith(tree1));
        // ... and leaves
        EdgeAttribute[] leaves1 = tree1.getLeafEdgeAttribs();
        EdgeAttribute[] leaves2 = tree2.getLeafEdgeAttribs(); // Assuming these are the same

        // Collect length differences for internal edges...
        for (PhyloTreeEdge pte : eic) {
            euc_value += Math.pow(pte.getLength(), 2);
        }

        for (PhyloTreeEdge pte : enic) {
            euc_value += Math.pow(pte.getLength(), 2);
        }

        // ... and leaves
        for (int i = 0; i < leaves1.length; i++) {
            euc_value += Math.pow(leaves1[i].norm() - leaves2[i].norm(), 2);
        }

        if (normalize) return Math.sqrt(euc_value) / (tree1.getDistanceFromOrigin() + tree2.getDistanceFromOrigin());
        return Math.sqrt(euc_value);
    }

    public static double getGeodesicDistance(PhyloTree tree1, PhyloTree tree2, boolean normalize) throws IOException {
        double distance = getGeodesic(tree1, tree2, null).getDist();
        if (normalize) return distance / (tree1.getDistanceFromOrigin() + tree2.getDistanceFromOrigin());
        return distance;
    }

    /**
     * Returns the distance between t1 and t2, accounting for any common edges and leaf edges.
     * Calls recursive getGeodesic
     * Does not assume t1 and t2 have the same number of edges.
     * Pass in null for geoFile to not write to a file.
     * XXX: how to deal with multifurcating trees
     */
    public static Geodesic getGeodesic(PhyloTree t1, PhyloTree t2, String geoFile) throws IOException {
        double leafContributionSquared = 0;
        EdgeAttribute[] t1LeafEdgeAttribs = t1.getLeafEdgeAttribs();
        EdgeAttribute[] t2LeafEdgeAttribs = t2.getLeafEdgeAttribs();
        Geodesic geo = new Geodesic(new RatioSequence());

        String verboseOutput = "";

        // get the leaf contributions
        for (int i = 0; i < t1.getLeaf2NumMap().size(); i++) {
            if (!(t1.getLeaf2NumMap().get(i).equals(t2.getLeaf2NumMap().get(i)))) {
                System.err.println("Error getting geodesic: trees do not have the same sets of leaves");
                System.err.println("Starting tree leaves: " + t1.getLeaf2NumMap());
                System.err.println("Target tree leaves: " + t2.getLeaf2NumMap());

                System.err.println("Starting tree: " + t1.getNewick(true));
                System.err.println("Target tree: " + t2.getNewick(true));

                throw new UnsupportedOperationException("Error getting geodesic: trees do not have the same sets of leaves");
            }
            leafContributionSquared = leafContributionSquared + Math.pow(EdgeAttribute.difference(t1LeafEdgeAttribs[i], t2LeafEdgeAttribs[i]).norm(), 2);
        }
        geo.setLeafContributionSquared(leafContributionSquared);

        aTreesNoCommonEdges = new Vector<PhyloTree>();
        bTreesNoCommonEdges = new Vector<PhyloTree>();

        // get the pairs of trees with no common edges put into aTreesNoCommonEdges and bTreesNoCommonEdges
        //  aTreesNoCommonEdges.get(i) goes with bTreesNoCommonEdges.get(i)
        splitOnCommonEdge(t1, t2);

        //set the common edges
        geo.setCommonEdges(PhyloTree.getCommonEdges(t1, t2));

        // find the geodesic between each pair of subtrees found by removing the common edges
        for (int i = 0; i < aTreesNoCommonEdges.size(); i++) {
            PhyloTree subTreeA = aTreesNoCommonEdges.get(i);
            PhyloTree subTreeB = bTreesNoCommonEdges.get(i);

            Geodesic newGeo = getGeodesicNoCommonEdges(subTreeA, subTreeB);

            geo.setRS(RatioSequence.interleave(geo.getRS(), newGeo.getRS()));
        }

        return geo;
    }

    /**
     * Returns the geodesic between t1 and t2, which are assumed to have no common edges.
     * Does not assume t1 and t2 have the same number of edges.
     * Does not take into account the leaf edges.
     * Uses polynomial algorithm.
     * XXX: how to deal with multifurcating trees
     * <p/>
     * Returns:  a Geodesic with just the ratio sequence set
     */
    public static Geodesic getGeodesicNoCommonEdges(PhyloTree t1, PhyloTree t2) {
        int numEdges1 = t1.getEdges().size(); // number of edges in tree 1
        int numEdges2 = t2.getEdges().size(); // number of edges in tree 2
        RatioSequence rs = new RatioSequence();
        int[] aVertices, bVertices;
        Vector<Ratio> queue = new Vector<Ratio>();
        Ratio ratio;
        int[][] cover;

        if (numEdges1 == 0 && numEdges2 == 0) {
            return new Geodesic(new RatioSequence());
        }

        // double-check no common edges
        Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1, t2);
        if (commonEdges.size() != 0) {
            System.out.println("Exiting: tried to compute geodesic between subtrees that should not have common edges, but do!  t1 = " + t1 + " and t2 = " + t2);
            System.exit(1);
        }

        // double-check that both trees have splits.  Otherwise didn't remove a common edge.
        if (numEdges1 == 0 || numEdges2 == 0) {
            System.out.println("Exiting: tried to compute geodesic between subtrees that should not have common/compatible edges, but do!  t1 = " + t1 + " and t2 = " + t2);
            System.exit(1);
        }

        // if we can't split the ratio because it has too few edges in either the numerator or denominator
        if ((numEdges1 == 1) || (numEdges2 == 1)) {
            rs.add(new Ratio(t1.getEdges(), t2.getEdges()));
            return new Geodesic(rs);
        }

        // initialize BipartiteGraph
        boolean[][] incidenceMatrix = Tools.getIncidenceMatrix(t1.getEdges(), t2.getEdges());


        BipartiteGraph bg = new BipartiteGraph(incidenceMatrix, t1.getIntEdgeAttribNorms(), t2.getIntEdgeAttribNorms());

        queue.add(new Ratio(t1.getEdges(), t2.getEdges()));

        while (queue.size() > 0) {
            ratio = queue.remove(0);

            aVertices = new int[ratio.getEEdges().size()];
            bVertices = new int[ratio.getFEdges().size()];

            // convert the ratio to what we pass to vertex cover
            for (int i = 0; i < ratio.getEEdges().size(); i++) {
                aVertices[i] = t1.getEdges().indexOf(ratio.getEEdges().get(i));
            }
            for (int i = 0; i < ratio.getFEdges().size(); i++) {
                bVertices[i] = t2.getEdges().indexOf(ratio.getFEdges().get(i));
            }

            // get the cover
            cover = bg.vertex_cover(aVertices, bVertices);

            // check if cover is trivial
            if ((cover[0][0] == 0) || (cover[0][0] == aVertices.length)) {
                // add ratio to geodesic
                rs.add(ratio);


            } else {  // cover not trivial
                // make two new ratios
                Ratio r1 = new Ratio();
                Ratio r2 = new Ratio();

                int j = 0;  // for index in cover array

                // split the ratio based on the cover
                for (int i = 0; i < aVertices.length; i++) {
                    if ((j < cover[2].length) && (aVertices[i] == cover[2][j])) {
                        r1.addEEdge(t1.getEdge(aVertices[i]));
                        j++;
                    } else { // the split is not in the cover, and hence dropped first
                        r2.addEEdge(t1.getEdge(aVertices[i]));
                    }
                }

                j = 0;   // reset index
                // split the ratio based on the cover
                for (int i = 0; i < bVertices.length; i++) {
                    if ((j < cover[3].length) && (bVertices[i] == cover[3][j])) {
                        r2.addFEdge(t2.getEdge(bVertices[i]));
                        j++;
                    } else { // the split is not in the cover, and hence dropped first
                        r1.addFEdge(t2.getEdge(bVertices[i]));
                    }
                }

                // add ratios to the queue
                queue.add(0, r2);
                queue.add(0, r1);
            }
        }

        return new Geodesic(rs);
    }

    /**
     * @param args
     */
    public static void main(String[] args) throws IOException {
        PhyloTree tree1;
        tree1 = new PhyloTree("(a:1, (b:2, c:4):8);", false);
        PhyloTree tree2 = new PhyloTree("(b:1, (c:2, b:4):8);", false);
        System.out.println(getGeodesicDistance(tree1, tree2, false));
        System.out.println(getRobinsonFouldsDistance(tree1, tree2, false));
        System.out.println(getWeightedRobinsonFouldsDistance(tree1, tree2, false));
        System.out.println(getEuclideanDistance(tree1, tree2, false));
        System.exit(0);
    }

}
