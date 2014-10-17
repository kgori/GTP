/** This file is part of GeodeMAPS, a program for computing the geodesic distance between phylogenetic trees.
 Copyright (C) 2008  Megan Owen

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

import java.io.*;
import java.util.*;


/**
 * This class computes all inter-tree geodesic distances for the trees given in the file.
 *
 * @author maowen
 */

public class TreeDistance {

    public static Vector<RatioSequence> finalRatioSeqs = new Vector<RatioSequence>();
    public static double minTreeDist = -1;
    public static RatioSequence minTreeDistRatioSeq = null;
    public static long numMaxPaths = 0;
    public static long numPrunes = 0;
    public static long numNodes = 0;

    public static int pathToSearch = 0;
    // stores pairs of trees with no common edges.  Should be reset at each new distance calculation
    public static Vector<PhyloTree> aTreesNoCommonEdges = new Vector<PhyloTree>();
    public static Vector<PhyloTree> bTreesNoCommonEdges = new Vector<PhyloTree>();
    public static double firstTreeDist = -1;
    public static boolean rooted = true;  //holds if the trees are rooted or not.
    public static boolean normalize = false;  // holds if we should normalize the tree split lengths
    public static int verbose = 0;
    public static Hashtable<String, Double> nodeHashtable = new Hashtable<String, Double>(5000);  // creates a new hash table with 5000 buckets
    public static Hashtable<String, Geodesic> subTreeHashtable = new Hashtable<String, Geodesic>(5000);  // creates a new hash table with 5000 buckets
    // key is the node, as represented by the Vector<Bipartition> m, and value is the shortest distance so far
    // to this node
    public static Hashtable<String, RatioSequence> subTreeRSHashtable = new Hashtable<String, RatioSequence>(5000);  // creates a new hash table with 5000 buckets
    // key is the subtree this problem represents, as represented by the string of leaf2NumMap, and value is the shortest distance between
    // the two subtrees
    public static String LEAF_CONTRIBUTION_SQUARED_DESCRIPTION = "(Leaf contribution squared = square of the length of the vector" +
            " whose i-th element is the absolute value of the difference between the length of the split ending in leaf i in the first tree" +
            " and the length of the split ending in leaf i in the second tree.)";
    // key is the subtree this problem represents, as represented by the string of leaf2NumMap, and value is the shortest distance between
    // the two subtrees
    private static long numBaseCase = 0;

    /**
     * To be called in between calculating min distances
     * Clears the static variables in TreeDistance
     */
    public static void resetTreeDistanceState() {
        finalRatioSeqs = new Vector<RatioSequence>();
        minTreeDist = -1;
        minTreeDistRatioSeq = null;
        nodeHashtable = new Hashtable<String, Double>(5000);
        numPrunes = 0;
        numNodes = 0;
        numMaxPaths = 0;
        firstTreeDist = -1;
    }

    /**
     * Changes all 1's in column col in the 01-vector representation of the Bipartitions in m to 0's,
     * and returns the resulting version of m.
     *
     * @param col
     * @return
     */
    public static Vector<Bipartition> zeroCol(int col, Vector<Bipartition> m) {
        Iterator<Bipartition> fEdgesIter = m.iterator();
        while (fEdgesIter.hasNext()) {
            Bipartition fEdge = fEdgesIter.next();
            if ((fEdge != null) && (fEdge.contains(col))) {
                fEdge.removeOne(col);
            }
        }
        return m;
    }

    /**
     * Checks for duplicate split in a vector of edges  (doesn't worry about split lengths)
     */
    public static boolean checkForDuplicateEdges(Vector<PhyloTreeEdge> edges) {
        boolean duplicate = false;
        for (int i = 0; i < edges.size() - 1; i++) {
            for (int j = i + 1; j < edges.size(); j++) {
                if (edges.get(i).sameBipartition(edges.get(j))) {
//					System.out.println("Duplicate edges i = " + i + " and j = " + j+ "; split " + edges.get(i) );
                    duplicate = true;
                }
            }
        }
        return duplicate;
    }


    /**
     * Returns the distance between t1 and t2, accounting for any common edges and leaf edges.
     * Calls recursive getGeodesic
     * Does not assume t1 and t2 have the same number of edges.
     * reset is true if this is the first call to
     * resets hash table for divide and conquer
     * XXX: how to deal with multifurcating trees
     */
    public static Geodesic getGeodesic2(PhyloTree t1, PhyloTree t2, String algorithm, String geoFile) throws IOException {
        double leafContributionSquared = 0;
        EdgeAttribute[] t1LeafEdgeAttribs = t1.getLeafEdgeAttribs();
        EdgeAttribute[] t2LeafEdgeAttribs = t2.getLeafEdgeAttribs();
        Geodesic geo = new Geodesic(new RatioSequence());

        String verboseOutput = "";

        // get the leaf contributions
        for (int i = 0; i < t1.getLeaf2NumMap().size(); i++) {
            if (!(t1.getLeaf2NumMap().get(i).equals(t2.getLeaf2NumMap().get(i)))) {
                System.out.println("Exiting: Leaves don't match for trees " + t1 + " and " + t2);
                throw new RuntimeException();
            }

            leafContributionSquared = leafContributionSquared + Math.pow(EdgeAttribute.difference(t1LeafEdgeAttribs[i], t2LeafEdgeAttribs[i]).norm(), 2);
// pre EdgeAttribute			leafContributionSquared = leafContributionSquared+ Math.pow(Math.abs(t1LeafEdgeAttribs[i] - t2LeafEdgeAttribs[i]), 2);
        }
        geo.setLeafContributionSquared(leafContributionSquared);

        if (TreeDistance.verbose > 0) {
            // print out leaf2NumMap
/*			System.out.println("leaf2NumMap for starting and target trees is ");
            for (int i = 0; i <t1.getLeaf2NumMap().size(); i++) {
				System.out.println("" + t1.getLeaf2NumMap().get(i) + "\t" + i);
			} */

            System.out.println("Starting tree: " + t1.getNewick(true));
            verboseOutput = verboseOutput + "Starting tree: " + t1.getNewick(true) + "\n";
            System.out.println("Target tree: " + t2.getNewick(true));
            verboseOutput = verboseOutput + "Target tree: " + t2.getNewick(true) + "\n";

            System.out.println("\nStarting tree edges:");
            verboseOutput = verboseOutput + "\nStarting tree edges:\n";
            verboseOutput = verboseOutput + PhyloTreeEdge.printEdgesVerbose(t1.getEdges(), t1.getLeaf2NumMap(), true);

            System.out.println("\nTarget tree edges:");
            verboseOutput = verboseOutput + "\nTarget tree edges:\n";
            verboseOutput = verboseOutput + PhyloTreeEdge.printEdgesVerbose(t2.getEdges(), t2.getLeaf2NumMap(), true);

            // leaf contributions
            System.out.println("\nLeaf contribution squared " + TreeDistance.truncate(leafContributionSquared, 6));
            verboseOutput = verboseOutput + "\nLeaf contribution squared " + TreeDistance.truncate(leafContributionSquared, 6) + "\n";
            System.out.println(LEAF_CONTRIBUTION_SQUARED_DESCRIPTION);
            verboseOutput = verboseOutput + LEAF_CONTRIBUTION_SQUARED_DESCRIPTION + "\n";
        }

        // reset hash table for divide and conquer
        subTreeHashtable = new Hashtable<String, Geodesic>(5000);
        aTreesNoCommonEdges = new Vector<PhyloTree>();
        bTreesNoCommonEdges = new Vector<PhyloTree>();

        // get the pairs of trees with no common edges put into aTreesNoCommonEdges and bTreesNoCommonEdges
        //  aTreesNoCommonEdges.get(i) goes with bTreesNoCommonEdges.get(i)
        splitOnCommonEdge(t1, t2);

        //set the common edges
        geo.setCommonEdges(PhyloTree.getCommonEdges(t1, t2));

        if (verbose > 0) {
            System.out.println("\nCommon edges are:  (Length = abs. value of difference in length between the two trees)");
            verboseOutput = verboseOutput + "\nCommon edges are:  (Length = abs. value of difference in length between the two trees)\n";

            Vector<PhyloTreeEdge> commonEdges = geo.getCommonEdges();

            verboseOutput = verboseOutput + PhyloTreeEdge.printEdgesVerbose(commonEdges, t1.getLeaf2NumMap(), true);
            double commonEdgeContributionSquared = 0;
            for (int i = 0; i < commonEdges.size(); i++) {
                commonEdgeContributionSquared = commonEdgeContributionSquared + Math.pow(commonEdges.get(i).getLength(), 2);
            }
            System.out.println("\nCommon edges contribution squared: " + TreeDistance.truncate(commonEdgeContributionSquared, 6));
            verboseOutput = verboseOutput + "\nCommon edges contribution squared: " + TreeDistance.truncate(commonEdgeContributionSquared, 6) + "\n";
            System.out.println("(sum of squares of above differences in length)");
            verboseOutput = verboseOutput + "(sum of squares of above differences in length)\n";
            System.out.println("=============================================================================================================================");
            verboseOutput = verboseOutput + "=============================================================================================================================\n";

            System.out.println("\nNow finding the geodesic between the following subtrees, which have no edges in common:");
            verboseOutput = verboseOutput + "\nNow finding the geodesic between the following subtrees, which have no edges in common:\n";
        }

        // find the geodesic between each pair of subtrees found by removing the common edges
        for (int i = 0; i < aTreesNoCommonEdges.size(); i++) {
            PhyloTree subTreeA = aTreesNoCommonEdges.get(i);
            PhyloTree subTreeB = bTreesNoCommonEdges.get(i);

            if (verbose > 0) {
                System.out.println("Leaves or subtree representatives in subtrees:");
                verboseOutput = verboseOutput + "Leaves or subtree representatives in subtrees:\n";
                for (int j = 0; j < subTreeA.getLeaf2NumMap().size(); j++) {
                    System.out.println("" + subTreeA.getLeaf2NumMap().get(j));
                    verboseOutput = verboseOutput + subTreeA.getLeaf2NumMap().get(j) + "\n";
                }

                System.out.println("\nStarting subtree edges:");
                verboseOutput = verboseOutput + "\nStarting subtree edges:\n";
                verboseOutput = verboseOutput + PhyloTreeEdge.printEdgesVerbose(subTreeA.getEdges(), subTreeA.getLeaf2NumMap(), false);

                System.out.println("\nTarget subtree edges:");
                verboseOutput = verboseOutput + "\nTarget subtree edges:\n";
                verboseOutput = verboseOutput + PhyloTreeEdge.printEdgesVerbose(subTreeB.getEdges(), subTreeB.getLeaf2NumMap(), false);
            }

            Geodesic newGeo = getGeodesicRecursive(subTreeA, subTreeB, algorithm);

            if (verbose > 0) {
                System.out.println("\nGeodesic distance between above subtrees, ignoring edges ending in leaves: " + TreeDistance.truncate(newGeo.getRS().getNonDesRSWithMinDist().getDistance(), 6));
                verboseOutput = verboseOutput + "\nGeodesic distance between above subtrees, ignoring edges ending in leaves: " + TreeDistance.truncate(newGeo.getRS().getNonDesRSWithMinDist().getDistance(), 6) + "\n";

                System.out.println("Ratio sequence corresponding to the geodesic:\nCombinatorial type: " + newGeo.getRS().getNonDesRSWithMinDist().toStringCombType());
                verboseOutput = verboseOutput + "Ratio sequence corresponding to the geodesic:\nCombinatorial type: " + newGeo.getRS().getNonDesRSWithMinDist().toStringCombType() + "\n";

                System.out.println(newGeo.getRS().getNonDesRSWithMinDist().toStringVerbose(subTreeA.getLeaf2NumMap()));
                verboseOutput = verboseOutput + newGeo.getRS().getNonDesRSWithMinDist().toStringVerbose(subTreeA.getLeaf2NumMap()) + "\n";

                System.out.println("------------------------------------------------------------------------------------------------------------");
                verboseOutput = verboseOutput + "------------------------------------------------------------------------------------------------------------\n";
            }

            geo.setRS(RatioSequence.interleave(geo.getRS(), newGeo.getRS()));
        }

        if (verbose > 0) {
            System.out.println("\nGeodesic distance between start and target tree is " + TreeDistance.truncate(geo.getDist(), 6));
            verboseOutput = verboseOutput + "\nGeodesic distance between start and target tree is " + TreeDistance.truncate(geo.getDist(), 6) + "\n";
        }

        // write verbose output to geofile, if in verbose mode
        if (verbose > 0) {
            PrintWriter outputStream = null;

            try {
                outputStream = new PrintWriter(new FileWriter(geoFile));

                outputStream.println(verboseOutput);

                if (outputStream != null) {
                    outputStream.close();
                }
            } catch (IOException e) {
                System.out.println("Error opening or writing to " + geoFile + ": " + e.getMessage());
                throw e;
            }
        }
        return geo;
    }

    /**
     * Stores subtrees with no common edges in the global variables aTreesNoCommonEdges (from t1)
     * and in bTreesNoCommonEdges (from t2).  Also returns a vector containing pairs of tree with no common edges.
     */
    public static Vector<PhyloTree> splitOnCommonEdge(PhyloTree t1, PhyloTree t2) {
        Vector<PhyloTree> disjointTreePairs = new Vector<PhyloTree>();
        int numEdges1 = t1.getEdges().size(); // number of edges in tree 1
        int numEdges2 = t2.getEdges().size(); /// number of edges in tree 2

        if (numEdges1 == 0 && numEdges2 == 0) {
            return disjointTreePairs;
        }
        // look for common edges
        Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1, t2);

        // if there are no common edges
        // XXX: need to check the following methods don't require the trees to have the same number of edges
        if (commonEdges.size() == 0) {
//			System.out.println("In splitOnCommonEdges, no common edges in " + t1 + " and " + t2);
            aTreesNoCommonEdges.add(t1);
            bTreesNoCommonEdges.add(t2);
            disjointTreePairs.add(t1);
            disjointTreePairs.add(t2);
            return disjointTreePairs;
        }
//		System.out.println("At least one common split; edges are " + commonEdges);

        // else if there exists a common split: split the trees along the first split in commonEdges
        // and recursively call getDistance for the two new pairs of trees.
        PhyloTreeEdge commonEdge = commonEdges.get(0);
//		System.out.println("Common edges is " + commonEdge);

        // A will be the tree with leaves corresponding to 1's in commonEdge
        Vector<String> leaf2NumMapA = new Vector<String>();
        Vector<String> leaf2NumMapB = new Vector<String>();

        Vector<PhyloTreeEdge> edgesA1 = new Vector<PhyloTreeEdge>();
        Vector<PhyloTreeEdge> edgesA2 = new Vector<PhyloTreeEdge>();
        Vector<PhyloTreeEdge> edgesB1 = new Vector<PhyloTreeEdge>();
        Vector<PhyloTreeEdge> edgesB2 = new Vector<PhyloTreeEdge>();

//		System.out.println("tree 1 is " + t1);

        for (int i = 0; i < numEdges1; i++) {
//			System.out.println("t1.getEdge(i).getOriginalID() is " + t1.getEdge(i).getOriginalID()  +" \n");
            edgesA1.add(new PhyloTreeEdge(t1.getEdge(i).getAttribute(), t1.getEdge(i).getOriginalEdge(), t1.getEdge(i).getOriginalID()));
            edgesB1.add(new PhyloTreeEdge(t1.getEdge(i).getAttribute(), t1.getEdge(i).getOriginalEdge(), t1.getEdge(i).getOriginalID()));
        }

        for (int i = 0; i < numEdges2; i++) {
            edgesA2.add(new PhyloTreeEdge(t2.getEdge(i).getAttribute(), t2.getEdge(i).getOriginalEdge(), t2.getEdge(i).getOriginalID()));
            edgesB2.add(new PhyloTreeEdge(t2.getEdge(i).getAttribute(), t2.getEdge(i).getOriginalEdge(), t2.getEdge(i).getOriginalID()));
        }

        Boolean aLeavesAdded = false;  // if we have added a leaf representing th
        int indexAleaves = 0;  // the index we are at in  the vectors holding the leaves in the A and B subtrees
        int indexBleaves = 0;
//		System.out.println("Edges are " + edgesB1);
        // step through the leafs represented in commonEdge
        // (there should be two more leaves than edges)
        for (int i = 0; i < t1.getLeaf2NumMap().size(); i++) {
//			System.out.println("leaf i = " +i);
            if (commonEdge.contains(i)) {
                // commonEdge contains leaf i

                leaf2NumMapA.add((String) t1.getLeaf2NumMap().get(i));
                // these leaves must be added as a group to the B trees
                if (!aLeavesAdded) {
                    leaf2NumMapB.add((String) t1.getLeaf2NumMap().get(i) + "*");    // add a one of the leaves of the A tree to represent all the A trees leaves
//					 add the column corresponding to this leaf to the B edges vector (for the corresponding trees)
                    for (int j = 0; j < numEdges1; j++) {
                        if (t1.getEdge(j).properlyContains(commonEdge)) {
                            edgesB1.get(j).addOne(indexBleaves);
//							System.out.println("Adding a one in position " + indexBleaves + " in !aLeavesAdded to split j = " +j+" of T1");
                        }
                    }
                    for (int j = 0; j < numEdges2; j++) {
                        if (t2.getEdge(j).properlyContains(commonEdge)) {
                            edgesB2.get(j).addOne(indexBleaves);
//							System.out.println("Adding a one in position " + indexBleaves + " in !aLeavesAdded to split j = " + j+ " of Tb");
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
                        //					System.out.println("t1 contains i = " + i + " in split j = " +j+" so added at 1 at position " + indexBleaves );
                    }
                }
                for (int j = 0; j < numEdges2; j++) {
                    if (t2.getEdges().get(j).contains(i)) {
                        edgesB2.get(j).addOne(indexBleaves);
//						System.out.println("t2 contains i = " + i + " in split j = " +j+ " so added at 1 at position " + indexBleaves );
                    }
                }
                indexBleaves++;
            }
        }
//		System.out.println("Edges are " + edgesA2);
//		System.out.println("B2 Edges are " + edgesB2);
        edgesA1 = TreeDistance.deleteZeroEdges(edgesA1);
        edgesA2 = TreeDistance.deleteZeroEdges(edgesA2);
        edgesB1 = TreeDistance.deleteZeroEdges(edgesB1);
        edgesB2 = TreeDistance.deleteZeroEdges(edgesB2);

//		System.out.println("B2 Edges are " + edgesB2);

        // make the 4 trees
        PhyloTree tA1 = new PhyloTree(TreeDistance.myVectorClonePhyloTreeEdge(edgesA1), TreeDistance.myVectorCloneString(leaf2NumMapA));
        PhyloTree tB1 = new PhyloTree(TreeDistance.myVectorClonePhyloTreeEdge(edgesB1), TreeDistance.myVectorCloneString(leaf2NumMapB));
        PhyloTree tA2 = new PhyloTree(TreeDistance.myVectorClonePhyloTreeEdge(edgesA2), TreeDistance.myVectorCloneString(leaf2NumMapA));
        PhyloTree tB2 = new PhyloTree(TreeDistance.myVectorClonePhyloTreeEdge(edgesB2), TreeDistance.myVectorCloneString(leaf2NumMapB));
//		System.out.println("About to find distance between A1 = " + tA1 + " and A2 = " + tA2 + "...");
//		System.out.println("Corresponding B's are B1 = " + tB1 + " and B2 = " + tB2);
        disjointTreePairs.addAll(splitOnCommonEdge(tA1, tA2));


//		System.out.println("geodesic " + geoA + " between A trees: A1 = " + tA1 + " and A2 = " + tA2);
//		System.out.println("About to find distance between B1 = " + tB1 + " and B2 = " + tB2 + "...");
        disjointTreePairs.addAll(splitOnCommonEdge(tB1, tB2));

        return disjointTreePairs;
    }


    /**
     * Returns the distance between t1 and t2, accounting for any common edges, but not leaf edges.
     * Called recursively.
     * Does not assume t1 and t2 have the same number of edges.
     * reset is true if this is the first call to
     * XXX: how to deal with multifurcating trees
     */
    public static Geodesic getGeodesicRecursive(PhyloTree t1, PhyloTree t2, String algorithm) {
        resetTreeDistanceState();  // reset the static variables in tree distance

        int numEdges1 = t1.getEdges().size(); // number of edges in tree 1
        int numEdges2 = t2.getEdges().size(); /// number of edges in tree 2

/*		if (numBaseCase%10000 == 0) {
            System.out.println("Called getGeodesic for trees " + t1 + " and " + t2);
		}*/
        if (numEdges1 == 0 && numEdges2 == 0) {
            return new Geodesic(new RatioSequence());
        }

        // look for common edges
        // XXX:  I should be able to clean up the common split part to match getGeodesic2
        Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1, t2);

        if (verbose > 1) {
            System.out.println("Common edges are " + commonEdges);
        }

        // if there are no common edges
        // XXX: need to check the following methods don't require the trees to have the same number of edges
        if (commonEdges.size() == 0) {
//			System.out.println("No common edges in " + t1 + " and " + t2);
            if (algorithm.equals("divide")) {
                // check for subproblem having been done
                Geodesic subTreeGeodesic = (Geodesic) subTreeHashtable.get(t1.getLeaf2NumMap().toString());
                if (subTreeGeodesic == null) {
                    subTreeGeodesic = TreeDistance.getDivideAndConquerGeodesicNoCommonEdges(t1, t2);
                    subTreeHashtable.put(t1.getLeaf2NumMap().toString(), subTreeGeodesic.clone());
                }
                return subTreeGeodesic;
            } else if (algorithm.equals("DivideAndConquerRS")) {
                // check for subproblem having been done
                RatioSequence subTreeRS = (RatioSequence) subTreeRSHashtable.get(t1.getLeaf2NumMap().toString());
                if (subTreeRS == null) {
                    subTreeRS = TreeDistance.getDivideAndConquerRSNoCommonEdges(t1, t2);
                    subTreeRSHashtable.put(t1.getLeaf2NumMap().toString(), subTreeRS.clone());
                }
//				System.out.println("returning " + subTreeRS);
                return new Geodesic(subTreeRS);
            } else if (algorithm.equals("dynamic")) {
                return getPruned2GeodesicNoCommonEdges(t1, t2);
            } else {
                System.out.println("" + algorithm + " is an invalid algorithm");
                throw new RuntimeException();
            }
        }
//		System.out.println("At least one common split; edges are " + commonEdges);

        // else if there exists a common split: split the trees along the first split in commonEdges
        // and recursively call getDistance for the two new pairs of trees.
        PhyloTreeEdge commonEdge = commonEdges.get(0);
//		System.out.println("Common edges is " + commonEdge);

        // A will be the tree with leaves corresponding to 1's in commonEdge
        Vector<String> leaf2NumMapA = new Vector<String>();
        Vector<String> leaf2NumMapB = new Vector<String>();

        Vector<PhyloTreeEdge> edgesA1 = new Vector<PhyloTreeEdge>();
        Vector<PhyloTreeEdge> edgesA2 = new Vector<PhyloTreeEdge>();
        Vector<PhyloTreeEdge> edgesB1 = new Vector<PhyloTreeEdge>();
        Vector<PhyloTreeEdge> edgesB2 = new Vector<PhyloTreeEdge>();

//		System.out.println("tree 1 is " + t1);

        for (int i = 0; i < numEdges1; i++) {
            edgesA1.add(new PhyloTreeEdge(t1.getEdge(i).getAttribute(), t1.getEdge(i).getOriginalEdge(), t1.getEdge(i).getOriginalID()));
            edgesB1.add(new PhyloTreeEdge(t1.getEdge(i).getAttribute(), t1.getEdge(i).getOriginalEdge(), t1.getEdge(i).getOriginalID()));
        }

        for (int i = 0; i < numEdges2; i++) {
            edgesA2.add(new PhyloTreeEdge(t2.getEdge(i).getAttribute(), t2.getEdge(i).getOriginalEdge(), t2.getEdge(i).getOriginalID()));
            edgesB2.add(new PhyloTreeEdge(t2.getEdge(i).getAttribute(), t2.getEdge(i).getOriginalEdge(), t2.getEdge(i).getOriginalID()));
        }

        Boolean aLeavesAdded = false;  // if we have added a leaf representing th
        int indexAleaves = 0;  // the index we are at in  the vectors holding the leaves in the A and B subtrees
        int indexBleaves = 0;
//		System.out.println("Edges are " + edgesB1);
        // step through the leafs represented in commonEdge
        // (there should be two more leaves than edges)
        for (int i = 0; i < t1.getLeaf2NumMap().size(); i++) {
//			System.out.println("leaf i = " +i);
            if (commonEdge.contains(i)) {
                // commonEdge contains leaf i

                leaf2NumMapA.add((String) t1.getLeaf2NumMap().get(i));
                // these leaves must be added as a group to the B trees
                if (!aLeavesAdded) {
                    leaf2NumMapB.add((String) t1.getLeaf2NumMap().get(i));    // add a one of the leaves of the A tree to represent all the A trees leaves
//					 add the column corresponding to this leaf to the B edges vector (for the corresponding trees)
                    for (int j = 0; j < numEdges1; j++) {
                        if (t1.getEdge(j).properlyContains(commonEdge)) {
                            edgesB1.get(j).addOne(indexBleaves);
//							System.out.println("Adding a one in position " + indexBleaves + " in !aLeavesAdded to split j = " +j+" of T1");
                        }
                    }
                    for (int j = 0; j < numEdges2; j++) {
                        if (t2.getEdge(j).properlyContains(commonEdge)) {
                            edgesB2.get(j).addOne(indexBleaves);
//							System.out.println("Adding a one in position " + indexBleaves + " in !aLeavesAdded to split j = " + j+ " of Tb");
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
                        //					System.out.println("t1 contains i = " + i + " in split j = " +j+" so added at 1 at position " + indexBleaves );
                    }
                }
                for (int j = 0; j < numEdges2; j++) {
                    if (t2.getEdges().get(j).contains(i)) {
                        edgesB2.get(j).addOne(indexBleaves);
//						System.out.println("t2 contains i = " + i + " in split j = " +j+ " so added at 1 at position " + indexBleaves );
                    }
                }
                indexBleaves++;
            }
        }
//		System.out.println("Edges are " + edgesB1);
        edgesA1 = TreeDistance.deleteZeroEdges(edgesA1);
        edgesA2 = TreeDistance.deleteZeroEdges(edgesA2);
        edgesB1 = TreeDistance.deleteZeroEdges(edgesB1);
        edgesB2 = TreeDistance.deleteZeroEdges(edgesB2);

//		System.out.println("Edges are " + edgesB1);

        // make the 4 trees
        PhyloTree tA1 = new PhyloTree(TreeDistance.myVectorClonePhyloTreeEdge(edgesA1), TreeDistance.myVectorCloneString(leaf2NumMapA));
        PhyloTree tB1 = new PhyloTree(TreeDistance.myVectorClonePhyloTreeEdge(edgesB1), TreeDistance.myVectorCloneString(leaf2NumMapB));
        PhyloTree tA2 = new PhyloTree(TreeDistance.myVectorClonePhyloTreeEdge(edgesA2), TreeDistance.myVectorCloneString(leaf2NumMapA));
        PhyloTree tB2 = new PhyloTree(TreeDistance.myVectorClonePhyloTreeEdge(edgesB2), TreeDistance.myVectorCloneString(leaf2NumMapB));
//		System.out.println("About to find distance between A1 = " + tA1 + " and A2 = " + tA2 + "...");
//		System.out.println("Corresponding B's are B1 = " + tB1 + " and B2 = " + tB2);
        Geodesic geoA = getGeodesicRecursive(tA1, tA2, algorithm);
//		System.out.println("geoA is " + geoA);
        //System.out.println("geoA is " + geoA.toStringVerbose(tA1, tA2));
//		System.out.println("geodesic " + geoA + " between A trees: A1 = " + tA1 + " and A2 = " + tA2);
//		System.out.println("About to find distance between B1 = " + tB1 + " and B2 = " + tB2 + "...");
        Geodesic geoB = getGeodesicRecursive(tB1, tB2, algorithm);
//		System.out.println("geoB is " + geoB);
//		System.out.println("geoB is " + geoB.toStringVerbose(tB1, tB2));
//		System.out.println("geodesic " + geoA + " between A trees: A1 = " + tA1 + " and A2 = " + tA2);
//		System.out.println("geodesic " + geoB + " between B trees, B1 = " + tB1 + " and B2 = " + tB2);
//XXX:  pre-EdgeAttribute		double distCommonEdge = Math.abs(t1.getLengthOfSplit(commonEdge) - t2.getLengthOfSplit(commonEdge));
//		System.out.println("Change in distance of common split " + commonEdge +  " is " + distCommonEdge);
//		System.out.println("Returning square root of " + distA + "^2 + " + distB + "^2 + " + distCommonEdge + "^2 = "  + Math.sqrt(Math.pow(distA, 2) + Math.pow(distB, 2) + Math.pow(distCommonEdge, 2)) );
//		System.out.println("Interleaved rs is " + RatioSequence.interleave(geoA.getRS(), geoB.getRS()));
//		return new Geodesic(Math.sqrt(Math.pow(geoA.getDist(), 2) + Math.pow(geoB.getDist(), 2) + Math.pow(distCommonEdge, 2)), RatioSequence.interleave(geoA.getRS(), geoB.getRS()) );
        Vector<PhyloTreeEdge> newCommonEdges = new Vector<PhyloTreeEdge>();
        newCommonEdges.addAll(geoA.getCommonEdges());
        newCommonEdges.addAll(geoB.getCommonEdges());
//XXX:  pre-EdgeAttribute		newCommonEdges.add(new PhyloTreeEdge(commonEdge.getSplit(), distCommonEdge, commonEdge.getOriginalID()));

        EdgeAttribute attribCommonEdge = EdgeAttribute.difference(t1.getAttribOfSplit(commonEdge), t2.getAttribOfSplit(commonEdge));
        attribCommonEdge.ensurePositive();
        newCommonEdges.add(new PhyloTreeEdge(commonEdge.asSplit(), attribCommonEdge, commonEdge.getOriginalID()));

        return new Geodesic(RatioSequence.interleave(geoA.getRS(), geoB.getRS()), newCommonEdges);
    }

    /**
     * Computes all the inter-tree distances between the trees in trees using algorithm, and returns them in a
     * matrix.  Prints the average time on the screen.  If doubleCheck is true, computes each distance both ways,
     * and displays a message if they differ.
     *
     * @param trees
     * @param numTrees
     * @param algorithm
     * @return
     */
    public static Geodesic[][] getAllInterTreeGeodesics(PhyloTree[] trees, int numTrees, String algorithm, boolean doubleCheck) throws IOException {
        Date startTime;
        Date endTime;
        long[][] compTimes = new long[numTrees][numTrees];
        long avgCompTime = 0;
        double[][] dists = new double[numTrees][numTrees];
        Geodesic[][] geos = new Geodesic[numTrees][numTrees];

        if (verbose >= 1) {
            System.out.println("Algorithm: " + algorithm + "\n");
        }

        for (int i = 0; i < numTrees; i++) {
            for (int j = i + 1; j < numTrees; j++) {
                if (algorithm.equals("divide") || algorithm.equals("DivideAndConquerRS") || algorithm.equals("ConjBothEnds") || algorithm.equals("dynamic")) {
                    startTime = new Date();
                    geos[i][j] = getGeodesic2(trees[i], trees[j], algorithm, "geo_" + i + "_" + j);
                    dists[i][j] = geos[i][j].getDist();
                    endTime = new Date();
                    compTimes[i][j] = endTime.getTime() - startTime.getTime();
                    // sum up all the times, then divide by number of computations
                    avgCompTime = avgCompTime + compTimes[i][j];
                    //	System.out.println("computed in time " + compTimes[i][j] + " geos[" + i + "][" + j + "].getRS().getMinAscRSDistance() is " + geos[i][j].getRS().getMinAscRSDistance() + "; commonEdges is " +geos[i][j].getCommonEdges() + "; and leafContributionSquared is " + geos[i][j].getLeafContributionSquared());
                }
                // algorithm just returns a distance
                else {
                    System.out.println("Unknown algorithm " + algorithm + "; exiting.");
                    throw new RuntimeException();
                }
            }
        }

        // compute average
        avgCompTime = avgCompTime / (numTrees * (numTrees - 1) / 2);
        System.out.println("Average dist. computation of " + algorithm + " was " + avgCompTime + " ms for " + numTrees * (numTrees - 1) / 2 + " trees.");

        // we want to doublecheck
        if (doubleCheck) {
            avgCompTime = 0;
            for (int i = 0; i < numTrees; i++) {
                for (int j = i + 1; j < numTrees; j++) {
                    if (algorithm.equals("divide") || algorithm.equals("DivideAndConquerRS") || algorithm.equals("ConjBothEnds") || algorithm.equals("dynamic")) {
                        startTime = new Date();
                        geos[j][i] = getGeodesic2(trees[j], trees[i], algorithm, "geo_" + j + "_" + "i");
                        dists[j][i] = geos[j][i].getDist();
                        endTime = new Date();
                        compTimes[j][i] = endTime.getTime() - startTime.getTime();
                        // sum up all the times, then divide by number of computations
                        avgCompTime = avgCompTime + compTimes[j][i];
                    }
                    // algorithm just returns a distance
                    else {
                        System.out.println("Unknown algorithm " + algorithm + "; exiting.");
                        throw new RuntimeException();
                    }
                    if (truncate(geos[i][j].getDist(), 10) != truncate(geos[j][i].getDist(), 10)) {
                        System.out.println("***" + algorithm + " distances don't match for trees " + i + " and " + j + "***");
                        System.out.println("Dist " + i + " -> " + j + " is " + geos[i][j].getDist() + " but dist " + j + " -> " + i + " is " + geos[j][i].getDist());
                        System.out.println("RS " + i + " -> " + j + "           : " + geos[i][j]);
                        System.out.println("geos[" + i + "][" + j + "].getRS().getAscRSWithMinDist().getDistance() is " + geos[i][j].getRS().getNonDesRSWithMinDist().getDistance() + "; commonEdges is " + geos[i][j].getCommonEdges() + "; and leafContributionSquared is " + geos[i][j].getLeafContributionSquared());


                        System.out.println("RS " + j + " -> " + i + " (reversed): " + geos[j][i].reverse());
                        System.out.println("geos[" + j + "][" + i + "].getRS().getAscRSWithMinDist().getDistance() is " + geos[j][i].getRS().getNonDesRSWithMinDist().getDistance() + "; commonEdges is " + geos[j][i].getCommonEdges() + "; and leafContributionSquared is " + geos[j][i].getLeafContributionSquared());

                    }
                }
            }
//			 compute average
            avgCompTime = avgCompTime / (numTrees * (numTrees - 1) / 2);
            System.out.println("In doubleCheck, average dist. computation of " + algorithm + " was " + avgCompTime + " ms for " + numTrees * (numTrees - 1) / 2 + " trees.");
        } else {
            // XXX: to avoid null pointers in outputting.  Fix this!
            for (int i = 0; i < numTrees; i++) {
                for (int j = i + 1; j < numTrees; j++) {
                    geos[j][i] = geos[i][j];
                }
            }
        }

        return geos;
    }

    /**
     * Returns the geodesic between t1 and t2 corresponding to the shortest distance, assuming that the two trees do not have any
     * edges in common.  Does a depth first search up the paths, choosing outgoing split with the lower ratio not
     * yet explored at each node.  Prunes by comparing distance so far along that path with the min distance found so far.
     *
     * @return
     */
    public static Geodesic getPruned1GeodesicNoCommonEdges(PhyloTree t1, PhyloTree t2) {
        resetTreeDistanceState();
        numMaxPaths = 0;
//		System.out.println("No common edges!");
//		double minDist = -1;
        int numEdges = t1.getEdges().size();

        // I was fixing this to return a geodesic instead of just the distance.


        if (numEdges == 0) {
            return new Geodesic(new RatioSequence());
        }
        if (numEdges == 1) {
            // return a geodesic containing the ratio
            Ratio r = new Ratio();
            r.addEEdge(t1.getEdge(0));
            r.addFEdge(t2.getEdge(0));

            RatioSequence rs = new RatioSequence();
            rs.add(r);

            return new Geodesic(rs);
        }

        Vector<Bipartition> m = t2.getCrossingsWith(t1);

//		System.out.println("In getDistance() crossings are: " + m);

        if (m == null) {
            System.out.println("The trees " + t1 + " and " + t2 + " do not have the same leaf labels.");
            throw new RuntimeException();
        }

        getMaxPathSpacesAsRatioSeqs(m, new RatioSequence(), t1.getEdges(), t2.getEdges());

//		System.out.println("real min ratio sequence was " + minTreeDistRatioSeq);
        return new Geodesic(minTreeDistRatioSeq);

    }

    /**
     * Returns the distance between t1 and t2, assuming that the two trees do not have any
     * edges in common.  Does a depth first search up the paths, choosing outgoing split with the lower ratio not
     * yet explored at each node.  Prunes by comparing the distance to the current node with the shortest distance
     * to that node found so far.  Stores the nodes in a hash table.
     * Doesn't assume t1 and t2 have the same number of edges.
     *
     * @return
     */
    public static Geodesic getPruned2GeodesicNoCommonEdges(PhyloTree t1, PhyloTree t2) {
        resetTreeDistanceState();
        numMaxPaths = 0;
//		System.out.println("No common edges!");
//		double minDist = -1;
        int numEdges1 = t1.getEdges().size();
        int numEdges2 = t2.getEdges().size();


//		System.out.println("t1 = " + t1);
//		System.out.println("t2 = " + t2);

        Vector<Bipartition> m = t2.getCrossingsWith(t1);

//		System.out.println("In getPruned2GeodesicNoCommonEdge() crossings are: " + m);

        if (m == null) {
            System.out.println("The trees " + t1 + " and " + t2 + " do not have the same leaf labels.");
            throw new RuntimeException();
        }

//		 the base case is when m contains all the same element
//		System.out.println("m is " + m);
        Bipartition mEl = m.get(0);
        Boolean inBaseCase = true;
        for (int i = 1; i < m.size(); i++) {
            if (!(m.get(i).equals(mEl))) {
                inBaseCase = false;
                break;
            }
        }
        if (inBaseCase) {
            numBaseCase++;

            RatioSequence rs = new RatioSequence();
            rs.add(new Ratio(myVectorClonePhyloTreeEdge(t1.getEdges()), myVectorClonePhyloTreeEdge(t2.getEdges())));
//			System.out.println("in base case: returning geo " + new Geodesic(rs));
//			System.out.println("trees were: " + t1 + " and " + t2);
//			return new Geodesic(t1.getDistanceFromOrigin() + t2.getDistanceFromOrigin(), rs);
            return new Geodesic(rs);
        }

        getPruned2MaxPathSpacesAsRatioSeqs(m, new RatioSequence(), t1.getEdges(), t2.getEdges());
        // the min ratio sequence found in the above method is stored in minTreeDistRatioSeq

//		System.out.println("real min ratio sequence was " + minTreeDistRatioSeq);
//		System.out.println("combined min ratio sequence was " + minTreeDistRatioSeq.getAscRSWithMinDist());

        if (numEdges1 > 2 && numEdges2 > 2) {
            System.out.println("Calculating distance between " + t1.getLeaf2NumMap().size() + "-leaved trees took " + numPrunes + " prunes, " + numMaxPaths + " shorter paths and " + numNodes + " nodes; min path was " + minTreeDist / firstTreeDist + " of first path");
//			if (numEdges1 ==2) {
//				System.out.println("ratio sequence is " + minTreeDistRatioSeq.getAscRSWithMinDist());
//			}
        }
        return new Geodesic(minTreeDistRatioSeq);
    }


    /**
     * Once an split has been brought in, uses the fact that the tree now has a common split.
     * Doesn't assume the two trees contain the same number of edges.
     * XXX: edges in ratio sequences will be a bit wonky - might not be original edges
     *
     * @return
     */
    public static Geodesic getDivideAndConquerGeodesicNoCommonEdges(PhyloTree t1, PhyloTree t2) {
        Bipartition minEl;
        PhyloTree newT1 = null;
        Geodesic geo = null;
        Geodesic potentialGeoToReturn = null;
        RatioSequence rsForPotentialGeoToReturn = null;

        Vector<Bipartition> m = t2.getCrossingsWith(t1);

//		System.out.println("In getDistance() crossings are: " + m);

        if (m == null) {
            System.out.println("The trees " + t1 + " and " + t2 + " do not have the same leaf labels.");
            throw new RuntimeException();
        }

        // the base case is when m contains all the same element
        Bipartition mEl = m.get(0);
        Boolean inBaseCase = true;
        for (int i = 1; i < m.size(); i++) {
            if (!(m.get(i).equals(mEl))) {
                inBaseCase = false;
                break;
            }
        }
        if (inBaseCase) {
            numBaseCase++;

            RatioSequence rs = new RatioSequence();
            rs.add(new Ratio(myVectorClonePhyloTreeEdge(t1.getEdges()), myVectorClonePhyloTreeEdge(t2.getEdges())));
//			System.out.println("in base case: returning geo " + new Geodesic(t1.getDistanceFromOrigin() + t2.getDistanceFromOrigin(), rs));
//			return new Geodesic(t1.getDistanceFromOrigin() + t2.getDistanceFromOrigin(), rs);
            return new Geodesic(rs);
        }

/*		if (t1.getLeaf2NumMap().size() > 25) {
			System.out.println("Called getDivideAndConquer... for trees t1 " + t1 + " and t2 " + t2);
		}*/


//		System.out.println("m is " + m);
//		get the minimal elements from m
        Vector<Bipartition> minEls = getMinElements(m);

        // for each min element, divide problem on the new edges
        for (int i = 0; i < minEls.size(); i++) {
            minEl = minEls.get(i);
//			System.out.println("minEl is " + minEl + " of " + minEls );
            Ratio ratio = calculateRatio(minEl, m, t1.getEdges(), t2.getEdges());

            // remove the edges from t1 corresponding to 1's in minEl
            newT1 = new PhyloTree(myVectorClonePhyloTreeEdge(t1.getEdges()), myVectorCloneString(t1.getLeaf2NumMap()));
            newT1.removeEdgesIndicatedByOnes(minEl);

            // add all edges equal to the minEl in m
            for (int j = 0; j < m.size(); j++) {
                if (m.get(j).equals(minEl)) {
//					System.out.println("Added split!");
//					System.out.println("newT1 before adding split " + t2.getEdge(j) + ": " + newT1);
                    newT1.addEdge(t2.getEdge(j).clone());
//					System.out.println("newT1 after adding: " + newT1);
                }
            }

            // find the shortest distance from newT1 to t2
//			System.out.println("Calling getGeodesic on newT1 " + newT1 + "; old t1 was " + t1);
//			System.out.println("num edges in newT1 = " + newT1.getEdges().size());

            Geodesic potentialGeo = getGeodesicRecursive(newT1, t2, "divide");
//			System.out.println("potential geo is " + potentialGeo);

            // combine the geodesic we found with minEl to get the geodesic to compare
            // with the other geodesics generated from other minEls
            rsForPotentialGeoToReturn = new RatioSequence();
            rsForPotentialGeoToReturn.add(ratio);
            rsForPotentialGeoToReturn.addAll(potentialGeo.getRS());
//			System.out.println("rsForGeoToReturn is " + rsForGeoToReturn);
            potentialGeoToReturn = new Geodesic(rsForPotentialGeoToReturn);

            if (geo == null) {
                geo = potentialGeoToReturn.clone();
            } else if (potentialGeoToReturn.getDist() < geo.getDist()) {
                geo = potentialGeoToReturn.clone();
            }

        }
//		System.out.println("Returning geo " + geo);
        return geo;
    }

    /**
     * Once an split has been brought in, uses the fact that the tree now has a common split.
     * Doesn't assume the two trees contain the same number of edges.
     * XXX: edges in ratio sequences will be a bit wonky - might not be original edges
     *
     * @return
     */
    public static RatioSequence getDivideAndConquerRSNoCommonEdges(PhyloTree t1, PhyloTree t2) {
        Bipartition minEl;
        PhyloTree newT1 = null;
        RatioSequence rs = null;
        RatioSequence rsForPotentialGeoToReturn = null;

        Vector<Bipartition> m = t2.getCrossingsWith(t1);

//		System.out.println("In getDistance() crossings are: " + m);

        if (m == null) {
            System.out.println("The trees " + t1 + " and " + t2 + " do not have the same leaf labels.");
            throw new RuntimeException();
        }

        // the base case is when m contains all the same element
        Bipartition mEl = m.get(0);
        Boolean inBaseCase = true;
        for (int i = 1; i < m.size(); i++) {
            if (!(m.get(i).equals(mEl))) {
                inBaseCase = false;
                break;
            }
        }
        if (inBaseCase) {
            numBaseCase++;

            rs = new RatioSequence();
            rs.add(new Ratio(myVectorClonePhyloTreeEdge(t1.getEdges()), myVectorClonePhyloTreeEdge(t2.getEdges())));
//			System.out.println("in base case: returning rs " + rs);
//			return new Geodesic(t1.getDistanceFromOrigin() + t2.getDistanceFromOrigin(), rs);
            return rs;
        }

/*		if (t1.getLeaf2NumMap().size() > 25) {
			System.out.println("Called getDivideAndConquer... for trees t1 " + t1 + " and t2 " + t2);
		}*/


//		System.out.println("m is " + m);
//		get the minimal elements from m
        Vector<Bipartition> minEls = getMinElements(m);

        // for each min element, divide problem on the new edges
        for (int i = 0; i < minEls.size(); i++) {
            minEl = minEls.get(i);
//			System.out.println("minEl is " + minEl + " of " + minEls );
            Ratio ratio = calculateRatio(minEl, m, t1.getEdges(), t2.getEdges());

            // remove the edges from t1 corresponding to 1's in minEl
            newT1 = new PhyloTree(myVectorClonePhyloTreeEdge(t1.getEdges()), myVectorCloneString(t1.getLeaf2NumMap()));
            newT1.removeEdgesIndicatedByOnes(minEl);

            // add all edges equal to the minEl in m
            for (int j = 0; j < m.size(); j++) {
                if (m.get(j).equals(minEl)) {
//					System.out.println("Added split!");
//					System.out.println("newT1 before adding split " + t2.getEdge(j) + ": " + newT1);
                    newT1.addEdge(t2.getEdge(j).clone());
//					System.out.println("newT1 after adding: " + newT1);
                }
            }

            // find the shortest distance from newT1 to t2
//			System.out.println("Calling getGeodesic on newT1 " + newT1 + "; old t1 was " + t1);
//			System.out.println("num edges in newT1 = " + newT1.getEdges().size());

//			System.out.println("Before splitOnCommonEdge aTreesNoCommonEdges is " + aTreesNoCommonEdges);
            int sizeOfATreesNoCommonEdges = aTreesNoCommonEdges.size();
            splitOnCommonEdge(newT1, t2);
//			System.out.println("After splitOnCommonEdge aTreesNoCommonEdges is " + aTreesNoCommonEdges);

//			System.out.println("Calling getDivideAndConquerRSNoCommonEdges for " + aTreesNoCommonEdges.get(aTreesNoCommonEdges.size() -1) + " and " + bTreesNoCommonEdges.get(bTreesNoCommonEdges.size() -1) );

            RatioSequence rs2 = getDivideAndConquerRSNoCommonEdges(aTreesNoCommonEdges.lastElement(), bTreesNoCommonEdges.lastElement());
//			System.out.println("rs2 is " + rs2);
            aTreesNoCommonEdges.remove(aTreesNoCommonEdges.lastElement());
            bTreesNoCommonEdges.remove(bTreesNoCommonEdges.lastElement());

            RatioSequence potentialRS;
            // check if splitting on the common edges gave us one tree (if only 1 split left in each tree after
            // removing common split) or two trees
            if (sizeOfATreesNoCommonEdges != aTreesNoCommonEdges.size()) {
                // two trees
                RatioSequence rs1 = getDivideAndConquerRSNoCommonEdges(aTreesNoCommonEdges.lastElement(), bTreesNoCommonEdges.lastElement());
                aTreesNoCommonEdges.remove(aTreesNoCommonEdges.lastElement());
                bTreesNoCommonEdges.remove(bTreesNoCommonEdges.lastElement());
                potentialRS = RatioSequence.interleave(rs1, rs2);
            } else {
                // one tree
                potentialRS = rs2;
            }

//			System.out.println("potential geo is " + potentialGeo);

            // combine the geodesic we found with minEl to get the geodesic to compare
            // with the other geodesics generated from other minEls
            rsForPotentialGeoToReturn = new RatioSequence();
            rsForPotentialGeoToReturn.add(ratio);
            rsForPotentialGeoToReturn.addAll(potentialRS);
//			System.out.println("rsForGeoToReturn is " + rsForPotentialGeoToReturn);

            if (rs == null) {
                rs = rsForPotentialGeoToReturn.clone();
            } else if (rsForPotentialGeoToReturn.getMinNonDesRSDistance() < rs.getMinNonDesRSDistance()) {
                rs = rsForPotentialGeoToReturn.clone();
            }

        }
//		System.out.println("Returning geo " + geo);
        return rs;
    }


    /**
     * Calculates the ratio if we remove the minimum element minEl from tree with partition poset given by m.
     * XXX:  fix description
     *
     * @param minEl
     * @param m
     * @return
     */
    public static Ratio calculateRatio(Bipartition minEl, Vector<Bipartition> m, Vector<PhyloTreeEdge> eEdges, Vector<PhyloTreeEdge> fEdges) {
        Ratio ratio = new Ratio();

        // add all edges in m equal to minEl to the ratio
        for (int i = 0; i < m.size(); i++) {
            if ((m.get(i) != null) && (m.get(i).equals(minEl))) {
                // add the f edges with this minEl to the f edges of the ratio
                ratio.addFEdge(fEdges.get(i).clone());
            }
        }

        // set 1's in columns corresponding to removed e's to be 0,
        // and add that e split to the ratio.
        Iterator<Integer> eEdgesIter = binaryToVector(minEl.getPartition()).iterator();
        while (eEdgesIter.hasNext()) {
            int eEdge = (Integer) eEdgesIter.next();
//			System.out.println("eEdges are " + eEdges);
//			System.out.println("and we want e-split " + eEdge + " as specified by minEl " + minEl);
            // following for debugging
            if (eEdge > eEdges.size()) {
                System.out.println("Trying to add split " + eEdge + " as indicated by minEl " + minEl);
                System.out.println("but edges are " + eEdges);
            }
            ratio.addEEdge(eEdges.get(eEdge));
        }
        return ratio;
    }


    /** Recursive method for finding all the ratio sequences.
     *
     * @param m
     * @param ratioSeq
     */
/*	public static void calculateFinalRatioSeqs(Vector<PhyloTreeEdge> m, RatioSequence ratioSeq, double minDist, RatioSequence rsWithMinDist, Vector<PhyloTreeEdge> eEdges) {
		Vector<PhyloTreeEdge> minEls = getMinElements(m);
		double d =0;
	    PhyloTreeEdge minEl = null;
	    Ratio minRatio = null;

		if (minEls.size() == 0) {
			// base case: we are at the end of the chain, so add to final ratioSeq
			System.out.println("Ratio sequence " + ratioSeq + " has distance " + ratioSeq.getDistance() );

			d = ratioSeq.getAscRSWithMinDist().getDistance();
			if (d < minDist) {
				minDist = d;
				rsWithMinDist = ratioSeq.clone();
			}
			return;
		}

//		 now find correct ratio of bipartition minEl  (numerator is based on 1's in bipartition;
		// denominator based on number of elements in m with the same bipartition)
		Iterator<PhyloTreeEdge> minElsIter = minEls.iterator();
		while (minElsIter.hasNext()) {
			minEl = minElsIter.next();

			Ratio ratio = calculateRatio(minEl, m, eEdges);
			if ( (minRatio == null) || (minRatio.getRatio() > ratio.getRatio() ) ) {
				minRatio = ratio.clone();
			}
		}

		// otherwise, call calculateFinalRatioSeqs for each min element in minEls
		Iterator<Bipartition> minElsIter = minEls.iterator();
		while (minElsIter.hasNext()) {
			Vector<Bipartition> newM = (Vector<Bipartition>) myVectorCloneBipartition(m);
			RatioSequence newRatioSeq = (RatioSequence) ratioSeq.clone();
			Bipartition minEl = minElsIter.next();

			// create a new Ratio for the edges being dropped and added
			Ratio ratio = new Ratio();

			// set all edges in m equal to minEl to be null
			// and add to the ratio.
			for (int i = 0; i < newM.size();i++) {
				if ((newM.get(i) != null) && (newM.get(i).equals(minEl))) {
					// add the f edges with this minEl to the f edges of the ratio
					ratio.addFEdge(fTree.getEdges().get(i));
					newM.set(i, null);
				}
			}

			// set 1's in columns corresponding to removed e's to be 0,
			// and add that e split to the ratio.
			Iterator eEdgesIter = binaryToVector(minEl.getEdge()).iterator();
			while (eEdgesIter.hasNext()) {
				int eEdge = (Integer) eEdgesIter.next();
				newM = zeroCol(eEdge, newM);
				ratio.addEEdge(eTree.getEdges().get(eEdge));
			}



			// push ratio onto ratio seq and call calculateFinalRatioSeqs recursively
			newRatioSeq.add(ratio);
			this.calculateFinalRatioSeqs(newM,newRatioSeq);
		}

	}*/

    /**
     * Recursive method for finding all the ratio sequences.
     *
     * @param m
     * @param ratioSeq
     */
    public static void getMaxPathSpacesAsRatioSeqs(Vector<Bipartition> m, RatioSequence ratioSeq, Vector<PhyloTreeEdge> eEdges, Vector<PhyloTreeEdge> fEdges) {
        // returns vector containing each minimal e-set.
        Vector<Bipartition> minEls = getMinElements(m);
        Bipartition minEl;
        Ratio ratio;

//lll	System.out.println("at start of getMaxPathSpacesAsRatioSeqsmin, minTreeDist is " + minTreeDist + " and minTreeDistRatioSeq is " + minTreeDistRatioSeq);


        if (minEls.size() == 0) {
            // base case: we are at the end of the chain, so add to final ratioSeq
            // first case to avoid memory overflow
            numMaxPaths++;
//			if (numMaxPaths%10== 0) {
//				System.out.println("Path " + numMaxPaths);
//			}
//			if (eEdges.size() > 5) {
            double dist = ratioSeq.getNonDesRSWithMinDist().getDistance();
            if (minTreeDist < 0 || dist < minTreeDist) {
//					System.out.println("New min dist of " + dist + " for rs " + ratioSeq);
                minTreeDist = dist;
                minTreeDistRatioSeq = ratioSeq.clone();
            }
            if (firstTreeDist == -1) {
                firstTreeDist = dist;
            }
/*			}
			else {
				finalRatioSeqs.add(ratioSeq.clone());
			} */
            return;
        }

        Vector<Bipartition> sortedMinEls = new Vector<Bipartition>();
        Vector<Ratio> sortedMinElRatios = new Vector<Ratio>();

        // reset minRatio
//		minRatio = null;

        // for each min element, find the one with the smallest corresponding ratio
        Iterator<Bipartition> minElsIter = minEls.iterator();
        // last index of non-null entry in sorted array of min elements
        while (minElsIter.hasNext()) {
            minEl = minElsIter.next();

            ratio = calculateRatio(minEl, m, eEdges, fEdges);

            Boolean inserted = false;
            // store both in sortedMinEls in order of ascending ratios
            for (int k = 0; k < sortedMinElRatios.size(); k++) {
                if (sortedMinElRatios.get(k).getRatio() > ratio.getRatio()) {
                    sortedMinElRatios.add(k, ratio.clone());
                    sortedMinEls.add(k, minEl.clone());
                    inserted = true;
                    break;
                }
            }
            if (!inserted) {
                // else bigger than all ratios currently in the vector
                sortedMinElRatios.add(ratio.clone());
                sortedMinEls.add(minEl.clone());
            }
//			System.out.println("sortedMinElRatios is " + sortedMinElRatios);
        }

        for (int j = 0; j < sortedMinElRatios.size(); j++) {
            minEl = sortedMinEls.get(j);
            ratio = sortedMinElRatios.get(j);
//			System.out.println("minEl is " + minEl);
//			System.out.println("ratio is " + ratio);

            Vector<Bipartition> newM = removeMinElFrom(m, minEl);

            // push ratio onto ratio seq and call getMaxPathSpacesAsRatioSeqs recursively
            ratioSeq.add(ratio);

            // this path is greater than the min path so far, so don't follow any further.
            double pruneDist = ratioSeq.getNonDesRSWithMinDist().getDistance();
            if (eEdges.size() > 5 && minTreeDist > 0 && pruneDist > minTreeDist) {
//				System.out.println("Pruning at distance " + pruneDist + " with minTreeDist " +  minTreeDist + " for rs: " + ratioSeq);
                numPrunes++;
                if (numPrunes % 1000 == 0) {
//					System.out.println("prunes: " + numPrunes);
                }
                // go on to next min element in the for loop
            } else {
                // this could still be a min path, so continue
                getMaxPathSpacesAsRatioSeqs(newM, ratioSeq, eEdges, fEdges);
            }
            // remove this ratio from the ratio sequence, as we will now consider a new path.
            ratioSeq.remove(ratio);
        }
        minEls = null;
    }

    /**
     * Recursive method for finding all the ratio sequences.
     *
     * @param m
     * @param ratioSeq
     */
    public static void getPruned2MaxPathSpacesAsRatioSeqs(Vector<Bipartition> m, RatioSequence ratioSeq, Vector<PhyloTreeEdge> eEdges, Vector<PhyloTreeEdge> fEdges) {
//		stepForDebugging++;

        // returns vector containing each minimal e-set.
        Vector<Bipartition> minEls = getMinElements(m);

//		System.out.println("m is " + m);
//		System.out.println("minEls is " + minEls);

        Bipartition minEl;
        Ratio ratio;
        double minNodeDist;

//lll	System.out.println("at start of getMaxPathSpacesAsRatioSeqsmin, minTreeDist is " + minTreeDist + " and minTreeDistRatioSeq is " + minTreeDistRatioSeq);

        if (minEls.size() == 0) {
            // base case: we are at the end of the chain, so add to final ratioSeq
            // first case to avoid memory overflow
            numMaxPaths++;
/*			if (numMaxPaths%5== 0) {
				System.out.println("Path " + numMaxPaths);
			} */

            double dist = ratioSeq.getNonDesRSWithMinDist().getDistance();
            if (pathToSearch == 1) {
                System.out.println("In base case with distance " + dist + " and rs " + ratioSeq.getNonDesRSWithMinDist());
            }
            if (minTreeDist < 0 || dist < minTreeDist) {
//				System.out.println("New min dist of " + dist + " for rs " + ratioSeq);
                minTreeDist = dist;
                minTreeDistRatioSeq = ratioSeq.clone();
            }
            if (firstTreeDist == -1) {
                firstTreeDist = dist;
            }
            return;
        }

        Vector<Bipartition> sortedMinEls = new Vector<Bipartition>();
        Vector<Ratio> sortedMinElRatios = new Vector<Ratio>();

        // order the min elements by their ratios
        Iterator<Bipartition> minElsIter = minEls.iterator();
        while (minElsIter.hasNext()) {
            minEl = minElsIter.next();

            ratio = calculateRatio(minEl, m, eEdges, fEdges);

            Boolean inserted = false;
            // store both in sortedMinEls in order of ascending ratios
            for (int k = 0; k < sortedMinElRatios.size(); k++) {
                if (sortedMinElRatios.get(k).getRatio() > ratio.getRatio()) {
                    sortedMinElRatios.add(k, ratio.clone());
                    sortedMinEls.add(k, minEl.clone());
                    inserted = true;
                    break;
                }
            }
            if (!inserted) {
                // else bigger than all ratios currently in the vector
                sortedMinElRatios.add(ratio.clone());
                sortedMinEls.add(minEl.clone());
            }

        }

        if (pathToSearch == 1 && minTreeDistRatioSeq == null) {
            System.out.println("sortedMinElRatios is " + sortedMinElRatios);
        }
        for (int j = 0; j < sortedMinElRatios.size(); j++) {
            minEl = sortedMinEls.get(j);
            ratio = sortedMinElRatios.get(j);
//			System.out.println("minEl is " + minEl);
//			System.out.println("ratio is " + ratio);

            Vector<Bipartition> newM = removeMinElFrom(m, minEl);


            // newM represents the node we are moving to.
            // Use the hashtable to see if we have already visited this node,
            // and only continue up this path if the node is unvisited or has a long distance to it so far.
            Double minNodeDistObject = (Double) nodeHashtable.get(newM.toString());
//			 add this ratio to the ratio sequence, and calculate the distance
            ratioSeq.add(ratio);

            double pruneDist = ratioSeq.getNonDesRSWithMinDist().getDistance();
/*			if (printOut > 0) {
				System.out.println("pruneDist is " + pruneDist);
			}*/
            if (minNodeDistObject != null) {
/*				if (printOut >0) {
					System.out.println("Found node " + newM + " in the hash table");	// 11
				}*/
                minNodeDist = minNodeDistObject.doubleValue();

                if (pruneDist < minNodeDist) {
//					 this distance is better than the best distance to this node
                    // and the best distance overall,
                    // so replace value in hashtable and keep going
                    nodeHashtable.put(newM.toString(), new Double(pruneDist));
                    if (pruneDist < minTreeDist) {
//					 this distance is better than the best distance overall, so keep going
                        if (pathToSearch == 1) {
//							System.out.println("Calling getPruned2MaxPathSpacesAsRatioSeqs for newM " + newM);
                        }

                        getPruned2MaxPathSpacesAsRatioSeqs(newM, ratioSeq, eEdges, fEdges);
                    } else {        // for debugging
                        numPrunes++;
                    }
                } else {        // for debugging
                    numPrunes++;
                }
            } else {
                numNodes++;
                // no entry for this node in the hash table yet,
                // so make one and continue
                nodeHashtable.put(newM.toString(), new Double(pruneDist));
                getPruned2MaxPathSpacesAsRatioSeqs(newM, ratioSeq, eEdges, fEdges);
            }
            ratioSeq.remove(ratio);
        }
        minEls = null;
    }


    /**
     * Converts a BitSet representing a 0-1 vector into a Vector containing the integers which
     * correspond to the posititions of the 1's in bin.
     *
     * @param bin
     * @return
     */
    public static Vector<Integer> binaryToVector(BitSet bin) {
        Vector<Integer> v = new Vector<Integer>();

        for (int i = 0; i < bin.length(); i++) {
            if (bin.get(i)) {
                v.add(i);
            }
        }
        return v;
    }

    /**
     * Given a set of bipartitions stored in vector m, returns a vector containing the
     * minimum bipartitions - namely every vector whose one's are not contained in the
     * set of some other vector's one's.
     * If two bipartitions are the same, just returns one of them.
     * Precondition: m is not null, but some elements may be null and should be skipped
     *
     * @param m
     * @return
     */
    public static Vector<Bipartition> getMinElements(Vector<Bipartition> m) {
//		System.out.println("m is " + m);
        Vector<Bipartition> minEls = new Vector<Bipartition>();
/*		boolean[] isMinEl = new boolean[m.size()];

		// run through m and set isMinEl to false if the corresponding entry in m is null
		for (int i = 0; i < m.size(); i++) {
			if (m.get(i) ==null) {
				isMinEl[i] = false;
			}
			else {
				isMinEl[i] = true;  // so that all elements are initialized to true or false
									// by default everything seems to initialize to false.
			}
		}

		// for each element in m, add to minEls unless:
		// 1) it contains some other element in m (then don't add to minEls)
		// 2) it equals some other element in m (then only add to minEls if has lowest index)
		// 3) it has already been marked as a non-minimal element by putting a false in the right place
		// in the array isMinEl
		for (int i = 0; i < m.size(); i++ ) {
			if (isMinEl[i] == false) {
				// then skip this element
				continue;
			}
			Bipartition lowerIndexedEdge = m.get(i);
			for (int j = i +1; j< m.size(); j++) {
				if (isMinEl[j] == false) {
					// then skip this element
					continue;
				}
				Bipartition higherIndexedEdge = m.get(j);
				if (lowerIndexedEdge.equals(higherIndexedEdge)) {
					isMinEl[j] = false;
				}
				else if (lowerIndexedEdge.properlyContains(higherIndexedEdge)) {
					isMinEl[i] = false;
				}
				else if (lowerIndexedEdge.equals(higherIndexedEdge)) {
					isMinEl[j] = false;
				}
			}
		}
		for (int i = 0; i < m.size(); i++) {
			if (isMinEl[i] != false) {
				minEls.add(m.get(i));
			}
		}
		return minEls;
	} */

/*	commented out Dec. 29, 2007 because seems more complicated than necessary,
 * and it might contain a bug.  However it is apparantly faster than the implementation above.*/
        boolean addToMinEls = true;
        // step through the f edges (elements of m)
        for (int i = 0; i < m.size(); i++) {
            if (m.get(i) == null) {
                // skip this element
                continue;
            }
            addToMinEls = true;

            // step through the elements in minEls (so none are null)
            int j = 0;
            while (j < minEls.size()) {
                Bipartition minEl = minEls.get(j);
//				 if minEl properly contains fi, then remove minEl and leave flag addToMinEls= true
                if (minEl.properlyContains(m.get(i))) {
                    minEls.remove(minEl);
                } else {
                    j++;
                }
//				 if fi contains minEl, change addToMinEls to false
                if (m.get(i).contains(minEl)) {
                    addToMinEls = false;
                }
//				 (above two cases can never both be true for different minEls.
            }
//			 if addToMinEls = true, add fi to minEls  (this covers the other case, where fi neither contains nor is contained in any minEls)
            if (addToMinEls) {
                minEls.add(m.get(i));
            }
        }
//		System.out.println("minEls are " + minEls);
        return minEls;
    }


    /**
     * Removes all the PhyloTreeEdge with split = 0 from the vector v.
     *
     * @param v
     * @return
     */
    public static Vector<PhyloTreeEdge> deleteZeroEdges(Vector<PhyloTreeEdge> v) {
        int k = 0;
        while (k < v.size()) {
//			System.out.println("v.get(k) is " + v.get(k));
            if (v.get(k).isZero()) {
                v.remove(k);
            } else {
                k++;
            }
        }
        return v;
    }

    public static Vector<PhyloTreeEdge> myVectorClonePhyloTreeEdge(Vector<PhyloTreeEdge> v) {
        if (v == null) {
            return null;
        }

        Vector<PhyloTreeEdge> newV = new Vector<PhyloTreeEdge>();
        for (int i = 0; i < v.size(); i++) {
            if (v.get(i) == null) {
                newV.add(null);
            } else {
                newV.add(v.get(i).clone());
            }
        }
        return newV;
    }

    public static Vector<RatioSequence> myVectorCloneRatioSequence(Vector<RatioSequence> v) {
        if (v == null) {
            return null;
        }

        Vector<RatioSequence> newV = new Vector<RatioSequence>();
        for (int i = 0; i < v.size(); i++) {
            if (v.get(i) == null) {
                newV.add(null);
            } else {
                newV.add(v.get(i).clone());
            }
        }
        return newV;
    }

    /** Returns the bipartition from subSetOfM of m, the partition poset, with the lowest ratio as a Ratio object.
     *
     * @param m
     * @param eEdges
     * @param fEdges
     * @return
     */
	/*	public static Ratio getMinRatio(Vector<Bipartition> subSetOfM, Vector<Bipartition> m, Vector<PhyloTreeEdge> eEdges, Vector<PhyloTreeEdge> fEdges) {
			Ratio minRatio = null;

			// for each min element, find the one with the smallest corresponding ratio
			Iterator<Bipartition> subSetIter = subSetOfM.iterator();
			while (subSetIter.hasNext()) {
				Bipartition b = subSetIter.next();

				Ratio ratio = calculateRatio(b, m, eEdges, fEdges );
				if ( (minRatio == null) || (minRatio.getRatio() > ratio.getRatio() ) ) {
					minRatio = ratio.clone();
					minMinEl = minEl.clone();
				}
			}
		}
		*/

    /**
     * Returns a new Vector<Bipartition> equal to m, but with minEl removed.
     *
     * @param m
     * @param minEl
     * @return
     */
    public static Vector<Bipartition> removeMinElFrom(Vector<Bipartition> m, Bipartition minEl) {

        Vector<Bipartition> newM = myVectorCloneBipartition(m);
        // remove the edges in minRatio from m and repeat
        //	 set all edges in m equal to minEl to be null
        for (int i = 0; i < newM.size(); i++) {
            if (newM.get(i) != null && newM.get(i).equals(minEl)) {
                newM.set(i, null);
            }
        }

        // set 1's in columns corresponding to removed e's to be 0,
        Iterator<Integer> eEdgesIter = binaryToVector(minEl.getPartition()).iterator();
        while (eEdgesIter.hasNext()) {
            int eEdge = (Integer) eEdgesIter.next();
            newM = zeroCol(eEdge, newM);
        }
        return newM;
    }

    public static Vector<String> myVectorCloneString(Vector<String> v) {
        if (v == null) {
            return null;
        }

        Vector<String> newV = new Vector<String>();
        for (int i = 0; i < v.size(); i++) {
            if (v.get(i) == null) {
                newV.add(null);
            } else {
                newV.add(new String(v.get(i)));
            }
        }
        return newV;
    }

    /**
     * Truncates the number d by p places.
     *
     * @param d
     * @param p
     * @return
     */
    public static double truncate(double d, int p) {
        return Math.floor(d * Math.pow(10, p)) / Math.pow(10, p);
    }

    public static double round(double d, int p) {
        return ((double) Math.round(d * Math.pow(10, p))) / Math.pow(10, p);
    }

    /**
     * Reads in all the phylogenetic trees from the file inFileName.
     * The trees should be one per line, in the Newick format.
     * There can also be a ";" at the end of each line.
     *
     * @param inFileName
     * @return
     */
    public static PhyloTree[] readInTreesFromFile(String inFileName) throws IOException {
        int numTrees = 0;  // count the number of trees read in
        Vector<String> stringTrees = new Vector<String>();

        BufferedReader inputStream = null;

        try {
            inputStream =
                    new BufferedReader(new FileReader(inFileName));

            String l;
            while ((l = inputStream.readLine()) != null) {
                // check for blank lines
                if (l != "") {
                    stringTrees.add(l);
                    numTrees++;
                }
            }
            if (inputStream != null) {
                inputStream.close();
            }
        } catch (IOException e) {
            System.out.println("Error opening or reading from " + inFileName + ": " + e.getMessage());
            throw e;
        }


        // convert the Vector of strings representing trees into an array of PhyloTrees
        PhyloTree[] trees = new PhyloTree[numTrees];
        for (int i = 0; i < numTrees; i++) {
            trees[i] = new PhyloTree(stringTrees.get(i), rooted);
            if (TreeDistance.normalize) {
                trees[i].normalize();
            }
            //       	System.out.println("Tree " + i + ": " + trees[i]);
        }

        return trees;
    }

    /**
     * Open file fileName, reads in the trees, and outputs the distances computed by algorithm.
     * Assumes first line of file is number of trees, and then one tree per line.
     */
    public static void computeAllInterTreeGeodesicsFromFile(String inFileName, String outFileName, String algorithm, boolean doubleCheck) throws IOException {


        PhyloTree[] trees = readInTreesFromFile(inFileName);
        int numTrees = trees.length;

        if (numTrees < 2) {
            System.out.println("Error:  tree file must contain at least 2 trees");
            throw new RuntimeException();
        }

        if (verbose >= 1) {
            System.out.println("" + numTrees + " trees read in from " + inFileName);
        }
        Geodesic[][] geos = getAllInterTreeGeodesics(trees, numTrees, algorithm, doubleCheck);

        // print distances to file
        PrintWriter outputStream = null;

        // Outputs the distances in a column, with the first two columns being the trees numbers and the third
        // number the geodesic distance between those trees
        try {
            outputStream = new PrintWriter(new FileWriter(outFileName));

            for (int i = 0; i < numTrees - 1; i++) {
                for (int j = i + 1; j < numTrees; j++) {
/*    				if (verbose >0) {
    					System.out.println("Geo " + i + " -> " + j + " is " + geos[i][j]);
    					System.out.println(geos[i][j].toStringVerbose(trees[i], trees[j]));
    				}*/
                    outputStream.println(i + "\t" + j + "\t" + round(geos[i][j].getDist(), 6));
                }
                outputStream.println();
            }
            if (outputStream != null) {
                outputStream.close();
            }
        } catch (IOException e) {
            System.out.println("Error opening or writing to " + outFileName + ": " + e.getMessage());
            throw e;
        }
    }

    /**
     * Compares algorithms alg1 and alg2 by checking that they give the same inter-distances for the input trees,
     * and by displaying the average time to compute these distances for each algorithm.
     *
     * @param trees
     * @param numTrees
     * @param alg1
     * @param alg2
     */
    public static void compareAlgorithms(PhyloTree[] trees, int numTrees, String alg1, String alg2) throws IOException {

        Date startTime, endTime;
        long alg1Time, alg2Time;
        long avgAlg1Time = 0;
        long avgAlg2Time = 0;
        double alg1Dist, alg2Dist;
        Geodesic alg1Geo, alg2Geo;

        // check input algorithms are valid
        if (!(alg1.equals("divide") || alg1.equals("DivideAndConquerRS") || alg1.equals("ConjBothEnds") || alg1.equals("dynamic")) || !(alg2.equals("divide") || alg2.equals("ConjBothEnds") || alg2.equals("dynamic"))) {
            System.out.println("Error:  either " + alg1 + " or " + alg2 + " is an invalid algorithm.");
            throw new RuntimeException();
        }

        for (int i = 0; i < numTrees; i++) {
            for (int j = 0; j < numTrees; j++) {
                if (i == j) {
                    continue;
                }
                startTime = new Date();
                alg1Geo = getGeodesic2(trees[i], trees[j], alg1, "geo_alg1_" + i + "_" + j);
                endTime = new Date();
                alg1Dist = alg1Geo.getDist();
                alg1Time = endTime.getTime() - startTime.getTime();
                // sum up all the times, then divide by number of computations
                avgAlg1Time = avgAlg1Time + alg1Time;

                startTime = new Date();
                alg2Geo = getGeodesic2(trees[i], trees[j], alg2, "geo_alg2_" + i + "_" + j);
                endTime = new Date();
                alg2Dist = alg2Geo.getDist();
                alg2Time = endTime.getTime() - startTime.getTime();
                // sum up all the times, then divide by number of computations
                avgAlg2Time = avgAlg2Time + alg2Time;

                if (truncate(alg1Dist, 10) != truncate(alg2Dist, 10)) {
                    System.out.println("***" + alg1 + " and " + alg2 + " distances don't match for trees " + i + " -> " + j + "***");
                    System.out.println(alg1 + " dist is " + alg1Dist + " but " + alg2 + " dist is " + alg2Dist);
                    System.out.println("RS " + i + " -> " + j + " for algorithm " + alg1 + ":");
                    ;
                    System.out.println("" + alg1Geo);
                    System.out.println("RS " + i + " -> " + j + " for algorithm " + alg2 + ":");
                    ;
                    System.out.println("" + alg2Geo);
                }
            }
        }

        // compute average
        avgAlg1Time = avgAlg1Time / (numTrees * (numTrees - 1));
        avgAlg2Time = avgAlg2Time / (numTrees * (numTrees - 1));
        System.out.println("Average dist. computation of " + alg1 + " was " + avgAlg1Time + " ms for " + numTrees * (numTrees - 1) + " trees.");
        System.out.println("Average dist. computation of " + alg2 + " was " + avgAlg2Time + " ms for " + numTrees * (numTrees - 1) + " trees.");
    }

    /**
     * Help message (ie. which arguments can be used, etc.)
     */
    public static void displayHelp() {
        System.out.println("Command line syntax:");
        System.out.println("geodeMAPS [options] treefile");
        System.out.println("Optional arguments:");
        System.out.println("\t -a <algorithm> \t uses <algorithm> to compute the geodesic distance.  Current options are divide to run GeodeMaps-Divide, and dynamic to run GeodeMaps-Dynamic.  The default is dynamic.");
        System.out.println("\t -d \t double check results, by computing each distance with the target tree as the starting tree and vice versa; default is false");
        System.out.println("\t -h || --help \t displays this message");
        System.out.println("\t -n \t normalize (vector of the lengths of all edges has length 1)");
        System.out.println("\t -o <outfile> \t store the output in the file <outfile>");
        System.out.println("\t -u \t unrooted trees (default is rooted trees)");
        System.out.println("\t -v || --verbose \t verbose output");
    }


    /**
     * @param args
     */
    public static void main(String[] args) {

        // variables maybe changed by program arguments
        String algorithm = "dynamic";   //XXX make the defaults be constants
        String treeFile = "";
        String outFile = "output.txt"; // default
        boolean doubleCheck = false;
//		boolean permutationTest = false;

        //"/Users/megan/research/from Anne/phyMLwoPoly.txt" output.txt Pruned2
        //"/Users/megan/research/data/75_taxa_from_LSU_data_from_RecIDcm_example.nex.run1.t_for_alg_just_trees.txt" output.txt Pruned2
        //"/Users/megan/research/data/43_taxa_trees.txt" -a Pruned2 -d -n -o distances_43_Pruned2_normalized_April20.txt

        // last argument is always the infile, and the rest are optional

        if (args.length < 1) {
            TreeDistance.displayHelp();
            System.exit(0);
        }
        treeFile = args[args.length - 1];
        for (int i = 0; i < args.length - 1; i++) {

            if (!args[i].startsWith("-")) {
                System.out.println("Invalid command line option");
                displayHelp();
                System.exit(0);
            }

            if (args[i].equals("--verbose")) {
                TreeDistance.verbose = 1;
            } else if (args[i].equals("--help")) {
                displayHelp();
                System.exit(0);
            }
            //algorithm
            else if (args[i].equals("-a")) {
                if (i < args.length - 2) {
                    algorithm = args[i + 1];
                    i++;
                } else {
                    displayHelp();
                    System.exit(0);
                }
            }
            // output file
            else if (args[i].equals("-o")) {
                if (i < args.length - 2) {
                    outFile = args[i + 1];
                    i++;
                } else {
                    displayHelp();
                    System.exit(0);
                }
            }

            // all other arguments.  Note we can have -vn
            else {
                for (int j = 1; j < args[i].length(); j++) {

                    switch (args[i].charAt(j)) {
                        // doublecheck distances
                        case 'd':
                            doubleCheck = true;
                            break;

                        // display help
                        case 'h':
                            displayHelp();
                            System.exit(0);
                            break;

                        // normalize trees?
                        case 'n':
                            TreeDistance.normalize = true;
                            break;

			/*			else if (args[i].equals("-p")) {
							permutationTest = true;
							i++;
						} */

                        // unrooted trees?
                        case 'u':
                            TreeDistance.rooted = false;
                            break;

                        // verbose output
                        case 'v':
                            TreeDistance.verbose = 1;
                            break;

                        default:
                            System.out.println("Illegal command line option.\n");
                            displayHelp();
                            System.exit(0);
                            break;
                    } // end switch
                } // end for j
            } // end parsing an individual argument
        }  // end for i (looping through arguments)


/*		if (permutationTest) {
			PermutationTest.basicTest(treeFile, outFile, algorithm);
			System.exit(0);
		}*/

        try {
            computeAllInterTreeGeodesicsFromFile(treeFile, outFile, algorithm, doubleCheck);
        } catch (IOException e) {
            System.err.println("IO Error detected");
            System.exit(1);
        }

        System.exit(0);

    }

    // XXX: fix the following three methods!!!  should be able to combine or something...
    public static Vector<Bipartition> myVectorCloneBipartition(Vector<Bipartition> v) {
        if (v == null) {
            return null;
        }

        Vector<Bipartition> newV = new Vector<Bipartition>();
        for (int i = 0; i < v.size(); i++) {
            if (v.get(i) == null) {
                newV.add(null);
            } else {
                newV.add(v.get(i).clone());
            }
        }
        return newV;
    }
}

/** Recursive method for finding all the ratio sequences.
 *
 * @param m
 * @param ratioSeq
 */
/* public static void getMaxPathSpacesAsRatioSeqs(Vector<Bipartition> m, RatioSequence ratioSeq, Vector<PhyloTreeEdge> eEdges, Vector<PhyloTreeEdge> fEdges) {
	// returns vector containing each minimal e-set.
	Vector<Bipartition> minEls = getMinElements(m);

	if (minEls.size() == 0) {
		// base case: we are at the end of the chain, so add to final ratioSeq
		// first case to avoid memory overflow
		numMaxPaths++;
//		if (numMaxPaths%10== 0) {
			System.out.println("Path " + numMaxPaths);
//		}
		if (eEdges.size() > 20) {
			double dist = ratioSeq.getAscRSWithMinDist().getDistance();
			if (minTreeDist < 0 || dist < minTreeDist) {
//		System.out.println("New min dist of " + dist + " for rs " + ratioSeq);
				minTreeDist = dist;
				minTreeDistRatioSeq = ratioSeq.clone();
			}
		}
		else {
			finalRatioSeqs.add(ratioSeq.clone());
		}
		return;
	}



	// otherwise, call calculateFinalRatioSeqs for each min element in minEls
	Iterator<Bipartition> minElsIter = minEls.iterator();
	while (minElsIter.hasNext()) {
		Vector<Bipartition> newM = (Vector<Bipartition>) myVectorCloneBipartition(m);
//		RatioSequence newRatioSeq = (RatioSequence) ratioSeq.clone();
		Bipartition minEl = minElsIter.next();

		// create a new Ratio for the edges being dropped and added
		Ratio ratio = new Ratio();

		// set all edges in m equal to minEl to be null
		// and add to the ratio.
		for (int i = 0; i < newM.size();i++) {
			if ((newM.get(i) != null) && (newM.get(i).equals(minEl))) {
				// add the f edges with this minEl to the f edges of the ratio
				ratio.addFEdge(fEdges.get(i));
				newM.set(i, null);
			}
		}

		// set 1's in columns corresponding to removed e's to be 0,
		// and add that e split to the ratio.
		Iterator eEdgesIter = binaryToVector(minEl.getEdge()).iterator();
		while (eEdgesIter.hasNext()) {
			int eEdge = (Integer) eEdgesIter.next();
			newM = zeroCol(eEdge, newM);
			ratio.addEEdge(eEdges.get(eEdge));
		}

		// push ratio onto ratio seq and call calculateFinalRatioSeqs recursively
//		System.out.println("Adding ratio " + ratio + " with bipartition " + minEl + " to ratio sequence");
//		newRatioSeq.add(ratio);
//		getMaxPathSpacesAsRatioSeqs(newM,newRatioSeq, eEdges, fEdges);
		ratioSeq.add(ratio);

		// this path is greater than the min path so far, so don't follow any further.
		double pruneDist = ratioSeq.getAscRSWithMinDist().getDistance();
		if ( eEdges.size() > 20 && minTreeDist > 0 && pruneDist > minTreeDist) {
//			System.out.println("Pruning at distance " + pruneDist + " with minTreeDist " +  minTreeDist + " for rs: " + ratioSeq);
			ratioSeq.remove(ratio);  // need to go back to the original ratio sequence we started this method with
			return;
		}
		getMaxPathSpacesAsRatioSeqs(newM,ratioSeq, eEdges, fEdges);

		ratioSeq.remove(ratio);
	}
	minEls = null;
} */

/** Recursive method for finding all the ratio sequences.
 *
 * @param m
 * @param ratioSeq
 */
// saved sept 15 2007
/*public static void getMaxPathSpacesAsRatioSeqs(Vector<Bipartition> m, RatioSequence ratioSeq, Vector<PhyloTreeEdge> eEdges, Vector<PhyloTreeEdge> fEdges) {
	// returns vector containing each minimal e-set.
	Vector<Bipartition> minEls = getMinElements(m);
	Bipartition minEl;
	Ratio ratio;

	if (minEls.size() == 0) {
		// base case: we are at the end of the chain, so add to final ratioSeq
		// first case to avoid memory overflow
		numMaxPaths++;
//		if (numMaxPaths%10== 0) {
			System.out.println("Path " + numMaxPaths);
//		}
		if (eEdges.size() > 5) {
			double dist = ratioSeq.getAscRSWithMinDist().getDistance();
			if (minTreeDist < 0 || dist < minTreeDist) {
		System.out.println("New min dist of " + dist + " for rs " + ratioSeq);
				minTreeDist = dist;
				minTreeDistRatioSeq = ratioSeq.clone();
			}
		}
		else {
			finalRatioSeqs.add(ratioSeq.clone());
		}
		return;
	}

	Vector<Bipartition> sortedMinEls = new Vector<Bipartition>();
	Vector<Ratio> sortedMinElRatios = new Vector<Ratio>();

	// reset minRatio
//	minRatio = null;

	// for each min element, find the one with the smallest corresponding ratio
	Iterator<Bipartition> minElsIter = minEls.iterator();
	// last index of non-null entry in sorted array of min elements
	while (minElsIter.hasNext()) {
		minEl = minElsIter.next();

		ratio = calculateRatio(minEl, m, eEdges, fEdges );

		Boolean inserted = false;
		// store both in sortedMinEls in order of ascending ratios
		for (int k = 0; k < sortedMinElRatios.size();k++ ) {
			if (sortedMinElRatios.get(k).getRatio() > ratio.getRatio()) {
				sortedMinElRatios.add(k, ratio.clone());
				sortedMinEls.add(k,minEl.clone());
				inserted = true;
				break;
			}
		}
		if (!inserted) {
			// else bigger than all ratios currently in the vector
			sortedMinElRatios.add(ratio.clone());
			sortedMinEls.add(minEl.clone());
		}
//		System.out.println("sortedMinElRatios is " + sortedMinElRatios);
	}

	for (int j = 0; j < sortedMinElRatios.size(); j ++) {
		minEl = sortedMinEls.get(j);
		ratio = sortedMinElRatios.get(j);
//		System.out.println("minEl is " + minEl);
//		System.out.println("ratio is " + ratio);
		// check if we can just use the following line?
		//newM = removeMinElFrom(m, minMinEl);

		Vector<Bipartition> newM = (Vector<Bipartition>) myVectorCloneBipartition(m);
//		System.out.println("newM is " + newM);
//		RatioSequence newRatioSeq = (RatioSequence) ratioSeq.clone();
//		Bipartition minEl = minElsIter.next();

		// create a new Ratio for the edges being dropped and added
//		Ratio ratio = new Ratio();

		// set all edges in m equal to minEl to be null
		// and add to the ratio.
		for (int i = 0; i < newM.size();i++) {
			if ((newM.get(i) != null) && (newM.get(i).equals(minEl))) {
				// add the f edges with this minEl to the f edges of the ratio
//				ratio.addFEdge(fEdges.get(i));
				newM.set(i, null);
			}
		}

		// set 1's in columns corresponding to removed e's to be 0,
		// and add that e split to the ratio.
		Iterator eEdgesIter = binaryToVector(minEl.getEdge()).iterator();
		while (eEdgesIter.hasNext()) {
			int eEdge = (Integer) eEdgesIter.next();
//			newM = zeroCol(eEdge, newM);
			ratio.addEEdge(eEdges.get(eEdge));
		}

		// push ratio onto ratio seq and call calculateFinalRatioSeqs recursively
//		System.out.println("Adding ratio " + ratio + " with bipartition " + minEl + " to ratio sequence");
//		newRatioSeq.add(ratio);
//		getMaxPathSpacesAsRatioSeqs(newM,newRatioSeq, eEdges, fEdges);
		ratioSeq.add(ratio);

		// this path is greater than the min path so far, so don't follow any further.
		double pruneDist = ratioSeq.getAscRSWithMinDist().getDistance();
		if ( eEdges.size() > 5 && minTreeDist > 0 && pruneDist > minTreeDist) {
//			System.out.println("Pruning at distance " + pruneDist + " with minTreeDist " +  minTreeDist + " for rs: " + ratioSeq);
			ratioSeq.remove(ratio);  // need to go back to the original ratio sequence we started this method with
			numPrunes++;
			if (numPrunes%10 == 0 ) {
				System.out.println("prunes: " + numPrunes);
			}

			return;
		}
		getMaxPathSpacesAsRatioSeqs(newM,ratioSeq, eEdges, fEdges);

		ratioSeq.remove(ratio);

		if (numPrunes == 50000) {
			break;
		}
	}
	minEls = null;
} */