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
import polyAlg.*;

/* XXX: (general problems)
1) we need a method to convert trees without the same leaf2NumMap to having the same one, if possible.

*/

public class PhyloTree {
	private Vector<PhyloTreeEdge> edges;	 // each element of edges represents an interior split of the tree.
	private Vector<String> leaf2NumMap;   // the position the leaf name occurs in this vector corresponds to its coordinate in an split vector.
							// ie. if leaf A is in position 5, then, a 1 in an split vector at position 5 means that edges contains leaf A
							// standardized to always be in the natural order for strings
	private EdgeAttribute[] leafEdgeAttribs;
	private String newick;  // the newick format of the tree, if available

	// constructor
	// assume these trees are already rooted according to convention;
	// this shouldn't be a problem, because we root any tree read in from a string according to convention
	public PhyloTree(Vector<PhyloTreeEdge> edges, Vector<String> leaf2NumMap, EdgeAttribute[] leafEdgeLengths) {
		this.edges = edges;
		this.leaf2NumMap = leaf2NumMap;
		this.leafEdgeAttribs = leafEdgeLengths;
	}

	// constructor
	// assume these trees are already rooted according to convention;
	// this shouldn't be a problem, because we root any tree read in from a string according to convention
	public PhyloTree(Vector<PhyloTreeEdge> edges, Vector<String> leaf2NumMap) {
		this.edges = edges;
		this.leaf2NumMap = leaf2NumMap;
	}

	// copy constructor
	public PhyloTree(PhyloTree t) {
		this.edges = Tools.myVectorClonePhyloTreeEdge(t.edges);
		this.leaf2NumMap = Tools.myVectorCloneString(t.leaf2NumMap);
		if (t.leafEdgeAttribs != null) {
			this.leafEdgeAttribs = Arrays.copyOf(t.leafEdgeAttribs, t.leafEdgeAttribs.length);
		}
		if (t.newick != null) {
			this.newick = new String(t.newick);
		}
	}


	private void setLeaf2NumMapFromNewick() {
		// go through the string and pull out all the leaf labels:  any string between '(' and ':' or ',' and ':'
		int i =0;
		while (i < newick.length()) {
			// however, the first character might be the beginning of a leaf label
			if ((newick.charAt(i) == '(' || newick.charAt(i) == ',')&&(newick.charAt(i+1) != '(')) {
				leaf2NumMap.add(newick.substring(i+1, newick.substring(i).indexOf(":")+i));
				i = nextIndex(newick, i, ",)");
			}
			else {
				i++;
			}
		}
		// sort the elements of leaf2NumMap
		Collections.sort(leaf2NumMap);
	}


	/** Constructor - turns t into a PhyloTree, by storing the split data in edges and creating a conversion chart between the leaf names
	 * and their position in split vectors.  Basically each name is assigned a different number from 1 to # of leaves.
	 * Example of Newick Standard:  (B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);
	 * If tree is really unrooted, as passed in as argument and stored in TreeDistance class variable,
	 * then reroot the tree so that the last leaf in leaf2NumMap is the root (and then everything automatically becomes one leaf smaller)
	 * Handles multifurcating trees.
	 * Removes all interior edges with length 0.  XXX:  should be up to tolerance, I think
	 * TODO:  insert some kind of check that this situation (((A:1,B:1):1):1,C:1); can't happen
	 * TODO:  insert check for same vertex (so can't have two of the same vertices)
	 * @param t representation of a tree in Newick standard.
	 */
	public PhyloTree(String t, boolean rooted) {
		int leafNum = 0;

		LinkedList<PhyloTreeEdge> queue = new LinkedList<PhyloTreeEdge>();   // queue for keeping track of edges we are still adding leaves to
		edges = new Vector<PhyloTreeEdge>();
		leaf2NumMap = new Vector<String>();
		int endOfLeafEdgeLength = 0;

		// remove anything before the first (
		t = t.substring(t.indexOf('('));

		newick = t;


		//pull off ';' if at end
		if (t.charAt(t.length()-1) == ';') {
			t = t.substring(0, t.length() -1);
		}

		// pull off the first and last brackets (and a root length, between the last bracket and ;, if there is one.
		t = t.substring(t.indexOf('(') + 1);
		t = t.substring(0, t.lastIndexOf(')') );

		// Sanity check brackets by adding 1 for each left bracket encountered
		// and subtracting 1 for each right bracket encountered.
		// Bracket "total" should always be >= 0.
		int bracketCount =0;
		for (int i =0; i < t.length(); i++) {
			if (t.charAt(i) == '(') {
				bracketCount++;
			}
			else if (t.charAt(i) == ')') {
				bracketCount--;
			}
			if (bracketCount < 0) {
				System.err.println("Error reading in tree:  bracket problem in Newick string " + newick);
				System.exit(1);
			}
		}
		if (bracketCount != 0) {
			System.err.println("Error reading in tree: uneven number of brackets in Newick string " + newick);
			System.exit(1);
		}



try{  // for stringIndexOutOfBoundsException
		this.setLeaf2NumMapFromNewick();

		leafEdgeAttribs = new EdgeAttribute[leaf2NumMap.size()];

		// now go back, and associate leaves with edges.
		int i = 0;
		while (i < t.length() && i > -1) {
			switch(t.charAt(i)) {
			case '(':
				queue.addFirst(new PhyloTreeEdge());
				i++;
				break;

			case ')':
				// Extract this split's attribute, which is between the : following the ), and either a , or a ).
				queue.peek().setAttribute(new EdgeAttribute(t.substring(i+2, nextIndex(t, i, ",)") ) ) );
				edges.add(queue.poll());
				// increment i, so it is the next index after the length
				i = nextIndex(t, i, ",)");  // while loop will end if this is -1, and hence we are at the end of the string
				break;

			case ',':
				i++;
				break;

				// this character is the beginning of a leaf name
			default:
				// get following leaf label, which ends with a :
				leafNum = leaf2NumMap.indexOf(t.substring(i, t.substring(i).indexOf(":")+i));

				// the following two lines get the leaf split length for this leaf, and store it
				if ((t.substring(i+1).indexOf(',') > -1) && (t.substring(i+1).indexOf(')') > -1) ) {
					endOfLeafEdgeLength = Math.min( t.substring(i+1).indexOf(','), t.substring(i+1).indexOf(')') ) + i + 1;
				}
				else if (t.substring(i+1).indexOf(')') > -1) {
					endOfLeafEdgeLength = t.substring(i+1).indexOf(')') + i+1;
				}
				else if (t.substring(i+1).indexOf(',') > -1) {
					endOfLeafEdgeLength = t.substring(i+1).indexOf(',') + i+1;
				}
				else {
					// we have removed the end bracket, so we are not finding that
					endOfLeafEdgeLength = t.length();
				}

				leafEdgeAttribs[leafNum] = new EdgeAttribute( t.substring( (t.substring(i+1).indexOf(':') + i + 2),endOfLeafEdgeLength  ) );

				// we now want to add this leaf to all the edges in our queue.
				for (PhyloTreeEdge e : queue) {
					e.addOne(leafNum);
				}

				i = nextIndex(t, i, ",)");  // while loop will end if this is -1, and hence we are at the end of the string
				break;

			}  //switch

		}	// while

	} 	//try
	catch(StringIndexOutOfBoundsException e) {
		System.err.println("Error reading in tree:  invalid Newick string: (" + t + ");");
		System.exit(1);
	}

	// if tree is really unrooted, reroot so that the last leaf in leaf2NumMap is the root
	if (!rooted) {
		for (PhyloTreeEdge e : edges) {
			int lastLeaf = leaf2NumMap.size()-1;
			if (e.contains(lastLeaf)) {
				e.complement(lastLeaf +1);
			}
		}
	}

	for (int k = 0; k < edges.size(); k++) {
		edges.get(k).setOriginalEdge(edges.get(k).asSplit());
		edges.get(k).setOriginalID(k);
	}

	} // constructor


	// Getters and Setters
	public Vector<PhyloTreeEdge> getEdges() {
		return edges;
	}

	/** Returns the split at position i in the split vector
	 *
	 * @param i
	 * @return
	 */
	public PhyloTreeEdge getEdge(int i) {
		return edges.get(i);
	}

	public void setEdges(Vector<PhyloTreeEdge> edges) {
		this.edges = edges;
	}


	public Vector<String> getLeaf2NumMap() {
		return leaf2NumMap;
	}

	public EdgeAttribute getAttribOfSplit(Bipartition edge) {
		Iterator<PhyloTreeEdge> edgesIter = edges.iterator();
		while (edgesIter.hasNext()){
			PhyloTreeEdge e = (PhyloTreeEdge) edgesIter.next();
			if (e.sameBipartition(edge)) {
				return e.getAttribute();
			}
		}
		return null;
	}

	/** Returns a vector containing the splits (as Bipartitions) corresponding to the edges of this tree.
	 *
	 * @return
	 */
	public Vector<Bipartition> getSplits() {
		Vector<Bipartition> splits = new Vector<Bipartition>();

		for(int i = 0; i < edges.size(); i++) {
			splits.add(getEdge(i).asSplit() );
		}

		return splits;
	}


	/* * Normalizes so that the vector of all edges (both internal and ones ending in leaves) has length 1.
	 *
	 */
	public void normalize() {
		double vecLength = getDistanceFromOrigin();

		// divide by the length of the split length vector to normalize
		for (int i =0; i < leafEdgeAttribs.length; i++) {
			leafEdgeAttribs[i].scaleBy(1.0/vecLength);
		}
		for (int i = 0; i < edges.size(); i++) {
			edges.get(i).getAttribute().scaleBy(1.0/vecLength);
		}
	}

	public int numLeaves() {
		return leaf2NumMap.size();
	}

	public void setLeaf2NumMap(Vector<String> leaf2NumMap) {
		this.leaf2NumMap = leaf2NumMap;
	}


	/**  Returns the distance from this tree to the origin of tree space (i.e. includes leaf edges).
	 *
	 * @return
	 */
	public double getDistanceFromOrigin() {
		double dist = 0;
		for (int i = 0; i< edges.size(); i++) {
			dist = dist + Math.pow(edges.get(i).getLength(), 2);
		}

		for (int i = 0; i < leafEdgeAttribs.length; i++) {
			dist = dist + Math.pow(leafEdgeAttribs[i].norm(),2);
		}

		return Math.sqrt(dist);
	}

	/**  Returns the distance from this tree to the origin of tree space \ leaf space.  (i.e. doesn't include leaves)
	 *
	 * @return
	 */
	public double getDistanceFromOriginNoLeaves() {
		double dist = 0;
		for (int i = 0; i< edges.size(); i++) {
			dist = dist + Math.pow(edges.get(i).getLength(), 2);
		}

		return Math.sqrt(dist);
	}



	/** Returns a vector containing edges with the same partition in t1 and t2.
	 *  If t1 and t2 do not have the same leaf2NumMap, then returns null.
	 *  If t1 and t2 have no common edges, then returns a vector of size 0.
	 *  We set the length of the returned edges to be the absolute value of the difference between the length
	 *  of the split in t1 and the length of the split in t2.
	 *  Doesn't assume t1 and t2 have same # of edges.
	 *  Edges from one tree that are compatible with all edges of the other tree are returned as common edges.
	 *  XXX:  can return 0 length common edges
	 * @param t1
	 * @param t2
	 * @return
	 */
	public static Vector<PhyloTreeEdge> getCommonEdges(PhyloTree t1, PhyloTree t2) {
		Vector<PhyloTreeEdge> commonEdges = new Vector<PhyloTreeEdge>();

		// if the two trees do not have the same leaf2NumMap
		if (!(t1.getLeaf2NumMap().equals(t2.getLeaf2NumMap()))){
			System.out.println("Error: the two trees do not have the same leaves!");
			System.out.println("First tree's leaves are " + t1.getLeaf2NumMap() );
			System.out.println("Second tree's leaves are " + t2.getLeaf2NumMap() );
			System.exit(1);
		}

		for (PhyloTreeEdge e1 : t1.edges) {
			// check that the split is not 0 (has no leaves and can be thought of as a null split; ignore such edges)
			if (!e1.isZero()) {
				if (t2.getSplits().contains(e1.asSplit() ) ){
					// then we have found the same split in each tree
					EdgeAttribute commonAttrib = EdgeAttribute.difference(e1.getAttribute(), t2.getAttribOfSplit(e1.asSplit()));
					commonEdges.add(new PhyloTreeEdge(e1.asSplit(), commonAttrib,e1.getOriginalID() ));
				}
				// otherwise check if the split is compatible with all splits in t2
				else if (e1.isCompatibleWith(t2.getSplits()) ) {
					EdgeAttribute commonAttrib = EdgeAttribute.difference(e1.getAttribute(), null);
					commonEdges.add(new PhyloTreeEdge(e1.asSplit(),commonAttrib,e1.getOriginalID() ));
				}
			}
		}
		// check for splits in t2 that are compatible with all splits in t1
		for (PhyloTreeEdge e2 : t2.getEdges()) {
			if (!e2.isZero()) {
				if (e2.isCompatibleWith(t1.getSplits()) && !(t1.getSplits().contains(e2.asSplit()))) {
					EdgeAttribute commonAttrib = EdgeAttribute.difference(null,e2.getAttribute());
					commonEdges.add(new PhyloTreeEdge(e2.asSplit(),commonAttrib,e2.getOriginalID() ));
				}
			}
		}
		return commonEdges;
	}

	/**  Returns the edges that are not in common with edges in t, excluding edges of length 0.
	 *
	 * @param t
	 * @return
	 */
	public Vector<PhyloTreeEdge> getEdgesNotInCommonWith(PhyloTree t) {
		Vector<PhyloTreeEdge> notCommonEdges = new Vector<PhyloTreeEdge>();

		// if the two trees do not have the same leaf2NumMap
		if (!(this.getLeaf2NumMap().equals(t.getLeaf2NumMap()))){
			System.out.println("Error: the two trees do not have the same leaves!");
			System.out.println("First tree's leaves are " + this.getLeaf2NumMap() );
			System.out.println("Second tree's leaves are " + t.getLeaf2NumMap() );
			System.exit(1);
		}

		for (int i = 0; i < this.edges.size(); i++) {
			Boolean common = false;
			for (int j = 0; j < t.edges.size(); j++) {

				// check that the split is not 0 (has no leaves and can be thought of as a null split; ignore such edges)
				if ( !this.getEdge(i).isZero() && !t.getEdge(j).isZero() && (this.getEdge(i).sameBipartition(t.getEdge(j))) ) {
					common = true;
					break;
				}
			}

			if (!common) {
				notCommonEdges.add(this.getEdge(i).clone());
			}
		}


		return notCommonEdges;
	}



	public void permuteLeaves() {
		Vector<PhyloTreeEdge> v = this.edges;
		int numEdges = v.size();
		int numLeaves = numEdges + 2;  // because assuming trees are rooted, or have been rooted

//		System.out.println("numEdges: " + numEdges + "; v is " + v);

		Vector<PhyloTreeEdge> permutedV = new Vector<PhyloTreeEdge>();
		Random r = new Random();

		// set all the original edges and lengths
		for (int i =0; i < numEdges; i++) {
			permutedV.add(new PhyloTreeEdge(v.get(i).getAttribute() ) );
			permutedV.get(i).setOriginalEdge(permutedV.get(i).asSplit() );
		}

		// initialize the vector to keep track of which leaves have already been permuted.
		Vector<Integer> leavesToBePermuted = new Vector<Integer>();
		for (int i = 0; i < numLeaves; i++ ) {
			leavesToBePermuted.add(new Integer(i));
		}

		// each time through the loop, one leaf is permuted.
		for (int i = numLeaves; i > 0; i--) {
			int leaf = leavesToBePermuted.remove(r.nextInt(i)).intValue();
//			System.out.println("Moving leaf " + leaf + " to " + (i-1));

			// make the v.size() - i column in permutedV be col in v
			for (int j = 0; j < numEdges; j++ ) {
				// if the bit in the i-th column of the j-th row is 1 ...
				if (v.get(j).getPartition().get(leaf) ) {
					permutedV.get(j).getPartition().set(i-1);
				}
			}
		}
//		System.out.println("permutedV is " + permutedV);
		edges = TreeDistance.myVectorClonePhyloTreeEdge(permutedV);
	}

	public void permuteLeaves(int[] permutation) {
		Vector<PhyloTreeEdge> v = this.edges;
		int numEdges = v.size();

		Vector<PhyloTreeEdge> permutedV = new Vector<PhyloTreeEdge>();
		// set all the original edges and lengths
		for (int i =0; i < numEdges; i++) {
			permutedV.add(new PhyloTreeEdge(v.get(i).getAttribute(), v.get(i).getOriginalID() ));
		}

		for (int j = 0; j < numEdges; j++ ) {
			for (int i = 0; i < permutation.length; i++) {
				//  changing leaf i to be in position permutation(i)
				if(v.get(j).getPartition().get(i)) {
					permutedV.get(j).getPartition().set(permutation[i]);
				}
			}

		}
//		System.out.println("permutedV is " + permutedV);
		edges = TreeDistance.myVectorClonePhyloTreeEdge(permutedV);
	}

	public void swapleaves(int l1,int l2) {
		Vector<PhyloTreeEdge> v = this.edges;
		int numEdges = v.size();

		int l1swap;   // saves the value of position l1 when swapping leaves

		for (int i = 0; i < numEdges; i++) {
		//  changing leaf l1 to be in position l2
			if(v.get(i).getPartition().get(l1)) {
				l1swap = 1;
			}
			else {
				l1swap = 0;
			}
			if (v.get(i).getPartition().get(l2) ) {
				v.get(i).getPartition().set(l1);
			}
			else {
				v.get(i).getPartition().clear(l1);
			}
			if (l1swap == 1) {
				v.get(i).getPartition().set(l2);
			}
			else {
				v.get(i).getPartition().clear(l2);
			}
		}
	}



	/** Assigns each interior split a random length between min and max.
	 *
	 * @param min
	 * @param max
	 */
//	public void randomIntEdgeLengths(double min, double max) {
//
//
//		for (int i = 0; i < edges.size();i++) {
//			edges.get(i).setLength(Math.random()*(max - min) + min);
//		}
//	}


	/** Deletes the split corresponding to bipartition e from the tree.
	 *
	 * @param e
	 */
	public boolean removeSplit(Bipartition e){
		boolean removed = false;

		for(int i = 0; i < edges.size(); i++) {
			if (edges.get(i).sameBipartition(e)) {
				edges.remove(i);
				removed = true;
				break;
			}
		}
		return removed;
	}

	/** Remove the edges with bipartitions corresponding to those in Vector<Bipartitions>
	 *
	 */
	public void removeSplits(Vector<Bipartition> splits) {
		for (int i = 0; i < splits.size(); i++) {
			this.removeSplit(splits.get(i));
		}
	}

	/**  Removes the edges corresponding to the ones in v.
	 *
	 * @param v
	 */
	public void removeEdgesIndicatedByOnes(Bipartition v) {
		int numRemovedSoFar = 0;
		for (int i = 0; i < leaf2NumMap.size() -2; i++) {
			if (v.contains(i)) {
				edges.remove(i - numRemovedSoFar);
				numRemovedSoFar++;
			}
		}
	}

	/** Adds the split e to the tree
	 *
	 * @param e
	 */
	public void addEdge(PhyloTreeEdge e) {
		edges.add(e);
	}

	/** Returns the smallest index in a string t after position i, at which one of the characters in s is located.
	 *
	 * @param i
	 * @param s
	 * @return t.length() if no character in s is located after position i in t; otherwise return the min index
	 */
	public static int nextIndex(String t, int i, String s) {
		int minIndex = t.length();
		int tempIndex = -1;

		for (int j = 0; j < s.length(); j++) {
			tempIndex = t.substring(i+1).indexOf(s.charAt(j));
			if ((tempIndex != -1)  && (tempIndex + i + 1 < minIndex)) {
				minIndex = tempIndex + i +1;
			}
		}
		return minIndex;
	}


	/** Returns a vector containing representations of 0-1 vectors which
	 * have a 1 in coordinate i if split i of tree t crosses the split of this tree, which this 0-1
	 * vector represents.  The entries of this vector correspond to the edges of this tree.
	 * Doesn't assume the trees have the same number of edges.
	 * XXX maybe should move exiting if trees don't have same leaf2numMap to this method??
	 * @param t
	 * @return
	 */
	public Vector<Bipartition> getCrossingsWith(PhyloTree t) {
		// if the two trees do not have the same leaf2NumMap
		if (!(this.getLeaf2NumMap().equals(t.getLeaf2NumMap()))){
			return null;
		}

/*		Vector<Bipartition> edges = new Vector<Bipartition>();
		Iterator<PhyloTreeEdge> thisIter = this.edges.iterator();
		int i = 0;  // keeps track of which vector coordinate the split
		// that we are comparing with will correspond to

		while (thisIter.hasNext()) {
			PhyloTreeEdge thisTreeEdge = thisIter.next();
			Bipartition newEdge = new Bipartition();

			Iterator<PhyloTreeEdge> tIter = t.edges.iterator();
			i = 0;  // reset i
			while (tIter.hasNext()) {
				PhyloTreeEdge e = (PhyloTreeEdge) tIter.next();
				if (e.getEdge() != 0 && e.crosses(thisTreeEdge)) {
					newEdge.addOne(i);
				}
				i++;
			}
			edges.add(newEdge);
		}
		*/
//		System.out.println("in getCrossings, this tree is " + this);
		Vector<Bipartition> edges = new Vector<Bipartition>();
		// for each vector in m
//		System.out.println("this.edges.size() is " + this.edges.size());
		for (int i = 0; i < this.edges.size(); i++) {

			PhyloTreeEdge thisTreeEdge = this.edges.get(i);;
			Bipartition newSplit = new Bipartition();

			// add the appropriate 1's
			for (int j = 0; j < t.edges.size(); j++) {
				PhyloTreeEdge e = t.edges.get(j);
				// e should not equal 0.
				if (e.crosses(thisTreeEdge)) {
					newSplit.addOne(j);
				}
			}
			edges.add((Bipartition) newSplit.clone());
//			System.out.println("Added split " + newEdge.clone());
		}
//		System.out.println("in getCrossingsWith returning " + edges);
		return edges;
	}

	public String toString() {
		return "Leaves: " + leaf2NumMap + "; edges: " + edges + "; leaf edges: " + Arrays.toString(leafEdgeAttribs);
	}


	// TODO:  not actually overriding the object clone method.  Also, clone should not call constructors.
	public PhyloTree clone() {
		return new PhyloTree(TreeDistance.myVectorClonePhyloTreeEdge(edges), TreeDistance.myVectorCloneString(leaf2NumMap), leafEdgeAttribs.clone() );
	}

	// TODO:  currently does not check equality of internal edges
	@Override public boolean equals(Object t) {
		if (t == null) {
			return false;
		}
		if (this == t) {
			return true;
		}

		if (!(t instanceof PhyloTree)) {
			return false;
		}

		//  compare the two vectors of edges without regard to order by making sure they
		//  are each contained in each other.
		return leaf2NumMap.equals( ((PhyloTree) t).leaf2NumMap) &&
			Arrays.equals(leafEdgeAttribs, ((PhyloTree)t).leafEdgeAttribs ) &&
				edges.containsAll( ((PhyloTree) t).edges) &&
					((PhyloTree) t).edges.containsAll(edges);
	}

	/** Returns true if the two trees are approximately equal.
	 *  Current definition of approximately equal is that
	 *  they are within epsilon geodesic distance of each other.
	 *
	 * @param t
	 * @return
	 */
	public boolean approxEquals(PhyloTree t, Double epsilon) {
		if (this == null || t == null) {
			return false;
		}
		return (PolyMain.getGeodesic(this,t,null).getDist() < epsilon);
	}


	// XXX for now this is the testing procedure
	public static void main(String args[]) {
		PhyloTree t1 = new PhyloTree("((A:0.1,B:0.2):1,(C:0.3,D:0.4):2)", true);
		// should produce:  [1.0 11, 2.0 1100]
		System.out.println("Edges of tree 1 (length, followed by leaves as 1's): " + t1.getEdges().toString());

		PhyloTree t2 = new PhyloTree("(((B:0.1,A:0.2):1,C:0.3):2,(D:0.4,E:0.5):3)", true);
		// should produce: [1.0 11, 2.0 111, 3.0 11000]
		System.out.println("Edges of tree 2: " + t2.getEdges().toString());

		// has the same leaf2NumMap as t2 and a common split
		PhyloTree t3 = new PhyloTree("((A:0.1,B:0.2):1,(C:0.3,(D:0.4,E:0.5):2):3)", true);
		// should produce [1.0 11, 2.0 11100, 3.0 11000]
		System.out.println("Edges of tree 3: " + t3.getEdges().toString());

		System.out.println("common edges of t2 and t3: " + getCommonEdges(t2,t3));

		// should have no common edges with t2
		PhyloTree t4 = new PhyloTree("((A:0,C:0):1,(D:0,(B:0,E:0):2):3)", true);
		System.out.println("Edges of tree 4:" + t4.getEdges().toString());
		System.out.println("1s indicate edges in tree 4 crossed by that split in tree2: " + t2.getCrossingsWith(t4));
		System.out.println("1s indicate edges in tree 1 crossed by that split in tree 4: " + t1.getCrossingsWith(t4));

		PhyloTree triTop = new PhyloTree("(((A:0,B:0):1,(C:0,D:0):2):3,E:0)",true);
		System.out.println("Edges of tree TriTop: " +triTop.getEdges());
		PhyloTree triBot = new PhyloTree("(((D:0,E:0):1,(C:0,B:0):2):3,A:0)",true);
		System.out.println("1s indicate edgdes in tree triBot crossed by that split in tree triTop:" + triBot.getCrossingsWith(triTop));
		PhyloTree noCommonSortingWithTriTop = new PhyloTree("(((A:0,C:0):1,E:0):2,(B:0,D:0):3)", true);
		System.out.println("Edges of tree noCommonSortingWithTriTop: " + noCommonSortingWithTriTop.getEdges());
		System.out.println("1s indicate edges in tree triTop crossed by that split in tree noCommonSortingWithTriTop: " +
				noCommonSortingWithTriTop.getCrossingsWith(triTop));
	}

	public EdgeAttribute[] getLeafEdgeAttribs() {
		return leafEdgeAttribs;
	}

	public EdgeAttribute[] getCopyLeafEdgeAttribs() {
		if (leafEdgeAttribs == null) {
			return null;
		}

		EdgeAttribute[] copy = new EdgeAttribute[leafEdgeAttribs.length];

		for (int i = 0; i < leafEdgeAttribs.length; i++) {
			copy[i] = leafEdgeAttribs[i].clone();
		}

		return copy;
	}

	/** Returns vector with the weights of the edges in the Vector of interior edges.
	 *
	 * @return
	 */
	public double[] getIntEdgeAttribNorms() {
		double[] norms = new double[edges.size()];

		for (int i = 0; i < edges.size(); i++) {
			norms[i] = edges.get(i).getAttribute().norm();
		}
		return norms;
	}

	public void setLeafEdgeAttribs(EdgeAttribute[] leafEdgeAttribs) {
		this.leafEdgeAttribs = leafEdgeAttribs;
	}

	public String getNewick(Boolean branchLengths) {
		String newNewick = "";
		Vector<String> strPieces = new Vector<String>();	// stores the pieces of the Newick format constructed so far
		Vector<PhyloTreeEdge> corrEdges = new Vector<PhyloTreeEdge>();	// stores the top split for each piece of Newick constructed so far

		if (newick != null) {
			return newick;
		}
		// otherwise we have to figure out the Newick representation for this tree
		if (edges.size() == 0) {
			// we have the star tree
			newNewick = "(";
			for (int i = 0;i < leaf2NumMap.size(); i++) {
				if (i >0 ) {
					newNewick = newNewick + ",";
				}
				newNewick = newNewick + leaf2NumMap.get(i);
				if (branchLengths) {
					newNewick = newNewick + ":" + leafEdgeAttribs[i];
				}
			}
			newNewick = newNewick  + ");";
			return newNewick;
		}

		Vector<PhyloTreeEdge> edgesLeft = TreeDistance.myVectorClonePhyloTreeEdge(edges);

		while (edgesLeft.size() > 0) {

			// pick the next minimal split to add
			PhyloTreeEdge minEdge = edgesLeft.get(0);

			for (int i = 1; i < edgesLeft.size(); i++) {
				if (minEdge.contains(edgesLeft.get(i))) {
					// minEdge was not the min elements
					minEdge = edgesLeft.get(i);
				}
			}
			// remove minEdge from edgesLeft
			edgesLeft.remove(minEdge);

			corrEdges.add(0,minEdge.clone());
			strPieces.add(0,"");

			// now we have a min element.
			// Start its Newick string.
			String str1 = "(";
			// Find out if it contains one of the min elements we have already processed.
			// Since we are allowing degenerate trees, there could be an arbitrary number of such edges.

			int k = 1;
			while (k < corrEdges.size()) {
				if (minEdge.contains(corrEdges.get(k))) {
					// add it to the string
					str1 = str1 + strPieces.get(k) + ",";

					// remove each leaf in this split from minEdge
					minEdge.getPartition().andNot(corrEdges.get(k).getPartition());

					// remove this split and its corresponding string from the vectors
					strPieces.remove(k);
					corrEdges.remove(k);
				}
				else {
					k++;
				}
//				System.out.println("in here");
			}

			// add all the elements still in minEdge (These are leaves that weren't already added as part of
			// a min split contained by minEdge.)
			if (!minEdge.getPartition().isEmpty()) {
				for (int i = 0; i < minEdge.getPartition().length(); i++) {
					if (minEdge.getPartition().get(i)) {
						str1 = str1 + leaf2NumMap.get(i);
						if (branchLengths) {
//							str1 = str1 + ":" + d6o.format(leafEdgeAttribs[i]);
							str1 = str1 + ":" + leafEdgeAttribs[i];
						}
						str1 = str1 +  ",";
					}
				}
			}
			// remove the last , and add the bracket and minEdge length
			str1 = str1.substring(0, str1.length()-1) + ")";
			if (branchLengths) {
//				str1 = str1 + ":" + d6o.format(minEdge.getLength());
				str1 = str1 + ":" + minEdge.getAttribute();
			}


			// store str1
			strPieces.set(0,new String(str1));
		}

		// now we need to combine all edges in corrEdges and all remaining leaves
		BitSet allLeaves = new BitSet();
		allLeaves.set(0,leaf2NumMap.size(),true);

		String newickString = "(";
		// add all the string pieces we've accumulated
		for (int i = 0; i < corrEdges.size(); i++) {
			newickString = newickString + strPieces.get(i) + ",";
			allLeaves.andNot(corrEdges.get(i).getPartition());
		}
		// add all remaining leaves
		if (!allLeaves.isEmpty()) {
			for (int i = 0; i < allLeaves.length(); i++) {
				if (allLeaves.get(i)) {
					newickString = newickString + leaf2NumMap.get(i);
					if (branchLengths) {
						newickString = newickString + ":" + leafEdgeAttribs[i];
//						newickString = newickString + ":" + d6o.format(leafEdgeAttribs[i]);
					}
					newickString = newickString + ",";
				}
			}
		}

		// remove the last ,
		newickString = newickString.substring(0, newickString.length()-1) + ");";

		return newickString;
	}

	public void setNewick(String newick) {
		this.newick = newick;
	}

	/** Returns number of edges in this tree, including any 0 length ones.
	 *
	 * @return
	 */
	public int numEdges() {
		return edges.size();
	}

	/**  Projects this tree onto the geodesic between trees t1 and t2;
	 * i.e. returns the tree on the geodesic that is closest (= has min geo distance) to this tree.
	 * Algorithm ends when distance between left and right is < epsilon
	 * @param geo
	 * @return
	 */
	public PhyloTree projectToGeo(PhyloTree t1, PhyloTree t2, double epsilon) {
		double rightIndex = 1;
		double leftIndex = 0;
		double middleIndex = 0.5;

		// test that the projection is not one of the ends of the geodesic
		// check left end
		if (PolyMain.getGeodesic(this,Geodesic.getTreeAt(t1,t2,0),null).getDist() < PolyMain.getGeodesic(this,Geodesic.getTreeAt(t1,t2,epsilon),null).getDist()) {
			return Geodesic.getTreeAt(t1,t2,0);
		}
		else if (PolyMain.getGeodesic(this,Geodesic.getTreeAt(t1,t2,1),null).getDist() < PolyMain.getGeodesic(this,Geodesic.getTreeAt(t1,t2,1-epsilon),null).getDist()) {
			return Geodesic.getTreeAt(t1,t2,1);
		}

		while ((rightIndex - leftIndex) > epsilon ) {
			double distToRight = PolyMain.getGeodesic(this,Geodesic.getTreeAt(t1,t2,rightIndex),null).getDist();
			double distToLeft = PolyMain.getGeodesic(this,Geodesic.getTreeAt(t1,t2,leftIndex),null).getDist();
			double distToMiddle = PolyMain.getGeodesic(this,Geodesic.getTreeAt(t1, t2, middleIndex),null).getDist();

//			System.out.println("Left index is " + leftIndex + "; distance to left tree is " + distToLeft);
//			System.out.println("Middle index is " + middleIndex + "; distance to middle tree is " + distToMiddle);
//			System.out.println("Right index is " + rightIndex + "; distance to right tree is " + distToRight);

			if ((distToLeft < distToMiddle) && (distToMiddle < distToRight)) {
				rightIndex = middleIndex;
				middleIndex = (rightIndex - leftIndex)/2 + leftIndex;
//				System.out.println("left < middle < right; changed to left = " + leftIndex + "; middle = " + middleIndex + "; right = " + rightIndex);
			}
			else if ((distToLeft > distToMiddle) && (distToMiddle > distToRight)) {
				leftIndex = middleIndex;
				middleIndex = (rightIndex - leftIndex)/2 + leftIndex;
//				System.out.println("left > middle > right; changed to left = " + leftIndex + "; middle = " + middleIndex + "; right = " + rightIndex);
			}
			else if ( ((distToMiddle < distToRight) && (distToRight < distToLeft)) || (( distToMiddle < distToLeft) && (distToLeft < distToRight)) ) {

				double distToTest, testIndex;

				// do a test on the bigger interval
				if ( (rightIndex - middleIndex) < (middleIndex - leftIndex) ) {
					// test between left and middle
					testIndex = (middleIndex - leftIndex)/2 + leftIndex;

					distToTest = PolyMain.getGeodesic(this,Geodesic.getTreeAt(t1,t2,testIndex),null).getDist();
					if (distToTest > distToMiddle) {
						leftIndex = testIndex;
					}
					else {  // distToTest < distToMiddle < distToRight < distToLeft
						rightIndex = middleIndex;
						middleIndex = testIndex;
					}
				}
				else { // test between middle and right
					testIndex = (rightIndex - middleIndex)/2 + middleIndex;
					distToTest = PolyMain.getGeodesic(this,Geodesic.getTreeAt(t1,t2,testIndex),null).getDist();

					if (distToMiddle < distToTest) {
						rightIndex = testIndex;
					}
					else { // distToTest < distToMiddle < distToRight < distToLeft
						leftIndex = middleIndex;
						middleIndex = testIndex;
					}
				}
			}
			else {
				// there is a problem
				System.err.println("Error projecting tree onto geodesic: illegal ordering in line search");
				System.exit(1);
			}
		} // end while
		middleIndex = (rightIndex - leftIndex)/2 + leftIndex;
		return Geodesic.getTreeAt(t1,t2,middleIndex);
	}
}
