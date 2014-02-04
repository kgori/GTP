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
import polyAlg.Tools;

/** A Ratio object represents some number of edges, which have been combined in the ratio sequence.
 * Made so you need the e and f lengths, but not necessary to associate edges with them.
 * If there are edges, the lengths must agree with them.
*/

//TODO:  TEST THIS CLASS!!!
// XXX: for now, everything from one tree labelled e, everything from the other tree labelled f
public class Ratio {
	private double eLength;  // always updated to correspond to the eEdges
	private double fLength;  // always updated to correspond to the fEdges
	private Vector<PhyloTreeEdge> eEdges;  // holds the numbers of the edges from the e-tree involved in this ratio
	private Vector<PhyloTreeEdge> fEdges;

	public Ratio() {
		eLength = 0;
		fLength = 0;
		eEdges = new Vector<PhyloTreeEdge>();
		fEdges = new Vector<PhyloTreeEdge>();
	}

	public Ratio(Vector<PhyloTreeEdge> eEdges, Vector<PhyloTreeEdge> fEdges ) {
		this.eLength = geoAvg(eEdges);
		this.fLength = geoAvg(fEdges);
		this.eEdges = eEdges;
		this.fEdges = fEdges;
	}

	public Ratio(double eLength, double fLength) {
		this.eLength = eLength;
		this.fLength = fLength;
		eEdges = new Vector<PhyloTreeEdge>();
		fEdges = new Vector<PhyloTreeEdge>();
	}

	public Vector<PhyloTreeEdge> getEEdges() {
		return eEdges;
	}


	public void addEEdge(PhyloTreeEdge edge) {
		eEdges.add(edge);
		eLength = geoAvg(eLength,edge.getLength());
	}

	public void addAllEEdges(Vector<PhyloTreeEdge> edges) {
		eEdges.addAll(edges);
		eLength = geoAvg(eEdges);
	}

	public double getELength() {
//		return eLength;
		return geoAvg(eEdges);
	}

	/** Only sets e length if no e edges are stored in the ratio.
	 * Otherwise does nothing.
	 *
	 * @param eLen
	 */
	public void setELength(double eLen) {
		if (eEdges.size() == 0) {
			eLength = eLen;
		}
	}

	/** Only sets f length if no f edges are stored in the ratio.
	 * Otherwise does nothing.
	 *
	 * @param fLen
	 */
	public void setFLength(double fLen) {
		if (fEdges.size() == 0) {
			fLength = fLen;
		}
	}


	public Vector<PhyloTreeEdge> getFEdges() {
		return fEdges;
	}

	public void addFEdge(PhyloTreeEdge edge) {
		fEdges.add(edge);
		fLength = geoAvg(fLength,edge.getLength());
	}

	public void addAllFEdges(Vector<PhyloTreeEdge> edges) {
		fEdges.addAll(edges);
		fLength = geoAvg(fEdges);
	}

	public double getFLength() {
//		return fLength;
		return geoAvg(fEdges);
	}


	/** Returns the actual ratio of the (combined) length of the e edges over the (combined) length of the f edges.
	 *
	 * @return
	 */
	public double getRatio() {
		return eLength/fLength;
	}

	/** Returns the "time" (between 0 and 1) that we cross the orthant boundary associated with this ratio.
	 *
	 * @return
	 */
	public double getTime() {
		return eLength/(eLength + fLength);
	}


	/**  Combines the two ratios.  If neither ratio has edges associated with it, just combine the lengths.
	 *
	 * @param r1
	 * @param r2
	 * @return
	 */
	public static Ratio combine(Ratio r1, Ratio r2) {
		Ratio r = new Ratio();
		if ((r1.eEdges.size() == 0) && (r2.eEdges.size() == 0)){
			r.setELength(geoAvg(r1.eLength, r2.eLength));
		}
		else {
			// changed from r.addAllEEdges( (Vector<PhyloTreeEdge>) r2.eEdges.clone() ); etc.  June 15; remove comment when code works after this date
			r.addAllEEdges( TreeDistance.myVectorClonePhyloTreeEdge(r1.eEdges) );
			r.addAllEEdges( TreeDistance.myVectorClonePhyloTreeEdge( r2.eEdges ) );
		}

		if ((r1.fEdges.size() == 0) && (r2.fEdges.size() == 0)){
			r.setFLength(geoAvg(r1.fLength, r2.fLength));
		}
		else {
			r.addAllFEdges( TreeDistance.myVectorClonePhyloTreeEdge( r1.fEdges ) );
			r.addAllFEdges( TreeDistance.myVectorClonePhyloTreeEdge( r2.fEdges ) );
		}

		return r;
	}

	/** Computes the geometric average of d1 and d2.
	 *  ie. returns sqrt(d1^2 + d2^2)
	 * @param d1
	 * @param d2
	 * @return
	 */
	public static double geoAvg(double d1, double d2) {
		return Math.sqrt(Math.pow(d1,2) + Math.pow(d2,2));
	}

	public static double geoAvg(Vector<PhyloTreeEdge> edges) {
		double gAvg = 0;

		for (PhyloTreeEdge e : edges) {
			gAvg = gAvg + Math.pow(e.getLength(),2);
		}

		return Math.sqrt(gAvg);
	}

	/** Returns the ratio with everything e switched with everything f.
	 *
	 * @return
	 */
	public Ratio reverse() {
		Ratio ratio = new Ratio();
		ratio.addAllEEdges(fEdges);
		ratio.addAllFEdges(eEdges);
		return ratio;
	}


	public Boolean containsOriginalEEdge(Bipartition edge) {
		PhyloTreeEdge ratioEdge;

		for(int i = 0; i < eEdges.size(); i++)  {
			ratioEdge = eEdges.get(i);
			if (ratioEdge.getOriginalEdge().equals(edge)) {
				return true;
			}
		}
		for(int i = 0; i < fEdges.size(); i++)  {
			ratioEdge = fEdges.get(i);
			if (ratioEdge.getOriginalEdge().equals(edge)) {
				return true;
			}
		}
		return false;
	}

/*	public Boolean containsOriginalEdge(long split) {
		PhyloTreeEdge ratioEdge;

		for(int i = 0; i < eEdges.size(); i++)  {
			ratioEdge = eEdges.get(i);
			if (ratioEdge.getOriginalEdge().equals(new Bipartition(split))) {
				return true;
			}
		}
		for(int i = 0; i < fEdges.size(); i++)  {
			ratioEdge = fEdges.get(i);
			if (ratioEdge.getOriginalEdge().equals(new Bipartition(split))) {
				return true;
			}
		}
		return false;
	} */

	/**
	 * XXX: Note that the toString representation of the Ratio is not unique - the edges could
	 * be ordered differently in the Vector, and hence represented differently as a String.
	 */
	public String toString() {
		DecimalFormat d4o = new DecimalFormat("#0.####");
		return eEdges.toString() + " " + d4o.format(eLength) + "/" + d4o.format(fLength) + " " + fEdges.toString();
	}

	/** Returns the edges, numbered by originalID, in order (?).
	 * {1,2}/{4,1}
	 *
	 * @return
	 */
	public String toStringCombTypeAndValue() {
		DecimalFormat d3o = new DecimalFormat("#0.###");
		String s = "{";

		for (int i = 0; i < eEdges.size(); i++) {
			s = s + eEdges.get(i).getOriginalID();
			if (i < eEdges.size() -1) {
				s = s + ",";
			}
		}
		s = s + "}/{";

		for (int i = 0; i < fEdges.size(); i++) {
			s = s + fEdges.get(i).getOriginalID();
			if (i < fEdges.size() -1) {
				s = s + ",";
			}
		}
		s = s + "} " + d3o.format(eLength) + "/" + d3o.format(fLength) + "=" + d3o.format(getRatio());
		return s;
	}

	/** Returns the edges, numbered by originalID, in order (?).
	 *
	 * @return
	 */
	public String toStringCombType() {
		String s = "{";

		for (int i = 0; i < eEdges.size(); i++) {
			s = s + eEdges.get(i).getOriginalID();
			if (i < eEdges.size() -1) {
				s = s + ",";
			}
		}
		s = s + "}/{";

		for (int i = 0; i < fEdges.size(); i++) {
			s = s + fEdges.get(i).getOriginalID();
			if (i < fEdges.size() -1) {
				s = s + ",";
			}
		}
		s = s + "}";
		return s;
	}


	/** Returns the ratio as just its value.
	 *
	 * @return
	 */
	public String toStringJustValue() {
		return "" + getRatio();
	}

	public String toStringVerbose(Vector<String> leaf2NumMap) {
		DecimalFormat d6o = new DecimalFormat("#0.########");
		String s = "" + d6o.format(getRatio()) + "\nTotal length and corresponding edges dropped:\n";

		s = s + d6o.format(eLength)  + "\t";

		// list the edges dropped
		for (int i = 0; i<eEdges.size(); i++) {
			if (i == 0) {   // nice formatting
				s = s + Bipartition.toStringVerbose(eEdges.get(i).getPartition(),leaf2NumMap) + "\n";
			}
			else {
				s = s + "\t\t" + Bipartition.toStringVerbose(eEdges.get(i).getPartition(),leaf2NumMap) + "\n";
			}
		}

		s = s + "\nTotal length and corresponding edges added:\n";

		s = s + d6o.format(fLength)  + "\t";

		// list the edges added
		for (int i = 0; i<fEdges.size(); i++) {
			if (i == 0) {   // nice formatting
				s = s + Bipartition.toStringVerbose(fEdges.get(i).getPartition(),leaf2NumMap) + "\n";
			}
			else {
				s = s + "\t\t" + Bipartition.toStringVerbose(fEdges.get(i).getPartition(),leaf2NumMap) + "\n";
			}
		}

		return s;

	}

	public Ratio clone() {
		Ratio r = new Ratio();

		// need if statements, as different behaviour depending on whether ratio contains edges or not
		if (eEdges.size() == 0) {
			r.eLength = this.eLength;
		}
		else {
			r.addAllEEdges(Tools.myVectorClonePhyloTreeEdge(eEdges));
		}

		if (fEdges.size() == 0) {
			r.fLength = this.fLength;
		}
		else {
			r.addAllFEdges(Tools.myVectorClonePhyloTreeEdge(fEdges));
		}

		return r;
	}

}
