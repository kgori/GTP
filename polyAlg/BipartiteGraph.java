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

import java.lang.reflect.*;

import static java.lang.System.*;
import static java.lang.Math.*;


/**
 * BipartiteGraph constructs and manipulates
 * a bipartite node-weighted graph, starting as
 * input with the node-node incidence matrix of
 * the graph plus the (squared) weights on the nodes sets.
 * NodeCover has as input the two subset of nodes
 * and returns the nodes in a min-weight vertex
 * cover on the induced subgraph, with weights
 * normalized by the sum of (squared) weights.
 */

public class BipartiteGraph {
    private boolean[][] edge; //node-node incidence matrix
    //	public double[] Aweight; //(non-squared) A-side weights
//	public double[] Bweight; //(non-squared) B-side weights
    private int nA, nB, n, i, j;  //nA=#A-side nodes, nB=#B-side nodes, n=max(nA,nB)
    private Vertex[] Avertex, Bvertex; //keeps information on the vertices
    private boolean debug = false;  // set to true to display more output messages

    public BipartiteGraph(boolean IncidenceMatrix[][], double Aweight[], double Bweight[]) {
        System.out.println("BipartiteGraph::constructor 1");
        nA = Array.getLength(Aweight);
        nB = Array.getLength(Bweight);
        n = max(nA, nB);
        edge = IncidenceMatrix;
        Avertex = new Vertex[n];
        Bvertex = new Vertex[n];
        for (i = 0; i <= nA - 1; i++) Avertex[i] = new Vertex(Aweight[i]);
        for (j = 0; j <= nB - 1; j++) Bvertex[j] = new Vertex(Bweight[j]);

        if (debug) {
            out.format("Size of Bipartite Graph = A-side: %d, B-side: %d\n", nA, nB);
            out.format("Data passed out of vertex_cover will be int [4][%d]\n\n", n);
            out.format("A (squared) weights:");
            for (i = 0; i <= nA - 1; i++) out.format(" %3.2f", Avertex[i].weight);
            out.format("\nB (squared) weights:");
            for (j = 0; j <= nB - 1; j++) out.format(" %3.2f", Bvertex[j].weight);
            out.format("\n\n");
            out.format("bipartite graph incidence matrix =\n");
            for (i = 0; i <= nA - 1; i++) {
                for (j = 0; j <= nB - 1; j++)
                    out.format("%5s ", edge[i][j]);
                out.format("\n");
            }
            out.format("\n");
        }
    }


    public static void main(String[] args) {
        boolean[][] inputBG = {{false, true, false, false}, {true, false, true, true}, {false, true, false, false}, {false, true, true, true}, {false, true, false, true}};
        int[][] vertex_cover;
        double[] inputAweight = {1, 1, 1, 1, 1};
        double[] inputBweight = {1, 1, 1, 1};
        int[] argA = {0, 1, 2, 3}, argB = {0, 1, 2, 3}; //arguments to VC call
        int i0, j0;


        BipartiteGraph myBG = new BipartiteGraph(inputBG, inputAweight, inputBweight);

        vertex_cover = myBG.vertex_cover(argA, argB);
        out.format("Vertex Cover = A-side: ");
        out.format("");
        for (i0 = 0; i0 <= vertex_cover[0][0] - 1; i0++) out.format(" %d", vertex_cover[2][i0]);
        out.format("\n               B-side: ");
        for (j0 = 0; j0 <= vertex_cover[1][0] - 1; j0++) out.format(" %d", vertex_cover[3][j0]);
        out.format("\n");
    }

    public int[][] vertex_cover(int[] Aindex, int[] Bindex) {
        System.out.println("BipartiteGraph::vertex_cover");
        /**indexA,indexB are the vectors of vertex indices of the subgraph.
         * nAVC, nVC are their cardinality
         * vertex_cover computes the min-normalized-square-weighted vertex cover
         *       on indexA x indexB subgraph
         * returns a 4xn matrix CD[0][]=#A-side cover elements
         *                      CD[1][]=#B-side cover elements
         *                      CD[2][]=list of A-side cover elements
         *                      CD[3][]=list of B-side cover elements
         * Unfortunately, I can't pass out the weight (double)
         */
        int nAVC = Array.getLength(Aindex), nBVC = Array.getLength(Bindex); //nAVC,nBVC=size of A and B
        double total;
        double[][] ABflow = new double[nA][nB];
        int i, j, k, AScanListSize, BScanListSize, augmentingPathEnd = -1, Apathnode, Bpathnode;
        int[][] CD = new int[4][n]; //output: incidence vectors of vertex covers, CD[0]=Aside; CD[1]=Bside;
        int[] AScanList = new int[nA], BScanList = new int[nB]; //list of newly scanned nodes

        if (debug) {
            out.format("A indices = ");
            for (i = 0; i <= nAVC - 1; i++) System.out.format("%3d ", Aindex[i]);
            out.format("\n");
            out.format("B indices = ");
            for (j = 0; j <= nBVC - 1; j++) System.out.format("%3d ", Bindex[j]);
            out.format("\n\n");
        }

		/* First set normalized weights */
        total = 0;
        for (i = 0; i <= nAVC - 1; i++) total = total + Avertex[Aindex[i]].weight;
        for (i = 0; i <= nAVC - 1; i++) Avertex[Aindex[i]].residual = Avertex[Aindex[i]].weight / total;
        total = 0;
        for (j = 0; j <= nBVC - 1; j++) total = total + Bvertex[Bindex[j]].weight;
        for (j = 0; j <= nBVC - 1; j++) Bvertex[Bindex[j]].residual = Bvertex[Bindex[j]].weight / total;

        if (debug) {
            out.format("normalized A weights = ");
            for (i = 0; i <= nAVC - 1; i++) System.out.format("%3.2f ", Avertex[Aindex[i]].residual);
            out.format("\n");
            out.format("normalized B weights = ");
            for (j = 0; j <= nBVC - 1; j++) System.out.format("%3.2f ", Bvertex[Bindex[j]].residual);
            out.format("\n\n");
        }


		/* Now comes the flow algorithm
         * Flow on outside arcs are represented by Vertex[i].residual
		 * Flow on inside arcs are represented by ABflow
		 * Initialize ABflow to 0, start scanlist
		 */

        for (i = 0; i <= nA - 1; i++) for (j = 0; j <= nB - 1; j++) ABflow[i][j] = 0;
        total = 1; //flow augmentation in last stage
        while (total > 0) {
//out.format("flow =\n");
//for (i=0;i<=nA-1;i++) {
//	for (j=0;j<=nB-1;j++) out.format("%3.2f ",ABflow[i][j]);
//	out.format("\n");
//	}
//out.format("\n");
            //Scan Phase
            //Set labels
            total = 0;
            for (i = 0; i <= nAVC - 1; i++) {
                Avertex[Aindex[i]].label = -1;
                Avertex[Aindex[i]].pred = -1;
            }
            for (j = 0; j <= nBVC - 1; j++) {
                Bvertex[Bindex[j]].label = -1;
                Bvertex[Bindex[j]].pred = -1;
            }
            AScanListSize = 0;
            for (i = 0; i <= nAVC - 1; i++) {
                if (Avertex[Aindex[i]].residual > 0) {
                    Avertex[Aindex[i]].label = Avertex[Aindex[i]].residual;
                    AScanList[AScanListSize] = Aindex[i];
                    AScanListSize++;
                } else Avertex[Aindex[i]].label = -1;
            }
            for (i = 0; i <= nBVC - 1; i++) Bvertex[i].label = -1;
//out.format("AScanList =");
//for (i=0;i<=AScanListSize-1;i++) out.format(" %d",AScanList[i]);
//out.format("\n\n");
            // scan for an augmenting path
            scanning:
            while (AScanListSize != 0) {
				/* Scan the A side nodes*/
                BScanListSize = 0;
                for (i = 0; i <= AScanListSize - 1; i++)
                    for (j = 0; j <= nBVC - 1; j++)
                        if (edge[AScanList[i]][Bindex[j]] && Bvertex[Bindex[j]].label == -1) {
                            Bvertex[Bindex[j]].label = Avertex[AScanList[i]].label;
                            Bvertex[Bindex[j]].pred = AScanList[i];
                            BScanList[BScanListSize] = Bindex[j];
                            BScanListSize++;
                        }
//out.format("BScanList =");
//for (j=0;j<=BScanListSize-1;j++) out.format(" %d",BScanList[j]);
//out.format("\n\n");
				/* Scan the B side nodes*/
                AScanListSize = 0;
                for (j = 0; j <= BScanListSize - 1; j++)
                    if (Bvertex[BScanList[j]].residual > 0) {
                        total = min(Bvertex[BScanList[j]].residual, Bvertex[BScanList[j]].label);
                        augmentingPathEnd = BScanList[j];
                        break scanning;
                    } else for (i = 0; i <= nAVC - 1; i++)
                        if (edge[Aindex[i]][BScanList[j]] && Avertex[Aindex[i]].label == -1 && ABflow[Aindex[i]][BScanList[j]] > 0) {
                            Avertex[Aindex[i]].label = min(Bvertex[BScanList[j]].label, ABflow[Aindex[i]][BScanList[j]]);
                            Avertex[Aindex[i]].pred = BScanList[j];
                            AScanList[AScanListSize] = Aindex[i];
                            AScanListSize++;
                        }
//out.format("AScanList =");
//for (i=0;i<=AScanListSize-1;i++) out.format(" %d",AScanList[i]);
//out.format("\n\n");
            }//scanning procedure

            if (total > 0) { // flow augmentation
                Bvertex[augmentingPathEnd].residual = Bvertex[augmentingPathEnd].residual - total;
                Bpathnode = augmentingPathEnd;
                Apathnode = Bvertex[Bpathnode].pred;
//out.format("flow augmenting path: B%d A%d",Bpathnode,Apathnode);
                ABflow[Apathnode][Bpathnode] = ABflow[Apathnode][Bpathnode] + total;
                while (Avertex[Apathnode].pred != -1) {
                    Bpathnode = Avertex[Apathnode].pred;
                    ABflow[Apathnode][Bpathnode] = ABflow[Apathnode][Bpathnode] - total;
                    Apathnode = Bvertex[Bpathnode].pred;
                    ABflow[Apathnode][Bpathnode] = ABflow[Apathnode][Bpathnode] + total;
//out.format(" B%d A%d",Bpathnode,Apathnode);
                }
                Avertex[Apathnode].residual = Avertex[Apathnode].residual - total;
//out.format(";  flow = %3.2f\n\n",total);
            } else { //min vertex cover found, unlabeled A's, labeled B's
                k = 0;
                for (i = 0; i <= nAVC - 1; i++)
                    if (Avertex[Aindex[i]].label == -1) {
                        CD[2][k] = Aindex[i];
                        k++;
                    }
                CD[0][0] = k;
                k = 0;
                for (j = 0; j <= nBVC - 1; j++)
                    if (Bvertex[Bindex[j]].label >= 0) {
                        CD[3][k] = Bindex[j];
                        k++;
                    }
                CD[1][0] = k;
            }
        }//flow algorithm
        return CD;

    } //vertex_cover

}
