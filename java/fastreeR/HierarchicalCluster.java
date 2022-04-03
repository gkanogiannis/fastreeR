/*
 *
 * BioInfoJava-Utils ciat.agrobio.hcluster.HierarchicalCluster
 *
 * Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */
package ciat.agrobio.hcluster;

import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;

import javax.swing.tree.DefaultTreeModel;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.core.JRITools_JavaUtils;

/**
 * Simple clustering algorithm, starting from distance matrix, 
 * single/complete/average linkage.
 * Creates a newick tree.
 * 
 */
public class HierarchicalCluster {
	public final static String AVERAGE = "Average";
	public final static String COMPLETE = "Complete";
	public final static String SINGLE = "Single";
	public final static String[] LINKAGE_METHODS = { AVERAGE, COMPLETE, SINGLE }; 
	
	private double[][] distanceMatrix; // modified during clustering
	private int linkageMethod;	
	private String[] leafLabels;  
	
	private int[] leafCount; // # leaves under each node
	private double joinDistance[]; // Inter-cluster distance; NOT branch length
	private Clade[] nodes;
	
	private GeneralTools gTools = GeneralTools.getInstance();

	/**
	 * Constructor using default linkage method.
	 */
	public HierarchicalCluster() {
		this(null, null, COMPLETE);
	}
	
	/**
	 * Constructor.
	 * @param distanceMatrix a square matrix with values in [0..1]
	 * @param labels leaf labels for each row/column in the distance matrix
	 */
	public HierarchicalCluster(double[][] distanceMatrix, String[] labels) {
		this(distanceMatrix, labels, COMPLETE);
	}
	
	/**
	 * Constructor.
	 * @param distanceMatrix a square matrix with values in [0..1]
	 * @param labels leaf labels for each row/column in the distance matrix
	 * @param linkageMethod either AVERAGE, COMPLETE or SINGLE
	 */
	public HierarchicalCluster(double[][] distanceMatrix, String[] labels, String linkageMethod) {
		setDistanceMatrix(distanceMatrix);
		setLeafLabels(labels);
		setLinkageMethod(linkageMethod);
	}
	 
	/**
	 * Cluster the data in the distance matrix.
	 * @return a tree consisting of ClusterTreeNodes containing node labels and
	 * branch lengths for each node
	 */
	public Clade cluster() {
		init();
		
		int rootIndex = -1; // catch bugs
		for (int i = 0; i < distanceMatrix.length - 1; i++)
			rootIndex = iterate();
		
		Clade root = nodes[rootIndex];
		restore(distanceMatrix);
		return root;
	}
	
	/**
	 * perform optional conversion of a similarity matrix to a distance form (it
	 * overwrites the matrix - hope that's OK!) Must do before iterations, after
	 * init !
	 * 
	 * PS. Could do it again to get back to the original sim form...
	 */
	public void similarityToDistance() {
		for (int i = 0; i < distanceMatrix.length; i++) {
			for (int j = 0; j < distanceMatrix[0].length; j++) {
				distanceMatrix[i][j] = 1.0 - distanceMatrix[i][j];
			}
		}
	}

	/**
	 * Perform necessary initialization. 
	 */
	private void init() {
		leafCount = new int[distanceMatrix.length];
		joinDistance = new double[distanceMatrix.length];
		nodes = new Clade[distanceMatrix.length];
		for (int i = 0; i < distanceMatrix.length; i++) {
			joinDistance[i] = 0.0; // not needed but it is clearer this way
			leafCount[i] = 1;
			nodes[i] = new Clade(leafLabels[i], 0);
		}
	}

	/**
	 * Compute distance between clusters.
	 * @param oneJoin index of first node
	 * @param twoJoin index of second node
	 * @param fixrow ???
	 * @return new distance between clusters
	 */
	private double newDistance(int oneJoin, int twoJoin, int fixrow) {
		double newD = -1111.; // catch bugs
		
		// Half of the distance matrix is overwritten to store 
		// the cluster-cluster distances.  The other half keeps
		// the original data.
		
		switch (linkageMethod) {
		// Average
		case 0:
			int oneN = leafCount[oneJoin];
			int twoN = leafCount[twoJoin];
			newD = (distanceMatrix[Math.max(oneJoin, fixrow)][Math.min(oneJoin, fixrow)]
					* oneN + distanceMatrix[Math.max(twoJoin, fixrow)][Math.min(twoJoin,
					fixrow)]
					* twoN)
					/ (oneN + twoN);
			break;
		// Complete
		case 1:
			newD = Math.max(distanceMatrix[Math.max(oneJoin, fixrow)][Math.min(oneJoin,
					fixrow)], distanceMatrix[Math.max(twoJoin, fixrow)][Math.min(twoJoin,
					fixrow)]);
			break;
		// Single
		case 2:
			newD = Math.min(distanceMatrix[Math.max(oneJoin, fixrow)][Math.min(oneJoin,
					fixrow)], distanceMatrix[Math.max(twoJoin, fixrow)][Math.min(twoJoin,
					fixrow)]);
			break;

		}
		return newD;
	}

	/**
	 * Iterate the clustering algorithm.
	 * @return index of the newest node.
	 */
	private int iterate() {
		double currentCloseness, newDistance = Double.MIN_VALUE; // catch bugs
		int one, two;
		one = two = -1; // catch bugs

		currentCloseness = Double.MAX_VALUE; // some big number!
		for (int i = 1; i < distanceMatrix.length; i++) {
			for (int j = 0; j < i; j++) {
				if (distanceMatrix[i][j] >= 0 && distanceMatrix[i][j] < currentCloseness) {
					currentCloseness = distanceMatrix[i][j];
					one = i;
					two = j;
				}
			}
		}
		// now do the merge and some bookkeeping.
		// above assures j<i ==> one > two
		// We write into d[big][little] only and retain the larger member
		// reference ...

		if (leafCount[one] > leafCount[two]) {
			int i = two;
			two = one;
			one = i;
		}
		
		// Update the distance matrix with new cluster-cluster distances
		for (int i = 0; i < distanceMatrix.length; i++) {
			if (one == i || two == i)
				continue;
			
			newDistance = newDistance(one, two, i);
			distanceMatrix[Math.max(two, i)][Math.min(two, i)] = newDistance; // merged distance
		}
		for (int i = 0; i < distanceMatrix.length; i++) {
			if (i != one)
				distanceMatrix[Math.max(one, i)][Math.min(one, i)] = -9999.; // it's gone now!
		}
			
		nodes[one].setBranchLength(currentCloseness - joinDistance[one]);
		nodes[two].setBranchLength(currentCloseness - joinDistance[two]);
		
		Clade joinedNode = new Clade();
		joinedNode.add(nodes[two]);
		joinedNode.add(nodes[one]);		
		nodes[two] = joinedNode;
		nodes[one] = null; 
		
		mergeMembers(leafCount, two, one); // Maybe only used for debugging
		joinDistance[two] = currentCloseness;
		return two;
	}

	private void mergeMembers(int[] numMember, int joinedNodeIndex, int oldNodeIndex) {
		numMember[joinedNodeIndex] += numMember[oldNodeIndex];
	}

    /**
     * Perform optional completion - repair the dm
     */
    private void restore(double[][] distanceMatrix) {
        for (int i=1;i<distanceMatrix.length;i++)
            for (int j=0; j<i;j++)
                  distanceMatrix[i][j] = distanceMatrix[j][i];
    }

	public double[][] getDistanceMatrix() {
		return distanceMatrix;
	}
	
	/**
	 * The distance matrix reordered based on the clustering order.
	 */
	public double[][] getReorderedMatrix(Clade root) {		
		return getReorderedMatrix(getReorderedLabels(root));
	}
	
	/**
	 * The distance matrix reordered based on the clustering order.
	 */
	public double[][] getReorderedMatrix(String[] reorderedLabels) {
		// find the index of the original label
		int[] indexes = new int[leafLabels.length];
		for(int i = 0; i < reorderedLabels.length; ++i) {
			for(int j = 0; j < reorderedLabels.length; ++j) {
				if(reorderedLabels[i] == leafLabels[j]) {
					indexes[i] = j;
					break;
				}
			}
		}
		
		// Create reordered matrix
		double[][] d = new double[distanceMatrix.length][distanceMatrix.length];
		for(int i = 0; i < d.length; ++i) {
			for(int j = 0; j < d.length; ++j) {
				d[i][j] = distanceMatrix[indexes[i]][indexes[j]];
			}
		}
		return d;
	}
	
	public String[] getReorderedLabels(Clade root) {
		String[] labels = new String[leafLabels.length];
		int k = 0;

		// iterate the tree, depth first
		@SuppressWarnings("rawtypes")
		Enumeration e = root.depthFirstEnumeration();
		while(e.hasMoreElements()) {
			Clade node = (Clade)e.nextElement();
			if(node.isLeaf() && node.getUserObject() != null)
				labels[k++] = node.toString();			
		}
		return labels;
	}

	public void setDistanceMatrix(double[][] distanceMatrix) {
		if(distanceMatrix!=null && distanceMatrix.length != distanceMatrix[0].length)
			throw new IllegalArgumentException("Distance matrix must be square");
		this.distanceMatrix = distanceMatrix;
	}

	public String[] getLeafLabels() {
		return leafLabels;
	}

	public void setLeafLabels(String[] labels) {
		this.leafLabels = labels;
	}

	public String getLinkageMethod() {
		return LINKAGE_METHODS[linkageMethod];
	}

	public void setLinkageMethod(String linkageMethod) {
		this.linkageMethod = -1;
		for (int i = 0; i < LINKAGE_METHODS.length; i++) {
			if (LINKAGE_METHODS[i].equals(linkageMethod)) {
				this.linkageMethod = i;
				return;
			}
		}
		if(this.linkageMethod == -1)
			throw new IllegalArgumentException("Invalid linkage method: " + linkageMethod);
	}
	
	public TreeMap<Integer, TreeSet<String>> hclusteringClusters(String[] sampleNames, double[][] distances, Integer minClusterSize, Double cutHeight, boolean extra){
		try {
			String treeString = this.hclusteringTree(sampleNames,distances);
			//System.err.println(treeString);
			//System.err.println(System.getProperty("java.library.path"));
			
			//Labels need to be in the order they appear in the tree string, before passed to dynamictreecut.
			//So we reorder the labels and then the distances matrix.
			String[] labelsReordered = gTools.reorderLabels(sampleNames, treeString);
			double[][] distancesReordered = gTools.reorderDistances(distances, sampleNames, labelsReordered);
		
			System.err.println("preJRI");
			JRITools_JavaUtils jritools = JRITools_JavaUtils.getInstance(null);
			System.err.println("precut");
			TreeMap<Integer, TreeSet<String>> clusters = jritools.dynamicTreeCut(treeString, distancesReordered, labelsReordered, minClusterSize, cutHeight, extra);
			System.err.println("afterCut");
			jritools.shutdown();
			
			System.err.println(GeneralTools.time()+" Clusters="+clusters.size()+"\n");
			for(Entry<Integer, TreeSet<String>> entry : clusters.entrySet()){
				int clusterId = entry.getKey();
				TreeSet<String> cluster = entry.getValue();
				//System.err.println("Cluster\t"+clusterId+"="+cluster.size());
				for(String name : cluster){
					//System.err.println("\t"+name);
					System.err.println(name+"\t"+clusterId+"\t"+cluster.size());
				}
			}
			System.err.println("\n\n");
			
			return clusters;
		}
		catch(Exception e) {
		    e.printStackTrace();
		    return null;
		}
	}
	
	public String[] hclusteringClustersNoJRI(TreeMap<Integer, TreeSet<String>> clusters){
		try {
			ArrayList<String> al = new ArrayList<String>();
			System.out.println(GeneralTools.time()+" Clusters="+clusters.size()+"\n");
			for(Entry<Integer, TreeSet<String>> entry : clusters.entrySet()){
				StringBuilder sb = new StringBuilder();
				int clusterId = entry.getKey();
				sb.append(clusterId);
				TreeSet<String> cluster = entry.getValue();
				sb.append(" "+cluster.size());
				//System.err.println("Cluster\t"+clusterId+"="+cluster.size());
				for(String name : cluster){
					//System.err.println("\t"+name);
					System.out.println(name+"\t"+clusterId+"\t"+cluster.size());
					sb.append(" "+name);
				}
				al.add(sb.toString());
			}
			System.out.println("\n\n");
			
			String[] ret = new String[al.size()];
	        ret = al.toArray(ret);
	        
	        return ret;
		}
		catch(Exception e) {
		    e.printStackTrace();
		    return null;
		}
	}
	
	public TreeMap<Integer, TreeSet<String>> findClusters(double[] result, String[] labels) {
		try {
			TreeMap<Integer, TreeSet<String>> clusters = new TreeMap<Integer, TreeSet<String>>();
			for(int i=0; i<result.length; i++){
				int clusterId = (int)result[i];
				TreeSet<String> cluster = clusters.get(clusterId);
				if(cluster==null){
					cluster = new TreeSet<String>();
					clusters.put(clusterId, cluster);
				}
				String label = labels[i];
				cluster.add(label);
			}
			return clusters;
		}
		catch(Exception e) {
		    e.printStackTrace();
		    return null;
		}
	}
	
	public TreeMap<Integer, TreeSet<String>> findClusters(int[] result, String[] labels) {
		return this.findClusters(Arrays.stream(result).asDoubleStream().toArray(), labels);
	}
	
	public String hclusteringTree(String[] sampleNames, double[][] distances){
		try {
			System.err.println(GeneralTools.time()+" Distances="+distances.length+"x"+distances[0].length);
			
			String method = HierarchicalCluster.COMPLETE;
			this.setDistanceMatrix(distances);
			this.setLeafLabels(sampleNames);
			this.setLinkageMethod(method);
			System.err.println("hierarchical method="+method);
			Clade root = this.cluster();
			
			PhylipWriter writer = new PhylipWriter();
			StringWriter sw = new StringWriter();
			writer.setOutput(sw);
			try {
				writer.write(new DefaultTreeModel(root));
			} 
			catch(IOException ex) {
				ex.printStackTrace();
			}
			
			String treeString = sw.toString().replace("\n", "");
			return treeString;
		}
		catch(Exception e) {
		    e.printStackTrace();
		    return null;
		}
	}
	
	/*
	public static void main(String[] args) {
		// Bad example: the order doesn't change...
		double[][] data5 = { 
				{0.00, 2.00, 2.83, 4.24, 5.00 },
                {2.00, 0.00, 2.00, 5.83, 6.40 },
                {2.83, 2.00, 0.00, 5.10, 5.38 },
                {4.24, 5.83, 5.10, 0.00, 1.00 },
                {5.00, 6.40, 5.38, 1.00, 0.00 }  
        };		
		double data[][] = data5; 
		boolean isPearson = (data[0][0] != 0.0);
		System.out.println("isPearson="+isPearson);

		System.out.print(String.format("\n%d\t%d\n", data.length, data[0].length));
		System.out.println("");

		String[] labels = new String[] {"one", "two", "three", "four", "five"};

		for (String method : LINKAGE_METHODS) {
			HierarchicalCluster hc = new HierarchicalCluster(data, labels);
			hc.setLinkageMethod(method);
			if (isPearson)
				hc.similarityToDistance();
			Clade root = hc.cluster();
			
			PhylipWriter writer = new PhylipWriter();
			StringWriter sw = new StringWriter();
			writer.setOutput(sw);
			try {
				writer.write(new DefaultTreeModel(root));
			} catch(IOException ex) {
				ex.printStackTrace();
			}
			System.out.println(method);
			System.out.println(sw.toString().replace("\n", ""));

			String[] relabeled = hc.getReorderedLabels(root);
			double[][] reordered = hc.getReorderedMatrix(relabeled);	
			for(int i = 0; i < relabeled.length; ++i) {
				System.out.print(relabeled[i] + "\t");
				for(int j = 0; j < relabeled.length; ++j) {
					System.out.print(reordered[i][j] + "\t");
				}
				System.out.println();
			}
			System.out.println();
			System.out.println();
			System.out.println();
		}
	}
	*/
	
}
