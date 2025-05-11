/*
 *
 * BioInfoJava-Utils ciat.agrobio.core.CalculateDistancesD2
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
package ciat.agrobio.core;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.atomic.AtomicInteger;

import ciat.agrobio.io.FastaManager;
import gnu.trove.iterator.TLongIntIterator;

public class CalculateDistancesD2 {

	public static boolean verbose = false;

	public CalculateDistancesD2(boolean verbose) {
		CalculateDistancesD2.verbose = verbose;
	}
	
	public static void resetCounters(){
		CalculateD2ChildTask.resetCounters();
	}
	
	public double[][] calculateDistances(int numOfThreads, ConcurrentHashMap<Integer, SequenceD2> seqVectors, 
										 List<String> seqNames, 
										 List<Integer> seqIds,
										 FastaManager frm) {
		try {
			int numOfSequences = seqNames.size();
			double[][] distances = new double[numOfSequences][numOfSequences]; 
			CalculateD2Task task = new CalculateD2Task(seqVectors, seqNames, seqIds, frm, distances);
			ForkJoinPool pool = new ForkJoinPool(numOfThreads);
			pool.execute(task);
			
			return task.get();
		} 
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
}

@SuppressWarnings("serial")
class CalculateD2Task extends RecursiveTask<double[][]> {
	private ConcurrentHashMap<Integer, SequenceD2> seqVectors;
	private List<String> seqNames;
	private List<Integer> seqIds;
	private FastaManager frm;
	
	private double[][] distances;
	
	public CalculateD2Task(ConcurrentHashMap<Integer, SequenceD2> seqVectors, 
					       List<String> seqNames, List<Integer> seqIds, FastaManager frm, double[][] distances) {
		this.seqVectors = seqVectors;
		this.seqNames = seqNames; 
		this.seqIds = seqIds;
		this.frm = frm;
		this.distances = distances;
	}

	@SuppressWarnings("rawtypes")
	protected double[][] compute() {
		List<ForkJoinTask> children = new ArrayList<ForkJoinTask>();
		for ( int row=seqNames.size()-1; row>=0; row--) {
			children.add(new CalculateD2ChildTask(seqVectors, seqNames, seqIds, frm, distances, row));
		}
		invokeAll(children);
		return distances;
	}
	
}

@SuppressWarnings("serial")
class CalculateD2ChildTask extends RecursiveAction {
	private static AtomicInteger sequenceCounter = new AtomicInteger(0);
	private static AtomicInteger taskCount = new AtomicInteger(0);
	private final int id = taskCount.getAndIncrement();
	
	private ConcurrentHashMap<Integer, SequenceD2> seqVectors;
	private List<String> seqNames;
	private List<Integer> seqIds;
	@SuppressWarnings("unused")
	private FastaManager frm;
	
	private double[][] distances;
	private int row;
	
	public static void resetCounters(){
		sequenceCounter = new AtomicInteger(0);
		taskCount = new AtomicInteger(0);
	}
	
	public CalculateD2ChildTask(ConcurrentHashMap<Integer, SequenceD2> seqVectors, 
							    List<String> seqNames, List<Integer> seqIds, FastaManager frm, double[][] distances, int row) {
		this.seqVectors = seqVectors;
		this.seqNames = seqNames; 
		this.seqIds = seqIds;
		this.frm = frm;
		this.distances = distances; 
		this.row = row;
	}

	@Override
	protected void compute() {
		try {
			SequenceD2 X = seqVectors.get(seqIds.get(row));
			for (int column=seqNames.size()-1;column>=row;column--) {
				SequenceD2 Y = seqVectors.get(seqIds.get(column));
				double d2_measure = DissimilarityMeasuresD2.d2_S_Dissimilarity(X, Y);
				if(Math.abs(d2_measure) < Double.valueOf("1E-15")) d2_measure=0.0;
				distances[row][column] = Double.parseDouble(GeneralTools.decimalFormat.format(d2_measure));
				distances[column][row] = distances[row][column];
				//System.out.println(Utils.time()+"\td2S for couple ["+seqNames.get(i)+" :: "+seqNames.get(j)+"]="+d2_measure);
			}
			int count = sequenceCounter.incrementAndGet();
			if(count % 50 == 0 && CalculateDistancesD2.verbose) System.err.println(GeneralTools.time()+" CalculateD2ChildTask("+id+"): "+"\t"+count);
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}

class DissimilarityMeasuresD2 {
	
	public static double d2_S_Dissimilarity(SequenceD2Interface X, SequenceD2Interface Y){
		double D2_S = 0.0;
		double tempX = 0.0;
		double tempY = 0.0;
		double cXi = 0.0;
		double cYi = 0.0;
		double cXi_bar = 0.0;
		double cYi_bar = 0.0;
		double temp3 = 0.0;
		
		HashSet<Long> set = new HashSet<Long>();
		
		for ( TLongIntIterator itX = X.iteratorCounts(); itX.hasNext(); ) {
			itX.advance();
			long kmerCodeX = itX.key();
			set.add(kmerCodeX);
		}
		for ( TLongIntIterator itY = Y.iteratorCounts(); itY.hasNext(); ) {
			itY.advance();
			long kmerCodeY = itY.key();
			set.add(kmerCodeY);
		}
		//set.addAll(X.getKmerCodes());
		//set.addAll(Y.getKmerCodes());
		
		//Iterator<Long> iter = X.iterator();
		//long kmerCode;
		for(long kmerCode : set){
		//while(iter.hasNext()){
			//kmerCode = iter.next();
			/*
			cXi = X.getCountForKmerCode(kmerCode);
			cYi = Y.getCountForKmerCode(kmerCode);
			if(cYi > 0){
				pXi = X.getProbForKmerCode(kmerCode);
				pYi = Y.getProbForKmerCode(kmerCode);
				cXi_bar = (double)cXi - (double)(X.getLength()-k+1)*pXi;
				cYi_bar = (double)cYi - (double)(Y.getLength()-k+1)*pYi;
				temp1 = Math.pow(cXi_bar, 2.0);
				temp2 = Math.pow(cYi_bar, 2.0);
				temp3 = Math.sqrt(temp1+temp2);
				if(temp3 == 0.0){
					temp3 = 1.0;
				}
				D2_S += cXi_bar*cYi_bar/temp3;
				tempX += cXi_bar*cXi_bar/temp3;
				tempY += cYi_bar*cYi_bar/temp3;
			}
			*/
			cYi = Y.getDoubleCountForKmerCode(kmerCode);
			cXi = X.getDoubleCountForKmerCode(kmerCode);
			//System.out.println("kmerCode="+kmerCode);
			//System.out.println("\tcountProbX="+X.getDoubleCountForKmerCode(kmerCode)+"  :  "+X.getDoubleProbForKmerCode(kmerCode));
			//System.out.println("\tcountProbY="+Y.getDoubleCountForKmerCode(kmerCode)+"  :  "+Y.getDoubleProbForKmerCode(kmerCode));
			//if(cYi > 0.0){
				cXi_bar = cXi - (double)X.getTotalCounts()*X.getDoubleProbForKmerCode(kmerCode);
				cYi_bar = cYi - (double)Y.getTotalCounts()*Y.getDoubleProbForKmerCode(kmerCode);
				//cXi_bar = X.getCountForKmerCodeDouble(kmerCode) - X.getProbForKmerCode(kmerCode);
				//cYi_bar = cYi - Y.getProbForKmerCode(kmerCode);
				//System.out.println("\tcXi_bar="+ cXi_bar);
				//System.out.println("\tcYi_bar="+ cYi_bar);
				temp3 = Math.sqrt(Math.pow(cXi_bar, 2.0) + Math.pow(cYi_bar, 2.0));
				if(temp3 == 0.0){
					temp3 = 1.0;
				}
				D2_S += (cXi_bar*cYi_bar)/temp3;
				tempX += (cXi_bar*cXi_bar)/temp3;
				tempY += (cYi_bar*cYi_bar)/temp3;
			//}
		}
		
		//System.out.println(tempY);
		
		tempX = Math.sqrt(tempX);
		tempY = Math.sqrt(tempY);
		double temp = D2_S/(tempX*tempY);
		return 0.5*(1.0 - temp);
	}
	
	public static double d2_Star_Dissimilarity(SequenceD2Interface X, SequenceD2Interface Y){
		double D2_Star = 0.0;
		double tempX = 0.0;
		double tempY = 0.0;
		double cYi = 0;
		double cXi_bar = 0.0;
		double cYi_bar = 0.0;
		double temp3 = 0.0;
		
		HashSet<Long> set = new HashSet<Long>();
		
		for ( TLongIntIterator itX = X.iteratorCounts(); itX.hasNext(); ) {
			itX.advance();
			long kmerCodeX = itX.key();
			set.add(kmerCodeX);
		}
		for ( TLongIntIterator itY = Y.iteratorCounts(); itY.hasNext(); ) {
			itY.advance();
			long kmerCodeY = itY.key();
			set.add(kmerCodeY);
		}
		//set.addAll(X.getKmerCodes());
		//set.addAll(Y.getKmerCodes());
		
		for(long kmerCode : set){
			cYi = Y.getDoubleCountForKmerCode(kmerCode);
			//if(cYi > 0.0){
				cXi_bar = X.getDoubleCountForKmerCode(kmerCode) - (double)X.getTotalCounts()*X.getDoubleProbForKmerCode(kmerCode);
				cYi_bar = cYi - (double)Y.getTotalCounts()*Y.getDoubleProbForKmerCode(kmerCode);
				temp3 = Math.sqrt((double)X.getTotalCounts()*X.getDoubleProbForKmerCode(kmerCode)) * Math.sqrt((double)Y.getTotalCounts()*Y.getDoubleProbForKmerCode(kmerCode));
				if(temp3 == 0.0){
					temp3 = 1.0;
				}
				D2_Star += (cXi_bar*cYi_bar)/temp3;
				tempX += (cXi_bar*cXi_bar)/temp3;
				tempY += (cYi_bar*cYi_bar)/temp3;
			//}
		}
		
		tempX = Math.sqrt(tempX);
		tempY = Math.sqrt(tempY);
		double temp = D2_Star/(tempX*tempY);
		return 0.5*(1.0 - temp);
	}
	
	/*
	public static double d2_Star_Dissimilarity(ReadD2 X, ReadD2 Y, int k){
		double D2_Star = 0.0;
		double d2_Star = 0.0;
		double tempX = 0.0;
		double tempY = 0.0;
		int cXi = 0;
		int cYi = 0;
		double pXi = 0.0;
		double pYi = 0.0;
		double cXi_bar = 0.0;
		double cYi_bar = 0.0;
		double temp1 = 0.0;
		double temp2 = 0.0;
		double temp3 = 0.0;
		
		Iterator<Long> iter = X.iterator();
		long kmerCode;
		while(iter.hasNext()){
			kmerCode = iter.next();
			cXi = X.getCountForKmerCode(kmerCode);
			cYi = Y.getCountForKmerCode(kmerCode);
			if(cYi > 0){
				pXi = X.getProbForKmerCode(kmerCode);
				pYi = Y.getProbForKmerCode(kmerCode);
				temp1 = (double)(X.getLength()-k+1)*pXi;
				temp2 = (double)(Y.getLength()-k+1)*pYi;
				cXi_bar = (double)cXi - temp1;
				cYi_bar = (double)cYi - temp2;
				temp3 = Math.sqrt(temp1*temp2);
				if(temp1 == 0.0){
					temp1 = 1.0;
					temp3 = 1.0;
				}
				if(temp2 == 0.0){
					temp2 = 1.0;
					temp3 = 1.0;
				}
				if(temp3 == 0.0){
					temp3 = 1.0;
				}
				
				D2_Star += cXi_bar*cYi_bar/temp3;
				tempX += cXi_bar*cXi_bar/temp1;
				tempY += cYi_bar*cYi_bar/temp2;
			}
		}
		
		tempX = Math.sqrt(tempX);
		tempY = Math.sqrt(tempY);
		double temp = D2_Star/(tempX*tempY);
		d2_Star = 0.5*(1-temp);
		return d2_Star;
	}
	*/
}
