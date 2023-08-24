/*
 *
 * BioInfoJava-Utils ciat.agrobio.core.CalculateDistancesCOSINE
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
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.atomic.AtomicInteger;

import ciat.agrobio.io.VCFManager;

public class CalculateDistancesCOSINE {

	public CalculateDistancesCOSINE() {
	}
	
	public static void resetCounters(){
		CalculateDistancesChildTask.resetCounters();
	}
	
	public double[][] calculateDistances(int numOfThreads, List<String> sampleNames, VariantManager vm, VCFManager vcfm, boolean ignoreHets, boolean onlyHets, boolean ignoremissing) {
		try {
			int numOfSamples = sampleNames.size();
			double[][] distances = new double[numOfSamples][numOfSamples]; 
			CalculateDistancesTask task = new CalculateDistancesTask(sampleNames, vm, vcfm, ignoreHets, onlyHets, ignoremissing, distances);
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
class CalculateDistancesTask extends RecursiveTask<double[][]> {
	private List<String> sampleNames;
	private VariantManager vm;
	private boolean ignoreHets, onlyHets, ignoremissing;
	
	private double[][] distances;
	
	public CalculateDistancesTask(List<String> sampleNames, VariantManager vm, VCFManager vcfm, boolean ignoreHets, boolean onlyHets, boolean ignoremissing, double[][] distances) {
		this.sampleNames = sampleNames; this.vm = vm; this.ignoreHets = ignoreHets;
		this.onlyHets = onlyHets; this.ignoremissing = ignoremissing; this.distances = distances;
	}

	@SuppressWarnings("rawtypes")
	protected double[][] compute() {
		List<ForkJoinTask> children = new ArrayList<ForkJoinTask>();
		for ( int row=sampleNames.size()-1; row>=0; row--) {
			children.add(new CalculateDistancesChildTask(sampleNames, vm, ignoreHets, onlyHets, ignoremissing, distances, row));
		}
		invokeAll(children);
		return distances;
	}
	
}

@SuppressWarnings("serial")
class CalculateDistancesChildTask extends RecursiveAction {
	private static AtomicInteger sampleCounter = new AtomicInteger(0);
	private static AtomicInteger taskCount = new AtomicInteger(0);
	private final int id = taskCount.getAndIncrement();
	
	private List<String> sampleNames;
	private VariantManager vm;
	private boolean ignoreHets, onlyHets, ignoremissing;
	
	private double[][] distances;
	private int row;
	
	public static void resetCounters(){
		sampleCounter = new AtomicInteger(0);
		taskCount = new AtomicInteger(0);
	}
	
	public CalculateDistancesChildTask(List<String> sampleNames, VariantManager vm,
									   boolean ignoreHets, boolean onlyHets, boolean ignoremissing, double[][] distances, int row) {
		this.sampleNames = sampleNames; this.vm = vm; this.ignoreHets = ignoreHets;
		this.onlyHets = onlyHets; this.ignoremissing = ignoremissing; this.distances = distances; this.row = row;
	}

	@Override
	protected void compute() {
		try {
			byte ploidy = vm.getPloidy();
			//String sampleName1 = sampleNames.get(row);
			//System.err.print(sampleName1);
			double dot;
			double notmissing;
			double norm1;
			double norm2;
			for (int column=sampleNames.size()-1;column>=row;column--) {
				//String sampleName2 = sampleNames.get(column);
				//System.err.print(sampleName2);
				dot = 0.0;
				notmissing = 0.0;
				norm1 = 0.0;
				norm2 = 0.0;
				int numVariants = vm.getNumVariants();
				if(ploidy==1) {
					for(int i=0; i<numVariants; i++) {
						byte dataS1P1 = vm.getSampleXvariantP1()[row][i];
						byte dataS2P1 = vm.getSampleXvariantP1()[column][i];

						if( GenotypeEncoder.isMis(dataS1P1, dataS1P1) && GenotypeEncoder.isMis(dataS2P1, dataS2P1) ) {
							continue;
						}
						else if( GenotypeEncoder.isMis(dataS1P1, dataS1P1) ) {
							norm2 += 1;
						}
						else if( GenotypeEncoder.isMis(dataS2P1, dataS2P1) ) {
							norm1 += 1;
						}
						else {
							dot += ((dataS1P1&dataS2P1)==(byte)0?-1.0:1.0);
							norm1 += 1;
							norm2 += 1;
							notmissing++;
						}
					}
				}
				else if(ploidy==2) {
					for(int i=0; i<numVariants; i++) {
						byte dataS1P1 = vm.getSampleXvariantP1()[row][i];
						byte dataS1P2 = vm.getSampleXvariantP2()[row][i];
						byte dataS2P1 = vm.getSampleXvariantP1()[column][i];
						byte dataS2P2 = vm.getSampleXvariantP2()[column][i];

						if( GenotypeEncoder.isMis(dataS1P1, dataS1P2) && GenotypeEncoder.isMis(dataS2P1, dataS2P2) ) {
							continue;
						}
						else if(ignoreHets &&
								(GenotypeEncoder.isHet(dataS1P1, dataS1P2)||GenotypeEncoder.isHet(dataS2P1, dataS2P2)) ) {
							continue;
						}
						else if(onlyHets &&
								(!GenotypeEncoder.isHet(dataS1P1, dataS1P2)&&!GenotypeEncoder.isHet(dataS2P1, dataS2P2)) ) {
							continue;
						}
						else if( GenotypeEncoder.isMis(dataS1P1, dataS1P2) ) {
							norm2 += 1;
						}
						else if( GenotypeEncoder.isMis(dataS2P1, dataS2P2) ) {
							norm1 += 1;
						}
						else {
							dot += ((dataS1P1&dataS2P1)==(byte)0?-0.5:0.5) + ((dataS1P2&dataS2P2)==(byte)0?-0.5:0.5);
							norm1 += 1;
							norm2 += 1;
							notmissing++;
						}
					}
				}

				double distance = 0.0;
				double cosine;
				if(ignoremissing) {	
					cosine = dot / notmissing;
				}
				else {
					cosine = dot / (Math.sqrt(norm1)*Math.sqrt(norm2));
				}
				if(Double.isInfinite(cosine) || Double.isNaN(cosine)) cosine = -1.0;
				distance = (1.0 - cosine)/2.0;
				if(distance<0.0) distance = 0.0;
				distances[row][column] = Double.parseDouble(GeneralTools.decimalFormat.format(distance));
				distances[column][row] = distances[row][column];
			}
			int count = sampleCounter.incrementAndGet();
			if(count % 50 == 0) System.err.println(GeneralTools.time()+" CalculateDistancesChildTask("+id+"): "+"\t"+count);
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
