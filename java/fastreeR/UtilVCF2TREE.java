/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilVCF2TREE
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
package ciat.agrobio.javautils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.CalculateDistancesCOSINE;
import ciat.agrobio.core.VariantManager;
import ciat.agrobio.core.VariantProcessor;
import ciat.agrobio.hcluster.HierarchicalCluster;
import ciat.agrobio.io.VCFManager;

@Parameters(commandDescription = "VCF2TREE")
public class UtilVCF2TREE {

	public UtilVCF2TREE() {
	}

	public static String getUtilName() {
		return "VCF2TREE";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(names = "--verbose")
	private boolean verbose = false;

	@Parameter(description = "Input_File", required = true)
	private String inputFileName;
	
	@Parameter(names = { "--numberOfThreads", "-t" })
	private int numOfThreads = 1;
	
	@Parameter(names = { "--ignoreMissing", "-m" })
	private boolean ignoreMissing = false;
	
	@Parameter(names={"--onlyHets", "-h"})
	private boolean onlyHets = false;
	
	@Parameter(names={"--ignoreHets", "-g"})
	private boolean ignoreHets = false;

	public void go() {
		try {
			int cpus = Runtime.getRuntime().availableProcessors();
			int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			if(verbose) System.err.println("cpus=" + cpus);
			if(verbose) System.err.println("using=" + usingThreads);

			CountDownLatch startSignal = new CountDownLatch(1);
			CountDownLatch doneSignal = new CountDownLatch(usingThreads + 1);

			//if(verbose) System.err.println(GeneralTools.time() + " START ");

			ExecutorService pool = Executors.newFixedThreadPool(usingThreads + 1);

			Map<Integer, VariantProcessor> variantProcessors = new HashMap<Integer, VariantProcessor>();
			VariantManager vm = new VariantManager();
			VCFManager vcfm = new VCFManager(vm, inputFileName, startSignal, doneSignal);
			pool.execute(vcfm);

			VariantProcessor.resetCounters();
			// Starting threads
			for (int i = 0; i < usingThreads; i++) {
				VariantProcessor vp = new VariantProcessor(vm, vcfm, startSignal, doneSignal, verbose);
				variantProcessors.put(vp.getId(), vp);
				pool.execute(vp);
			}

			doneSignal.await();
			pool.shutdown();
			
			vm.populateSampleVariant();

			List<String> sampleNames = new ArrayList<String>();
			for (int i = 9; i < vcfm.getHeaderData().length; i++) {
				sampleNames.add(new String(vcfm.getHeaderData()[i]));
			}

			if(verbose) System.err.printf("\rProcessed variants : \t%8d\n", VariantProcessor.getVariantCount().get());

			// Calculate distances
			//double[][] distances = SampleVariantTools.calculateDistances(samplesToVariantsData, sampleNames, vcfm, ignoreHets, onlyHets, false);
			CalculateDistancesCOSINE.resetCounters();
			CalculateDistancesCOSINE fj = new CalculateDistancesCOSINE(verbose);
			double[][] distances = fj.calculateDistances(usingThreads, sampleNames, vm, vcfm, ignoreHets, onlyHets, ignoreMissing);
			
			//HCluster tree
			HierarchicalCluster hc = new HierarchicalCluster();
			String treeString = hc.hclusteringTree(sampleNames.toArray(new String[sampleNames.size()]), distances);
			System.out.println(treeString);
			
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
