/*
 *
 * BioInfoJava-Utils ciat.agrobio.core.VariantProcessor
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

import java.util.Arrays;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import ciat.agrobio.io.VCFManager;

public class VariantProcessor implements Runnable {

	private static AtomicInteger variantCount = new AtomicInteger(0);
	
	private static AtomicInteger taskCount = new AtomicInteger(0);
	
	private final int id = taskCount.getAndIncrement();
	private VariantManager vm = null;
	private VCFManager vcfm = null;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	
	//private ConcurrentSkipListSet<String> variantNames;
	
	public static void resetCounters(){
		variantCount = new AtomicInteger(0);
		taskCount = new AtomicInteger(0);
	}
	
	public VariantProcessor(VariantManager vm, VCFManager vcfm, CountDownLatch startSignal, CountDownLatch doneSignal) {
		this.vm = vm;
		this.vcfm = vcfm;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
	}
	
	public static AtomicInteger getVariantCount() {
		return variantCount;
	}

	public void setStartSignal(CountDownLatch startSignal) {
		this.startSignal = startSignal;
	}

	public void setDoneSignal(CountDownLatch doneSignal) {
		this.doneSignal = doneSignal;
	}

	public int getId() {
		return id;
	}

	public void run() {
		try{
			startSignal.await();
			boolean done = false;
			while(!done){
				Variant variant = vm.getNextVariantRaw();
				int numSamples = vm.getNumSamples();
				byte ploidy = vm.getPloidy();
				if(variant==null){
					if(!vm.hasMoreRaw() && vcfm.isDone()){
						done = true;
						break;
					}
					continue;
				}
				//Process variant data
				int count = variantCount.incrementAndGet(); 
				if(count % 50000 == 0){
					System.err.println(GeneralTools.time()+" VariantProcessor: "+id+"\t"+count);
				}
				
				byte numAlleles=0;
				int indexGT = Arrays.asList(variant.getFormat().split(":")).indexOf("GT");
				if(ploidy==2) {
					byte[][] variantDataRaw = variant.getDataRaw();
					byte[] variantDataSamplesP1 = new byte[numSamples];
					byte[] variantDataSamplesP2 = new byte[numSamples];
					variant.setDataSamplesP1(variantDataSamplesP1);
					variant.setDataSamplesP2(variantDataSamplesP2);
					String GT;
					byte[] GTCode;
					for(int i=0; i<numSamples; i++) {
						GT = new String(variantDataRaw[i+9]).split(":")[indexGT];
						GTCode = GenotypeEncoder.encodeGT(GT,ploidy);
						variantDataSamplesP1[i] = GTCode[0];
						variantDataSamplesP2[i] = GTCode[1];
						if(GTCode[ploidy]>numAlleles) {
							numAlleles = GTCode[ploidy];
						}
					}
				}
				else {
					System.err.println(" Ploidy "+ploidy+" is not yet supported.");
					System.exit(ploidy);
				}
				variant.setNumAlleles(numAlleles);
				if(vm.isCleanVariantData()) {
					variant.cleanDataRaw();
				}
				vm.putVariantQueue(variant);
			}
			doneSignal.countDown();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
}
