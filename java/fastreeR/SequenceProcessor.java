/*
 *
 * BioInfoJava-Utils ciat.agrobio.core.SequenceProcessor
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

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicInteger;

import ciat.agrobio.io.FastaManager;

public class SequenceProcessor implements Runnable {

	private static AtomicInteger sequenceCount = new AtomicInteger(0);
	
	private static AtomicInteger taskCount = new AtomicInteger(0);
	
	private final int id = taskCount.getAndIncrement();
	private FastaManager frm = null;
	private int k;
	private boolean normalize;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	
	private ConcurrentHashMap<Integer, SequenceD2> seqVectors;
	
	public static void resetCounters(){
		sequenceCount = new AtomicInteger(0);
		taskCount = new AtomicInteger(0);
	}
	
	public SequenceProcessor(ConcurrentHashMap<Integer, SequenceD2> seqVectors, FastaManager frm, int k, boolean normalize, CountDownLatch startSignal, CountDownLatch doneSignal) {
		this.seqVectors = seqVectors;
		this.frm = frm;
		this.k = k;
		this.normalize = normalize;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
	}
	
	public static AtomicInteger getSequenceCount() {
		return sequenceCount;
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

	@SuppressWarnings("unused")
	public void run() {
		try{
			startSignal.await();
			boolean done = false;
			while(!done){
				Sequence sequence = frm.getNextSequence();
				SequenceD2 sequenceD2 = null;
				if(sequence!=null){
					sequenceD2 = new SequenceD2(sequence);
					sequenceD2.getName(); sequenceD2.getShortName();
				}
				if(sequenceD2==null){
					if(!frm.hasMore()){
						done = true;
						break;
					}
					continue;
				}
				sequenceCount.incrementAndGet();
				
				long kmerCode;
				for(int i=0; i+k<=sequenceD2.getLength(); i++){	
					int oldAs = sequenceD2.as;
					int oldTs = sequenceD2.ts;
					int oldCs = sequenceD2.cs;
					int oldGs = sequenceD2.gs;
					kmerCode = SequenceEncoder.encodeToLong(sequenceD2, i, i+k, false, true);
					if(kmerCode>=0L){
						sequenceD2.insertKmerCount(kmerCode, 1);
						sequenceD2.insertKmerProb(kmerCode, (short)(sequenceD2.as-oldAs), (short)(sequenceD2.ts-oldTs), (short)(sequenceD2.cs-oldCs), (short)(sequenceD2.gs-oldGs));
					}
					//reverse
					oldAs = sequenceD2.as;
					oldTs = sequenceD2.ts;
					oldCs = sequenceD2.cs;
					oldGs = sequenceD2.gs;
					kmerCode = SequenceEncoder.encodeToLong(sequenceD2, i, i+k, true, true);
					if(kmerCode>=0L){
						sequenceD2.insertKmerCount(kmerCode, 1);
						sequenceD2.insertKmerProb(kmerCode, (short)(sequenceD2.as-oldAs), (short)(sequenceD2.ts-oldTs), (short)(sequenceD2.cs-oldCs), (short)(sequenceD2.gs-oldGs));
					}
				}
				//calculate per sequence probs
				double sumprob = sequenceD2.calculateProbs(k);
				if(normalize) {
					sequenceD2.normalizeProbs(sequenceD2.getNorm());
				}
				SequenceD2 seqVector = sequenceD2;
				//ReadD2Centroid seqVector = new ReadD2Centroid(read);
				seqVectors.put(seqVector.getSequenceId(), seqVector);
								
				System.err.println(sequenceCount.get()+"\t"+seqVector.getShortName());
				
				//Clean sequence
				sequenceD2.clearHeadSeq();

				/*
				System.err.print("\tkmer_Count="+seqVector.getTotalCounts());
				System.err.print("\tATCG="+seqVector.getTotalATCG());
				System.err.print("\tAs="+seqVector.getAs());
				System.err.print("\tTs="+seqVector.getTs());
				System.err.print("\tCs="+seqVector.getCs());
				System.err.print("\tGs="+seqVector.getGs());
				System.err.print("\tnorm="+seqVector.getNorm());
				System.err.println("\tsumprob="+sumprob);
				*/
				//System.err.println(GeneralTools.time()+" END of Vector Creation for Sequence="+seqName);
			}
			doneSignal.countDown();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	
}
