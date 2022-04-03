/*
 *
 * BioInfoJava-Utils ciat.agrobio.io.FastaManager
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
package ciat.agrobio.io;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
//import java.util.concurrent.ConcurrentLinkedDeque;
//import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import ciat.agrobio.core.GeneralTools;
import ciat.agrobio.core.Sequence;

public class FastaManager implements Runnable{
	private BlockingQueue<Sequence> sequences = null;
	private ConcurrentHashMap<Integer,Sequence> staticSequences = null;
	private List<Integer> sequenceIds = null;
	private List<String> sequenceNames = null;
	
	//private ConcurrentLinkedDeque<Integer> sequencesIds = null;
	//private ConcurrentSkipListSet<Integer> sequencesIds = null;
	
	private AtomicInteger currSequenceId = new AtomicInteger(0);
	
	public boolean isFastq = false;
	private boolean keepQualities = false;
	private List<String> inputFileNames = null;
	private boolean done = false;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	
	private GeneralTools gTools = GeneralTools.getInstance();
	
	//private BTreeMap<Integer, Sequence> sequences = null;
	
	public FastaManager(boolean keepQualities, List<String> inputFileNames, CountDownLatch startSignal, CountDownLatch doneSignal) {
		this.keepQualities = keepQualities;
		this.inputFileNames = inputFileNames;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
		
		this.sequences = new LinkedBlockingQueue<Sequence>();
		this.staticSequences = new ConcurrentHashMap<Integer,Sequence>();
		this.sequenceIds = new ArrayList<Integer>();
		this.sequenceNames = new ArrayList<String>();
		
		//this.sequencesIds = new ConcurrentLinkedDeque<Integer>();
		//this.sequencesIds = new ConcurrentSkipListSet<Integer>();
		
		//this.sequences = Utils.getDB().getTreeMap("sequences");// getHashMap("sequences");
	}
	
	public void clear(){
		if(sequences!=null){
			sequences.clear();
		}
		sequences = new LinkedBlockingQueue<Sequence>();
		
		if(staticSequences!=null){
			staticSequences.clear();
		}
		staticSequences = new ConcurrentHashMap<Integer,Sequence>();
		
		if(sequenceIds!=null){
			sequenceIds.clear();
		}
		sequenceIds = new ArrayList<Integer>();
		
		if(sequenceNames!=null){
			sequenceNames.clear();
		}
		sequenceNames = new ArrayList<String>();
	}
	
	public final ConcurrentHashMap<Integer,Sequence> getStaticSequences(){
		return staticSequences;
	}
	
	public final List<Integer> getSequenceIds() {
		return sequenceIds;
	}
	
	public final List<String> getSequenceNames() {
		return sequenceNames;
	}
	
	public boolean hasMore() {
		if(!done){
			return true;
		}
		else{
			return !this.sequences.isEmpty();
		}
	}
	
	private boolean putSequence(Sequence sequence) {
		try {
			sequences.put(sequence);
			staticSequences.put(sequence.getSequenceId(), sequence);
			sequenceIds.add(sequence.getSequenceId());
			sequenceNames.add(sequence.getShortName());
			return true;
		} 
		catch (Exception e) {
			sequences.remove(sequence);
			staticSequences.remove(sequence.getSequenceId());
			sequenceIds.remove(sequence.getSequenceId());
			sequenceNames.remove(sequence.getShortName());
			return false;
		}
	}
	
	public Sequence getNextSequence() {
		try {
			return sequences.poll(100, TimeUnit.MILLISECONDS);
		} 
		catch (Exception e) {
			return null;
		}
	}
	
	/*
	public List<Sequence> getSequences() {
		return Collections.unmodifiableList(sequences);
	}
	
	public boolean hasMore() {
		if(!done){
			return true;
		}
		else{
			if(this.sequenceIds.isEmpty()){
				return false;
			}
			else{
				return true;
			}
		}
	}
	
	private void putSequence(Sequence sequence) {
		*
		System.err.println("id="+read.getReadId());
		System.err.println("length="+read.getLength());
		if(read.getHeader()!=null)
		System.err.println("header="+new String(read.getHeader()));
		if(read.getSeq()!=null)
		System.err.println(" seq="+new String(read.getSeq()));
		if(read.getQual()!=null)
		System.err.println("qual="+new String(read.getQual()));
		//System.exit(0);
		*
		
		sequences.add(sequence);
		sequencesIds.add(sequence.getSequenceId());
	}
	
	public Sequence getNextSequence() {
		Integer sequenceId = sequencesIds.pollFirst();
		if(sequenceId == null){
			return null;
		}
		else{
			return sequences.get(sequenceId-1);
		}
	}
	*/
	
	public void run() {
		try{
			done = false;
			
		    System.err.println(GeneralTools.time()+" FastaManager: START READ");
		    startSignal.countDown();
		    
		    List<File> inputFiles = new ArrayList<File>();
		    for(String inputFileName : inputFileNames){
		    	inputFileName = inputFileName.trim();
		    	if( !isFastq && (inputFileName.endsWith(".fastq") || inputFileName.endsWith(".fq"))){
		    		isFastq = true;
		    	}
		    	inputFiles.add(new File(inputFileName).getCanonicalFile());
		    }
		    
		    //Check if files exist
		    for(File f : inputFiles){
		    	if(!f.exists() || !f.canRead()){
		    		System.err.println("\tERROR : File "+f+"\n\tdoes not exist ot cannot be read. Exiting.");
		    		System.exit(1);
		    	}
		    }
		    
		    FastaDecoder decoder = new FastaDecoder();
		    FastaIterator<byte[]> iterator = FastaIterator.create(decoder, inputFiles);
		    byte[] lastHeader = null;
		    ArrayList<byte[]> lastSeq = null;
		    byte[] lastQualSeq = null;
		    for (List<byte[]> chunk : iterator) {
		    	System.err.println(GeneralTools.time()+" FastaManager: Chunk with "+chunk.size()+" lines.");
		    	int i = 0;
		    	
		    	//From previous chunk
		    	if(lastHeader!=null){
		    		
		    		//FASTA
		    		if(!isFastq){
		    			while( i < chunk.size()){
		    				byte[] line = chunk.get(i);
		    				//Read sequence
							if(line[0] != '>'){
								lastSeq.add(line);
							}
							else{
								break;
							}
							
							i++;
		    			}
		    			//System.err.println(numOfReads+" From previous chunk A");
		    			//putRead(new Read(numOfReads++, lastHeader, Utils.concat(lastSeq)));
		    			int seqId = gTools.getRandomInteger();
		    			while(sequenceIds.contains(seqId)) {
		    				seqId = gTools.getRandomInteger();
		    			}
		    			Sequence sequence = new Sequence(seqId, lastHeader, gTools.concat(lastSeq));
		    			currSequenceId.incrementAndGet();
		    			if(!putSequence(sequence)) {
		    				currSequenceId.decrementAndGet();
		    			}
		    			lastHeader = null;
		    		}
		    		
		    		//FASTQ
		    		else{
		    			byte[] seq = null;
		    			byte[] qual = null;
						if(lastQualSeq!=null){
							seq = lastQualSeq;
							
							if(chunk.get(i)[0]=='+'){
								i++;
							}
							//Read qualities
							qual = chunk.get(i);
						}
						else{
							//Read sequence
							seq = chunk.get(i);
								
							i++;
							i++;
							//Read qualities
							qual = chunk.get(i);
						}
						//System.err.println(numOfReads+" From previous chunk Q");
						if(keepQualities){
			    			//putRead(new Read(numOfReads++, lastHeader, seq, qual));
							int seqId = gTools.getRandomInteger();
			    			while(sequenceIds.contains(seqId)) {
			    				seqId = gTools.getRandomInteger();
			    			}
							Sequence sequence = new Sequence(seqId, lastHeader, seq, qual);
			    			currSequenceId.incrementAndGet();
			    			if(!putSequence(sequence)) {
			    				currSequenceId.decrementAndGet();
			    			}
						}
			    		else{
			    			//putRead(new Read(numOfReads++, lastHeader, seq));
			    			int seqId = gTools.getRandomInteger();
			    			while(sequenceIds.contains(seqId)) {
			    				seqId = gTools.getRandomInteger();
			    			}
			    			Sequence sequence = new Sequence(seqId, lastHeader, seq);
			    			currSequenceId.incrementAndGet();
			    			if(!putSequence(sequence)) {
			    				currSequenceId.decrementAndGet();
			    			}
			    		}
						lastHeader = null;
						lastQualSeq = null;
						
						i++;
		    		}	
		    		
		    	}
		    	
		    	
		    	while( i < chunk.size() ){
		    		byte[] line = chunk.get(i);
		    		//System.err.println(i+" "+chunk.size()+" "+new String(line));
		    		
		    		//FASTA
		    		if(!isFastq && line[0] == '>'){
		    			ArrayList<byte[]> seq = new ArrayList<byte[]>();
		    			byte[] header = line;
						while( ++i < chunk.size()){
							line = chunk.get(i);
							//Read sequence
							if(line[0] != '>'){
								seq.add(line);
							}
							else{
								break;
							}
						}
						
						if(i == chunk.size()){
							lastHeader = header;
							lastSeq = seq;
						}
						else{
							//System.err.println(numOfReads+" normalA");
							//putRead(new Read(numOfReads++, header, Utils.concat(seq)));
							int seqId = gTools.getRandomInteger();
			    			while(sequenceIds.contains(seqId)) {
			    				seqId = gTools.getRandomInteger();
			    			}
							Sequence sequence = new Sequence(seqId, header, gTools.concat(seq));
			    			currSequenceId.incrementAndGet();
			    			if(!putSequence(sequence)) {
			    				currSequenceId.decrementAndGet();
			    			}
						}
					}
		    			
		    		//FASTQ
		    		else if(isFastq && line[0] == '@'){
		    			byte[] seq = null;
		    			byte[] qual = null;
		    			byte[] header = line;
						if( ++i < chunk.size()){
							//Read sequence
							seq = chunk.get(i);
							
							if( ++i < chunk.size()){
								if( ++i < chunk.size()){
									//Read qualities
									qual = chunk.get(i);
								}
							}
						}
						
						if(i == chunk.size()){
							lastHeader = header;
							lastQualSeq = seq;
						}
						else{
							//System.err.println(numOfReads+" normalQ");
				    		if(keepQualities){
				    			//putRead(new Read(numOfReads++, header, seq, qual));
				    			int seqId = gTools.getRandomInteger();
				    			while(sequenceIds.contains(seqId)) {
				    				seqId = gTools.getRandomInteger();
				    			}
				    			Sequence sequence = new Sequence(seqId, header, seq, qual);
				    			currSequenceId.incrementAndGet();
				    			if(!putSequence(sequence)) {
				    				currSequenceId.decrementAndGet();
				    			}
				    		}
				    		else{
				    			//putRead(new Read(numOfReads++, header, seq));
				    			int seqId = gTools.getRandomInteger();
				    			while(sequenceIds.contains(seqId)) {
				    				seqId = gTools.getRandomInteger();
				    			}
				    			Sequence sequence = new Sequence(seqId, header, seq);
				    			currSequenceId.incrementAndGet();
				    			if(!putSequence(sequence)) {
				    				currSequenceId.decrementAndGet();
				    			}
				    		}
						}
						
						i++;
					}
		    		
		    		else{
		    			i++;
		    		}
		    	
		    	}
		    	chunk.clear();
		    	chunk = null;
		    }
		    
		    
		    //From last chunk
	    	if(lastHeader!=null){
	 		
	    		//FASTA
	    		if(!isFastq){	
	    			//System.err.println(numOfReads+" From last chunk A");
	    			//putRead(new Read(numOfReads++, lastHeader, Utils.concat(lastSeq)));
	    			int seqId = gTools.getRandomInteger();
	    			while(sequenceIds.contains(seqId)) {
	    				seqId = gTools.getRandomInteger();
	    			}
	    			Sequence sequence = new Sequence(seqId, lastHeader, gTools.concat(lastSeq));
	    			currSequenceId.incrementAndGet();
	    			if(!putSequence(sequence)) {
	    				currSequenceId.decrementAndGet();
	    			}
	    			lastHeader = null;
	    		}
	    		
	    		//FASTQ
	    		else{
					if(!keepQualities){
						//System.err.println(numOfReads+" From last chunk Q");
		    			//putRead(new Read(numOfReads++, lastHeader, lastQualSeq));
						int seqId = gTools.getRandomInteger();
		    			while(sequenceIds.contains(seqId)) {
		    				seqId = gTools.getRandomInteger();
		    			}
						Sequence sequence = new Sequence(seqId, lastHeader, lastQualSeq);
		    			currSequenceId.incrementAndGet();
		    			if(!putSequence(sequence)) {
		    				currSequenceId.decrementAndGet();
		    			}
					}
					lastHeader = null;
					lastQualSeq = null;
	    		}	
	    	}
		    
			System.err.println(GeneralTools.time()+" FastaManager: END READ");
			System.err.println(GeneralTools.time()+" FastaManager: "+(isFastq?"FASTQ":"FASTA"));
			done = true;
			//startSignal.countDown();
			doneSignal.countDown();
		}
		catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}

}
