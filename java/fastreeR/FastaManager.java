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
	
	private AtomicInteger currSequenceId = new AtomicInteger(0);
	
	public boolean isFastq = false;
	private boolean keepQualities = false;
	private List<String> inputFileNames = null;
	private boolean done = false;
	private CountDownLatch startSignal = null;
	private CountDownLatch doneSignal = null;
	
	public FastaManager(boolean keepQualities, List<String> inputFileNames, CountDownLatch startSignal, CountDownLatch doneSignal) {
		this.keepQualities = keepQualities;
		this.inputFileNames = inputFileNames;
		this.startSignal = startSignal;
		this.doneSignal = doneSignal;
		
		this.sequences = new LinkedBlockingQueue<Sequence>();
		this.staticSequences = new ConcurrentHashMap<Integer,Sequence>();
		this.sequenceIds = new ArrayList<Integer>();
		this.sequenceNames = new ArrayList<String>();
	}
	
	public void clear(){
		if(sequences!=null) sequences.clear();
		sequences = new LinkedBlockingQueue<Sequence>();
		if(staticSequences!=null) staticSequences.clear();
		staticSequences = new ConcurrentHashMap<Integer,Sequence>();
		if(sequenceIds!=null) sequenceIds.clear();
		sequenceIds = new ArrayList<Integer>();
		if(sequenceNames!=null) sequenceNames.clear();
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
		if(!done) return true;
		else return !this.sequences.isEmpty();
	}
	
	private boolean putSequence(Sequence sequence) {
		try {
			if(sequence==null || sequence.getSeq()==null || sequence.getHeader()==null) return false;
			
			sequences.put(sequence);
			//staticSequences.put(sequence.getSequenceId(), sequence);
			sequenceIds.add(sequence.getSequenceId());
			sequenceNames.add(sequence.getShortName());
			return true;
		} 
		catch (Exception e) {
			sequences.remove(sequence);
			//staticSequences.remove(sequence.getSequenceId());
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
	
	private boolean addSequence(byte[] header, byte[] seq) {
		return addSequence(header, seq, null);
	}
	
	private boolean addSequence(byte[] header, byte[] seq, byte[] qual) {
		try {
			boolean ret = true;
			int seqId = currSequenceId.incrementAndGet();
			Sequence sequence;
			if(qual==null) sequence = new Sequence(seqId, header, seq);
			else sequence = new Sequence(seqId, header, seq, qual);
			//System.err.println(sequence.toString());
			if(!putSequence(sequence)) {
				currSequenceId.decrementAndGet();
				ret = false;
			}
			return ret;
		} 
		catch (Exception e) {
			return false;
		}
	}
	
	public void run() {
		try{
			done = false;
			
		    System.err.println(GeneralTools.time()+" FastaManager: START READ");
		    startSignal.countDown();
		    
		    List<File> inputFiles = new ArrayList<File>();
		    for(String inputFileName : inputFileNames){
		    	inputFileName = inputFileName.trim();
		    	if( inputFileName.contains(".fastq") || inputFileName.contains(".fq")) isFastq = true;
		    	else if( inputFileName.contains(".fasta") || inputFileName.contains(".fa")) isFastq = false;
		    	File f = new File(inputFileName);
		    	if(!f.exists() || !f.canRead()){
		    		System.err.println("\tERROR : File "+f+"\n\tdoes not exist ot cannot be read. Exiting.");
		    		System.exit(1);
		    	}
		    	inputFiles.add(f.getCanonicalFile());
		    }
		    
		    //Get list of lines of full reads (reads_chunk) from the fasta/fastq file (possibly multiline)
		    FastaIterator<byte[]> iterator = FastaIterator.create(inputFiles);
		    for (List<byte[]> reads_chunk : iterator) {
		    	if(reads_chunk!=null) {
		    		System.err.println(GeneralTools.time()+" FastaManager: Chunk with "+reads_chunk.size()+" lines.");
			    	boolean on_qual = false;
		    		byte[] header = null;
		    		ArrayList<byte[]> seq = new ArrayList<byte[]>();
		    		ArrayList<byte[]> qual = new ArrayList<byte[]>();
			    	
			    	//FASTA
			    	if(!isFastq) {
			    		for(byte[] line : reads_chunk) {
			    			//skip empty lines
			    			if(line.length==0) continue;
			    			//line is fasta header
			    			else if(line[0]=='>') {
			    				//Add the previous read
			    				if(!seq.isEmpty()) {
			    					addSequence(header, GeneralTools.concat(seq));
			    					seq.clear();
			    				}
			    				header = line;
			    			}
			    			//line is fasta sequence part
			    			else seq.add(line);
			    		}
			    		//Add last read
			    		if(!seq.isEmpty()) {
			    			addSequence(header, GeneralTools.concat(seq));
	    					seq.clear();
	    				}
			    	}
			    	
			    	//FASTQ
			    	else {
			    		for(byte[] line : reads_chunk) {
			    			//skip empty lines
			    			if(line.length==0) continue;
			    			//line is fastq header
			    			else if(line[0]=='@') {
			    				//Add the previous read
			    				if(!seq.isEmpty() && !qual.isEmpty()) {
			    					if(keepQualities) addSequence(header, GeneralTools.concat(seq), GeneralTools.concat(qual));
			    					else addSequence(header, GeneralTools.concat(seq));
			    					seq.clear();
			    					qual.clear();
			    					on_qual = false;
			    				}
			    				header = line;
			    			}
			    			//we are on the qualities
			    			else if(line[0]=='+') on_qual = true;
			    			//line is fasta quality part
			    			else if(on_qual) qual.add(line);
			    			//line is fasta sequence part
			    			else seq.add(line);
			    		}
			    		//Add last read
			    		if(!seq.isEmpty() && !qual.isEmpty()) {
	    					if(keepQualities) addSequence(header, GeneralTools.concat(seq), GeneralTools.concat(qual));
	    					else addSequence(header, GeneralTools.concat(seq));
	    					seq.clear();
	    					qual.clear();
	    				}
			    	}
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
