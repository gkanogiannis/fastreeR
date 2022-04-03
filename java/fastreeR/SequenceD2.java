/*
 *
 * BioInfoJava-Utils ciat.agrobio.io.SequenceD2
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

import java.nio.ByteBuffer;

import gnu.trove.iterator.TLongDoubleIterator;
import gnu.trove.iterator.TLongIntIterator;
import gnu.trove.map.hash.TLongDoubleHashMap;

public class SequenceD2 extends Sequence implements SequenceD2Interface{

	public int as = 0;
	public int ts = 0;
	public int cs = 0;
	public int gs = 0;
	
	private long totalCounts = 0L;
	private double norm = 0.0;
	
	private TLongDoubleHashMap kmerProbs = null;
	
	public SequenceD2(Sequence read) {
		super(read.getSequenceId(), read.getHeader(), read.getSeq(), read.getQual());
		this.kmerProbs = new TLongDoubleHashMap(256);
	}

	public long getAs() {return (long)as;}
	public long getTs() {return (long)ts;}
	public long getCs() {return (long)cs;}
	public long getGs() {return (long)gs;}

	public long getNumOfElements() {
		return 1;
	}
	
	public long getTotalATCG() {
		return (long)as+(long)ts+(long)cs+(long)gs;
	}
	
	public double getDoubleCountForKmerCode(long kmerCode){
		return (double)super.getCountForKmerCode(kmerCode);
	}
	
	public double getDoubleProbForKmerCode(long kmerCode){
		return kmerProbs.get(kmerCode);
	}
	
	public long getTotalCounts(){
		if(totalCounts==0L){
			for ( TLongIntIterator it = iteratorCounts(); it.hasNext(); ) {
				it.advance();
				totalCounts += it.value();
			}
		}
		return totalCounts;
	}
	
	public double getNorm(){
		norm=0.0;
		if(norm==0.0){
			for ( TLongDoubleIterator it = iteratorProbs(); it.hasNext(); ) {
				it.advance();
				norm += Math.pow(it.value(), 2.0);
			}
			norm = Math.sqrt(norm);
		}
		return norm;
	}
	
	public void normalizeProbs(double norm){
		for ( TLongDoubleIterator it = iteratorProbs(); it.hasNext(); ) {
			it.advance();
			kmerProbs.put(it.key(), it.value()/norm);
		}
	}
	
	public void insertKmerProb(long kmerCode, short A, short T, short C, short G){
		if(kmerCode<0L){
			return;
		}
		//kmerProbs.putIfAbsent(kmerCode, ((ByteBuffer)ByteBuffer.wrap(new byte[8]).putShort(A).putShort(T).putShort(C).putShort(G).position(0)).getDouble());
		kmerProbs.put(kmerCode, ((ByteBuffer)ByteBuffer.wrap(new byte[8]).putShort(A).putShort(T).putShort(C).putShort(G).position(0)).getDouble());
		
		/*
		double previousProb = kmerProbs.get(kmerCode);
		if(previousProb ==  0.0){
			kmerProbs.put(kmerCode, ((ByteBuffer)ByteBuffer.wrap(new byte[8]).putShort(A).putShort(T).putShort(C).putShort(G).position(0)).getDouble());
		}
		*/
	}
	
	public void insertKmerProb(long kmerCode, double prob){
		if(kmerCode<0L){
			return;
		}
		kmerProbs.putIfAbsent(kmerCode, prob);
		
		/*
		double previousProb = kmerProbs.get(kmerCode);
		if(previousProb ==  0.0){
			kmerProbs.put(kmerCode, prob);
		}
		*/
	}
	
	public void insertKmerProbForce(long kmerCode, double prob){
		if(kmerCode<0L){
			return;
		}
		kmerProbs.put(kmerCode, prob);
	}
	
	public void adjustKmerProb(long kmerCode, double prob){
		if(kmerCode<0L){
			return;
		}
		kmerProbs.adjustOrPutValue(kmerCode, prob, prob);
	}
	
	public double calculateProbs(int k){	
		double sum = 0.0;
		long total = getTotalATCG();
		//System.out.println("id:"+getReadId()+"\ttotal="+total);
		//System.out.println("id:"+getReadId()+"\tas="+getAs()+"\tts="+getTs()+"\tcs="+getCs()+"\tgs="+getGs());
		ByteBuffer bb = ByteBuffer.wrap(new byte[8]);
		for ( TLongDoubleIterator it = iteratorProbs(); it.hasNext(); ) {
			it.advance();
			long kmerCode = it.key();
			
			((ByteBuffer) bb.position(0)).putDouble(kmerProbs.get(kmerCode)).position(0);	
			
			short a = bb.getShort();
			short t = bb.getShort();
			short c = bb.getShort();
			short g = bb.getShort();
			
			double prob  =  Math.pow((double)as/(double)total, (double)a);
			prob *=  		Math.pow((double)ts/(double)total, (double)t);
			prob *=  		Math.pow((double)cs/(double)total, (double)c);
			prob *=  		Math.pow((double)gs/(double)total, (double)g);
							
			insertKmerProbForce(kmerCode, prob);
			sum += prob;
		}
		//System.err.print("\tsumprob="+sum+"\t");
		//System.err.println("distinct kmers="+kmerProbs.size());
		return sum;
	}
	
	public void clear(){
		super.clearFull();
		if(kmerProbs!=null){
			kmerProbs.clear();
			kmerProbs = null;
		}
	}
	
	public TLongDoubleIterator iteratorProbs() {
		return new MyIterator();
	}
	
	private class MyIterator implements TLongDoubleIterator {

        private final TLongDoubleIterator iterator;
        private int current;
        private boolean hasNext;

        MyIterator() {
        	iterator = kmerProbs.iterator();
        	current = kmerProbs.size();
        }

        public boolean hasNext() {
        	hasNext = (current > 0);
            return hasNext;
        }

        public void advance() {
			if(!hasNext){
    			return;
    		}
        	iterator.advance();
            current--;
		}
        
        public long key() {
        	return iterator.key();
		}
        
        public double value() {
        	return iterator.value();
		}
        
        
        //Unused
        public void remove() {
      	}
      	public double setValue(double arg0) {
      		return 0.0;
      	}	
      	
	}

}
