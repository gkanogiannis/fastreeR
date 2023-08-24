package ciat.agrobio.core;

import java.util.concurrent.ConcurrentLinkedDeque;
import java.util.concurrent.ConcurrentLinkedQueue;

public class VariantManager {
	private int numVariants;
	private int numSamples;
	private byte ploidy = 0;
	private byte[][] sampleXvariantP1 = null;
	private byte[][] sampleXvariantP2 = null;
	private boolean cleanVariantData = true;
	
	private ConcurrentLinkedDeque<Variant> variantsRaw = null;
	private ConcurrentLinkedQueue<Variant> variantsQueue = null;

	public VariantManager() {
		this.variantsRaw = new ConcurrentLinkedDeque<Variant>();
		this.variantsQueue = new ConcurrentLinkedQueue<Variant>();
	}
	
	public byte[][] getSampleXvariantP1() {
		return this.sampleXvariantP1;
	}
	
	public byte[][] getSampleXvariantP2() {
		return this.sampleXvariantP2;
	}
	
	public void setNumSamples(int numSamples) {
		this.numSamples = numSamples;
	}
	
	public final int getNumSamples() {
		return numSamples;
	}
	
	public void setNumVariants(int numVariants) {
		this.numVariants = numVariants;
	}
	
	public final int getNumVariants() {
		return numVariants;
	}
	
	public void setPloidy(byte ploidy) {
		this.ploidy = ploidy;
	}
	
	public byte getPloidy() {
		return ploidy;
	}
	
	@SuppressWarnings("unused")
	private void clear(){
		if(variantsRaw!=null){
			variantsRaw.clear();
		}
		variantsRaw = new ConcurrentLinkedDeque<Variant>();
		
		if(variantsQueue!=null){
			variantsQueue.clear();
		}
		variantsQueue = new ConcurrentLinkedQueue<Variant>();
	}
	
	public void setCleanVariantData(boolean cleanVariantData) {
		this.cleanVariantData = cleanVariantData;
	}
	
	public boolean isCleanVariantData() {
		return cleanVariantData;
	}
	
	public void populateSampleVariant() {
		if(ploidy==1){
			this.sampleXvariantP1 = new byte[numSamples][numVariants];
			Variant variant;
			while((variant = variantsQueue.poll())!=null) {
				for(int i=0; i<numSamples; i++) {
					this.sampleXvariantP1[i][variant.getVariantId()-1] = variant.getDataSamplesP1()[i];
				}
				variant.cleanDataSamples();
			}
		}
		else if(ploidy==2){
			this.sampleXvariantP1 = new byte[numSamples][numVariants];
			this.sampleXvariantP2 = new byte[numSamples][numVariants];
			Variant variant;
			while((variant = variantsQueue.poll())!=null) {
				for(int i=0; i<numSamples; i++) {
					this.sampleXvariantP1[i][variant.getVariantId()-1] = variant.getDataSamplesP1()[i];
					this.sampleXvariantP2[i][variant.getVariantId()-1] = variant.getDataSamplesP2()[i];
				}
				variant.cleanDataSamples();
			}
		}
		else{
			return;
		}
	}
	
	public boolean hasMoreRaw() {
		return !this.variantsRaw.isEmpty();
	}
	
	public boolean putVariantRaw(Variant variant) {
		try {
			long free = Runtime.getRuntime().freeMemory();
			//long max = Runtime.getRuntime().maxMemory();
			//long total = Runtime.getRuntime().totalMemory();
			if(free < 48*1024*1024) {
				return false;
			}
			this.variantsRaw.offer(variant);
			return true;
		} 
		catch (Exception e) {
			return false;
		}
	}
	
	public Variant getNextVariantRaw() {
		return this.variantsRaw.poll();
	}
	
	public void putVariantQueue(Variant variant) {
		this.variantsQueue.offer(variant);
	}
	
}
