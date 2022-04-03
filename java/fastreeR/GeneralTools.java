/*
 *
 * BioInfoJava-Utils ciat.agrobio.core.GeneralTools
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
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.security.SecureRandom;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;

//import org.mapdb.DB;
//import org.mapdb.DBMaker;

public class GeneralTools {
	
	public static final DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	public static final DecimalFormat decimalFormat = new DecimalFormat("#.############", new DecimalFormatSymbols(Locale.US));
	
	private static SecureRandom prng;
	static {
		try {
			prng = SecureRandom.getInstance("SHA1PRNG");
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private static GeneralTools instance = new GeneralTools();
	
	private GeneralTools() {
	}

	public static GeneralTools getInstance() {
		return instance;
	}
	
	public int getRandomInteger() {
		return prng.nextInt();
	}
	
	/*
	private static DB db = null;
	
	static{
		db = DBMaker
				.newHeapDB()
			    .transactionDisable()
			    .asyncWriteFlushDelay(100)
			    .make();
		//db = DBMaker.newMemoryDB().asyncWriteEnable(). /*cacheHardRefEnable()./transactionDisable().make();
		//db = DBMaker.newAppendFileDB(new File("/env/cns/bigtmp1/agkanogi/mapdb.temp")).closeOnJvmShutdown().transactionDisable().make();
		//db = DBMaker.newFileDB(new File("/env/cns/bigtmp1/agkanogi/mapdb.temp")).closeOnJvmShutdown().randomAccessFileEnable().transactionDisable().make();
	}
	
	public static DB getDB(){
		return db;
	}
	*/
		
	public static String time() {
		Calendar cal = Calendar.getInstance();
		return dateFormat.format(cal.getTime());
	}
	
	public String RAMInfo(Runtime runtime){
		double gb = (double)1024*1024*1024;
        NumberFormat format = NumberFormat.getInstance(); 
        StringBuilder sb = new StringBuilder();
        
        sb.append("\n##### Heap utilization statistics [GB] #####\n");
        sb.append("Max Memory=" + format.format(runtime.maxMemory() / gb) + "\n");
        sb.append("Current Total Memory=" + format.format(runtime.totalMemory() / gb) + "\n");
        sb.append("Current Used Memory="  + format.format((runtime.totalMemory() - runtime.freeMemory()) / gb) + "\n");
        sb.append("Current Free Memory="  + format.format(runtime.freeMemory() / gb) + "\n");
        sb.append("############################################\n");
	
        return sb.toString();
	}

	final public byte[][] concatSeqQual(final ArrayList<byte[]> arrays) {
		byte[] all = concat(arrays);
		int size = all.length;
		int sizeSeq = 0;
		int sizeQual = 0;
		for(byte b : all){
			if('+' == (char)(b & 0xff)){
				sizeQual = size - sizeSeq - 1;
				break;
			}
			else{
				sizeSeq++;
			}
		}
		
		byte[][] ret = new byte[2][];
		if(sizeSeq>0){
			byte[] destSeq = new byte[sizeSeq];
			System.arraycopy(all, 0, destSeq, 0, sizeSeq);
			ret[0] = destSeq;
		}
		if(sizeQual>0){
			byte[] destQual = new byte[sizeQual];
			System.arraycopy(all, sizeSeq+1, destQual, 0, sizeQual);
			ret[1] = destQual;
		}
		
		return ret;
	}
	
	final public byte[] concat(final ArrayList<byte[]> arrays) {
		int size = 0;
		for (byte[] a : arrays)
			size += a.length;

		byte[] dest = new byte[size];

		int destPos = 0;
		for (int i = 0; i < arrays.size(); i++) {
			if (i > 0)
				destPos += arrays.get(i - 1).length;
			int length = arrays.get(i).length;
			System.arraycopy(arrays.get(i), 0, dest, destPos, length);
		}

		return dest;
	}
	
	public int round(double d){
	    double dAbs = Math.abs(d);
	    int i = (int) dAbs;
	    double result = dAbs - (double) i;
	    if(result<0.5){
	        return d<0 ? -i : i;            
	    }else{
	        return d<0 ? -(i+1) : i+1;          
	    }
	}

	public char[] bytesToStringUTFCustom(byte[] bytes) {
		char[] buffer = new char[bytes.length >> 1];
		for (int i = 0; i < buffer.length; i++) {
			int bpos = i << 1;
			char c = (char) (((bytes[bpos] & 0x00FF) << 8) + (bytes[bpos + 1] & 0x00FF));
			buffer[i] = c;
		}
		return buffer;
		//return new String(buffer);
	}
	
	public byte[] stringToBytesUTFCustom(char[] str) {
		byte[] b = new byte[str.length << 1];
		for (int i = 0; i < str.length; i++) {
			char strChar = str[i];
			int bpos = i << 1;
			b[bpos] = (byte) ((strChar & 0xFF00) >> 8);
			b[bpos + 1] = (byte) (strChar & 0x00FF);
		}
		return b;
	}
	
	public List<String> directoriesToFiles(List<String> inputFastaFileNames){
		ArrayList<String> outputFastaFileNames = new ArrayList<String>();
		try {
			for(String s : inputFastaFileNames){
				File f = new File(s);
				if(!f.isDirectory()){
					outputFastaFileNames.add(f.getCanonicalPath());
				}
				else{
					for(File ff : f.listFiles()){
						outputFastaFileNames.add(ff.getCanonicalPath());
					}
				}
			}
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
		return outputFastaFileNames;
	}
	
	public <K, V extends Comparable<? super V>> Map<K, V> sortByValue( Map<K, V> map ){
		List<Map.Entry<K, V>> list = new LinkedList<Map.Entry<K, V>>( map.entrySet() );
		Collections.sort( list, new Comparator<Map.Entry<K, V>>(){
			public int compare( Map.Entry<K, V> o1, Map.Entry<K, V> o2 ){
				int compare = (o1.getValue()).compareTo(o2.getValue());
				if(compare == 0){
					Long key1 = (Long)o1.getKey();
					Long key2 = (Long)o2.getKey();
					return +1 * (key1).compareTo(key2);
				}
				else{
					return -1 * compare;
				}
			}
		} );

		Map<K, V> result = new LinkedHashMap<K, V>();
		for (Map.Entry<K, V> entry : list){
			result.put( entry.getKey(), entry.getValue() );
		}
		return result;
	}
	
	public ArrayList<String> generateKmerVocabulary(int k){
		ArrayList<String> vocabulary = new ArrayList<String>();
		
		char[] chars = "ATCG".toCharArray();
        int len = k;
        iterate(chars, len, new char[len], 0, vocabulary);
        
		return vocabulary;
	}
	
	private void iterate(char[] chars, int len, char[] build, int pos, ArrayList<String> vocabulary) {
        if (pos == len) {
            String word = new String(build);
            vocabulary.add(word);
            return;
        }

        for (int i = 0; i < chars.length; i++) {
            build[pos] = chars[i];
            iterate(chars, len, build, pos + 1, vocabulary);
        }
    }
	
	/*
	public static String getRandomSequence(int length){
		return RandomStringUtils.random(length, "ATCG");
	}
	*/
	
	public final String reverseComplement(String sequence){
		StringBuilder sb = new StringBuilder();
		for(int i=sequence.length()-1; i>=0; i--){
			switch (sequence.charAt(i)){
				case 'A':
				case 'a':
					sb.append("T");
				break;

				case 'T':
				case 't':
					sb.append("A");
				break;

				case 'C':
				case 'c':
					sb.append("G");
				break;

				case 'G':
				case 'g':
					sb.append("C");
				break;

				default:
					;
				break;
			}
		}
		return sb.toString();
	}

	@SuppressWarnings("rawtypes")
	public long classBuildTimeMillis(Class c){
	    URL resource = c.getResource(c.getSimpleName() + ".class");
	    if (resource == null) {
	        System.out.println("Failed to find class file for class: " + c.getName());
	        return 0L;
	    }
	    
	    if (resource.getProtocol().equals("file")) {
	    	try {
				return new File(resource.toURI()).lastModified();
			} 
	    	catch (URISyntaxException e) {
				e.printStackTrace();
				return 0L;
	    	}

	    } 
	    else if (resource.getProtocol().equals("jar")) {
	        String path = resource.getPath();
	        return new File(path.substring(5, path.indexOf("!"))).lastModified();
	    } 
	    else {
	        System.out.println("Unhandled url protocol: " + resource.getProtocol() + " for class: " + c.getName() + " resource: " + resource.toString());
	        return 0L;
	    }
	}

	public <T extends Comparable<? super T>> List<T> asSortedList(Collection<T> c) {
		List<T> list = new ArrayList<T>(c);
		Collections.sort(list);
		return list;
	}
	
	public String[] reorderLabels(String[] original, String treeString){
		String[] labelsReordered = new String[original.length];
		
		String[] data = treeString.split(",");
		for(int i=0; i<data.length; i++){
			String label = data[i].trim();
			
			if(label.contains("(")){
				label = label.substring(label.lastIndexOf("(")+1);
			}
			if(label.contains(")")){
				label = label.substring(0, label.indexOf(")"));
			}
			if(label.contains(":")){
				label = label.substring(0, label.indexOf(":"));
			}
			
			//System.out.println("i="+i+"\tlabel="+label);
			labelsReordered[i] = label;
		}
		
		return labelsReordered;
	}
	
	public double[][] reorderDistances(double[][] distancesOriginal, String[] labelsOriginal, String[] labelsReordered){
		int size = distancesOriginal.length;
		double[][] distancesReordered = new double[size][size];
		
		int[] indexToReorder = new int[size];
		for(int i=0; i<size; i++){
			String labelOriginal = labelsOriginal[i].trim();
			for(int j=0; j<size; j++){
				String labelReordered = labelsReordered[j].trim();
				if(labelOriginal.equalsIgnoreCase(labelReordered)){
					indexToReorder[i]=j;
					break;
				}
			}
		}
		
		//reorder
		for(int i=0; i<size; i++){
			for(int j=0; j<size; j++){
				//distancesReordered[i][j] = distancesOriginal[indexToReorder[i]][indexToReorder[j]];
				distancesReordered[indexToReorder[i]][indexToReorder[j]] = distancesOriginal[i][j];
			}
		}
		
		return distancesReordered;
	}
	
	public Object[] readDistancesSamples(String input) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(input));
			String line = br.readLine();
			System.err.println(line);
			int numOfSamples = Integer.parseInt(line.split("[\\s,\\t]+")[0]);
			double[][] distances = new double[numOfSamples][numOfSamples];
			String[] sampleNames = new String[numOfSamples];
			int i = 0;
			while((line=br.readLine())!=null) {
				String[] data = line.split("[\\s,\\t]+");
				sampleNames[i] = data[0].trim();
				for(int j=0; j<numOfSamples; j++) {
					try {
						distances[i][j] = Double.parseDouble(data[j+1]);
					} 
					catch (NumberFormatException e) {
						distances[i][j] = 1.0;
					}
				}
				i++;
			}
			br.close();
			
			Object[] ret = new Object[2];
			ret[0]=distances;
			ret[1]=sampleNames;
			
			return ret;
		} 
		catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
}
