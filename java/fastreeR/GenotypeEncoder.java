/*
 *
 * BioInfoJava-Utils ciat.agrobio.core.GenotypeEncoder
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

import java.util.HashMap;
import java.util.Map;

public class GenotypeEncoder {
	
	//For haploid
	//[0]: chromosome1
	//[1]: numAlleles
	private static Map<String, byte[]> map1;

	//For diploid
	//[0]: chromosome1
	//[1]: chromosome2
	//[2]: numAlleles
	private static Map<String, byte[]> map2;
	static {
		map1 = new HashMap<String, byte[]>();
		map2 = new HashMap<String, byte[]>();
		
		for(int i=0; i<=7; i++) {
			map1.put(i+"", new byte[] {(byte) (1<<i), (byte)(i+1)});
		}
		map1.put(".", new byte[]{0,0});

		for(int i=0; i<=7; i++) {
			for(int j=0; j<=7; j++) {
				map2.put(i+"/"+j, new byte[] {(byte) (1<<i),(byte) (1<<j), (i>j?(byte)(i+1):(byte)(j+1))});
				map2.put(i+"|"+j, new byte[] {(byte) (1<<i),(byte) (1<<j), (i>j?(byte)(i+1):(byte)(j+1))});
			}
		}
		map2.put("./.", new byte[]{0,0,0});
		map2.put(".|.", new byte[]{0,0,0});
	}
	
	public static byte[] encodeGT(String GT, byte ploidy) {
		if(ploidy==1) {
			byte[] ret =  map1.get(GT);
			if(ret!=null)
				return ret;
			else {
				return new byte[]{0,0};
			}	
		}
		else if(ploidy==2) {
			byte[] ret =  map2.get(GT);
			if(ret!=null)
				return ret;
			else {
				return new byte[]{0,0,0};
			}	
		}
		else {
			System.err.println("Ploidy "+ploidy+" is not yet implemented. Exciting");
			System.exit(ploidy);
			return null;
		}
	}
	
	public static boolean isHet(byte[] encoded, byte ploidy) {
		if(ploidy==1) {
			return false;
		}
		else if(ploidy==2) {
			return encoded[0]!=encoded[1];
		}
		else {
			System.err.println("Ploidy "+ploidy+" is not yet implemented. Exciting");
			System.exit(ploidy);
		}
		return false;
	}
	
	public static boolean isHet(byte dataP1, byte dataP2) {
		return dataP1!=dataP2;
	}
	
	public static boolean isMis(byte[] encoded, byte ploidy) {
		if(ploidy==1) {
			return encoded[0]==0;
		}
		else if(ploidy==2) {
			return encoded[0]==0 && encoded[1]==0;
		}
		else {
			System.err.println("Ploidy "+ploidy+" is not yet implemented. Exciting");
			System.exit(ploidy);
		}
		return false;
	}

	public static boolean isMis(byte dataP1, byte dataP2) {
		return dataP1==0 && dataP2==0;
	}
	
}
