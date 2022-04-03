/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilVCF2ISTATS
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

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.zip.GZIPInputStream;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.GeneralTools;

@Parameters(commandDescription = "VCF2ISTATS")
public class UtilVCF2ISTATS {
	
	public UtilVCF2ISTATS() {}
	
	public static String getUtilName() {
		return "VCF2ISTATS";
	}
	
	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(description = "Input_File", required = true)
	private List<String> inputFileNames = new ArrayList<String>();
	
	@SuppressWarnings("unchecked")
	public void go() {
		try {
			InputStream fis = Files.newInputStream(Paths.get(inputFileNames.get(0)));
			BufferedReader br;
			if(inputFileNames.get(0).endsWith(".gz")) {
				br = new BufferedReader(new InputStreamReader(new GZIPInputStream(fis),"UTF-8"));
			}
			else {
				br = new BufferedReader(new InputStreamReader(fis, "UTF-8"));
			}
			
			int varCounter = 0;
			String line; 
			List<String> headerData = null;
			int numOfInd = 0;
			ArrayList<Object> indivStats = new ArrayList<Object>();
			while((line=br.readLine()) != null) {
				if(line.startsWith("##")) continue;
				else if (line.startsWith("#")) {
					headerData = Arrays.asList(line.split("\\t"));
					numOfInd = headerData.size() - 9;
					System.err.println("Num of Ind = "+numOfInd);
					for(int i=9; i<headerData.size(); i++) {
						ArrayList<Object> indiv = new ArrayList<Object>();
						indiv.add(headerData.get(i));//name
						indiv.add(0);//sites
						indiv.add(0);//het
						indiv.add(0);//alt
						indiv.add(0);//ref
						indiv.add(0);//miss
						indivStats.add(indiv);
					}
					continue;
				}
				else {
					List<String> variantData = Arrays.asList(line.split("\\t"));
					if(variantData.size() < 10) {
						continue;
					}
					String format = variantData.get(8);
					int indexGT = Arrays.asList(format.split(":")).indexOf("GT");
					for(int i=9; i<headerData.size(); i++) {
						ArrayList<Object> indiv = (ArrayList<Object>)indivStats.get(i-9);
						indiv.set(1, (int)indiv.get(1)+1);
						String GT = variantData.get(i).split(":")[indexGT];
						int allele1=-1;
						int allele2=-1;
						try {
							allele1 = Integer.parseInt(GT.split("[/|]")[0]);
							allele2 = Integer.parseInt(GT.split("[/|]")[1]);
						} 
						catch (Exception e) {
							allele1=-1;
							allele2=-1;
						}
						
						if(allele1<0)
							indiv.set(5, (int)indiv.get(5)+1);
						else if(allele1!=allele2)
							indiv.set(2, (int)indiv.get(2)+1);
						else if(allele1==0)
							indiv.set(4, (int)indiv.get(4)+1);
						else
							indiv.set(3, (int)indiv.get(3)+1);
	
					}
					if(++varCounter % 50000 == 0) System.err.println(GeneralTools.time()+" Variants Processed : \t"+varCounter);
				}
			}
			br.close();
			
			NumberFormat formatter = new DecimalFormat("#0.0000", new DecimalFormatSymbols(Locale.US));
			System.out.println("INDIV\tN_SITES\tN_HET\tN_ALT\tN_REF\tN_MISS\tP_HET\tP_ALT\tP_REF\tP_MISS");
			for(int i=0; i<numOfInd;i++) {
				String indivName = (String)((ArrayList<Object>)indivStats.get(i)).get(0);
				int n_sites = (int)((ArrayList<Object>)indivStats.get(i)).get(1);
				int n_het = (int)((ArrayList<Object>)indivStats.get(i)).get(2);
				int n_alt = (int)((ArrayList<Object>)indivStats.get(i)).get(3);
				int n_ref = (int)((ArrayList<Object>)indivStats.get(i)).get(4);
				int n_miss = (int)((ArrayList<Object>)indivStats.get(i)).get(5);
				double p_het = (double)n_het/((double)n_het+n_alt+n_ref);
				double p_alt = (double)n_alt/((double)n_het+n_alt+n_ref);
				double p_ref = (double)n_ref/((double)n_het+n_alt+n_ref);
				double p_miss = (double)n_miss/((double)n_het+n_alt+n_ref+n_miss);
				System.out.println(indivName+"\t"+n_sites+"\t"+n_het+"\t"+n_alt+"\t"+n_ref+"\t"+n_miss+"\t"+
				                   formatter.format(p_het)+"\t"+formatter.format(p_alt)+"\t"+formatter.format(p_ref)+"\t"+formatter.format(p_miss));
			}
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
}
