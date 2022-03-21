/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.UtilDIST2TREE
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
package fastreeR;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

@Parameters(commandDescription = "DIST2TREE")
public class UtilDIST2TREE {

	private static UtilDIST2TREE instance = new UtilDIST2TREE();

	private UtilDIST2TREE() {
	}

	public static UtilDIST2TREE getInstance() {
		return instance;
	}

	public static String getUtilName() {
		return "DIST2TREE";
	}

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(description = "Input_File", required = true)
	private String inputFileName;
	
	//@Parameter(names = { "--numberOfThreads", "-t" })
	//private int numOfThreads = 1;

	@SuppressWarnings("unused")
	public void go() {
		try {
			//int cpus = Runtime.getRuntime().availableProcessors();
			//int usingThreads = (cpus < numOfThreads ? cpus : numOfThreads);
			//System.err.println("cpus=" + cpus);
			//System.err.println("using=" + usingThreads);

			//Read distances matrix and sample names
			Object[] data = GeneralTools.readDistancesSamples(inputFileName);
			
			//HCluster tree
			String treeString = HierarchicalCluster.hclusteringTree((String[])data[1], (double[][])data[0]);
			System.out.println(treeString);
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}