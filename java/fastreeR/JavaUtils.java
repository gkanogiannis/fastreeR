/*
 *
 * BioInfoJava-Utils ciat.agrobio.javautils.JavaUtils
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

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.Date;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.MissingCommandException;
import com.beust.jcommander.ParameterException;

import ciat.agrobio.core.GeneralTools;

public class JavaUtils {
	
	public static void main(String[] args) {
		JCommander jc = null;
		try {
			jc = JCommander.newBuilder().build();
			
			//jc.addCommand(UtilVCF2ARFF.getUtilName(), UtilVCF2ARFF.getInstance());
			//jc.addCommand(UtilVCF2CSV.getUtilName(),  UtilVCF2CSV.getInstance());
			//jc.addCommand(UtilVCF2TABLE.getUtilName(),  UtilVCF2TABLE.getInstance());
			jc.addCommand(UtilVCF2ISTATS.getUtilName(),  UtilVCF2ISTATS.getInstance()); // in R
			//jc.addCommand(UtilVCF2SVM.getUtilName(), UtilVCF2SVM.getInstance());
			jc.addCommand(UtilVCF2DIST.getUtilName(), UtilVCF2DIST.getInstance()); // in R
			jc.addCommand(UtilVCF2TREE.getUtilName(), UtilVCF2TREE.getInstance()); // in R
			//jc.addCommand(UtilVCFRemoveClones.getUtilName(), UtilVCFRemoveClones.getInstance()); //
			//jc.addCommand(UtilVCFNonrelated.getUtilName(), UtilVCFNonrelated.getInstance()); //
			//jc.addCommand(UtilVCFKeepVariants.getUtilName(), UtilVCFKeepVariants.getInstance());
			//jc.addCommand(UtilVCFKeepSamples.getUtilName(), UtilVCFKeepSamples.getInstance());
			//jc.addCommand(UtilVCFRemoveSamples.getUtilName(), UtilVCFRemoveSamples.getInstance());
			//jc.addCommand(UtilVCFClassifyFromTxtModel.getUtilName(), UtilVCFClassifyFromTxtModel.getInstance());
			//jc.addCommand(UtilVCFsIntersection.getUtilName(), UtilVCFsIntersection.getInstance()); //
			//jc.addCommand(UtilVCFFilter.getUtilName(), UtilVCFFilter.getInstance());
			
			jc.addCommand(UtilFASTA2DIST.getUtilName(), UtilFASTA2DIST.getInstance()); //
			
			jc.addCommand(UtilDIST2TREE.getUtilName(), UtilDIST2TREE.getInstance()); // in R
			jc.addCommand(UtilDIST2Clusters.getUtilName(), UtilDIST2Clusters.getInstance()); // in R
			jc.addCommand(UtilDIST2Hist.getUtilName(), UtilDIST2Hist.getInstance()); // in R
			
			//jc.addCommand(UtilARFFClassifyFromWeka.getUtilName(), UtilARFFClassifyFromWeka.getInstance());
			//jc.addCommand(UtilARFFaddTrait.getUtilName(), UtilARFFaddTrait.getInstance());
			
			jc.parse(args);
			
			Class<?> utilClass = Class.forName(JavaUtils.class.getPackage().getName()+".Util"+jc.getParsedCommand());
			Method getUtilNameMethod = utilClass.getMethod("getUtilName");
			Method getUtilInstanceMethod = utilClass.getMethod("getInstance");
			Method goMethod = utilClass.getMethod("go");
			Object util = getUtilInstanceMethod.invoke(null);
			String selectedUtilName = (String)getUtilNameMethod.invoke(util);
			Field f = utilClass.getDeclaredField("help");
			f.setAccessible(true);
			boolean help = f.getBoolean(util);
			
			if(help) {
				StringBuilder sb = new StringBuilder();
				jc.usage(jc.getParsedCommand(),sb);
				System.err.println(sb.toString());
				System.exit(0);
			}
			
			System.err.println( selectedUtilName + " : " + new Date(GeneralTools.classBuildTimeMillis(utilClass)).toString());
			goMethod.invoke(util);	
		} 
		catch (MissingCommandException e) {
			System.err.println("Invalid JavaUtil selection!");
		    System.err.println("Use one of : " + jc.getCommands().keySet());
		}
		catch (ParameterException e) {
			System.err.println(e.getMessage());
			StringBuilder sb = new StringBuilder();
			jc.usage(jc.getParsedCommand(),sb);
			System.err.println(sb.toString());
		}
		catch (Exception e) {
			//e.printStackTrace();
			System.err.println("Invalid JavaUtil selection!");
		    System.err.println("Use one of : " + jc.getCommands().keySet());
		    System.exit(0);
		}
	}
}
