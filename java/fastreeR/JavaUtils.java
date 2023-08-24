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
	
	private GeneralTools gTools = GeneralTools.getInstance();
	
	public void go(String[] args) {
		JCommander jc = null;
		try {
			jc = JCommander.newBuilder()
			
			//jc.addCommand(UtilVCF2ARFF.getUtilName(), new UtilVCF2ARFF());
			//jc.addCommand(UtilVCF2CSV.getUtilName(),  new UtilVCF2CSV());
			//jc.addCommand(UtilVCF2TABLE.getUtilName(),  new UtilVCF2TABLE());
			.addCommand(UtilVCF2ISTATS.getUtilName(),  new UtilVCF2ISTATS()) // in R
			//jc.addCommand(UtilVCF2SVM.getUtilName(), new UtilVCF2SVM());
			.addCommand(UtilVCF2DIST.getUtilName(), new UtilVCF2DIST()) // in R
			.addCommand(UtilVCF2TREE.getUtilName(), new UtilVCF2TREE()) // in R
			//jc.addCommand(UtilVCFRemoveClones.getUtilName(), new UtilVCFRemoveClones()); //
			//jc.addCommand(UtilVCFNonrelated.getUtilName(), new UtilVCFNonrelated()); //
			//jc.addCommand(UtilVCFKeepVariants.getUtilName(), new UtilVCFKeepVariants());
			//jc.addCommand(UtilVCFKeepSamples.getUtilName(), new UtilVCFKeepSamples());
			//jc.addCommand(UtilVCFRemoveSamples.getUtilName(), new UtilVCFRemoveSamples());
			//jc.addCommand(UtilVCFClassifyFromTxtModel.getUtilName(), new UtilVCFClassifyFromTxtModel());
			//jc.addCommand(UtilVCFsIntersection.getUtilName(), new UtilVCFsIntersection()); //
			//jc.addCommand(UtilVCFFilter.getUtilName(), new UtilVCFFilter());
			
			.addCommand(UtilFASTA2DIST.getUtilName(), new UtilFASTA2DIST()) //
			
			.addCommand(UtilDIST2TREE.getUtilName(), new UtilDIST2TREE()) // in R
			.addCommand(UtilDIST2Clusters.getUtilName(), new UtilDIST2Clusters()) // in R
			//jc.addCommand(UtilDIST2Hist.getUtilName(), new UtilDIST2Hist()); // in R
			
			//jc.addCommand(UtilARFFClassifyFromWeka.getUtilName(), new UtilARFFClassifyFromWeka());
			//jc.addCommand(UtilARFFaddTrait.getUtilName(), new UtilARFFaddTrait());
			
			.build();
			jc.parse(args);
			
			Class<?> utilClass = Class.forName(JavaUtils.class.getPackage().getName()+".Util"+jc.getParsedCommand());
			//Method getUtilInstanceMethod = utilClass.getMethod("getInstance");
			Method goMethod = utilClass.getMethod("go");
			Object util = jc.getCommands().get(jc.getParsedCommand()).getObjects().get(0); // getUtilInstanceMethod.invoke(null);

			Field f = utilClass.getDeclaredField("help");
			f.setAccessible(true);
			boolean help = f.getBoolean(util);
			
			if(help) {
				StringBuilder sb = new StringBuilder();
				jc.usage(jc.getParsedCommand(),sb);
				System.err.println(sb.toString());
				System.exit(0);
			}
			
			System.err.println( jc.getParsedCommand() + " : " + new Date(GeneralTools.classBuildTimeMillis(utilClass)).toString());
			goMethod.invoke(util);	
		} 
		catch (MissingCommandException|ClassNotFoundException e) {
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
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		JavaUtils ju = new JavaUtils();
		ju.go(args);
	}
}
