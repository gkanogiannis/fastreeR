/*
 *
 * BioInfoJava-Utils ciat.agrobio.hcluster.Clade
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
package ciat.agrobio.hcluster;

import javax.swing.tree.DefaultMutableTreeNode;

public class Clade extends DefaultMutableTreeNode {
	private static final long serialVersionUID = -8392803888823319176L;
	private double branchLength;

    public Clade() {
        this(null, 1);
    }

    public Clade(String nodeLabel, double branchLength) {
    	super(nodeLabel);
    	this.branchLength = branchLength;    
    }

    public double getBranchLength() {
        return branchLength;
    }

    public void setBranchLength(double distance) {
    	if(distance < 0)
    		throw new IllegalArgumentException("Branch length cannot be negative.");
        branchLength = distance;
    }

}
