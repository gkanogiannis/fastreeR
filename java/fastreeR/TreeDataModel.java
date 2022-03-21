/*
 *
 * BioInfoJava-Utils ciat.agrobio.hcluster.TreeDataModel
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

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;

import javax.swing.tree.DefaultMutableTreeNode;

public class TreeDataModel {
	@SuppressWarnings("rawtypes")
	protected List data = new ArrayList();
    protected DefaultMutableTreeNode root;

    public DefaultMutableTreeNode getRoot() {
        return root;
    }

    @SuppressWarnings({ "rawtypes", "unchecked" })
	public void setRoot(DefaultMutableTreeNode root) { //, boolean isBreadthFirst)
        this.root = root;
        data.clear();

        Enumeration en = root.breadthFirstEnumeration();
        //Enumeration en = isBreadthFirst ? root.breadthFirstEnumeration()
        //                   : root.depthFirstEnumeration();
        while(en.hasMoreElements()) {
            data.add(en.nextElement());
        }
    }

    @SuppressWarnings("unchecked")
	public void add(Object item) {
        data.add(item);
    }
    
    public int size() {
        return data.size();
    }

    public Object get(int i) {
        return data.get(i);
    }

    public int indexOf(Object obj) {
        return data.indexOf(obj);
    }

    @SuppressWarnings("rawtypes")
	public Iterator iterator() {
        return data.iterator();
    }
}
