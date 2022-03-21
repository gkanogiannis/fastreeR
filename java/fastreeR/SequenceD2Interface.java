/*
 *
 * BioInfoJava-Utils ciat.agrobio.io.SequenceD2Interface
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
package ciat.agrobio.io;

import gnu.trove.iterator.TLongDoubleIterator;
import gnu.trove.iterator.TLongIntIterator;

public interface SequenceD2Interface {
	
	public double  getDoubleCountForKmerCode(long kmerCode);
	public double  getDoubleProbForKmerCode(long kmerCode);
	
	public TLongIntIterator    iteratorCounts();
	public TLongDoubleIterator iteratorProbs();
	
	public long getTotalATCG();
	public long getTotalCounts();
	
	public long getNumOfElements();
	
	public long getAs();
	public long getTs();
	public long getCs();
	public long getGs();
}
