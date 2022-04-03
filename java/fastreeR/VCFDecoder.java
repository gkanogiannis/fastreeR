/*
 *
 * BioInfoJava-Utils ciat.agrobio.io.VCFDecoder
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

import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;

public class VCFDecoder implements VCFDecoderInterface<byte[][]> {
	private static final byte LF = 10;
	private static final byte TAB = 9;
	private static final byte CR = 13;
	//private boolean skip = false;

	public VCFDecoder() {
	}

	public byte[][] decode(ByteBuffer buffer) {
		int lineStartPos = buffer.position();
		int limit = buffer.limit();
		int crs = 0;
		while (buffer.hasRemaining()) {
			byte b = buffer.get();
			if(b == CR) {
				crs++;
			}
			else if (b == LF) { // reached line feed so parse line
				int lineEndPos = buffer.position();
				// set positions for one row duplication
				if (buffer.limit() < lineEndPos + 1) {
					buffer.position(lineStartPos).limit(lineEndPos);
				} else {
					buffer.position(lineStartPos).limit(lineEndPos + 1);
				}
				
				byte[] line = parseBytes(buffer, lineEndPos - lineStartPos - 1, crs);
				crs = 0;
				
				if (line != null) {
					// reset main buffer
					buffer.position(lineEndPos);
					buffer.limit(limit);
					// set start after LF
					lineStartPos = lineEndPos;
				}
				if(line.length==0){
					continue;
				}
				
				/*
				if(line[0] == '>' || line[0] == '@'){
					skip = false;
				}
				else if(line[0] == '+'){
					skip = true;
				}
				*/
				
				/*
				if(skip){
					continue;
				}
				else{
					return line;
				}
				*/
				
				byte[][] split = splitInTabs(line);
				line = null;
				return split;
			}
		}
		buffer.position(lineStartPos);
		return null;
	}

	private byte[] parseBytes(ByteBuffer buffer, int length, int crs) {
		byte[] bytes = new byte[length-crs];
		for (int i = 0; i < length; i++) {
			byte b = buffer.get();
			if(b != CR)
				bytes[i] = b;
		}
		return bytes;
	}
	
	private byte[][] splitInTabs(byte[] bytes) {
		int numTabs = 0;
		ArrayList<Integer> tabs = new ArrayList<Integer>();
		for(int i = 0; i < bytes.length; i++) {
			if(bytes[i]==TAB) {
				numTabs++;
				tabs.add(i);
			}
		}
		
		byte[][] ret;
		if(numTabs==0 || (bytes[0]=='#' && bytes[1]=='#')) {
			ret = new byte[1][bytes.length];
			ret[0] = Arrays.copyOf(bytes, bytes.length);
		}
		else {
			ret = new byte[numTabs+1][];
			//first
			ret[0] = Arrays.copyOfRange(bytes, 0, tabs.get(0));
			for(int i = 1; i < numTabs; i++) {
				ret[i] = Arrays.copyOfRange(bytes, tabs.get(i-1)+1, tabs.get(i));
			}
			//last
			ret[numTabs] = Arrays.copyOfRange(bytes, tabs.get(numTabs-1)+1, bytes.length);
		}
		tabs.clear();
		tabs = null;
		bytes = null;
		return ret;
	}
}

interface VCFDecoderInterface<T> {

	public T decode(ByteBuffer buffer);

}