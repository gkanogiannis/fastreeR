/*
 *
 * BioInfoJava-Utils ciat.agrobio.io.FastaDecoder
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

public class FastaDecoder implements FastaFileDecoderInterface<byte[]> {
	private static final byte LF = 10;
	//private boolean skip = false;

	public FastaDecoder() {
	}

	public byte[] decode(ByteBuffer buffer) {
		int lineStartPos = buffer.position();
		int limit = buffer.limit();
		while (buffer.hasRemaining()) {
			byte b = buffer.get();
			if (b == LF) { // reached line feed so parse line
				int lineEndPos = buffer.position();
				// set positions for one row duplication
				if (buffer.limit() < lineEndPos + 1) {
					buffer.position(lineStartPos).limit(lineEndPos);
				} else {
					buffer.position(lineStartPos).limit(lineEndPos + 1);
				}
				
				byte[] line = parseBytes(buffer, lineEndPos - lineStartPos - 1);
				
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
				
				return line;
			}
		}
		buffer.position(lineStartPos);
		return null;
	}

	private byte[] parseBytes(ByteBuffer buffer, int length) {
		byte[] bytes = new byte[length];
		for (int i = 0; i < bytes.length; i++) {
			bytes[i] = buffer.get();
		}
		return bytes;
	}
}

interface FastaFileDecoderInterface<T> {

	public T decode(ByteBuffer buffer);

}