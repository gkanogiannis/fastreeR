/*
 *
 * BioInfoJava-Utils ciat.agrobio.io.FastaIterator
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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;

import ciat.agrobio.core.GeneralTools;

public class FastaIterator<T> implements Iterable<List<byte[]>> {
	private static final byte LF = 10;
	private static int CHUNK_SIZE;
	private static ByteBuffer buffer_static;
	static {
		CHUNK_SIZE = (int)Runtime.getRuntime().maxMemory()/32;
		CHUNK_SIZE = (CHUNK_SIZE/1024)*1024;
		if(CHUNK_SIZE < 8*1024*1024) CHUNK_SIZE = 8*1024*1024;
		if(CHUNK_SIZE > Integer.MAX_VALUE) CHUNK_SIZE = Integer.MAX_VALUE;
		System.err.println("Max RAM: " + Runtime.getRuntime().maxMemory());
		System.err.println("IO Chunk: " + CHUNK_SIZE);
		buffer_static = ByteBuffer.allocate(CHUNK_SIZE);
	}
	
	private Iterator<File> files;

	private FastaIterator(File... files) {
		this(Arrays.asList(files));
	}

	private FastaIterator(List<File> files) {
		this.files = files.iterator();
	}

	public static <T> FastaIterator<T> create(List<File> files) {
		return new FastaIterator<T>(files);
	}

	public static <T> FastaIterator<T> create(File... files) {
		return new FastaIterator<T>(files);
	}

	public Iterator<List<byte[]>> iterator() {
		return new Iterator<List<byte[]>>() {
			private int chunkPos = -1;
			private ByteBuffer buffer;
			private BufferedInputStream bis;

			//Return a list of lines of full reads from the fasta/fastq file (possibly multiline)
			//If buffer is finished without reaching the end of the read (haven't seen yet the start of the next read),
			//Then a new buffer is requested until the end of the current read is reached.
			public List<byte[]> next() {
				try {
					//Check that current buffer is active or there is a new buffer.
					if(!hasNext()) {
						return null;
					}
					//List of full reads to return. Each entry is a line for fasta/fastq.
					//In cases that buffer ends 
					List<byte[]> entries = new ArrayList<byte[]>();
					
					boolean last_round = false;
					boolean ready = false;
					int prev_line_start = 0;
					int pos = 0;
					while(hasNext()) {
						while(buffer.hasRemaining()) {
							byte b = buffer.get();
							if (b == LF) { // reached line feed so parse line
								entries.add(Arrays.copyOfRange(buffer.array(), prev_line_start, pos));
								prev_line_start = pos+1;
								chunkPos = -1;
								ready = false;
							}
							else if(!entries.isEmpty() && (b == '>' || b == '@')) { //start of new fasta or fastq read
								chunkPos = pos;
								ready = true;
								if(last_round) break;
							}
							pos++;
						}
						last_round = true;
						if(entries.isEmpty() && pos>0) {//fasta header too long, increase chunk
							CHUNK_SIZE = 2*CHUNK_SIZE;
							if(CHUNK_SIZE > Integer.MAX_VALUE) CHUNK_SIZE = Integer.MAX_VALUE;
							System.err.println("IO Chunk: " + CHUNK_SIZE);
							ByteBuffer temp = ByteBuffer.allocate(CHUNK_SIZE);
							buffer_static.position(0);
							temp.position(CHUNK_SIZE/2);
							temp.put(buffer_static);
							buffer_static = temp;
							buffer = null;
							chunkPos = CHUNK_SIZE/2;
							prev_line_start = 0;
							pos = 0;
							continue;
						}
						if(ready) {
							buffer = null;
							break;
						}
						else {
							if(pos>prev_line_start)
								entries.add(Arrays.copyOfRange(buffer.array(), prev_line_start, pos));
							buffer = null;
							chunkPos = -1;
							prev_line_start = 0;
							pos = 0;
							continue;
						}
					}
					if(entries.isEmpty()) return null;
					return entries;
				} 
				catch (Exception e) {
					e.printStackTrace();
					return null;
				}
			}
			
			//Returns true if current buffer is active(not null) or if the current is inactive and there is a new buffer.
			//Returns false if the current buffer is inactive (null) and there is not anymore new buffers
			public boolean hasNext() {
				return buffer==null? ((buffer=nextBuffer())==null?false:true):true;
			}
			
			//Get a new buffer
			//If there are residual byte not read from the previous one (marked by chunkPos), 
			//copy those old bytes (unused and not returned as a line) to the new buffer and fill the rest with new data
			private ByteBuffer nextBuffer() {
				try {
					if (bis == null) {
						if (files.hasNext()) {
							File f = files.next();
							if(GeneralTools.isGZipped(f)) bis = new BufferedInputStream(new GZIPInputStream(new FileInputStream(f),CHUNK_SIZE));
							else bis = new BufferedInputStream(new FileInputStream(f),CHUNK_SIZE);
							chunkPos = -1;
						} 
						else return null;
					}
					//Recover residual bytes, if any
					int residual = 0;
					if(chunkPos >= 0 && chunkPos<buffer_static.limit()) {
						residual = buffer_static.limit() - chunkPos;
						byte[] b = new byte[residual];
						buffer_static.position(chunkPos);
						buffer_static.get(b);
						buffer_static.position(0);
						buffer_static.put(b);
						chunkPos = -1;
					}
					//Fill rest of buffer
					int read = bis.read(buffer_static.array(), residual, CHUNK_SIZE-residual);
					if(read>=0) buffer_static.limit(residual+read);
					else {
						buffer_static.limit(buffer_static.arrayOffset());
						bis.close();
						bis = null;
					}
					buffer_static.position(0);
					return buffer_static;
				} 
				catch (Exception e) {
					e.printStackTrace();
					return null;
				}
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}
}
