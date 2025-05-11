/*
 *
 * BioInfoJava-Utils ciat.agrobio.io.VCFIterator
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

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class VCFIterator<T> implements Iterable<List<T>> {
	private static long CHUNK_SIZE;
	static {
		CHUNK_SIZE = Runtime.getRuntime().maxMemory()/32;
		CHUNK_SIZE = (CHUNK_SIZE/1024)*1024;
		if(CHUNK_SIZE < 8*1024*1024) CHUNK_SIZE = 8*1024*1024;
		if(CHUNK_SIZE > Integer.MAX_VALUE) CHUNK_SIZE = Integer.MAX_VALUE;
		System.err.println("Max RAM: " + Runtime.getRuntime().maxMemory());
		System.err.println("IO Chunk: " + CHUNK_SIZE);
	}
	private final VCFDecoderInterface<T> decoder;
	private Iterator<File> files;

	private VCFIterator(VCFDecoderInterface<T> decoder, File... files) {
		this(decoder, Arrays.asList(files));
	}

	private VCFIterator(VCFDecoderInterface<T> decoder, List<File> files) {
		this.files = files.iterator();
		this.decoder = decoder;
	}

	public static <T> VCFIterator<T> create(VCFDecoderInterface<T> decoder, List<File> files) {
		return new VCFIterator<T>(decoder, files);
	}

	public static <T> VCFIterator<T> create(VCFDecoderInterface<T> decoder, File... files) {
		return new VCFIterator<T>(decoder, files);
	}

	public Iterator<List<T>> iterator() {
		return new Iterator<List<T>>() {
			private List<T> entries;
			private long chunkPos = 0;
			private MappedByteBuffer buffer;
			private FileChannel channel;
			private final ByteOrder byteOrder = java.nio.ByteOrder.nativeOrder();

			public boolean hasNext() {
				if (buffer == null || !buffer.hasRemaining()) {
					buffer = nextBuffer(chunkPos);
					if (buffer == null) {
						return false;
					} else {
						buffer.order(byteOrder);
					}
				}
				T result = null;
				while ((result = decoder.decode(buffer)) != null) {
					if (entries == null) {
						entries = new ArrayList<T>();
					}
					entries.add(result);
				}
				// set next MappedByteBuffer chunk
				chunkPos += buffer.position();
				//try {
					//sun.misc.Cleaner cleaner = ((sun.nio.ch.DirectBuffer) buffer).cleaner();
					//if (cleaner != null) {
		                //cleaner.clean();
		            //}
				//} 
				//catch (Exception e) {System.out.println("sun.misc.Cleaner error");}
				
				buffer = null;
				if (entries != null) {
					return true;
				} else {
					try {
						channel.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
					return false;
				}
			}

			@SuppressWarnings("resource")
			private MappedByteBuffer nextBuffer(long position) {
				try {
					if (channel == null || channel.size() == position) {
						if (channel != null) {
							channel.close();
							channel = null;
						}
						if (files.hasNext()) {
							File file = files.next();
							channel = new RandomAccessFile(file, "r").getChannel();
							chunkPos = 0;
							position = 0;
						} else {
							return null;
						}
					}
					long chunkSize = CHUNK_SIZE;
					if (channel.size() - position < chunkSize) {
						chunkSize = channel.size() - position;
					}
					return channel.map(FileChannel.MapMode.READ_ONLY, chunkPos, chunkSize);
				} catch (IOException e) {
					try {
						channel.close();
					} catch (IOException e1) {
						e1.printStackTrace();
					}
					throw new RuntimeException(e);
				}
			}

			public List<T> next() {
				List<T> res = entries;
				entries = null;
				return res;
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}
}
