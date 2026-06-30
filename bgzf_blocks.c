/* The MIT License

   Copyright (c) 2026 The bamdst Authors

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

/*
 * bgzf_blocks.c — Block-level API for BGZF files.
 *
 * Provides block-boundary scanning, raw compressed block reading, and
 * stateless single-block decompression — for multi-threaded pipeline
 * processing where a producer reads compressed blocks and workers
 * decompress them independently.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "bgzf.h"
#include "bgzf_blocks.h"

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8

/* ---- little-endian helpers ---- */

static inline int unpackInt16(const uint8_t *buffer)
{
	return buffer[0] | buffer[1] << 8;
}

static inline int unpackInt32(const uint8_t *buffer)
{
	return buffer[0] | buffer[1] << 8 | buffer[2] << 16 | buffer[3] << 24;
}

/* ---- BGZF header check (mirrors bgzf.c check_header) ---- */

static int check_bgzf_header(const uint8_t *header)
{
	return (header[0] == 31 && header[1] == 139 && header[2] == 8
		&& (header[3] & 4) != 0
		&& unpackInt16((uint8_t*)&header[10]) == 6
		&& header[12] == 'B' && header[13] == 'C'
		&& unpackInt16((uint8_t*)&header[14]) == 2);
}

/* ---- public API ---- */

int64_t bgzf_scan_blocks(const char *path, int64_t **offsets)
{
	FILE *fp;
	uint8_t header[BLOCK_HEADER_LENGTH];
	int64_t pos, n = 0, m = 0;
	int64_t *off = NULL;

	fp = fopen(path, "rb");
	if (!fp) return -1;

	while (1) {
		pos = ftello(fp);
		if (fread(header, 1, BLOCK_HEADER_LENGTH, fp) != BLOCK_HEADER_LENGTH)
			break;
		if (!check_bgzf_header(header)) break;

		if (n == m) {
			m = m == 0 ? 4096 : m * 2;
			off = realloc(off, m * sizeof(int64_t));
			if (!off) { fclose(fp); return -1; }
		}
		off[n++] = pos;

		/* skip to next block using BSIZE */
		int block_len = unpackInt16((uint8_t*)&header[16]) + 1;
		if (block_len <= BLOCK_HEADER_LENGTH) break;
		if (fseeko(fp, pos + block_len, SEEK_SET) != 0) break;
	}
	fclose(fp);
	*offsets = off;
	return n;
}

int bgzf_read_raw_block(BGZF *fp, uint8_t **compressed, int *compressed_len)
{
	uint8_t header[BLOCK_HEADER_LENGTH];
	int block_len;

	if (bgzf_read(fp, header, BLOCK_HEADER_LENGTH) != BLOCK_HEADER_LENGTH)
		return -1;
	if (!check_bgzf_header(header)) {
		fp->errcode |= BGZF_ERR_HEADER;
		return -1;
	}
	block_len = unpackInt16((uint8_t*)&header[16]) + 1;

	uint8_t *buf = malloc(block_len);
	if (!buf) return -1;
	memcpy(buf, header, BLOCK_HEADER_LENGTH);

	int remaining = block_len - BLOCK_HEADER_LENGTH;
	if (bgzf_read(fp, buf + BLOCK_HEADER_LENGTH, remaining) != remaining) {
		fp->errcode |= BGZF_ERR_IO;
		free(buf);
		return -1;
	}

	*compressed = buf;
	*compressed_len = block_len;
	return block_len;
}

int bgzf_inflate_raw(const uint8_t *compressed, int compressed_len,
                      uint8_t *uncompressed, int uncompressed_size)
{
	z_stream zs;
	int block_len;

	if (compressed_len < BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH)
		return -1;
	if (!check_bgzf_header(compressed))
		return -1;

	block_len = unpackInt16((uint8_t*)&compressed[16]) + 1;
	if (block_len != compressed_len) return -1;

	memset(&zs, 0, sizeof(z_stream));
	zs.next_in   = (uint8_t*)compressed + BLOCK_HEADER_LENGTH;
	zs.avail_in  = compressed_len - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
	zs.next_out  = uncompressed;
	zs.avail_out = uncompressed_size;

	if (inflateInit2(&zs, -15) != Z_OK) return -1;
	if (inflate(&zs, Z_FINISH) != Z_STREAM_END) {
		inflateEnd(&zs);
		return -1;
	}
	if (inflateEnd(&zs) != Z_OK) return -1;

	return (int)zs.total_out;
}
