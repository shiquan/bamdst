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
 * bgzf_blocks.h — Block-level API for BGZF files.
 *
 * These functions operate on raw BGZF blocks, enabling multi-threaded
 * pipeline processing where a producer thread reads compressed blocks
 * and worker threads decompress and parse them independently.
 *
 * No BGZF handle state is required for decompression (bgzf_inflate_raw
 * is stateless), and bgzf_scan_blocks locates all block boundaries
 * without needing a BAI index.
 */

#ifndef BGZF_BLOCKS_H
#define BGZF_BLOCKS_H

#include <stdint.h>
#include "bgzf.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Scan a BGZF file and return the byte offsets of all blocks.
 * Reads only 18-byte headers; each block is skipped via the BSIZE
 * field in the extra subfield.
 *
 * @param path    file path
 * @param offsets output: array of block byte offsets (caller must free)
 * @return        number of blocks, or -1 on error
 */
int64_t bgzf_scan_blocks(const char *path, int64_t **offsets);

/**
 * Read one raw (compressed) BGZF block from the current file position.
 * The fp must be positioned at a valid BGZF block header.
 * Allocates *compressed; caller must free.
 *
 * @param fp              BGZF file handle (read mode)
 * @param compressed      output: compressed block data (caller frees)
 * @param compressed_len  output: length of compressed data
 * @return                compressed length, or -1 on error
 */
int bgzf_read_raw_block(BGZF *fp, uint8_t **compressed, int *compressed_len);

/**
 * Inflate a raw BGZF block.  Stateless — no BGZF handle required.
 *
 * @param compressed        full block including 18-byte header + 8-byte footer
 * @param compressed_len    total compressed block length
 * @param uncompressed      output buffer (≥ 64KB recommended)
 * @param uncompressed_size size of output buffer
 * @return                  uncompressed data length, or -1 on error
 */
int bgzf_inflate_raw(const uint8_t *compressed, int compressed_len,
                      uint8_t *uncompressed, int uncompressed_size);

#ifdef __cplusplus
}
#endif

#endif
