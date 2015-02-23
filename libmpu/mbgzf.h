/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology

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


#ifndef _M_BGZF_H
#define _M_BGZF_H
//#define _mpu_bgzf _mpu__mpu_bgzf
#ifdef _WIN32
#include <windows.h>
#define	ftello	_ftelli64
#define	fseeko	_fseeki64
#endif	//	#ifdef _WIN32

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif

//typedef int8_t bool;

typedef struct {
    int file_descriptor;
    char open_mode;  // 'r' or 'w'
    int16_t owned_file, compress_level;
#ifdef _USE_KNETFILE
	union {
		knetFile *fpr;
		FILE *fpw;
	} x;
#else
    FILE* file;
#endif
    int uncompressed_block_size;
    int compressed_block_size;
    void* uncompressed_block;
    void* compressed_block;
    int64_t block_address;
    int block_length;
    int block_offset;
	int cache_size;
    const char* error;
	void *cache; // a pointer to a hash table
} _mpu_bgzf;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Open an existing file descriptor for reading or writing.
 * Mode must be either "r" or "w".
 * A subsequent _mpu_bgzf_close will not close the file descriptor.
 * Returns null on error.
 */
_mpu_bgzf* _mpu_bgzf_fdopen(int fd, const char* __restrict mode);

/*
 * Open the specified file for reading or writing.
 * Mode must be either "r" or "w".
 * Returns null on error.
 */
_mpu_bgzf* _mpu_bgzf_open(const char* path, const char* __restrict mode);

/*
 * Close the BGZ file and free all associated resources.
 * Does not close the underlying file descriptor if created with _mpu_bgzf_fdopen.
 * Returns zero on success, -1 on error.
 */
int _mpu_bgzf_close(_mpu_bgzf* fp);

/*
 * Read up to length bytes from the file storing into data.
 * Returns the number of bytes actually read.
 * Returns zero on end of file.
 * Returns -1 on error.
 */
int _mpu_bgzf_read(_mpu_bgzf* fp, void* data, int length);

/*
 * Write length bytes from data to the file.
 * Returns the number of bytes written.
 * Returns -1 on error.
 */
int _mpu_bgzf_write(_mpu_bgzf* fp, const void* data, int length);

/*
 * Return a virtual file pointer to the current location in the file.
 * No interpetation of the value should be made, other than a subsequent
 * call to _mpu_bgzf_seek can be used to position the file at the same point.
 * Return value is non-negative on success.
 * Returns -1 on error.
 */
#define _mpu_bgzf_tell(fp) ((fp->block_address << 16) | (fp->block_offset & 0xFFFF))

/*
 * Set the file to read from the location specified by pos, which must
 * be a value previously returned by _mpu_bgzf_tell for this file (but not
 * necessarily one returned by this file handle).
 * The where argument must be SEEK_SET.
 * Seeking on a file opened for write is not supported.
 * Returns zero on success, -1 on error.
 */
int64_t _mpu_bgzf_seek(_mpu_bgzf* fp, int64_t pos, int where);

/*
 * Set the cache size. Zero to disable. By default, caching is
 * disabled. The recommended cache size for frequent random access is
 * about 8M bytes.
 */
void _mpu_bgzf_set_cache_size(_mpu_bgzf *fp, int cache_size);

int _mpu_bgzf_check_EOF(_mpu_bgzf *fp);
int _mpu_bgzf_read_block(_mpu_bgzf* fp);
int _mpu_bgzf_flush(_mpu_bgzf* fp);
int _mpu_bgzf_flush_try(_mpu_bgzf *fp, int size);
int _mpu_bgzf_check__mpu_bgzf(const char *fn);

#ifdef __cplusplus
}
#endif

static inline int _mpu_bgzf_getc(_mpu_bgzf *fp)
{
	int c;
	if (fp->block_offset >= fp->block_length) {
		if (_mpu_bgzf_read_block(fp) != 0) return -2; /* error */
		if (fp->block_length == 0) return -1; /* end-of-file */
	}
	c = ((unsigned char*)fp->uncompressed_block)[fp->block_offset++];
    if (fp->block_offset == fp->block_length) {
#ifdef _USE_KNETFILE
        fp->block_address = knet_tell(fp->x.fpr);
#else
        fp->block_address = ftello(fp->file);
#endif
        fp->block_offset = 0;
        fp->block_length = 0;
    }
	return c;
}

#endif
