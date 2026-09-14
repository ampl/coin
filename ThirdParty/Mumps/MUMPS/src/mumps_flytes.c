/*
 *
 *  This file is part of MUMPS 5.9.1, released
 *  on Mon Jul 20 09:00:43 UTC 2026
 *
 *
 *  Copyright 1991-2026 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
#include <stdio.h>
#include <stdint.h>
#if defined(USE_AVX512_VBMI)
#include <immintrin.h>
#endif
// NB : mumps_flytes undef __AVX512{F__/VBMI__} flags if USE_AVX512_VBMI is not defined
#include "mumps_flytes.h"
/* this implementation exists to avoid depending on a c++ compiler
 * this is inspired from
 * https://gitlab.com/AntJego/adaptative-precision-blr */
#if defined(NOMIXHPC)
void MUMPS_CALL mumps_flyte_return() {};
#else
#if defined(__x86_64__) || defined(_M_X64)
#define ARCH_X86 1
#elif defined(__i386__) || defined(_M_IX86)
#define ARCH_X86 1
#endif
/*
  we do not want misalignement reports on X86 machines because it works
 */
#if defined(ARCH_X86) && defined(__GNUC__) &&  (__GNUC__ > 3) && defined(__has_attribute)
#if (__has_attribute(no_sanitize))
#define NO_SAN_ALIGN __attribute__((no_sanitize("alignment")))
#else
#define NO_SAN_ALIGN
#endif
#else
#define NO_SAN_ALIGN
#endif
/* 
 * TODO:
 * endian-ness-dependent portability should be achievable with memcpy 
 *
 * the permutation arrays should be updtaed accordingly
 *
 * e.g. current version does not work on big-endian architecture
 */
void MUMPS_CALL avx512_available(MUMPS_INT* f, MUMPS_INT* vbmi) {
	*f = 0;
	*vbmi = 0;
#if defined(__AVX512F__)
	*f = 1;
#endif
#if defined(__AVX512VBMI__)
	*vbmi = 1;
#endif
}
/* The following routines assume column-major format */
/* Real domain compression */
void MUMPS_CALL double_fp56(const double* const in, uint8_t* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_double(7, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL double_fp48(const double* const in, uint8_t* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_double(6, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL double_fp40(const double* const in, uint8_t* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_double(5, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL double_fp24s(const double* const in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_double_float(3, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL double_fp16s(const double* const in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_double_float(2, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL float_fp24s(const float* const in, uint8_t* const out,   const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_float(3, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL float_fp16s(const float* const in, uint8_t* const out,   const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_float(2, in, out, *nrows, *ncols, *ldi, *ldo);
}
/* Real domain decompression */
void MUMPS_CALL fp56_double(const uint8_t* const in, double* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) { 
	decompress_double(7, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL fp48_double(const uint8_t* const in, double* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_double(6, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL fp40_double(const uint8_t* const in, double* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_double(5, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL fp24s_double(const uint8_t* const in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_double_float(3, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL fp16s_double(const uint8_t* const in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_double_float(2, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL fp24s_float(const uint8_t* const in, float* const out,   const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_float(3, in, out, *nrows, *ncols, *ldi, *ldo);
}
void MUMPS_CALL fp16s_float(const uint8_t* const in, float* const out,   const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_float(2, in, out, *nrows, *ncols, *ldi, *ldo);
}
/* Complex domain compression */
void MUMPS_CALL zdouble_fp56z(const double* const in, uint8_t* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_double(7, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL zdouble_fp48z(const double* const in, uint8_t* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_double(6, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL zdouble_fp40z(const double* const in, uint8_t* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_double(5, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL zdouble_fp24c(const double* const in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_double_float(3, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL zdouble_fp16c(const double* const in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_double_float(2, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL cfloat_fp24c(const float* const in, uint8_t* const out,   const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_float(3, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL cfloat_fp16c(const float* const in, uint8_t* const out,   const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	compress_float(2, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
/* Complex domain decompression */
void MUMPS_CALL fp56z_zdouble(const uint8_t* const in, double* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) { 
	decompress_double(7, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL fp48z_zdouble(const uint8_t* const in, double* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_double(6, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL fp40z_zdouble(const uint8_t* const in, double* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_double(5, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL fp24c_zdouble(const uint8_t* const in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_double_float(3, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL fp16c_zdouble(const uint8_t* const in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_double_float(2, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL fp24c_cfloat(const uint8_t* const in, float* const out,   const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_float(3, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
void MUMPS_CALL fp16c_cfloat(const uint8_t* const in, float* const out,   const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo) {
	decompress_float(2, in, out, 2*(*nrows), *ncols, 2*(*ldi), 2*(*ldo));
}
/* CONSIDER whether or not column-major format would be useful */
MUMPS_INLINE void MUMPS_CALL compress_double(uint8_t nb_bytes, const double* const in, uint8_t* const out, const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo) {
	int i;
	if (ldi == ldo && nrows == ldi) compress_vector_double(nb_bytes, in, out, nrows*ncols);
	else
		for (i = 0; i < ncols; i++)
			compress_vector_double(nb_bytes, in + i*ldi, out + i*nb_bytes*ldo, nrows);
}
MUMPS_INLINE void MUMPS_CALL decompress_double(uint8_t nb_bytes, const uint8_t* const in, double* const out, const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo) {
	int i;
	if (ldi == ldo && nrows == ldi) decompress_vector_double(nb_bytes, in, out, nrows*ncols);
	else
		for (i = 0; i < ncols; i++)
			decompress_vector_double(nb_bytes, in + i*nb_bytes*ldi, out + i*ldo, nrows);
}
MUMPS_INLINE void MUMPS_CALL compress_double_float(uint8_t nb_bytes, const double* const in, uint8_t* const out, const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo) {
	int i;
	if (ldi == ldo && nrows == ldi) compress_vector_double_float(nb_bytes, in, out, nrows*ncols);
	else
		for (i = 0; i < ncols; i++)
			compress_vector_double_float(nb_bytes, in + i*ldi, out + i*nb_bytes*ldo, nrows);
}
MUMPS_INLINE void MUMPS_CALL decompress_double_float(uint8_t nb_bytes, const uint8_t* const in, double* const out, const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo) {
	int i;
	if (ldi == ldo && nrows == ldi) decompress_vector_double_float(nb_bytes, in, out, nrows*ncols);
	else
		for (i = 0; i < ncols; i++)
			decompress_vector_double_float(nb_bytes, in + i*nb_bytes*ldi, out + i*ldo, nrows);
}
MUMPS_INLINE void MUMPS_CALL compress_float(uint8_t nb_bytes, const float* const in, uint8_t* const out, const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo) {
	int i;
	if (ldi == ldo && nrows == ldi) compress_vector_float(nb_bytes, in, out, nrows*ncols);
	else
		for (i = 0; i < ncols; i++)
			compress_vector_float(nb_bytes, in + i*ldi, out + i*nb_bytes*ldo, nrows);
}
MUMPS_INLINE void MUMPS_CALL decompress_float(uint8_t nb_bytes, const uint8_t* const in, float* const out, const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo) {
	int i;
	if (ldi == ldo && nrows == ldi) decompress_vector_float(nb_bytes, in, out, nrows*ncols);
	else
		for (i = 0; i < ncols; i++)
			decompress_vector_float(nb_bytes, in + i*nb_bytes*ldi, out + i*ldo, nrows);
}
MUMPS_INLINE void MUMPS_CALL compress_vector_double(const uint8_t nb_bytes, const double* const in, 
  		              uint8_t* const out, const MUMPS_INT length) {
	int start_base = 0;
	int j;
#if defined(__AVX512VBMI__)
	const char* p = get_compress_permutation(nb_bytes);
	const __m512i shuffle = _mm512_set_epi8(
			p[63],p[62],p[61],p[60], p[59],p[58],p[57],p[56],
			p[55],p[54],p[53],p[52], p[51],p[50],p[49],p[48],
			p[47],p[46],p[45],p[44], p[43],p[42],p[41],p[40],
			p[39],p[38],p[37],p[36], p[35],p[34],p[33],p[32],
			p[31],p[30],p[29],p[28], p[27],p[26],p[25],p[24],
			p[23],p[22],p[21],p[20], p[19],p[18],p[17],p[16],
			p[15],p[14],p[13],p[12], p[11],p[10],p[9], p[8],
			p[7], p[6], p[5], p[4],  p[3], p[2], p[1], p[0]
			);
	int offset = 512/(8*sizeof(double));
	start_base = length - length % offset;
	for (j = 0; j < start_base; j+=offset) {
		__m512i vec = _mm512_loadu_si512(in + j);
		vec = _mm512_permutexvar_epi8(shuffle, vec);
		_mm512_mask_storeu_epi8(out + j * nb_bytes, (1ULL << (nb_bytes*offset)) - 1, vec);
	}
#endif
	for (j = start_base; j < length; j++) {
		compress_scalar_double(nb_bytes, &out[j * nb_bytes], in[j]);
	}
}
MUMPS_INLINE void MUMPS_CALL decompress_vector_double(const uint8_t nb_bytes, 
				       const uint8_t* const in, double* const out, const MUMPS_INT length) {
	int start_base = 0;
	int j;
#if defined(__AVX512VBMI__)
	const __mmask64 zero_mask = build_zero_mask(nb_bytes);
	const __mmask32 zero_mask_256 = zero_mask >> 32;
	const __mmask16 zero_mask_128 = zero_mask_256 >> 16;
	const __mmask64 half_mask = build_half_mask(nb_bytes);
	const char* p = get_decompress_permutation(nb_bytes);
	const __m512i shuffle = _mm512_set_epi8(
			p[63],p[62],p[61],p[60], p[59],p[58],p[57],p[56],
			p[55],p[54],p[53],p[52], p[51],p[50],p[49],p[48],
			p[47],p[46],p[45],p[44], p[43],p[42],p[41],p[40],
			p[39],p[38],p[37],p[36], p[35],p[34],p[33],p[32],
			p[31],p[30],p[29],p[28], p[27],p[26],p[25],p[24],
			p[23],p[22],p[21],p[20], p[19],p[18],p[17],p[16],
			p[15],p[14],p[13],p[12], p[11],p[10],p[9], p[8],
			p[7], p[6], p[5], p[4],  p[3], p[2], p[1], p[0]
			);
	const __m512i half = _mm512_maskz_set1_epi8(half_mask, 0x80);
	const __m256i shuffle_256 = _mm512_castsi512_si256(shuffle);
	const __m256i half_256 = _mm256_maskz_set1_epi8(half_mask >> 32, 0x80);
	const __m128i shuffle_128 = _mm256_castsi256_si128(shuffle_256);
	const __m128i half_128 = _mm_maskz_set1_epi8(half_mask >> (32+16), 0x80);
	int offset = 512/(8*sizeof(double));
	int end_base = length - length % offset;
	for (j = start_base; j < end_base; j+=offset) {
		__m512i vec = _mm512_loadu_si512(in + j * nb_bytes);
		vec = _mm512_maskz_permutexvar_epi8(zero_mask, shuffle, vec);
		vec = _mm512_add_epi64(vec, half);
		_mm512_storeu_epi32(out + j, vec);
	}
	start_base = end_base;
	/* this ~complicated stuff is done to squeeze avx as much as possible */
	offset /= 2;
	end_base = length - length % offset;
	for (j = start_base; j < end_base; j+=offset) {
		__m256i vec = _mm256_loadu_si256(in + j * nb_bytes);
		vec = _mm256_maskz_permutexvar_epi8(zero_mask_256, shuffle_256, vec);
		vec = _mm256_add_epi64(vec, half_256);
		_mm256_storeu_epi32(out + j, vec);
	}
	start_base = end_base;
	offset /= 2;
	end_base = length - length % offset;
	for (j = start_base; j < end_base; j+=offset) {
		__m128i vec = _mm_loadu_si128(in + j * nb_bytes);
		vec = _mm_maskz_permutexvar_epi8(zero_mask_128, shuffle_128, vec);
		vec = _mm_add_epi64(vec, half_128);
		_mm_storeu_epi32(out + j, vec);
	}
	start_base = end_base;
#endif
	for (j = start_base; j < length; j++) {
		out[j] = decompress_scalar_double(nb_bytes, &in[j * nb_bytes]);
	}
}
MUMPS_INLINE void MUMPS_CALL compress_vector_double_float(const uint8_t nb_bytes, 
                                    const double* const in, uint8_t* const out, 
                                    const MUMPS_INT length) {
	int start_base = 0;
	int j;
#if defined(__AVX512VBMI__)
	const char* p = get_compress_permutation(nb_bytes);
	const __m256i shuffle = _mm256_set_epi8(
			p[31],p[30],p[29],p[28], p[27],p[26],p[25],p[24],
			p[23],p[22],p[21],p[20], p[19],p[18],p[17],p[16],
			p[15],p[14],p[13],p[12], p[11],p[10],p[9], p[8],
			p[7], p[6], p[5], p[4],  p[3], p[2], p[1], p[0]
			);
	int offset = 512/(8*sizeof(double));
	start_base = length - length % offset;
	for (int j = 0; j < start_base; j+=offset) {
		__m512d vecd = _mm512_loadu_pd(in + j);
		__m256 vecfl = _mm512_cvtpd_ps(vecd);
		__m256i vec = _mm256_permutexvar_epi8(shuffle, _mm256_castps_si256(vecfl));
		_mm256_mask_storeu_epi8(out + j * nb_bytes, (1 << (nb_bytes*offset)) - 1, vec);
	}
#endif
	for (j = start_base; j < length; j++) {
		compress_scalar_double(nb_bytes, &out[j * nb_bytes], in[j]);
	}
}
#if defined(__AVX512F__)
void print_mask64(__mmask64 mask) {
	uint64_t i;
	for (i = 0; i < 64; i++)
		printf("%d", ((mask >> i) & 1) ? 1 : 0);
	printf("\n");
}
void print_mask32(__mmask32 mask) {
	uint64_t i;
	for (i = 0; i < 32; i++)
		printf("%d", ((mask >> i) & 1) ? 1 : 0);
	printf("\n");
}
#endif
MUMPS_INLINE void MUMPS_CALL decompress_vector_double_float(const uint8_t nb_bytes, 
                                      const uint8_t* const in, double* const out, 
                                      const MUMPS_INT length) {
	int start_base = 0;
	int j;
#if defined(__AVX512F__)
#if defined(__AVX512VBMI__)
	const __mmask64 flyte_zero_mask = build_zero_mask_float(nb_bytes);
	const __mmask32 zero_mask = flyte_zero_mask >> 32;
	const __mmask64 half_mask = build_half_mask_float(nb_bytes);
	const char* p = get_decompress_permutation(nb_bytes);
	const __m256i shuffle = _mm256_set_epi8(
			p[31],p[30],p[29],p[28], p[27],p[26],p[25],p[24],
			p[23],p[22],p[21],p[20], p[19],p[18],p[17],p[16],
			p[15],p[14],p[13],p[12], p[11],p[10],p[9], p[8],
			p[7], p[6], p[5], p[4],  p[3], p[2], p[1], p[0]
			);
	const __m256i half = _mm256_maskz_set1_epi8(half_mask >> 32, 0x80);
	int offset = 512/(8*sizeof(double));
	int end_base = length - length % offset;
	for (j = start_base; j < end_base; j+=offset) {
		__m256i vec;
		__m256 vec_ps;
		vec = _mm256_loadu_si256(in + j*nb_bytes);
		vec = _mm256_maskz_permutexvar_epi8(zero_mask, shuffle, vec);
		vec = _mm256_add_epi64(vec, half);
		vec_ps = _mm256_castsi256_ps(vec);
		_mm512_storeu_pd(out + j, _mm512_cvtps_pd(vec_ps));
	}
	start_base = end_base;
	const __mmask32 zero_mask_128 = flyte_zero_mask >> (32+16);
	const __m128i shuffle_128 = _mm256_castsi256_si128(shuffle);
	const __m128i half_128 = _mm_maskz_set1_epi8(half_mask >> (32+16), 0x80);
	offset /= 2;
	end_base = length - length % offset;
	for (j = start_base; j < end_base; j+=offset) {
		__m128i vec;
		__m128 vec_ps;
		vec = _mm_loadu_si128(in + j*nb_bytes);
		vec = _mm_maskz_permutexvar_epi8(zero_mask_128, shuffle_128, vec);
		vec = _mm_add_epi64(vec, half_128);
		vec_ps = _mm_castsi128_ps(vec);
		_mm256_storeu_pd(out + j, _mm256_cvtps_pd(vec_ps));
	}
	start_base = end_base;
#endif
#endif
	for (j = start_base; j < length; j++) {
		out[j] = decompress_scalar_double(nb_bytes, &in[j*nb_bytes]);
	}
}
MUMPS_INLINE void MUMPS_CALL compress_vector_float(const uint8_t nb_bytes, const float* const in, 
  		              uint8_t* const out, const MUMPS_INT length) {
	int start_base = 0;
	int j;
#if defined(__AVX512VBMI__)
	const char* p = get_compress_permutation(nb_bytes);
	const __m512i shuffle = _mm512_set_epi8(
			p[63],p[62],p[61],p[60], p[59],p[58],p[57],p[56],
			p[55],p[54],p[53],p[52], p[51],p[50],p[49],p[48],
			p[47],p[46],p[45],p[44], p[43],p[42],p[41],p[40],
			p[39],p[38],p[37],p[36], p[35],p[34],p[33],p[32],
			p[31],p[30],p[29],p[28], p[27],p[26],p[25],p[24],
			p[23],p[22],p[21],p[20], p[19],p[18],p[17],p[16],
			p[15],p[14],p[13],p[12], p[11],p[10],p[9], p[8],
			p[7], p[6], p[5], p[4],  p[3], p[2], p[1], p[0]
			);
	int offset = 512/(8*sizeof(float));
	start_base = length - length % offset;
	for (j = 0; j < start_base; j+=offset) {
		__m512i vec = _mm512_loadu_si512(in + j);
		vec = _mm512_permutexvar_epi8(shuffle, vec);
		_mm512_mask_storeu_epi8(out + j * nb_bytes, (1ULL << (nb_bytes*offset)) - 1, vec);
	}
#endif
	for (j = start_base; j < length; j++) {
		compress_scalar_float(nb_bytes, &out[j * nb_bytes], in[j]);
	}
}
MUMPS_INLINE void MUMPS_CALL decompress_vector_float(const uint8_t nb_bytes, 
				       const uint8_t* const in, float* const out, const MUMPS_INT length) {
	int start_base = 0;
	int j;
#if defined(__AVX512F__)
#if defined(__AVX512VBMI__)
	const __mmask64 zero_mask = build_zero_mask_float(nb_bytes);
	const __mmask32 zero_mask_256 = zero_mask >> 32;
	const __mmask16 zero_mask_128 = zero_mask_256 >> 16;
	const __mmask64 half_mask = build_half_mask_float(nb_bytes);
	const char* p = get_decompress_permutation(nb_bytes);
	const __m512i shuffle = _mm512_set_epi8(
			p[63],p[62],p[61],p[60], p[59],p[58],p[57],p[56],
			p[55],p[54],p[53],p[52], p[51],p[50],p[49],p[48],
			p[47],p[46],p[45],p[44], p[43],p[42],p[41],p[40],
			p[39],p[38],p[37],p[36], p[35],p[34],p[33],p[32],
			p[31],p[30],p[29],p[28], p[27],p[26],p[25],p[24],
			p[23],p[22],p[21],p[20], p[19],p[18],p[17],p[16],
			p[15],p[14],p[13],p[12], p[11],p[10],p[9], p[8],
			p[7], p[6], p[5], p[4],  p[3], p[2], p[1], p[0]
			);
	const __m512i half = _mm512_maskz_set1_epi8(half_mask, 0x80);
	const __m256i shuffle_256 = _mm512_castsi512_si256(shuffle);
	const __m256i half_256 = _mm256_maskz_set1_epi8(half_mask >> 32, 0x80);
	const __m128i shuffle_128 = _mm256_castsi256_si128(shuffle_256);
	const __m128i half_128 = _mm_maskz_set1_epi8(half_mask >> (32+16), 0x80);
	int offset = 512/(8*sizeof(float));
	int end_base = length - length % offset;
	for (j = start_base; j < end_base; j+=offset) {
		__m512i vec = _mm512_loadu_si512(in + j * nb_bytes);
		vec = _mm512_maskz_permutexvar_epi8(zero_mask, shuffle, vec);
		vec = _mm512_add_epi64(vec, half);
		_mm512_storeu_epi32(out + j, vec);
	}
	start_base = end_base;
	offset /= 2;
	end_base = length - length % offset;
	for (j = start_base; j < end_base; j+=offset) {
		__m256i vec = _mm256_loadu_si256(in + j * nb_bytes);
		vec = _mm256_maskz_permutexvar_epi8(zero_mask_256, shuffle_256, vec);
		vec = _mm256_add_epi64(vec, half_256);
		_mm256_storeu_epi32(out + j, vec);
	}
	start_base = end_base;
	offset /= 2;
	end_base = length - length % offset;
	for (j = start_base; j < end_base; j+=offset) {
		__m128i vec = _mm_loadu_si128(in + j * nb_bytes);
		vec = _mm_maskz_permutexvar_epi8(zero_mask_128, shuffle_128, vec);
		vec = _mm_add_epi64(vec, half_128);
		_mm_storeu_epi32(out + j, vec);
	}
	start_base = end_base;
#endif
#endif
	for (j = start_base; j < length; j++) {
		out[j] = decompress_scalar_float(nb_bytes, &in[j * nb_bytes]);
	}
}
#if defined(__AVX512VBMI__)
MUMPS_INLINE const char* MUMPS_CALL get_compress_permutation(const uint8_t nb_bytes) {
	const char* p = NULL;
	switch (nb_bytes) {
		case 7:
			p = (const char*) FP56_COMPRESS_PERMUTATION;
			break;
		case 6:
			p = (const char*) FP48_COMPRESS_PERMUTATION;
			break;
		case 5:
			p = (const char*) FP40_COMPRESS_PERMUTATION;
			break;
		case 3:
			p = (const char*) FP24_COMPRESS_PERMUTATION;
			break;
		case 2:
			p = (const char*) FP16_COMPRESS_PERMUTATION;
			break;
	}
	return p;
} 
MUMPS_INLINE const char* MUMPS_CALL get_decompress_permutation(const uint8_t nb_bytes) {
	const char* p = NULL;
	switch (nb_bytes) {
		case 7:
			p = (const char*) FP56_DECOMPRESS_PERMUTATION;
			break;
		case 6:
			p = (const char*) FP48_DECOMPRESS_PERMUTATION;
			break;
		case 5:
			p = (const char*) FP40_DECOMPRESS_PERMUTATION;
			break;
		case 3:
			p = (const char*) FP24_DECOMPRESS_PERMUTATION;
			break;
		case 2:
			p = (const char*) FP16_DECOMPRESS_PERMUTATION;
			break;
	}
	return p;
} 
#endif
#if defined(__AVX512F__)
MUMPS_INLINE __mmask64 build_zero_mask(const uint8_t nb_bytes) {
	int i;
	__mmask64 mask = 0;
	uint8_t base = 0;
	base += ((1 << nb_bytes) - 1);
	base <<= sizeof(double) - nb_bytes;
	for (i = 0; i < 512/(8*sizeof(double)); i++) {
		mask += base;
		mask <<= sizeof(double); 
	}
	mask += base;
	return mask;
};
MUMPS_INLINE __mmask64 build_half_mask(const uint8_t nb_bytes) {
	int i;
	__mmask64 mask = 0;
	uint8_t base = 0;
	base += ( 1 << (sizeof(double) - nb_bytes - 1) );
	for (i = 0; i < 512/(8*sizeof(double)); i++) {
		mask += base;
		mask <<= sizeof(double); 
	}
	mask += base;
	return mask;
};
MUMPS_INLINE __mmask64 build_zero_mask_float(const uint8_t nb_bytes) {
	int i;
	__mmask64 mask = 0;
	uint8_t base = 0;
	base += ((1 << nb_bytes) - 1);
	base <<= sizeof(float) - nb_bytes;
	for (i = 0; i < 512/(8*sizeof(float)); i++) {
		mask += base;
		mask <<= sizeof(float); 
	}
	mask += base;
	return mask;
};
MUMPS_INLINE __mmask64 build_half_mask_float(const uint8_t nb_bytes) {
	int i;
	__mmask64 mask = 0;
	uint8_t base = 0;
	base += ( 1 << (sizeof(float) - nb_bytes - 1) );
	for (i = 0; i < 512/(8*sizeof(float)); i++) {
		mask += base;
		mask <<= sizeof(float); 
	}
	mask += base;
	return mask;
};
#endif
MUMPS_INLINE void MUMPS_CALL compress_scalar_double(const uint8_t nb_bytes, uint8_t* const c, const double f) {
	switch (nb_bytes) {
		case 7:
			compress_fp56_scalar(c, f);
			break;
		case 6:
			compress_fp48_scalar(c, f);
			break;
		case 5:
			compress_fp40_scalar(c, f);
			break;
		case 3:
			compress_fp24_scalar(c, (float) f);
			break;
		case 2:
			compress_fp16_scalar(c, (float) f);
			break;
	}
}
MUMPS_INLINE double MUMPS_CALL decompress_scalar_double(const uint8_t nb_bytes, const uint8_t* const c) {
	double f;
	switch (nb_bytes) {
		case 7:
			f = decompress_fp56_scalar(c);
			break;
		case 6:
			f = decompress_fp48_scalar(c);
			break;
		case 5:
			f = decompress_fp40_scalar(c);
			break;
		case 3:
			f = (double) decompress_fp24_scalar(c);
			break;
		case 2:
			f = (double) decompress_fp16_scalar(c);
			break;
	}
	return f;
}
MUMPS_INLINE void MUMPS_CALL compress_scalar_float(const uint8_t nb_bytes, uint8_t* c, const float f) {
	switch (nb_bytes) {
		case 3:
			compress_fp24_scalar(c, f);
			break;
		case 2:
			compress_fp16_scalar(c, f);
			break;
	}
}
MUMPS_INLINE float MUMPS_CALL decompress_scalar_float(const uint8_t nb_bytes, const uint8_t* const c) {
	float f;
	switch (nb_bytes) {
		case 3:
			f = decompress_fp24_scalar(c);
			break;
		case 2:
			f = decompress_fp16_scalar(c);
			break;
	}
	return f;
}
NO_SAN_ALIGN MUMPS_INLINE void MUMPS_CALL compress_fp56_scalar(uint8_t* const c, const double f) {
	union fl_fp64 u64;
	u64.f = f;
	uint32_t *i32 = (uint32_t*) &c[0];
	uint16_t *i16 = (uint16_t*) &c[4];
	*i32 = (uint32_t) (u64.i >> 32);
	*i16 = (uint16_t) (u64.i >> 16);
	c[6] = (uint8_t) (u64.i >> 8);
}
NO_SAN_ALIGN MUMPS_INLINE void MUMPS_CALL compress_fp48_scalar(uint8_t* const c, const double f) {
	union fl_fp64 u64;
	u64.f = f;
	uint32_t *i32 = (uint32_t*) (&c[0]);
	uint16_t *i16 = (uint16_t*) (&c[4]);
	*i32 = (uint32_t) (u64.i >> 32);
	*i16 = (uint16_t) (u64.i >> 16);
}
NO_SAN_ALIGN MUMPS_INLINE void MUMPS_CALL compress_fp40_scalar(uint8_t* const c, const double f) {
	union fl_fp64 u64;
	u64.f = f;
	uint32_t *i32 = (uint32_t*) (&c[0]);
	*i32 = (uint32_t) (u64.i >> 32);
	c[4] = (uint8_t) (u64.i >> 24);
}
NO_SAN_ALIGN MUMPS_INLINE void MUMPS_CALL compress_fp24_scalar(uint8_t* const c, const float f) {
	union fl_fp32 u32;
	u32.f = f;
	uint16_t *i16 = (uint16_t*)(&c[0]);
	*i16 = (uint16_t) (u32.i >> 16);
	c[2] = (uint8_t) (u32.i >> 8);
}
NO_SAN_ALIGN MUMPS_INLINE void MUMPS_CALL compress_fp16_scalar(uint8_t* const c, const float f) {
	union fl_fp32 u32;
	u32.f = f;
	uint16_t *i16 = (uint16_t*) (&c[0]);
	*i16 = (uint16_t) (u32.i >> 16);
}
NO_SAN_ALIGN MUMPS_INLINE double MUMPS_CALL decompress_fp56_scalar(const uint8_t* const c) {
	uint32_t *i32 = (uint32_t*) &c[0];
	uint16_t *i16 = (uint16_t*) &c[4];
	union fl_fp64 u64;
	u64.i = 0x80ULL + (((uint64_t) *i32) << 32) + (((uint64_t) *i16) << 16) + (((uint64_t) c[6]) << 8) ;
	return u64.f;
}
NO_SAN_ALIGN MUMPS_INLINE double MUMPS_CALL decompress_fp48_scalar(const uint8_t* const c) {
	uint32_t *i32 = (uint32_t*) &c[0];
	uint16_t *i16 = (uint16_t*) &c[4];
	union fl_fp64 u64;
	u64.i = 0x8000ULL + (((uint64_t) *i32) << 32) + (((uint64_t) *i16) << 16) ;
	return u64.f;
}
NO_SAN_ALIGN MUMPS_INLINE double MUMPS_CALL decompress_fp40_scalar(const uint8_t* const c) {
	uint32_t *i32 = (uint32_t*) &c[0];
	union fl_fp64 u64;
	u64.i = 0x800000ULL + (((uint64_t) *i32) << 32) + (((uint64_t) c[4]) << 24) ;
	return u64.f;
}
NO_SAN_ALIGN MUMPS_INLINE float MUMPS_CALL decompress_fp24_scalar(const uint8_t* const c) {
	uint16_t *i16 = (uint16_t*)(&c[0]);
	union fl_fp32 u32;
	u32.i = 0x80 +(((uint32_t) *i16) << 16) + (((uint32_t) c[2]) << 8) ;
	return u32.f;
}
NO_SAN_ALIGN MUMPS_INLINE float MUMPS_CALL decompress_fp16_scalar(const uint8_t* const c) {
	uint16_t *i16 = (uint16_t*)(&c[0]);
	union fl_fp32 u32;
	u32.i = 0x8000 +(((uint32_t) *i16) << 16);
	return u32.f;
}
#endif
