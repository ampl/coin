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
/* $Id */
#ifndef MUMPS_FLYTES_H
#define MUMPS_FLYTES_H
#include <stdint.h>
#include "mumps_common.h"
#if !defined(USE_AVX512_VBMI)
#undef __AVX512F__
#undef __AVX512VBMI__
#endif
#if defined(NOMIXHPC)
void MUMPS_CALL mumps_flyte_return();
#else
union fl_fp64 {
	uint64_t i;
	double f;
};
union fl_fp32 {
	uint32_t i;
	float f;
};
#define avx512_available F_SYMBOL(avx512_available,AVX512_AVAILABLE)
#define double_fp56 F_SYMBOL(double_fp56,DOUBLE_FP56)
#define double_fp48 F_SYMBOL(double_fp48,DOUBLE_FP48)
#define double_fp40 F_SYMBOL(double_fp40,DOUBLE_FP40)
#define double_fp24s F_SYMBOL(double_fp24s,DOUBLE_FP24S)
#define double_fp16s F_SYMBOL(double_fp16s,DOUBLE_FP16S)
#define float_fp24s F_SYMBOL(float_fp24s,FLOAT_FP24S)
#define float_fp16s F_SYMBOL(float_fp16s,FLOAT_FP16S)
#define fp56_double F_SYMBOL(fp56_double,FP56_DOUBLE)
#define fp48_double F_SYMBOL(fp48_double,FP48_DOUBLE)
#define fp40_double F_SYMBOL(fp40_double,FP40_DOUBLE)
#define fp24s_double F_SYMBOL(fp24s_double,FP24S_DOUBLE)
#define fp16s_double F_SYMBOL(fp16s_double,FP16S_DOUBLE)
#define fp24s_float F_SYMBOL(fp24s_float,FP24S_FLOAT)
#define fp16s_float F_SYMBOL(fp16s_float,FP16S_FLOAT)
#define zdouble_fp56z F_SYMBOL(zdouble_fp56z,ZDOUBLE_FP56Z)
#define zdouble_fp48z F_SYMBOL(zdouble_fp48z,ZDOUBLE_FP48Z)
#define zdouble_fp40z F_SYMBOL(zdouble_fp40z,ZDOUBLE_FP40Z)
#define zdouble_fp24c F_SYMBOL(zdouble_fp24c,ZDOUBLE_FP24C)
#define zdouble_fp16c F_SYMBOL(zdouble_fp16c,ZDOUBLE_FP16C)
#define cfloat_fp24c F_SYMBOL(cfloat_fp24c,CFLOAT_FP24C)
#define cfloat_fp16c F_SYMBOL(cfloat_fp16c,CFLOAT_FP16C)
#define fp56z_zdouble F_SYMBOL(fp56z_zdouble,FP56Z_ZDOUBLE)
#define fp48z_zdouble F_SYMBOL(fp48z_zdouble,FP48Z_ZDOUBLE)
#define fp40z_zdouble F_SYMBOL(fp40z_zdouble,FP40Z_ZDOUBLE)
#define fp24c_zdouble F_SYMBOL(fp24c_zdouble,FP24C_ZDOUBLE)
#define fp16c_zdouble F_SYMBOL(fp16c_zdouble,FP16C_ZDOUBLE)
#define fp24c_cfloat F_SYMBOL(fp24c_cfloat,FP24C_CFLOAT)
#define fp16c_cfloat F_SYMBOL(fp16c_cfloat,FP16C_CFLOAT)
/* utils */
void MUMPS_CALL avx512_available(MUMPS_INT* f, MUMPS_INT* vbmi);
/* Real domain compression */
void MUMPS_CALL double_fp56(const double* const in,  uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL double_fp48(const double* const in,  uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL double_fp40(const double* const in,  uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL double_fp24s(const double* const in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL double_fp16s(const double* const in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL float_fp24s(const float* const in,   uint8_t* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL float_fp16s(const float* const in,   uint8_t* const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* ldi, const MUMPS_INT* const ldo);
/* Real domain decompression */
void MUMPS_CALL fp56_double(const uint8_t* const in,  double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp48_double(const uint8_t* const in,  double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp40_double(const uint8_t* const in,  double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp24s_double(const uint8_t* const in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp16s_double(const uint8_t* const in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp24s_float(const uint8_t*  const in, float*  const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp16s_float(const uint8_t* const in,  float*  const out,  const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
/* Complex domain compression */
void MUMPS_CALL zdouble_fp56z(const double* in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL zdouble_fp48z(const double* in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL zdouble_fp40z(const double* in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL zdouble_fp24c(const double* in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL zdouble_fp16c(const double* in, uint8_t* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL cfloat_fp24c(const float* in, uint8_t*   const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL cfloat_fp16c(const float* in, uint8_t*   const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
/* Complex domain decompression */
void MUMPS_CALL fp56z_zdouble(const uint8_t* in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp48z_zdouble(const uint8_t* in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp40z_zdouble(const uint8_t* in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp24c_zdouble(const uint8_t* in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp16c_zdouble(const uint8_t* in, double* const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp24c_cfloat(const uint8_t* in, float*   const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
void MUMPS_CALL fp16c_cfloat(const uint8_t* in, float*   const out, const MUMPS_INT* const nrows, const MUMPS_INT* const ncols, const MUMPS_INT* const ldi, const MUMPS_INT* const ldo);
const char FP56_COMPRESS_PERMUTATION[64] = {
		   4, 5, 6, 7, 2, 3, 1, 	12,13,14,15,10,11, 9, 
		  20,21,22,23,18,19,17, 	28,29,30,31,26,27,25, 
		  36,37,38,39,34,35,33, 	44,45,46,47,42,43,41, 
		  52,53,54,55,50,51,49, 	60,61,62,63,58,59,57, 
		 -1,-1,-1,-1,	-1,-1,-1,-1
	};
const char FP56_DECOMPRESS_PERMUTATION[64] = {
		-1, 6, 4, 5,   0, 1, 2, 3,
		-1,13,11,12,   7, 8, 9,10,
		-1,20,18,19,  14,15,16,17,
		-1,27,25,26,  21,22,23,24,
		-1,34,32,33,  28,29,30,31,
		-1,41,39,40,  35,36,37,38,
		-1,48,46,47,  42,43,44,45,
		-1,55,53,54,  49,50,51,52
	};
const char FP48_COMPRESS_PERMUTATION[64] = {
		   4, 5, 6, 7, 2, 3, 	12,13,14,15,10,11, 
		  20,21,22,23,18,19, 	28,29,30,31,26,27, 
		  36,37,38,39,34,35, 	44,45,46,47,42,43, 
		  52,53,54,55,50,51, 	60,61,62,63,58,59, 
		 -1,-1,-1,-1,	-1,-1,-1,-1,
		 -1,-1,-1,-1,	-1,-1,-1,-1
	};
const char FP48_DECOMPRESS_PERMUTATION[64] = { 
		-1,-1, 4, 5,  0, 1, 2, 3, 
		-1,-1,10,11,  6, 7, 8, 9,
		-1,-1,16,17, 12,13,14,15,
		-1,-1,22,23, 18,19,20,21,
		-1,-1,28,29, 24,25,26,27,
		-1,-1,34,35, 30,31,32,33,
		-1,-1,40,41, 36,37,38,39,
		-1,-1,46,47, 42,43,44,45
	};
const char FP40_COMPRESS_PERMUTATION[64] = {
		   4, 5, 6, 7, 3, 	12,13,14,15,11,
		  20,21,22,23,19, 	28,29,30,31,27,
		  36,37,38,39,35, 	44,45,46,47,43,
		  52,53,54,55,51, 	60,61,62,63,59,
		 -1,-1,-1,-1,	-1,-1,-1,-1,
		 -1,-1,-1,-1,	-1,-1,-1,-1,
		 -1,-1,-1,-1,	-1,-1,-1,-1
	};
const char FP40_DECOMPRESS_PERMUTATION[64] = {
		-1,-1,-1, 4, 0, 1, 2, 3, 
		-1,-1,-1, 9, 5, 6, 7, 8, 
		-1,-1,-1,14,10,	11,12,13,
		-1,-1,-1,19,15,	16,17,18,
		-1,-1,-1,24,20,	21,22,23,
		-1,-1,-1,29,25,	26,27,28,
		-1,-1,-1,34,30,	31,32,33,
		-1,-1,-1,39,35,	36,37,38
	};
const char FP24_COMPRESS_PERMUTATION[64] = {
		 2, 3, 1,   6, 7, 5,
		10,11, 9,  14,15,13,
		18,19,17,  22,23,21,
		26,27,25,  30,31,29,
		34,35,33,  38,39,37,
		42,43,41,  46,47,45,
		50,51,49,  54,55,53,
		58,59,57,  62,63,61,
		/* 64 - 9x6 = 10 */	
		-1,-1,-1,-1,-1, -1,-1,-1,-1,-1
	};
const char FP24_DECOMPRESS_PERMUTATION[64] = {
		-1, 2, 0, 1,	-1, 5, 3, 4, 
		-1, 8, 6, 7,	-1,11, 9,10, 
		-1,14,12,13,	-1,17,15,16, 
		-1,20,18,19,	-1,23,21,22, 
		-1,26,24,25,	-1,29,27,28, 
		-1,32,30,31,	-1,35,33,34, 
		-1,38,36,37,	-1,41,39,40, 
		-1,44,42,43, 	-1,47,45,46
	};
const char FP16_COMPRESS_PERMUTATION[64] = {
		 2, 3, 6, 7,	10,11,14,15,
		18,19,22,23,	26,27,30,31,
		34,35,38,39,	42,43,46,47,
		50,51,54,55,	58,59,62,63,
		-1,-1,-1,-1,	-1,-1,-1,-1,
		-1,-1,-1,-1, 	-1,-1,-1,-1,
		-1,-1,-1,-1, 	-1,-1,-1,-1,
		-1,-1,-1,-1, 	-1,-1,-1,-1
	};
const char FP16_DECOMPRESS_PERMUTATION[64] = {
		-1,-1, 0, 1, 	-1,-1, 2, 3,  
		-1,-1, 4, 5, 	-1,-1, 6, 7,  
		-1,-1, 8, 9, 	-1,-1,10,11,  
		-1,-1,12,13, 	-1,-1,14,15,  
		-1,-1,16,17, 	-1,-1,18,19,  
		-1,-1,20,21, 	-1,-1,22,23,  
		-1,-1,24,25, 	-1,-1,26,27,  
		-1,-1,28,29,  	-1,-1,30,31
	};
MUMPS_INLINE void MUMPS_CALL compress_double(          const uint8_t nb_bytes, const double* const in, uint8_t* const out, const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo);
MUMPS_INLINE void MUMPS_CALL decompress_double(        const uint8_t nb_bytes, const uint8_t* const in, double* const out, const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo);
MUMPS_INLINE void MUMPS_CALL compress_double_float(    const uint8_t nb_bytes, const double* const in, uint8_t* const out, const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo);
MUMPS_INLINE void MUMPS_CALL decompress_double_float(  const uint8_t nb_bytes, const uint8_t* const in, double* const out, const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo);
MUMPS_INLINE void MUMPS_CALL compress_float(           const uint8_t nb_bytes, const float* const in, uint8_t*  const out,  const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo);
MUMPS_INLINE void MUMPS_CALL decompress_float(         const uint8_t nb_bytes, const uint8_t* const in, float*  const out,  const MUMPS_INT nrows, const MUMPS_INT ncols, const MUMPS_INT ldi, const MUMPS_INT ldo);
MUMPS_INLINE void MUMPS_CALL compress_vector_double(        const uint8_t nb_bytes, const double* const in, uint8_t* const out, const MUMPS_INT length);
MUMPS_INLINE void MUMPS_CALL decompress_vector_double(      const uint8_t nb_bytes, const uint8_t* const in, double* const out, const MUMPS_INT length);
MUMPS_INLINE void MUMPS_CALL compress_vector_double_float(  const uint8_t nb_bytes, const double* const in, uint8_t* const out, const MUMPS_INT length);
MUMPS_INLINE void MUMPS_CALL decompress_vector_double_float(const uint8_t nb_bytes, const uint8_t* const in, double* const out, const MUMPS_INT length);
MUMPS_INLINE void MUMPS_CALL compress_vector_float(         const uint8_t nb_bytes, const float* const in, uint8_t*  const out, const MUMPS_INT length);
MUMPS_INLINE void MUMPS_CALL decompress_vector_float(       const uint8_t nb_bytes, const uint8_t* const in, float*  const out, const MUMPS_INT length);
MUMPS_INLINE const char* MUMPS_CALL get_compress_permutation(   const uint8_t nb_bytes);
MUMPS_INLINE const char* MUMPS_CALL get_decompress_permutation( const uint8_t nb_bytes);
#if defined(__AVX512F__)
MUMPS_INLINE __mmask64 MUMPS_CALL build_zero_mask(      const uint8_t nb_bytes);
MUMPS_INLINE __mmask64 MUMPS_CALL build_zero_mask_float(const uint8_t nb_bytes);
MUMPS_INLINE __mmask64 MUMPS_CALL build_half_mask(      const uint8_t nb_bytes);
MUMPS_INLINE __mmask64 MUMPS_CALL build_half_mask_float(const uint8_t nb_bytes);
#endif
MUMPS_INLINE void MUMPS_CALL compress_scalar_double(    const uint8_t nb_bytes, uint8_t* const c, const double f);
MUMPS_INLINE double MUMPS_CALL decompress_scalar_double(const uint8_t nb_bytes, const uint8_t* const c);
MUMPS_INLINE void MUMPS_CALL compress_scalar_float(     const uint8_t nb_bytes, uint8_t* const c, const float f);
MUMPS_INLINE float MUMPS_CALL decompress_scalar_float(  const uint8_t nb_bytes, const uint8_t* const c);
MUMPS_INLINE void MUMPS_CALL compress_fp56_scalar(      uint8_t* const c, const double f);
MUMPS_INLINE double MUMPS_CALL decompress_fp56_scalar(  const uint8_t* const c);
MUMPS_INLINE void MUMPS_CALL compress_fp48_scalar(      uint8_t* const c, const double f);
MUMPS_INLINE double MUMPS_CALL decompress_fp48_scalar(  const uint8_t* const c);
MUMPS_INLINE void MUMPS_CALL compress_fp40_scalar(      uint8_t* const c, const double f);
MUMPS_INLINE double MUMPS_CALL decompress_fp40_scalar(  const uint8_t* const c);
MUMPS_INLINE void MUMPS_CALL compress_fp24_scalar(      uint8_t* const c, const float f);
MUMPS_INLINE float MUMPS_CALL decompress_fp24_scalar(   const uint8_t* const c);
MUMPS_INLINE void MUMPS_CALL compress_fp16_scalar(      uint8_t* const c, const float f);
MUMPS_INLINE float MUMPS_CALL decompress_fp16_scalar(   const uint8_t* const c);
#endif
#endif
