#include "tic.h"

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

// ====== DEBUG SWITCHES ======
#define TDBG_FORCE_444        0   // force no subsampling
#define TDBG_ONLY_LUMA        0   // encode Y only; set U=V=128 on decode
#define TDBG_DISABLE_FILL     0   // remove FILL; all blocks go DCT
#define TDBG_FORCE_DCT        0   // same as above; belt and braces
#define TDBG_NO_QUANT         0   // bypass quantization (lossless-ish float<->int16)
#define TDBG_NO_ZIGZAG        0   // use raster order
#define TDBG_NO_RLE           0   // write 64 coeffs raw, no RLE/EOB
#define TDBG_VALIDATE_BLOCKS  0   // check block counts per plane
#define TDBG_TRACE_FIRST_N    0   // set >0 to print first N blocks per plane
#define TDBG_SHOW_MODES       0   // 1 = paint blocks: FILL=blue, GRAD=green, GRAD_0=pink, DCT=red
#define USE_RICE              1   // 0 = varint, 1 = Rice-Golomb
#define USE_DEADZONE_QUANT    1   // more aggressive quantisation
#define USE_NEW_CODEC         1   // 0 = old DCT-only codec, 1 = new codec, slower but better quality/compression


// ---- RD tuning ----
#define LAMBDA_MUL        0.003f   // was 0.02f*q; reduce bits penalty
#define GRAD_ERR_GAIN     0.90f    // GRAD must beat DCT by 10%
#define GRAD_ERR_VS_FILL  0.98f    // and beat FILL a bit too
#define GRAD_SLOPE_MAX    12       // clamp |gx|,|gy|

#define MAGIC 0x54444330u /* "TDC0" */
#define BLK 8
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

typedef enum { BM_FILL = 0, BM_DCT = 1, BM_GRAD = 2, BM_GRAD_0 = 3 } BlockMode;
typedef enum { PK_Y = 0, PK_U = 1, PK_V = 2 } PlaneKind;

typedef struct {
	uint32_t magic;
	uint32_t w, h;
	uint32_t subsample; // 0=444,1=420
	uint32_t q;         // quality factor
} Header;

/* ---------- plane storage helpers ---------- */
typedef struct { int w, h, stride; float* p; } Plane;

static void plane_alloc(Plane* pl, int w, int h) {
	pl->w = w; pl->h = h; pl->stride = w;
	pl->p = (float*)malloc(sizeof(float) * w * h);
}
static void plane_free(Plane* pl) { free(pl->p); pl->p = NULL; }


/* ---------- bitstream helpers ---------- */
#if 1
typedef struct {
	uint8_t* data;
	size_t cap, len;
	size_t rd;
} Buf;

static void bw_init(Buf* b) { b->data = NULL; b->cap = 0; b->len = 0; b->rd = 0; }
static void bw_reserve(Buf* b, size_t add) {
	if (b->len + add > b->cap) {
		size_t ncap = MAX(b->cap * 2, b->len + add + 1024);
		b->data = (uint8_t*)realloc(b->data, ncap);
		b->cap = ncap;
	}
}
static void bw_put_u8(Buf* b, uint8_t v) { bw_reserve(b, 1); b->data[b->len++] = v; }
static void bw_put_i16(Buf* b, int16_t v) { bw_reserve(b, 2); b->data[b->len++] = (uint8_t)(v & 0xFF); b->data[b->len++] = (uint8_t)((v >> 8) & 0xFF); }
static void bw_put_u32(Buf* b, uint32_t v) { bw_reserve(b, 4); for (int i = 0; i < 4; i++) bw_put_u8(b, (uint8_t)((v >> (8 * i)) & 0xFF)); }

static uint8_t br_get_u8(Buf* b) { return b->data[b->rd++]; }
static int16_t br_get_i16(Buf* b) { int16_t v = (int16_t)(b->data[b->rd] | (b->data[b->rd + 1] << 8)); b->rd += 2; return v; }
static uint32_t br_get_u32(Buf* b) { uint32_t v = 0; for (int i = 0; i < 4; i++) v |= ((uint32_t)b->data[b->rd++]) << (8 * i); return v; }

// ----- signed LEB128 (zigzag) -----
static void bw_put_varint(Buf* b, int v) {
	uint32_t zz = ((uint32_t)(v << 1)) ^ (uint32_t)((int32_t)v >> 31);
	while (zz >= 0x80) { bw_put_u8(b, (uint8_t)(0x80 | (zz & 0x7F))); zz >>= 7; }
	bw_put_u8(b, (uint8_t)zz);
}
static int br_get_varint(Buf* b) {
	uint32_t x = 0, s = 0;
	for (;;) {
		uint8_t byte = br_get_u8(b);
		x |= (uint32_t)(byte & 0x7F) << s; s += 7;
		if (!(byte & 0x80)) break;
	}
	int v = (int)((x >> 1) ^ (~(x & 1) + 1));  // unzigzag
	return v;
}

typedef struct { Buf* b; uint32_t acc; int n; } BitW;
static void bwb_init(BitW* w, Buf* b) { w->b = b; w->acc = 0; w->n = 0; }
static void bwb_put1(BitW* w, int bit) { w->acc = (w->acc << 1) | (bit & 1); if (++w->n == 8) { bw_put_u8(w->b, (uint8_t)w->acc); w->acc = 0; w->n = 0; } }
static void bwb_putn(BitW* w, uint32_t v, int n) { for (int i = n - 1; i >= 0; i--) bwb_put1(w, (v >> i) & 1); }
static void bwb_flush(BitW* w) { if (w->n) { w->acc <<= (8 - w->n); bw_put_u8(w->b, (uint8_t)w->acc); w->acc = 0; w->n = 0; } }

typedef struct { Buf* b; uint8_t cur; int n; } BitR;
static void brb_init(BitR* r, Buf* b) { r->b = b; r->cur = 0; r->n = 0; }
static int brb_get1(BitR* r) { if (r->n == 0) { r->cur = br_get_u8(r->b); r->n = 8; } int bit = (r->cur >> 7) & 1; r->cur <<= 1; r->n--; return bit; }
static uint32_t brb_getn(BitR* r, int n) { uint32_t v = 0; for (int i = 0; i < n; i++) v = (v << 1) | brb_get1(r); return v; }
static void brb_align(BitR* r) { if (r->n) r->n = 0; }  // drop partial byte

#if USE_RICE

static inline uint32_t zz32(int v) { return ((uint32_t)v << 1) ^ (uint32_t)((int32_t)v >> 31); }
static inline int unzz32(uint32_t z) { return (int)((z >> 1) ^ (~(z & 1) + 1)); }

static void rice_put(BitW* w, int v, int k) {
	uint32_t z = zz32(v);
	uint32_t q = (k ? (z >> k) : z);
	for (uint32_t i = 0; i < q; i++) bwb_put1(w, 0);  // q zeros
	bwb_put1(w, 1);                            // stop bit
	if (k) bwb_putn(w, z & ((1u << k) - 1), k);    // remainder
}
static int rice_get(BitR* r, int k) {
	uint32_t q = 0; while (brb_get1(r) == 0) q++;  // zeros until 1
	uint32_t rem = k ? brb_getn(r, k) : 0;
	uint32_t z = (q << k) | rem;
	return unzz32(z);
}

// simple per-coeff k
static inline int rice_k_dc(void) { return 2; }          // DC delta
static inline int rice_k_ac(int idx) {
	if (idx <= 3) return 0;        // very small
	if (idx <= 15) return 1;
	return 2;
}

#endif

#endif // bitstream region

/* ---------- colour space ---------- */
#if 1
static inline void rgb_to_yuv(uint8_t r, uint8_t g, uint8_t b, float* y, float* u, float* v) {
	// BT.601 full-range-ish
	float R = r, G = g, B = b;
	*y = 0.299f * R + 0.587f * G + 0.114f * B;          // 0..255
	*u = -0.168736f * R - 0.331264f * G + 0.5f * B + 128.0f;
	*v = 0.5f * R - 0.418688f * G - 0.081312f * B + 128.0f;
}
static inline void yuv_to_rgb(float y, float u, float v, uint8_t* r, uint8_t* g, uint8_t* b) {
	float U = u - 128.0f, V = v - 128.0f;
	float R = y + 1.402f * V;
	float G = y - 0.344136f * U - 0.714136f * V;
	float B = y + 1.772f * U;
	int Ri = (int)lrintf(R), Gi = (int)lrintf(G), Bi = (int)lrintf(B);
	*r = (uint8_t)MIN(255, MAX(0, Ri));
	*g = (uint8_t)MIN(255, MAX(0, Gi));
	*b = (uint8_t)MIN(255, MAX(0, Bi));
}

static void rgb_to_yuv_planes(uint8_t* rgb, int w, int h, int ss, Plane* Y, Plane* U, Plane* V) {
	int cw = (ss == TIC_SUBSAMPLING_420) ? (w + 1) / 2 : w;
	int ch = (ss == TIC_SUBSAMPLING_420) ? (h + 1) / 2 : h;
	plane_alloc(Y, w, h); plane_alloc(U, cw, ch); plane_alloc(V, cw, ch);
	// init zeros
	memset(Y->p, 0, sizeof(float) * w * h);
	memset(U->p, 0, sizeof(float) * cw * ch);
	memset(V->p, 0, sizeof(float) * cw * ch);

	if (ss == TIC_SUBSAMPLING_444) {
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				uint8_t r = rgb[3 * (y * w + x) + 0], g = rgb[3 * (y * w + x) + 1], b = rgb[3 * (y * w + x) + 2];
				float fy, fu, fv; rgb_to_yuv(r, g, b, &fy, &fu, &fv);
				Y->p[y * w + x] = fy; U->p[y * w + x] = fu; V->p[y * w + x] = fv;
			}
		}
	}
	else { // 4:2:0 average box 2x2
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				uint8_t r = rgb[3 * (y * w + x) + 0], g = rgb[3 * (y * w + x) + 1], b = rgb[3 * (y * w + x) + 2];
				float fy, fu, fv; rgb_to_yuv(r, g, b, &fy, &fu, &fv);
				Y->p[y * w + x] = fy;
			}
		}
		for (int cy = 0; cy < ch; cy++) {
			for (int cx = 0; cx < cw; cx++) {
				int sx = cx * 2, sy = cy * 2;
				float su = 0, sv = 0; int cnt = 0;
				for (int dy = 0; dy < 2; dy++) for (int dx = 0; dx < 2; dx++) {
					int ix = sx + dx, iy = sy + dy;
					if (ix < w && iy < h) {
						uint8_t r = rgb[3 * (iy * w + ix) + 0], g = rgb[3 * (iy * w + ix) + 1], b = rgb[3 * (iy * w + ix) + 2];
						float fy, fu, fv; rgb_to_yuv(r, g, b, &fy, &fu, &fv);
						su += fu; sv += fv; cnt++;
					}
				}
				U->p[cy * cw + cx] = su / cnt;
				V->p[cy * cw + cx] = sv / cnt;
			}
		}
	}
}

static void yuv_to_rgb_from_planes(uint8_t* rgb, int w, int h, int ss, Plane* Y, Plane* U, Plane* V) {
	if (ss == TIC_SUBSAMPLING_444) {
		for (int y = 0; y < h; y++)for (int x = 0; x < w; x++) {
			float fy = Y->p[y * w + x], fu = U->p[y * w + x], fv = V->p[y * w + x];
			uint8_t r, g, b; yuv_to_rgb(fy, fu, fv, &r, &g, &b);
			rgb[3 * (y * w + x) + 0] = r; rgb[3 * (y * w + x) + 1] = g; rgb[3 * (y * w + x) + 2] = b;
		}
	}
	else { // 4:2:0 nearest upsample
		int cw = (w + 1) / 2, chh = (h + 1) / 2;
		for (int y = 0; y < h; y++)for (int x = 0; x < w; x++) {
			int cx = x / 2, cy = y / 2;
			float fy = Y->p[y * w + x], fu = U->p[cy * cw + cx], fv = V->p[cy * cw + cx];
			uint8_t r, g, b; yuv_to_rgb(fy, fu, fv, &r, &g, &b);
			rgb[3 * (y * w + x) + 0] = r; rgb[3 * (y * w + x) + 1] = g; rgb[3 * (y * w + x) + 2] = b;
		}
	}
}

#endif // colour space region

/* ---------- DCT/IDCT (float) ---------- */
#if 1
static float ctab[BLK][BLK];
static void dct_init() {
	for (int x = 0; x < BLK; x++) for (int k = 0; k < BLK; k++)
		ctab[x][k] = cosf(((2 * x + 1) * k * M_PI) / (2.0f * BLK));
}
static void fdct8x8(const float in[64], float out[64]) {
	// out[u,v] = alpha(u)*alpha(v)/4 * sum_x sum_y in[x,y]*c(x,u)*c(y,v)
	const float a0 = 1.0f / sqrtf(2.0f);
	for (int v = 0; v < BLK; v++) {
		for (int u = 0; u < BLK; u++) {
			float sum = 0.f;
			for (int y = 0; y < BLK; y++)
				for (int x = 0; x < BLK; x++)
					sum += in[y * BLK + x] * ctab[x][u] * ctab[y][v];
			float au = (u == 0) ? a0 : 1.0f;
			float av = (v == 0) ? a0 : 1.0f;
			out[v * BLK + u] = 0.25f * au * av * sum;
		}
	}
}
static void idct8x8(const float in[64], float out[64]) {
	const float a0 = 1.0f / sqrtf(2.0f);
	for (int y = 0; y < BLK; y++) {
		for (int x = 0; x < BLK; x++) {
			float sum = 0.f;
			for (int v = 0; v < BLK; v++) {
				for (int u = 0; u < BLK; u++) {
					float au = (u == 0) ? a0 : 1.0f, av = (v == 0) ? a0 : 1.0f;
					sum += au * av * in[v * BLK + u] * ctab[x][u] * ctab[y][v];
				}
			}
			out[y * BLK + x] = 0.25f * sum;
		}
	}
}
#endif // DCT region

/* ---------- zigzag ---------- */
#if 1
static const uint8_t zigzag[64] = {
	 0, 1, 8,16, 9, 2, 3,10,
	17,24,32,25,18,11, 4, 5,
	12,19,26,33,40,48,41,34,
	27,20,13, 6, 7,14,21,28,
	35,42,49,56,57,50,43,36,
	29,22,15,23,30,37,44,51,
	58,59,52,45,38,31,39,46,
	53,60,61,54,47,55,62,63
};
static void zigzag_f(const float in[64], float out[64]) {
	for (int i = 0; i < 64; i++) out[i] = in[zigzag[i]];
}
static void inv_zigzag_f(const float in[64], float out[64]) {
	for (int i = 0; i < 64; i++) out[zigzag[i]] = in[i];
}
#endif // zigzag region

/* ---------- block analysis ---------- */
#if 1
static int is_flat_block(const float b[64], float thr) {
	// variance threshold
	float mean = 0.f; for (int i = 0; i < 64; i++) mean += b[i]; mean /= 64.f;
	float var = 0.f; for (int i = 0; i < 64; i++) { float d = b[i] - mean; var += d * d; }
	var /= 64.f;
	return var <= thr;
}
#endif // block analysis region

/* ---------- quantisation ---------- */
#if 1
static const uint8_t qtbl[64] = {
	8,6,6,8,10,12,15,16,
	6,6,7,8,12,15,16,16,
	6,7,8,10,12,16,16,16,
	8,8,10,12,15,16,16,16,
	10,12,12,15,16,16,16,16,
	12,15,16,16,16,16,16,16,
	15,16,16,16,16,16,16,16,
	16,16,16,16,16,16,16,16
};

static void quantize(const float in[64], int16_t qout[64], int q) {
	// scalar quantisation by position-weighted step (mildly JPEG-like)
#if USE_DEADZONE_QUANT
	// deadzone quant
	for (int i = 0; i < 64; i++) {
		float step = (qtbl[i] * q) / 16.0f; if (step < 1.0f) step = 1.0f;
		float v = in[i] / step;
		float dz = 0.5f; // dead-zone in steps
		int s = (v > 0) - (v < 0);
		qout[i] = s ? (int16_t)lrintf(fmaxf(0.f, fabsf(v) - dz) * s) : 0;
	}
#else
	for (int i = 0; i < 64; i++) {
		float step = (qtbl[i] * q) / 16.0f;
		if (step < 1.0f) step = 1.0f;
		qout[i] = (int16_t)lrintf(in[i] / step);
	}
#endif
}
static void dequantize(const int16_t in[64], float out[64], int q) {
	for (int i = 0; i < 64; i++) {
		float step = (qtbl[i] * q) / 16.0f;
		if (step < 1.0f) step = 1.0f;
		out[i] = in[i] * step;
	}
}
#endif // quantisation region

/* ---------- gradient analysis ---------- */
#if 1
static void encode_grad_params(const float blk[64], uint8_t* c, int8_t* gx, int8_t* gy) {
	float mean = 0; for (int i = 0; i < 64; i++) mean += blk[i]; mean /= 64.f;

	float sx = 0, sy = 0, sxx = 0, syy = 0;
	for (int y = 0; y < 8; y++) {
		for (int x = 0; x < 8; x++) {
			float dx = x - 3.5f, dy = y - 3.5f;
			float d = blk[y * 8 + x] - mean;
			sx += d * dx;  sy += d * dy;
			sxx += dx * dx; syy += dy * dy;
		}
	}
	// Least-squares slopes (per-pixel change)
	float gx_f = (sxx > 0 ? sx / sxx : 0.f);
	float gy_f = (syy > 0 ? sy / syy : 0.f);

	// Clamp and quantize to int8 with 1.0 step
	if (gx_f >  GRAD_SLOPE_MAX) gx_f =  GRAD_SLOPE_MAX;
	if (gx_f < -GRAD_SLOPE_MAX) gx_f = -GRAD_SLOPE_MAX;
	if (gy_f >  GRAD_SLOPE_MAX) gy_f =  GRAD_SLOPE_MAX;
	if (gy_f < -GRAD_SLOPE_MAX) gy_f = -GRAD_SLOPE_MAX;

	int C  = (int)lrintf(mean); if (C < 0)C = 0; if (C > 255)C = 255;
	int GX = (int)lrintf(gx_f); if (GX < -127)GX = -127; if (GX > 127)GX = 127;
	int GY = (int)lrintf(gy_f); if (GY < -127)GY = -127; if (GY > 127)GY = 127;

	*c = (uint8_t)C; *gx = (int8_t)GX; *gy = (int8_t)GY;
}
static void apply_grad_block(float outblk[64], uint8_t c, int8_t gx, int8_t gy) {
	for (int y = 0; y < 8; y++) for (int x = 0; x < 8; x++)
		outblk[y * 8 + x] = (float)c + (float)gx * (x - 3.5f) + (float)gy * (y - 3.5f);
}

static int varint_len(int v) {
	uint32_t zz = ((uint32_t)v << 1) ^ (uint32_t)((int32_t)v >> 31);
	int n = 1; while (zz >= 0x80) { zz >>= 7; n++; } return n;
}
#endif // gradient analysis region

/* ---------- entropy coding ---------- */
#if 1
	#if USE_RICE
// DC delta in qz[0]; mask = 8 bytes as before
static void write_coeffs(Buf* b, const int16_t qz[64]) {
	// DC
	BitW bw; bwb_init(&bw, b);
	rice_put(&bw, qz[0], rice_k_dc()); bwb_flush(&bw);

	// mask
	uint64_t mask = 0; for (int i = 1; i < 64; i++) if (qz[i] != 0) mask |= 1ull << (i - 1);
	for (int k = 0; k < 8; k++) bw_put_u8(b, (uint8_t)((mask >> (8 * k)) & 0xFF));

	// values
	bwb_init(&bw, b);
	for (int i = 1; i < 64; i++) if (mask & (1ull << (i - 1))) rice_put(&bw, qz[i], rice_k_ac(i));
	bwb_flush(&bw);
}

static int read_coeffs(Buf* b, int16_t qz[64]) {
	memset(qz, 0, 64 * sizeof(int16_t));
	BitR br; brb_init(&br, b);
	qz[0] = (int16_t)rice_get(&br, rice_k_dc()); brb_align(&br);
	uint64_t mask = 0; for (int k = 0; k < 8; k++) mask |= (uint64_t)br_get_u8(b) << (8 * k);

	brb_init(&br, b);
	for (int i = 1; i < 64; i++) if (mask & (1ull << (i - 1))) qz[i] = (int16_t)rice_get(&br, rice_k_ac(i));
	brb_align(&br);
	return 0;
}

	#else // USE_RICE
// eob = highest nonzero AC index in zigzag (1..63). 0 => only DC present.

// DC delta already stored in qz[0]
static void write_coeffs(Buf* b, const int16_t qz[64]) {
	bw_put_varint(b, qz[0]);                 // DC delta

	int eob = 0;
	for (int i = 63; i >= 1; i--) if (qz[i] != 0) { eob = i; break; }
	bw_put_u8(b, (uint8_t)eob);              // 0..63

	int nbytes = (eob + 7) >> 3;                 // ceil(eob/8)
	uint8_t mb[8] = { 0 };
	for (int i = 1; i <= eob; i++) if (qz[i] != 0) mb[(i - 1) >> 3] |= (uint8_t)(1u << ((i - 1) & 7));
	for (int k = 0; k < nbytes; k++) bw_put_u8(b, mb[k]);

	for (int i = 1; i <= eob; i++) if (qz[i] != 0) bw_put_varint(b, qz[i]);
}

static int read_coeffs(Buf* b, int16_t qz[64]) {
	memset(qz, 0, 64 * sizeof(int16_t));
	qz[0] = (int16_t)br_get_varint(b);
	int eob = br_get_u8(b);                  // 0..63
	if (eob < 0 || eob>63) return -1;

	int nbytes = (eob + 7) >> 3;
	uint8_t mb[8] = { 0 };
	for (int k = 0; k < nbytes; k++) mb[k] = br_get_u8(b);

	for (int i = 1; i <= eob; i++) {
		if (mb[(i - 1) >> 3] & (1u << ((i - 1) & 7)))
			qz[i] = (int16_t)br_get_varint(b);
	}
	return 0;
}

	#endif // USE_RICE
#endif // entropy coding region

// old varint encoding
#if 0
// DC delta already stored in qz[0]
static void write_coeffs(Buf* b, const int16_t qz[64]) {
	bw_put_varint(b, qz[0]);                // DC delta
	uint64_t mask = 0; for (int i = 1; i < 64; i++) if (qz[i] != 0) mask |= 1ull << (i - 1);
	for (int k = 0; k < 8; k++) bw_put_u8(b, (uint8_t)((mask >> (8 * k)) & 0xFF));
	for (int i = 1; i < 64; i++) if (mask & (1ull << (i - 1))) bw_put_varint(b, qz[i]);
}
static int read_coeffs(Buf* b, int16_t qz[64]) {
	memset(qz, 0, 64 * sizeof(int16_t));
	qz[0] = (int16_t)br_get_varint(b);      // DC delta
	uint64_t mask = 0; for (int k = 0; k < 8; k++) mask |= (uint64_t)br_get_u8(b) << (8 * k);
	for (int i = 1; i < 64; i++) if (mask & (1ull << (i - 1))) qz[i] = (int16_t)br_get_varint(b);
	return 0;
}
#endif

// old bitmap block coefficients compression
#if 0
// ---- Bitmap coder: DC:int16, MASK:64 bits for AC1..AC63, then only present ACs:int16 ----
static void write_coeffs2(Buf* b, const int16_t qz[64]) {
	// DC
	bw_put_i16(b, qz[0]);
	// MASK: bit (i-1) corresponds to AC index i in zigzag order
	uint64_t mask = 0;
	for (int i = 1; i < 64; i++) if (qz[i] != 0) mask |= (1ull << (i - 1));
	// little-endian 8 bytes
	for (int k = 0; k < 8; k++) bw_put_u8(b, (uint8_t)((mask >> (8 * k)) & 0xFF));
	// payload
	for (int i = 1; i < 64; i++) if (mask & (1ull << (i - 1))) bw_put_i16(b, qz[i]);
}

static int read_coeffs2(Buf* b, int16_t qz[64]) {
	memset(qz, 0, 64 * sizeof(int16_t));
	if (b->rd + 2 > b->len) return -1;
	qz[0] = br_get_i16(b);
	if (b->rd + 8 > b->len) return -2;
	uint64_t mask = 0; for (int k = 0; k < 8; k++) mask |= ((uint64_t)br_get_u8(b)) << (8 * k);
	for (int i = 1; i < 64; i++) {
		if (mask & (1ull << (i - 1))) {
			if (b->rd + 2 > b->len) return -3;
			qz[i] = br_get_i16(b);
		}
	}
	return 0;
}
#endif


/* ---------- format encoding ---------- */
#if 1
	#define NEW_LAMBDA 1
	#define NEW_ENCODER 1

// TODO: remove - used for test
int g_fill_blocks = 0;
int g_grad_blocks = 0;
int g_grad0_blocks = 0;
int g_dct_blocks = 0;

#if NEW_ENCODER
static void encode_plane(Buf* out, Plane* P, int q, float flat_thr, PlaneKind pk) {
	int W = P->w, H = P->h, BW = (W + 7) / 8;
	int16_t* prev_dc_row = (int16_t*)calloc(BW, sizeof(int16_t));
	for (int by = 0, byi = 0; by < H; by += BLK, byi++) {
		int16_t left_dc = 0;
		for (int bx = 0, bxi = 0; bx < W; bx += BLK, bxi++) {
			// gather
			float blk[64];
			for (int y = 0; y < 8; y++) {
				for (int x = 0; x < 8; x++) {
					int ix = bx + x, iy = by + y;
					int cx = ix < 0 ? 0 : (ix >= W ? W - 1 : ix);
					int cy = iy < 0 ? 0 : (iy >= H ? H - 1 : iy);
					blk[y * 8 + x] = P->p[cy * W + cx];
				}
			}

#if NEW_LAMBDA
			// FILL
			float mean = 0; for (int i = 0; i < 64; i++) mean += blk[i]; mean /= 64.f;
			float err_fill = 0; for (int i = 0; i < 64; i++) { float d = blk[i] - mean; err_fill += d * d; }

			// GRAD params and model
			uint8_t gc; int8_t ggx, ggy; encode_grad_params(blk, &gc, &ggx, &ggy);
			float gblk[64]; apply_grad_block(gblk, gc, ggx, ggy);

			// GRAD residual -> DCT path
			float rblk[64]; for (int i = 0; i < 64; i++) rblk[i] = blk[i] - gblk[i];
			float gdct[64]; fdct8x8(rblk, gdct);
			float gz[64];   zigzag_f(gdct, gz);
			int16_t gqz[64]; quantize(gz, gqz, q);

			// DCT baseline
			float dct[64]; fdct8x8(blk, dct);
			float z[64];  zigzag_f(dct, z);
			int16_t qz[64]; quantize(z, qz, q);

			// --- DC prediction (for baseline DCT only) ---
			int pred = (bxi > 0 && byi > 0) ? ((int)left_dc + (int)prev_dc_row[bxi]) >> 1
				: (bxi > 0 ? (int)left_dc : (byi > 0 ? (int)prev_dc_row[bxi] : 0));
			int16_t dc_actual = qz[0];
			qz[0] = (int16_t)(dc_actual - pred);

			// cost: DCT
			int nz = 0, bits_vals = 0;          // exact varint length
			bits_vals += varint_len(qz[0]);
			uint64_t mask = 0;
			for (int i = 1; i < 64; i++) { if (qz[i]) { nz++; mask |= 1ull << (i - 1); bits_vals += varint_len(qz[i]); } }
			float bits_dct = 1.0f + 8.0f + (float)bits_vals;
			float err_dct = 0; {
				float dq[64]; dequantize(qz, dq, q); float idin[64]; inv_zigzag_f(dq, idin);
				float rec[64]; idct8x8(idin, rec);
				for (int i = 0; i < 64; i++) { float d = blk[i] - rec[i]; err_dct += d * d; }
			}

			// cost: GRAD + residual (no DC predictor on residual; it’s near zero anyway)
			int gr_nz = 0, gr_bits_vals = 0;
			gr_bits_vals += varint_len(gqz[0]);
			uint64_t gr_mask = 0;
			for (int i = 1; i < 64; i++) { if (gqz[i]) { gr_nz++; gr_mask |= 1ull << (i - 1); gr_bits_vals += varint_len(gqz[i]); } }
			float bits_grad = 1.0f /*mode*/ + 3.0f /*c,gx,gy*/ + 8.0f /*mask*/ + (float)gr_bits_vals;
			float err_grad = 0; {
				float gdq[64]; dequantize(gqz, gdq, q); float gidin[64]; inv_zigzag_f(gdq, gidin);
				float grec[64]; idct8x8(gidin, grec);    // residual recon
				for (int i = 0; i < 64; i++) { float v = gblk[i] + grec[i]; float d = v - blk[i]; err_grad += d * d; }
			}
			// detect perfect predictor (all-zero quantized residual)
			int gr_all_zero = 1; for (int i = 0; i < 64; i++) if (gqz[i] != 0) { gr_all_zero = 0; break; }

			// FILL bits
			float bits_fill = 2.0f;

			// RD
			float lambda = 0.003f * (float)q;
			float Jf = err_fill + lambda * bits_fill;
			float Jg = err_grad + lambda * bits_grad;
			float Jd = err_dct + lambda * bits_dct;

			BlockMode pick;
#if TDBG_DISABLE_FILL || TDBG_FORCE_DCT
			pick = BM_DCT;
#else
			// require GRAD to actually win against DCT
			if (Jg < Jd * 0.95f && Jg <= Jf) pick = BM_GRAD;
			else if (Jf <= Jd) pick = BM_FILL;
			else pick = BM_DCT;
#endif      

			// --- emit ---
#if TDBG_SHOW_MODES
			{
				// Encode a solid-color block via FILL, colored by chosen mode for THIS plane.
				uint8_t c; dbg_color(pick, pk, &c);
				bw_put_u8(out, (uint8_t)BM_FILL);
				bw_put_u8(out, c);
				// predictor update uses an approximate DC from this solid value
				left_dc = (int16_t)c; prev_dc_row[bxi] = left_dc;
			}
#else

		#define ENABLE_FILL 1
		#define DISABLE_GRAD0 0

#if DISABLE_GRAD0
			pick = BM_GRAD;
			gr_all_zero = 0;
#endif

			if (pick == BM_FILL && ENABLE_FILL) {
				bw_put_u8(out, (uint8_t)BM_FILL);
				bw_put_u8(out, (uint8_t)MIN(255, MAX(0, (int)lrintf(mean))));
				left_dc = (int16_t)lrintf(mean); prev_dc_row[bxi] = left_dc;

				g_fill_blocks = g_fill_blocks + 1;
			}
			else if (pick == BM_GRAD && gr_all_zero) {
				// GRAD block with no residuals - perfect predictor
				bw_put_u8(out, (uint8_t)BM_GRAD_0);
				bw_put_u8(out, gc); bw_put_u8(out, (uint8_t)ggx); bw_put_u8(out, (uint8_t)ggy);
				left_dc = (int16_t)gc; prev_dc_row[bxi] = left_dc;

				g_grad0_blocks = g_grad0_blocks + 1;
			}
			else if (pick == BM_GRAD) {
				bw_put_u8(out, (uint8_t)BM_GRAD);
				bw_put_u8(out, gc); bw_put_u8(out, (uint8_t)ggx); bw_put_u8(out, (uint8_t)ggy);
				// residual payload for GRAD:
				// store as varint-bitmap using existing write_coeffs()
				// but do not poison DC predictor; encode residual with qz[0] as-is (no DPCM)
				write_coeffs(out, gqz);
				left_dc = (int16_t)gc; prev_dc_row[bxi] = left_dc;   // predictor uses mean

				g_grad_blocks = g_grad_blocks + 1;
			}
			else {
				int16_t dc_actual = qz[0];

				bw_put_u8(out, (uint8_t)BM_DCT);
				write_coeffs(out, qz);

				left_dc = dc_actual;
				prev_dc_row[bxi] = left_dc;

				g_dct_blocks = g_dct_blocks + 1;
			}
#endif // TDBG_SHOW_MODES

#else // NEW_LAMBDA

			// FILL
			float mean = 0; for (int i = 0; i < 64; i++) mean += blk[i]; mean /= 64.f;
			float err_fill = 0; for (int i = 0; i < 64; i++) { float d = blk[i] - mean; err_fill += d * d; }

			// GRAD
			uint8_t gc; int8_t ggx, ggy; encode_grad_params(blk, &gc, &ggx, &ggy);
			float gblk[64]; apply_grad_block(gblk, gc, ggx, ggy);
			float err_grad = 0; for (int i = 0; i < 64; i++) { float d = blk[i] - gblk[i]; err_grad += d * d; }

			// DCT
			float dct[64]; fdct8x8(blk, dct);
			float z[64];  zigzag_f(dct, z);
			int16_t qz[64]; quantize(z, qz, q);

			// DC DPCM
			int pred = (bxi > 0 && byi > 0) ? ((int)left_dc + (int)prev_dc_row[bxi]) >> 1
				: (bxi > 0 ? (int)left_dc : (byi > 0 ? (int)prev_dc_row[bxi] : 0));
			int16_t dc_actual = qz[0];
			qz[0] = (int16_t)(dc_actual - pred);

			// cheap RD
			// compute nz and approx bit cost using real varint lengths
			int nz = 0, bits_vals = 0;
			bits_vals += varint_len(qz[0]); // DC delta
			uint64_t mask = 0;
			for (int i = 1; i < 64; i++) {
				if (qz[i]) { nz++; mask |= 1ull << (i - 1); bits_vals += varint_len(qz[i]); }
			}
			float bits_fill = 2.0f;                    // mode + byte
			float bits_grad = 4.0f;                    // mode + c,gx,gy
			float bits_dct = 1.0f + 8.0f + (float)bits_vals; // mode + mask + varints

			float err_dct = 0;
			{
				float dq[64]; dequantize(qz, dq, q); float idin[64]; inv_zigzag_f(dq, idin);
				float rec[64]; idct8x8(idin, rec); for (int i = 0; i < 64; i++) { float d = blk[i] - rec[i]; err_dct += d * d; }
			}

			float lambda = LAMBDA_MUL * (float)q;

			// base decisions
			float Jf = err_fill + lambda * bits_fill;
			float Jg = err_grad + lambda * bits_grad;
			float Jd = err_dct + lambda * bits_dct;

			BlockMode pick = BM_DCT;
#if TDBG_DISABLE_FILL || TDBG_FORCE_DCT
			pick = BM_DCT;
#else
			// require GRAD to clearly beat both
			// if (Jg <= Jd * GRAD_ERR_GAIN && Jg <= Jf * GRAD_ERR_VS_FILL) pick = BM_GRAD; else
			if (Jf <= Jd && Jf <= Jg) pick = BM_FILL;
			else pick = BM_DCT;
#endif


#if TDBG_SHOW_MODES
			{
				// Encode a solid-color block via FILL, colored by chosen mode for THIS plane.
				uint8_t c; dbg_color(pick, pk, &c);
				bw_put_u8(out, (uint8_t)BM_FILL);
				bw_put_u8(out, c);
				// predictor update uses an approximate DC from this solid value
				left_dc = (int16_t)c; prev_dc_row[bxi] = left_dc;
			}
#else
			// emit
			if (pick == BM_FILL) {
				bw_put_u8(out, (uint8_t)BM_FILL);
				bw_put_u8(out, (uint8_t)MIN(255, MAX(0, (int)lrintf(mean))));
				left_dc = (int16_t)lrintf(mean); prev_dc_row[bxi] = left_dc;
			}
			else if (pick == BM_GRAD) {
				bw_put_u8(out, (uint8_t)BM_GRAD);
				bw_put_u8(out, gc); bw_put_u8(out, (uint8_t)ggx); bw_put_u8(out, (uint8_t)ggy);
				left_dc = (int16_t)gc; prev_dc_row[bxi] = left_dc;
			}
			else { // DCT
				bw_put_u8(out, (uint8_t)BM_DCT);
				write_coeffs(out, qz);
				left_dc = dc_actual; prev_dc_row[bxi] = left_dc;
			}
#endif

#endif
		}
	}
	free(prev_dc_row);
}
static void decode_plane(Buf* inb, Plane* P, int q) {
	int W = P->w, H = P->h, BW = (W + 7) / 8;
	int16_t* prev_dc_row = (int16_t*)calloc(BW, sizeof(int16_t));
	// YUV colors for debug: FILL=blue, GRAD=green, DCT=red

	for (int by = 0, byi = 0; by < H; by += BLK, byi++) {
		int16_t left_dc = 0;
		for (int bx = 0, bxi = 0; bx < W; bx += BLK, bxi++) {
			uint8_t mode = br_get_u8(inb);
			float outblk[64] = { 0 };
			if (mode == BM_FILL) {
				uint8_t m = br_get_u8(inb);
				for (int i = 0; i < 64; i++) outblk[i] = (float)m;
				left_dc = (int16_t)m; prev_dc_row[bxi] = left_dc;
			}
#if NEW_LAMBDA
			else if (mode == BM_GRAD_0) {
				uint8_t c = br_get_u8(inb); int8_t gx = (int8_t)br_get_u8(inb); int8_t gy = (int8_t)br_get_u8(inb);
				apply_grad_block(outblk, c, gx, gy);              // no residual to read
				left_dc = (int16_t)c; prev_dc_row[bxi] = left_dc;
			}
			else if (mode == BM_GRAD) {
				uint8_t c = br_get_u8(inb); int8_t gx = (int8_t)br_get_u8(inb); int8_t gy = (int8_t)br_get_u8(inb);
				int16_t gqz[64]; read_coeffs(inb, gqz);                   // residual coeffs
				float gdq[64]; dequantize(gqz, gdq, q); float gidin[64]; inv_zigzag_f(gdq, gidin);
				float grec[64]; idct8x8(gidin, grec);                      // residual recon
				float base[64]; apply_grad_block(base, c, gx, gy);
				for (int i = 0; i < 64; i++) outblk[i] = base[i] + grec[i];      // grad + residual
				left_dc = (int16_t)c; prev_dc_row[bxi] = left_dc;
			}
#else
			else if (mode == BM_GRAD) {
				uint8_t c = br_get_u8(inb); int8_t gx = (int8_t)br_get_u8(inb); int8_t gy = (int8_t)br_get_u8(inb);
				apply_grad_block(outblk, c, gx, gy);
				left_dc = (int16_t)c; prev_dc_row[bxi] = left_dc;
			}
#endif
			else { // BM_DCT
				int16_t qz[64]; read_coeffs(inb, qz);
				int pred = (bxi > 0 && byi > 0) ? ((int)left_dc + (int)prev_dc_row[bxi]) >> 1
					: (bxi > 0 ? (int)left_dc : (byi > 0 ? (int)prev_dc_row[bxi] : 0));
				int16_t dc_actual = (int16_t)(qz[0] + pred);
				qz[0] = dc_actual;

				float zf[64]; dequantize(qz, zf, q);
				float idct_in[64]; inv_zigzag_f(zf, idct_in);
				float sp[64]; idct8x8(idct_in, sp);
				for (int i = 0; i < 64; i++) outblk[i] = sp[i];

				left_dc = dc_actual; prev_dc_row[bxi] = left_dc;
			}
			for (int y = 0; y < 8; y++)for (int x = 0; x < 8; x++) {
				int ix = bx + x, iy = by + y; if (ix < W && iy < H) P->p[iy * W + ix] = outblk[y * 8 + x];
			}
		}
	}
	free(prev_dc_row);
}
#else // NEW_ENCODER

// old original experiment: FILL+DCT only
// faster but worse compression than the new one

static void encode_plane(Buf* out, Plane* P, int q, float flat_thr, PlaneKind pk) {
	int W = P->w, H = P->h;
	// store plane by 8x8 blocks: for each block -> mode byte, payload
	for (int by = 0; by < H; by += BLK) {
		for (int bx = 0; bx < W; bx += BLK) {
			float blk[64] = { 0 };
			for (int y = 0; y < BLK; y++) for (int x = 0; x < BLK; x++) {
				int ix = bx + x, iy = by + y;
				int cx = ix < 0 ? 0 : (ix >= W ? W - 1 : ix);
				int cy = iy < 0 ? 0 : (iy >= H ? H - 1 : iy);
				blk[y * BLK + x] = P->p[cy * W + cx];
			}

#if TDBG_DISABLE_FILL || TDBG_FORCE_DCT
			const int force_dct = 1;
#else
			const int force_dct = 0;
#endif

			if (!force_dct && is_flat_block(blk, flat_thr)) {
				// FILL
				float mean = 0.f; for (int i = 0; i < 64; i++) mean += blk[i]; mean /= 64.f;
				bw_put_u8(out, (uint8_t)BM_FILL);
				bw_put_u8(out, (uint8_t)MIN(255, MAX(0, (int)lrintf(mean))));
			}
			else {
				// DCT
				float dct[64]; fdct8x8(blk, dct);
				float z[64]; zigzag_f(dct, z);
				int16_t qz[64]; quantize(z, qz, q);
				bw_put_u8(out, (uint8_t)BM_DCT);
				write_coeffs(out, qz);
			}
		}
	}
}

static void decode_plane(Buf* inb, Plane* P, int q) {
	int W = P->w, H = P->h;
	for (int by = 0; by < H; by += BLK) {
		for (int bx = 0; bx < W; bx += BLK) {
			uint8_t mode = br_get_u8(inb);
			float outblk[64] = { 0 };
			if (mode == BM_FILL) {
				uint8_t m = br_get_u8(inb);
				for (int i = 0; i < 64; i++) outblk[i] = (float)m;
			}
			else {
				int16_t qz[64]; read_coeffs(inb, qz);
				float zf[64]; dequantize(qz, zf, q);
				float idct_in[64]; inv_zigzag_f(zf, idct_in);
				float sp[64]; idct8x8(idct_in, sp);
				for (int i = 0; i < 64; i++) outblk[i] = sp[i];
			}
			for (int y = 0; y < BLK; y++) for (int x = 0; x < BLK; x++) {
				int ix = bx + x, iy = by + y; if (ix < W && iy < H) P->p[iy * W + ix] = outblk[y * BLK + x];
			}
		}
	}
}

#endif // NEW_ENCODER


#endif 

// ----------- file api ----------- //
#if 1
tic_result_t tic_read_image(const char* filename, void** out_data, int* out_w, int* out_h, tic_format_t* fmt) {
	// read whole file
	FILE* f = fopen(filename, "rb");
	if (!f) { return TIC_FILE_ERR; }

	fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET); // file length
	uint8_t* buf = (uint8_t*)malloc(sz);
	if ((long)fread(buf, 1, sz, f) != sz) { fclose(f); free(buf); return TIC_FILE_ERR; }
	fclose(f);

	Buf bs; bw_init(&bs); bs.data = buf; bs.cap = bs.len = sz; bs.rd = 0;

	// read header
	Header hdr;
	hdr.magic = br_get_u32(&bs); hdr.w = br_get_u32(&bs); hdr.h = br_get_u32(&bs);
	hdr.subsample = br_get_u32(&bs); hdr.q = br_get_u32(&bs);
	if (hdr.magic != MAGIC) { free(buf); return TIC_READ_BAD_MAGIC; }
	if (hdr.subsample > 1) { free(buf); return TIC_READ_BAD_SUBSAMPLING; }

	Plane Y, U, V;
	int w = hdr.w, h = hdr.h; int ss = (int)hdr.subsample;
	int cw = (ss == TIC_SUBSAMPLING_420) ? (w + 1) / 2 : w;
	int ch = (ss == TIC_SUBSAMPLING_420) ? (h + 1) / 2 : h;

	*out_w = w;
	*out_h = h;

	plane_alloc(&Y, w, h); plane_alloc(&U, cw, ch); plane_alloc(&V, cw, ch);
	dct_init();

	decode_plane(&bs, &Y, (int)hdr.q);
#if TDBG_ONLY_LUMA
	// synth chroma
	for (int i = 0; i < cw * ch; i++) { U.p[i] = 128.f; V.p[i] = 128.f; }
#else
	decode_plane(&bs, &U, (int)hdr.q);
	decode_plane(&bs, &V, (int)hdr.q);
#endif

	*out_data = malloc(3 * w * h);
	yuv_to_rgb_from_planes(*out_data, w, h, ss, &Y, &U, &V);

	return TIC_OK;
}

tic_result_t tic_write_image(const char* filename, void* data, int w, int h, int quality, tic_format_t fmt, tic_subsampling_t subsampling) {
	Plane Y, U, V; rgb_to_yuv_planes(data, w, h, subsampling, &Y, &U, &V);

	Buf bs; bw_init(&bs);

	// header
	Header hdr = { MAGIC,(uint32_t)w,(uint32_t)h,(uint32_t)subsampling,(uint32_t)quality };
	bw_put_u32(&bs, hdr.magic); bw_put_u32(&bs, hdr.w); bw_put_u32(&bs, hdr.h);
	bw_put_u32(&bs, hdr.subsample); bw_put_u32(&bs, hdr.q);

	// planes
	dct_init();
	float flatY = 3.0f; // flat thresholds
	float flatC = 2.0f;

	// TODO: R, RGB, RGBA modes

	encode_plane(&bs, &Y, quality, flatY, PK_Y);
#if !TDBG_ONLY_LUMA
	encode_plane(&bs, &U, quality, flatC, PK_U);
	encode_plane(&bs, &V, quality, flatC, PK_V);
#endif

	// write data
	FILE* f = fopen(filename, "wb");
	if (!f) { plane_free(&Y); plane_free(&U); plane_free(&V); free(bs.data); return TIC_FILE_ERR; }
	fwrite(bs.data, 1, bs.len, f);
	fclose(f);

	plane_free(&Y); plane_free(&U); plane_free(&V); free(bs.data);

	return TIC_OK;
}
#endif // file api block