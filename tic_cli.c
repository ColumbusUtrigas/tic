#include "tic.h"

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

// TODO: better residuals
// TODO: quadtree, superblocks
// TODO: V/H predictors
// TODO: perceptual weighting
// TODO: auto-monochrome mode

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

static void usage(void) {
	fprintf(stderr,
		"Usage:\n"
		"  encode: tic encode in.png out.tic [--q Q] [--ss 444|420]\n"
		"  decode: tic decode in.tic out.png\n");
}

int main(int argc, char** argv) {
	if (argc < 2) { usage(); return 1; }
	if (strcmp(argv[1], "encode") == 0) {
		if (argc < 4) { usage(); return 1; }
		const char* in_png = argv[2]; const char* out_tdc = argv[3];
		int q = 28; tic_subsampling_t ss = TIC_SUBSAMPLING_420;
		for (int i = 4; i < argc; i++) {
			if (strcmp(argv[i], "--q") == 0 && i + 1 < argc) { q = atoi(argv[++i]); q = MAX(1, MIN(99, q)); }
			else if (strcmp(argv[i], "--ss") == 0 && i + 1 < argc) {
				i++; if (strcmp(argv[i], "444") == 0) ss = TIC_SUBSAMPLING_444; else if (strcmp(argv[i], "420") == 0) ss = TIC_SUBSAMPLING_420; else { fprintf(stderr, "bad --ss\n"); return 1; }
			}
		}

		int w, h, n;
		uint8_t* rgb = stbi_load(in_png, &w, &h, &n, 3);
		if (!rgb) { fprintf(stderr, "load fail: %s\n", in_png); return 1; }

		tic_result_t res = tic_write_image(out_tdc, rgb, w, h, q, TIC_FMT_RGB, ss);
		if (res != TIC_OK) {
			fprintf(stderr, "tic write fail: %i\n", res);
			return 1;
		}
	}
	else if (strcmp(argv[1], "decode") == 0) {
		if (argc < 4) { usage(); return 1; }

		int w, h;
		void* data;
		tic_format_t fmt;
		tic_result_t res = tic_read_image(argv[2], &data, &w, &h, &fmt);
		if (res != TIC_OK) {
			fprintf(stderr, "tic read fail: %i\n", res);
			return 1;
		}

		printf("%ix%i\n", w, h);

		int ok = stbi_write_png(argv[3], w, h, 3, data, w * 3);
		if (!ok) fprintf(stderr, "png write fail\n");
	}
	else {
		usage(); return 1;
	}
}