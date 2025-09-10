#ifndef TIC_H
#define TIC_H

// TIC - Tiny Image Codec
// 
// an experimental image codec with multiple compression modes
// 
// supports chroma subsampling, solid FILL blocks, gradient blocks with residuals (GRAD),
// gradient blocks without residuals (GRAD0), discrete cosine blocks (DCT)
//
// the algorithm will try to pick the smallest possible type of block based on rate-distortion metric

typedef enum {
	TIC_OK,
	TIC_FILE_ERR,
	TIC_READ_BAD_MAGIC,
	TIC_READ_BAD_SUBSAMPLING
} tic_result_t;

typedef enum { TIC_FMT_MONO, TIC_FMT_RGB, TIC_FMT_RGBA } tic_format_t;
typedef enum { TIC_SUBSAMPLING_444, TIC_SUBSAMPLING_420 } tic_subsampling_t;

typedef struct {
	int w, h;
	tic_format_t fmt;
	tic_subsampling_t subsampling;
} tic_info_t;

// ----------- file api ----------- //

tic_result_t tic_read_image(const char* filename, void** out_data, int* w, int* h, tic_format_t* fmt);

tic_result_t tic_write_image(const char* filename, void* data, int w, int h, int quality, tic_format_t fmt, tic_subsampling_t subsampling);


#endif // TIC_H
