# Tiny Image Codec (TIC)

Tiny Image Codec (TIC) is an experimental image compression library designed for high efficiency and flexibility. TIC supports multiple compression modes and chroma subsampling, aiming to deliver excellent compression ratios while maintaining image quality.

## Features

- **Multiple Compression Modes**
  - **FILL blocks:** Solid colour blocks for flat regions
  - **GRAD blocks:** Gradient blocks with (GRAD) or without (GRAD0) residuals for smooth transitions
  - **DCT blocks:** Discrete Cosine Transform blocks for complex textures
- **Adaptive Block Selection**
  - Automatically selects the optimal block type per region using a rate-distortion metric
- **Chroma Subsampling**
  - Supports 4:4:4 (no subsampling) and 4:2:0 (reduced chroma resolution)
- **Colour Formats**
  - Handles grayscale (MONO), RGB, and RGBA images
- **Efficient Entropy Coding**
  - Uses Rice-Golomb or varint encoding for coefficients
- **Gradient Prediction**
  - Advanced gradient modeling for improved compression of smooth areas
- **Simple API**
  - Read and write images with a minimal interface:
    - `tic_read_image()`
    - `tic_write_image()`

## Feature Checklist

- [x] Basic 8x8 DCT blocks
- [x] FILL blocks
- [x] Chroma subsampling
- [x] RLE
- [x] Golomb-Rice entropy coding
- [x] Bit-stream
- [x] Quantisation
- [x] Dead zone
- [x] GRAD blocks
- [x] GRAD block residuals
- [x] Rate-Distortion metric
- [x] DPCM
- [x] Library API
- [ ] Quadtree subdivided superblocks (AVIF- and HEIC-like)
- [ ] Optional palette mode
- [ ] V/H/DIAG predictors
- [ ] Perceptual weighting
- [ ] Auto-monochrome mode
- [ ] Adaptive chroma subsampling
- [ ] Intra-blocks from neighbours
- [ ] Wavelet
- [ ] Deblocking filters
- [ ] CDEF

## Efficiency

Comprehensive tests TBD

Current **estimate is 30-40% size reduction** compared to JPEG, but still losing to AVIF

## API Example
```c
#include "tic.h"
void* image_data; int w, h; tic_format_t fmt;
// Read a TIC image
tic_result_t res = tic_read_image("input.tic", &image_data, &w, &h, &fmt);
if (res == TIC_OK) {
    // error
}

// Write a TIC image
tic_result_t res = tic_write_image("output.tic", image_data, w, h, 90, TIC_FMT_RGB, TIC_SUBSAMPLING_420);
if (res == TIC_OK) {
    // error
}
```


## How It Works

TIC divides images into 8x8 blocks and analyzes each block to choose the most space-efficient representation:
- Flat blocks use a single colour.
- Gradient blocks model smooth colour transitions.
- DCT blocks capture detailed textures.
The codec applies quantization and entropy coding to further reduce file size.

## Status

TIC is experimental and intended for research, prototyping, or educational use. The codec is implemented in portable C89.

## License

[Specify your license here]