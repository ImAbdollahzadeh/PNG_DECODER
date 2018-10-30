#include <iostream>
#include <windows.h>
#include "zlib.h"

#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996)

#define INIT_WINDOW_SIZE 10000000
#define WINDOW_EDGE      20
#define ANCILLARY
#define BUTTON_ID        0x321
#define TEXT_ID          0x322

typedef unsigned char  BYTE;
typedef unsigned short WORD;
typedef unsigned int   DOUBLEWORD;

BYTE  pre_decoded[INIT_WINDOW_SIZE];
BYTE  pre_encoded[INIT_WINDOW_SIZE]; 
BYTE* fin_decoded    = nullptr;
int   chunks_counter = 0;
int   idat_counter   = 0;
int   tEXt_counter   = 0;
int   zTXt_counter   = 0;

BYTE IHDR_MARKER[4] = { 'I', 'H', 'D', 'R' };
BYTE IEND_MARKER[4] = { 'I', 'E', 'N', 'D' };
BYTE IDAT_MARKER[4] = { 'I', 'D', 'A', 'T' };
BYTE PLTE_MARKER[4] = { 'P', 'L', 'T', 'E' };
BYTE sRGB_MARKER[4] = { 's', 'R', 'G', 'B' };
BYTE bKGD_MARKER[4] = { 'b', 'K', 'G', 'D' };
BYTE cHRM_MARKER[4] = { 'c', 'H', 'R', 'M' };
BYTE gAMA_MARKER[4] = { 'g', 'A', 'M', 'A' };
BYTE sBIT_MARKER[4] = { 's', 'B', 'I', 'T' };
BYTE zTXt_MARKER[4] = { 'z', 'T', 'X', 'T' };
BYTE tEXt_MARKER[4] = { 't', 'E', 'X', 't' };
BYTE tIME_MARKER[4] = { 't', 'I', 'M', 'E' };
BYTE pHYs_MARKER[4] = { 'p', 'H', 'Y', 's' };
BYTE hIST_MARKER[4] = { 'h', 'I', 'S', 'T' };

namespace Mathematics {
	int abs(int num) {
		return num < 0 ? -num : num;
	}
}

typedef enum {
	IHDR_ID = 0x00,
	IEND_ID = 0x01,
	IDAT_ID = 0x02,
	PLTE_ID = 0x03,
	sRGB_ID = 0x04,
	bKGD_ID = 0x05,
	cHRM_ID = 0x06,
	gAMA_ID = 0x07,
	sBIT_ID = 0x08,
	zTXt_ID = 0x09,
	tEXt_ID = 0x0A,
	tIME_ID = 0x0B,
	pHYs_ID = 0x0C,
	hIST_ID = 0x0D
} CHUNK_ID;

typedef enum {
	RGB_PURE    = 0x0002,
	RGB_ALPHA   = 0x0006,
	RGB_UNKNOWN = 0x00A0
};

typedef struct __Signature {
	BYTE signature[8];
} Signature;

#pragma pack(push, 1)
typedef struct __IHDR {
	DOUBLEWORD Width;
	DOUBLEWORD Height;
	BYTE       Bit_depth;
	BYTE       Color_type;
	BYTE       Compression_method;
	BYTE       Filter_method;
	BYTE       Interlace_method;
} IHDR;
#pragma pack(pop)

typedef struct __IEND {
	bool if_ended;
} IEND;

#pragma pack(push, 1)
typedef struct __IDAT {
	DOUBLEWORD total_pixel_size;
	BYTE       Compression_method_from_ihdr;
	BYTE       Filter_method_from_ihdr;
	BYTE       Bit_depth_from_ihdr;
	BYTE       Color_type_from_ihdr;
	BYTE       Interlace_method_from_ihdr;
	BYTE*      data;
} IDAT;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct __PLTE {
	BYTE red;
	BYTE green;
	BYTE blue;
} PLTE;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct ANCILLARY __sRGB {

} sRGB;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct ANCILLARY __bKGD {
	BYTE Palette_index;
	WORD Gray;
	WORD Red;
	WORD Green;
	WORD Blue;
} bKGD;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct ANCILLARY __cHRM {
	DOUBLEWORD White_Point_x;
	DOUBLEWORD White_Point_y;
	DOUBLEWORD Red_x;
	DOUBLEWORD Red_y;
	DOUBLEWORD Green_x;
	DOUBLEWORD Green_y;
	DOUBLEWORD Blue_x;
	DOUBLEWORD Blue_y;
} cHRM;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct ANCILLARY __gAMA {
	DOUBLEWORD Image_gamma;
} gAMA;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct ANCILLARY __sBIT {
	BYTE  color_type;
	BYTE* significant_bits;
} sBIT;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct ANCILLARY __zTXt {
	BYTE* Keyword;
	BYTE  Null_separator;
	BYTE  Compression_method;
	BYTE* Compressed_text;
} zTXt;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct ANCILLARY __tEXt {
	BYTE* Keyword;
	BYTE  Null_separator;
	BYTE* Text;
} tEXt;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct ANCILLARY __tIME {
	WORD Year;
	BYTE Month;
	BYTE Day;
	BYTE Hour;
	BYTE Minute;
	BYTE Second;
} tIME;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct ANCILLARY __pHYs {
	DOUBLEWORD Pixels_per_unit_X;
	DOUBLEWORD Pixels_per_unit_Y;
	BYTE       Unit_specifier;
} pHYs;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct ANCILLARY __tRNS {
	BYTE* alphas;
} tRNS;
#pragma pack(pop)

typedef struct ANCILLARY __hIST {
	WORD* serie;
} hIST;

#pragma pack(push, 1)
typedef struct __Chunk {
	DOUBLEWORD Length;
	DOUBLEWORD Chunk_Type;
	BYTE*      Chunk_Data;
	DOUBLEWORD CRC;
	CHUNK_ID   id;
} Chunk;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct __IPng {
	const char* address;
	Signature   signature;
	IHDR		ihdr;
	PLTE		plte;
	IDAT*		idat;
	sRGB		srgb;
	bKGD		bkgd;
	cHRM		chrm;
	gAMA		gama;
	sBIT		sbit;
	tIME		time;
	pHYs		phys;
	hIST		hist;
	zTXt*		ztxt;
	tEXt*       text;
	BYTE		number_of_chunks;
	Chunk*		chunks;
	IEND        iend;
} IPng;
#pragma pack(pop)

int inflate(const void* src,
	        int         srcLen,
	        void*       dst,
	        int         dstLen)
{
	z_stream strm  = { nullptr };
	strm.total_in  = strm.avail_in = srcLen;
	strm.total_out = strm.avail_out = dstLen;
	strm.next_in   = (Bytef*)src;
	strm.next_out  = (Bytef*)dst;
	strm.zalloc    = Z_NULL;
	strm.zfree     = Z_NULL;
	strm.opaque    = Z_NULL;
	int err        = -1;
	int ret        = -1;
	err            = inflateInit2(&strm, (15 + 32));
	if (err == Z_OK) {
		err = inflate(&strm, Z_FINISH);
		if (err == Z_STREAM_END) ret = strm.total_out;
		else {
			inflateEnd(&strm);
			return err;
		}
	}
	else {
		inflateEnd(&strm);
		return err;
	}
	inflateEnd(&strm);
	return ret;
}

BYTE PaethPredictor(BYTE a, BYTE b, BYTE c) {
	BYTE p  = a + b - c;
	BYTE pa = Mathematics::abs(p - a);
	BYTE pb = Mathematics::abs(p - b);
	BYTE pc = Mathematics::abs(p - c);
	return (pa <= pb && pa <= pc) ? a : (pb <= pc) ? b : c;
}

void _ImageReconstruction(BYTE* data, int w, int h, int mode) {
	//MessageBox(nullptr, L"JUST__RGB!!!", L"ERROR", MB_OK);
	fin_decoded = (BYTE*)malloc(w * h * 3);
	if (mode == RGB_PURE) {
		int chunk_len = 1 + (w * 3);
		BYTE* chunk = (BYTE*)malloc(chunk_len);
		for (size_t k = 0; k < h; k++) {
			for (size_t i = 0; i != chunk_len; ++i) chunk[i] = data[k * chunk_len + i];
			fin_decoded[k * (w * 3) + 0] = chunk[1];
			fin_decoded[k * (w * 3) + 1] = chunk[2];
			fin_decoded[k * (w * 3) + 2] = chunk[3];
			if (chunk[0] == 0) {
				for (size_t i = 3; i != chunk_len - 1; ++i) fin_decoded[k * (w * 3) + i] = chunk[i + 1];
			}
			else if (chunk[0] == 1) {
				for (size_t i = 3; i != chunk_len - 1; ++i) {
					int tt = (int)(chunk[i + 1] + fin_decoded[k * (w * 3) + i - 3]);
					fin_decoded[k * (w * 3) + i] = (BYTE)tt;
				}
			}
			else if (chunk[0] == 2) {
				for (size_t i = 3; i != chunk_len - 1; ++i) {
					int tt = (int)(chunk[i + 1] + fin_decoded[(k - 1) * (w * 3) + i]);
					fin_decoded[k * (w * 3) + i] = (BYTE)tt;
				}
			}
			else if (chunk[0] == 3) {
				for (size_t i = 3; i != chunk_len - 1; ++i) {
					int tt = (int)(chunk[i + 1] + (int)((fin_decoded[(k - 1) * (w * 3) + i] + fin_decoded[k * (w * 3) + i - 3]) / 2));
					fin_decoded[k * (w * 3) + i] = (BYTE)tt;
				}
			}
			else if (chunk[0] == 4) {
				for (size_t i = 3; i != chunk_len - 1; ++i) {
					int tt = (int)(chunk[i + 1] + (int)PaethPredictor(fin_decoded[k * (w * 3) + i - 3], fin_decoded[(k - 1) * (w * 3) + i], fin_decoded[(k - 1) * (w * 3) + i - 3]));
					fin_decoded[k * (w * 3) + i] = (BYTE)tt;
				} 
			}
		}
		if (chunk != nullptr) {
			free(chunk);
			chunk = nullptr;
		}
	}	
	else if (mode == RGB_ALPHA) {
		//MessageBox(nullptr, L"WITH__ALPHA!!!", L"ERROR", MB_OK);
		BYTE* __fin_decoded = (BYTE*)malloc(w * h * 4);
		int chunk_len = 1 + (w * 4);
		BYTE* chunk = (BYTE*)malloc(chunk_len);
		for (size_t k = 0; k < h; k++) {
			for (size_t i = 0; i != chunk_len; ++i) chunk[i] = data[k * chunk_len + i];
			__fin_decoded[k * (w * 4) + 0] = chunk[1];
			__fin_decoded[k * (w * 4) + 1] = chunk[2];
			__fin_decoded[k * (w * 4) + 2] = chunk[3];
			__fin_decoded[k * (w * 4) + 3] = chunk[4];
			if (chunk[0] == 0) {
				for (size_t i = 4; i != chunk_len - 1; ++i) __fin_decoded[k * (w * 4) + i] = chunk[i + 1];
			}
			else if (chunk[0] == 1) {
				for (size_t i = 4; i != chunk_len - 1; ++i) {
					int tt = (int)(chunk[i + 1] + __fin_decoded[k * (w * 4) + i - 4]);
					__fin_decoded[k * (w * 4) + i] = (BYTE)tt;
				}
			}
			else if (chunk[0] == 2) {
				for (size_t i = 4; i != chunk_len - 1; ++i) {
					int tt = (int)(chunk[i + 1] + __fin_decoded[(k - 1) * (w * 4) + i]);
					__fin_decoded[k * (w * 4) + i] = (BYTE)tt;
				}
			}
			else if (chunk[0] == 3) {
				for (size_t i = 4; i != chunk_len - 1; ++i) {
					int tt = (int)(chunk[i + 1] + (int)((__fin_decoded[(k - 1) * (w * 4) + i] + __fin_decoded[k * (w * 4) + i - 4]) / 2));
					__fin_decoded[k * (w * 4) + i] = (BYTE)tt;
				}
			}
			else if (chunk[0] == 4) {
				for (size_t i = 4; i != chunk_len - 1; ++i) {
					int tt = (int)(chunk[i + 1] + (int)PaethPredictor(__fin_decoded[k * (w * 4) + i - 4], __fin_decoded[(k - 1) * (w * 4) + i], __fin_decoded[(k - 1) * (w * 4) + i - 4]));
					__fin_decoded[k * (w * 4) + i] = (BYTE)tt;
				}
			}
		}
		if (chunk != nullptr) {
			free(chunk);
			chunk = nullptr;
		}
		int gap = 0;
		for (size_t i = 0; i != w * h * 3; i+=3) {
			fin_decoded[i]   = __fin_decoded[i + gap];
			fin_decoded[i+1] = __fin_decoded[i + 1 + gap];
			fin_decoded[i+2] = __fin_decoded[i + 2 + gap];
			++gap;
		}
	}
}

DOUBLEWORD reverse_netwok_Byte_order(DOUBLEWORD num) {
	return ((num & 0xFF) << 24)      |
		   ((num & 0xFF00) << 8)     |
		   ((num & 0xFF0000) >> 8)   |
		   ((num & 0xFF000000) >> 24);
}

BYTE* byte_by_byteD(DOUBLEWORD num) {
	void* ptr = &num;
	BYTE TMP[4];
	for (size_t i = 0; i < 4; i++) 
		TMP[i] = *(char*)((char*)ptr + i);
	return TMP;
}

void parse_IHDR(Chunk* cnk, IHDR* ihdr, DOUBLEWORD code_length, FILE* file) {
	if (file == nullptr) return;
	fread(ihdr, code_length, 1, file);
	void* ptr = ihdr;
	for (size_t i = 0; i != sizeof(IHDR); ++i) cnk->Chunk_Data[i] = *(BYTE*)((char*)ptr + i);
	return;
}

int merged = 0;
void merge_idat_data(BYTE* data, DOUBLEWORD len) {
	for (size_t i = 0; i < len; i++) pre_encoded[i + merged] = data[i];
	merged += len;
}

void parse_IDAT(Chunk* cnk, IDAT* idat, IHDR* ihdr, DOUBLEWORD code_length, FILE* file) {
	if (file == nullptr) 
		return;
	idat->Bit_depth_from_ihdr = ihdr->Bit_depth;
	idat->Color_type_from_ihdr = ihdr->Color_type;
	idat->Compression_method_from_ihdr = ihdr->Compression_method;
	idat->Filter_method_from_ihdr = ihdr->Filter_method;
	idat->Interlace_method_from_ihdr = ihdr->Interlace_method;
	idat->total_pixel_size = reverse_netwok_Byte_order(ihdr->Width) * reverse_netwok_Byte_order(ihdr->Height);
	idat->data = (BYTE*)malloc(code_length);
	fread(idat->data, code_length, 1, file);
	merge_idat_data(idat->data, code_length);
	void* ptr = idat->data;
	for (size_t i = 0; i != code_length; ++i) 
		cnk->Chunk_Data[i] = *(BYTE*)((char*)ptr + i);
	return;
}

void parse_IEND(IEND* iend) {
	if (iend == nullptr)
		return;
	iend->if_ended = true;
	return;
}

void ImageReconstruction(IPng* png) {
	inflate(pre_encoded, INIT_WINDOW_SIZE, pre_decoded, INIT_WINDOW_SIZE);
	_ImageReconstruction(pre_decoded,
		                 reverse_netwok_Byte_order(png->ihdr.Width),
		                 reverse_netwok_Byte_order(png->ihdr.Height),
		                 (png->idat->Color_type_from_ihdr == 2) ? RGB_PURE : (png->idat->Color_type_from_ihdr == 6) ? RGB_ALPHA : RGB_UNKNOWN);
}

void parse_PLTE(Chunk* cnk, PLTE* plte, DOUBLEWORD code_length, FILE* file) {
	if (file == nullptr) 
		return;
	fread(plte, code_length, 1, file);
	void* ptr = plte;
	for (size_t i = 0; i != sizeof(PLTE); ++i) 
		cnk->Chunk_Data[i] = *(BYTE*)((char*)ptr + i);
	return;
}

BYTE data_rubbish[1024];
void Null_chunk_read(DOUBLEWORD code_length, FILE* file) {
	DOUBLEWORD crc_rubbish = 0;
	fread(data_rubbish, code_length, 1, file);
	fread(&crc_rubbish, 4, 1, file);
	return;
}

int read_chunks(IPng* png, FILE* file) {
	if (file == nullptr) 
		return -1;
	DOUBLEWORD tmp_length = 0;
	DOUBLEWORD tmp_Chunk_Type = 0;
	DOUBLEWORD tmp_crc = 0;
	BYTE       code[4];
	png->number_of_chunks = chunks_counter;
	fread(&tmp_length, 4, 1, file);
	fread(&tmp_Chunk_Type, 4, 1, file);
	BYTE* code_type = byte_by_byteD(tmp_Chunk_Type);
	for (size_t i = 0; i != 4; ++i) code[i] = code_type[i];
	code_type = nullptr;
	DOUBLEWORD code_length = reverse_netwok_Byte_order(tmp_length);
	while (code[1] != IEND_MARKER[1] && code[2] != IEND_MARKER[2] && code[3] != IEND_MARKER[3]) {
		if (png->number_of_chunks == 0) 
			png->chunks = (Chunk*)malloc(sizeof(Chunk));
		else 
			png->chunks = (Chunk*)realloc(png->chunks, sizeof(Chunk) * (chunks_counter + 1));
		if (code[0] == IHDR_MARKER[0] && code[1] == IHDR_MARKER[1] && code[2] == IHDR_MARKER[2] && code[3] == IHDR_MARKER[3]) {
			png->chunks[chunks_counter].Chunk_Data = (BYTE*)malloc(sizeof(IHDR));
			png->chunks[chunks_counter].Length = tmp_length;
			png->chunks[chunks_counter].Chunk_Type = tmp_Chunk_Type;
			png->chunks[chunks_counter].id = IHDR_ID;
			parse_IHDR(&png->chunks[chunks_counter], &png->ihdr, code_length, file);
			png->chunks[chunks_counter].CRC = fread(&tmp_crc, 4, 1, file);
			++chunks_counter;
			++png->number_of_chunks;
			fread(&tmp_length, 4, 1, file);
			fread(&tmp_Chunk_Type, 4, 1, file);
			code_type = byte_by_byteD(tmp_Chunk_Type);
			for (size_t i = 0; i != 4; ++i) 
				code[i] = code_type[i];
			code_type = nullptr;
			code_length = reverse_netwok_Byte_order(tmp_length);
		}
		else if (code[0] == IDAT_MARKER[0] && code[1] == IDAT_MARKER[1] && code[2] == IDAT_MARKER[2] && code[3] == IDAT_MARKER[3]) {
			png->chunks[chunks_counter].Chunk_Data = (BYTE*)malloc(code_length);
			png->chunks[chunks_counter].Length = tmp_length;
			png->chunks[chunks_counter].Chunk_Type = tmp_Chunk_Type;
			png->chunks[chunks_counter].id = IDAT_ID;
			if (idat_counter == 0) 
				png->idat = (IDAT*)malloc(sizeof(IDAT));
			else 
				png->idat = (IDAT*)realloc(png->idat, sizeof(IDAT) * (idat_counter + 1));
			parse_IDAT(&(png->chunks[chunks_counter]), &(png->idat[idat_counter]), &png->ihdr, code_length, file);
			png->chunks[chunks_counter].CRC = fread(&tmp_crc, 4, 1, file);
			++chunks_counter;
			++idat_counter;
			++png->number_of_chunks;
			fread(&tmp_length, 4, 1, file);
			fread(&tmp_Chunk_Type, 4, 1, file);
			code_type = byte_by_byteD(tmp_Chunk_Type);
			for (size_t i = 0; i != 4; ++i) 
				code[i] = code_type[i];
			code_type = nullptr;
			code_length = reverse_netwok_Byte_order(tmp_length);
		}
		else if (code[0] == PLTE_MARKER[0] && code[1] == PLTE_MARKER[1] && code[2] == PLTE_MARKER[2] && code[3] == PLTE_MARKER[3]) {
			png->chunks[chunks_counter].Chunk_Data = (BYTE*)malloc(sizeof(PLTE));
			png->chunks[chunks_counter].Length = tmp_length;
			png->chunks[chunks_counter].Chunk_Type = tmp_Chunk_Type;
			png->chunks[chunks_counter].id = PLTE_ID;
			parse_PLTE(&png->chunks[chunks_counter], &png->plte, code_length, file);
			png->chunks[chunks_counter].CRC = fread(&tmp_crc, 4, 1, file);
			++chunks_counter;
			++png->number_of_chunks;
			fread(&tmp_length, 4, 1, file);
			fread(&tmp_Chunk_Type, 4, 1, file);
			code_type = byte_by_byteD(tmp_Chunk_Type);
			for (size_t i = 0; i != 4; ++i) 
				code[i] = code_type[i];
			code_type = nullptr;
			code_length = reverse_netwok_Byte_order(tmp_length);
		}
		else ANCILLARY {
			Null_chunk_read(code_length, file);
			++chunks_counter;
			++png->number_of_chunks;
			fread(&tmp_length, 4, 1, file);
			fread(&tmp_Chunk_Type, 4, 1, file);
			code_type = byte_by_byteD(tmp_Chunk_Type);
			for (size_t i = 0; i != 4; ++i)
				code[i] = code_type[i];
			code_type = nullptr;
			code_length = reverse_netwok_Byte_order(tmp_length);
		}
	}
	parse_IEND(&png->iend);
	if (png->iend.if_ended) 
		ImageReconstruction(png);
	return 0;
}

static IPng* png = nullptr;

void read_(IPng* png) {
	FILE* file = fopen(png->address, "rb");
	if (file == nullptr) 
		return;
	BYTE* tmp = png->signature.signature;
	if (tmp == nullptr)
		return;
	fread(tmp, 1, 8, file);
	read_chunks(png, file);
	fclose(file);
	return;
}

BYTE* ARRAY = nullptr;
void __SYSTEMCALL_draw(HWND hWnd, IPng* png) {
	int w = reverse_netwok_Byte_order(png->ihdr.Width);
	int h = reverse_netwok_Byte_order(png->ihdr.Height);
	int tot = 3 * w * h;
	HDC dc = GetDC(hWnd);
	BITMAPINFO info;
	info.bmiHeader.biBitCount = 32;
	info.bmiHeader.biWidth  = w;
	info.bmiHeader.biHeight = -h;
	info.bmiHeader.biPlanes = 1;
	info.bmiHeader.biSizeImage = (((tot + 31) & ~31) / 8) * h;
	info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	info.bmiHeader.biCompression = BI_RGB;
	int result = StretchDIBits(dc, 0, 0, w, h, 0, 0, w, h, ARRAY, &info, DIB_RGB_COLORS, SRCCOPY);
}

void DrawOnWindow(HWND hWnd, IPng* png) {
	int w = reverse_netwok_Byte_order(png->ihdr.Width);
	int h = reverse_netwok_Byte_order(png->ihdr.Height);
	int tot = w * h * 3;
	ARRAY = (BYTE*)malloc(tot + (w*h));
	int space = 0;
	for (size_t i = 0; i < tot; i += 3) {
		ARRAY[i+space ]      = fin_decoded[i + 2];
		ARRAY[i + 1 + space] = fin_decoded[i + 1];
		ARRAY[i + 2 + space] = fin_decoded[i];
		++space;
	} 
	__SYSTEMCALL_draw(hWnd, png);
	free(ARRAY);
}

void __Clean_All() {
	// to do
}

void CleanPng(IPng* png) {
	void* __idat = png->idat;
	int CHNK = chunks_counter;
	for (size_t i = 0; i < CHNK; i++) {
		void* __chunks = png->chunks[i].Chunk_Data;
		if (__chunks != nullptr) {
			free(__chunks);
			__chunks = nullptr;
		}
	}
	void* __chunks = png->chunks;
	if (__chunks != nullptr) {
		free(__chunks);
		__chunks = nullptr;
	}
	if (__idat != nullptr) {
		free(__idat);
		__idat = nullptr;
	}
	__Clean_All();
}

int SafeEdge(int ref) {
	if (ref <= 50) return 50;
	else           return 30;
}

LPWSTR nname = L" ";
LRESULT CALLBACK WndProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam) {
	HWND hButton;
	HWND hText;
	switch (message) {
		case WM_PAINT:
			break;
		case WM_CLOSE: break;
		case WM_DESTROY:
			PostQuitMessage(0);
			return 0;
		case WM_CREATE:
			hButton = CreateWindow(L"button", 
				                   L"Check The Provided Image Address", 
				                   WS_CHILD | WS_VISIBLE | BS_DEFPUSHBUTTON, 
				                   0, 
				                   0, 
				                   250,
				                   20, 
				                   hwnd, 
				                   (HMENU)BUTTON_ID, 
				                   (HINSTANCE)GetWindowLong(hwnd, 0),
				                   NULL);

			hText   = CreateWindow(L"EDIT",
				                   nname,
				                   WS_CHILD | WS_VISIBLE,
				                   0, 
				                   25,
				                   1000,
				                   20,
				                   hwnd, 
				                   (HMENU)TEXT_ID, 
				                   (HINSTANCE)GetWindowLong(hwnd, 0),
				                   NULL);
			break;
		case WM_COMMAND:
			switch (wParam) {
			case BUTTON_ID:
				if (png->signature.signature[0] == 137 &&
					png->signature.signature[1] == 80  &&
					png->signature.signature[2] == 78  &&
					png->signature.signature[3] == 71  &&
					png->signature.signature[4] == 13  &&
					png->signature.signature[5] == 10  &&
					png->signature.signature[6] == 26  &&
					png->signature.signature[7] == 10) {
					MessageBox(hwnd, L".PNG image found\nPress OK to show it ...", L" ", MB_OK);
				}
				DrawOnWindow(hwnd, png);
				break;
			case TEXT_ID:
				break;
			}
			break;
		default: break;
	}
	return DefWindowProc(hwnd, message, wParam, lParam);
}

int __stdcall WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int iCmdShow) {
	png                    = (IPng*)malloc(sizeof(IPng));
	png->address           = "__Mont";
	nname                  = L"__Mont";

	/* png call */ read_(png);

	TCHAR szAppName[]      = L"Test";
	WNDCLASS wndclass;
	wndclass.style         = CS_HREDRAW | CS_VREDRAW;
	wndclass.lpfnWndProc   = WndProc;
	wndclass.cbClsExtra    = 0;
	wndclass.cbWndExtra    = 0;
	wndclass.hInstance     = hInstance;
	wndclass.hIcon         = LoadIcon(nullptr, IDI_APPLICATION);
	wndclass.hCursor       = LoadCursor(nullptr, IDC_ARROW);
	wndclass.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
	wndclass.lpszMenuName  = NULL;
	wndclass.lpszClassName = szAppName;
	if (!RegisterClass(&wndclass)) exit(0);
	HWND hwnd = CreateWindow(szAppName, 
		                       L"IMAN", 
		                       WS_OVERLAPPEDWINDOW, 
		                       0, 
		                       0, 
		                       reverse_netwok_Byte_order(png->ihdr.Width),
		                       reverse_netwok_Byte_order(png->ihdr.Height) + SafeEdge(png->ihdr.Height),
		                       nullptr, 
		                       nullptr,
		                       hInstance, 
		                       nullptr);

	if (!hwnd) 
		exit(0);
	ShowWindow(hwnd, iCmdShow);
	UpdateWindow(hwnd);
	MSG msg;
	while (GetMessage(&msg, nullptr, 0, 0)) {
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	//CleanPng(png);
	return 0;
}
