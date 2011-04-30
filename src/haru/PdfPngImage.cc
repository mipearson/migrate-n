/*
 * << H a r u -- Free PDF Library >> -- PdfPngImage.cpp
 *
 * Copyright (c) 1999-2003 Takeshi Kanno <takeshi_kanno@est.hi-ho.ne.jp>
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 * It is provided "as is" without express or implied warranty.
 *
 */

#ifndef NOPNG

#include <assert.h>
#include <stdio.h>
#include <errno.h>
#include "libharu.h"

void pdf_png_read_data(png_structp png_ptr, png_bytep data, png_uint_32 length);

void pdf_error_func(png_structp png_ptr, const char* msg);

void pdf_png_read_data(png_structp png_ptr, png_bytep data, png_uint_32 length)
{
    FILE* f = (FILE*)png_get_io_ptr(png_ptr);

    if (fread(data, length, 1, f) < 1)
        throw PdfException(errno, "ERROR: PdfPngImage fread failed.");
}

void pdf_error_func(png_structp png_ptr, const char* msg)
{
    throw PdfException(PDF_RUNTIME_ERROR, "ERROR: PdfPngImage[%s].", msg);
}

/*----------------------------------------------------------------------------*/
/*----- PdfPngImage class ----------------------------------------------------*/

PdfPngImage::PdfPngImage(PdfDoc* doc)
        : PdfImage(doc)
{
    fWidth = 0;
    fHeight = 0;
    fBitsPerComponent = 0;
    fHasImage = false;
    fNumPallet = 0;
    fPallet = NULL;
}

PdfPngImage::~PdfPngImage()
{
    PDF_DEBUG_PRINT(("PdfPngImage::~PdfPngImage()\n"));
    FreeImage();
}

void
PdfPngImage::CreatePallet(png_structp png_ptr, png_infop info_ptr)
{
    /* make a pallet array for indexed image. */
    int num_pl = 0;
    png_color* src_pl = NULL;
	
	if (fPallet != NULL)
		delete[] fPallet;
	
    png_get_PLTE(png_ptr, info_ptr, (png_color**)&src_pl, &num_pl);

    fPallet = new unsigned char[num_pl * 3];
	PDF_DEBUG_PRINT(("++ [%x] fPallet new.\n", (int)fPallet));

	unsigned char* p = fPallet;
    for (int j = 0; j < num_pl; j++, src_pl++) {
        *p++ = src_pl->red;
        *p++ = src_pl->green;
        *p++ = src_pl->blue;
    }

    fNumPallet = num_pl;
    fColorSpace = PDF_CS_INDEXED;

	PdfArray* cs = new PdfArray(GetXref());
	AddElement("ColorSpace", cs);

	cs->Add(new PdfName("Indexed"));
	cs->Add(new PdfName("DeviceRGB"));
	cs->Add(new PdfNumber(fNumPallet - 1));

	PdfBinary* pl = new PdfBinary();
	pl->SetData((void*)fPallet, fNumPallet * 3);
	cs->Add(pl);
}

#define PNG_BYTES_TO_CHECK   8

void
PdfPngImage::LoadFromFile(const char* filename)
{
    FreeImage();

    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    FILE* infile = NULL;

    /* create read_struct. */
    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
        NULL, pdf_error_func, pdf_error_func);

    if (png_ptr == NULL)
        throw PdfException(PDF_ERR_MALLOC,
            "ERROR: PdfPdfImage::LoadFromFile --png_create_read_struct.\n");

    /* create info-struct */
    info_ptr = png_create_info_struct(png_ptr);

    if (info_ptr == NULL) {
        png_destroy_read_struct(&png_ptr, NULL, NULL);
        throw PdfException(PDF_ERR_MALLOC,
            "ERROR: PdfPdfImage::LoadFromFile --png_create_info_struct.\n");
    }

    try {
        /* io setting */
        infile = OpenImageFile(filename);

        png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK);
        png_set_read_fn(png_ptr, (void*)infile, (png_rw_ptr)&pdf_png_read_data);

        /* reading info structure. */
        png_read_info(png_ptr, info_ptr);

        /* 16bit image and alpha color type are not supported. */
        if (info_ptr->bit_depth == 16)
            png_set_strip_16(png_ptr);

        if (PNG_COLOR_MASK_ALPHA & info_ptr->color_type) {
            /* png image with alpha chanel is not supported. */
            png_set_strip_alpha(png_ptr);
        }

        png_read_update_info(png_ptr, info_ptr);

        /* if the image has color palette, copy the pallet of the image to
         * create color map.
         */
        if (info_ptr->color_type == PNG_COLOR_TYPE_PALETTE)
            CreatePallet(png_ptr, info_ptr);
        else if (info_ptr->color_type == PNG_COLOR_TYPE_GRAY)
            fColorSpace = PDF_CS_DEVICE_GRAY;
        else
            fColorSpace = PDF_CS_DEVICE_RGB;

        png_uint_32 len = png_get_rowbytes(png_ptr, info_ptr);

        /* if the image is interlaced, read whole image at once. */
        if (png_get_interlace_type(png_ptr, info_ptr) != PNG_INTERLACE_NONE) {
            png_bytep* row_pointers = new png_bytep[info_ptr->height];

            memset(row_pointers, 0x00, info_ptr->height);
            try {
                for (int i = 0; i < (int)info_ptr->height; i++) {
                    row_pointers[i] = new png_byte[len];
                    PDF_DEBUG_PRINT(("++ [%x] png_byte new.\n", 
                                (int)row_pointers[i]));

                    if (row_pointers[i] == NULL)
                        throw PdfException(PDF_ERR_MALLOC,
                        "ERROR: PdfPdfImage::LoadFromFile --malloc failed.");
                }
                png_read_image(png_ptr, row_pointers);
                for (int j = 0; j < (int)info_ptr->height; j++)
                    GetStream()->Write(row_pointers[j], len);
            } catch (PdfException& e) {
                for (int i = 0; i < (int)info_ptr->height; i++)
                    if (row_pointers[i] != NULL)
                        delete[] row_pointers[i];
                delete[] row_pointers;
                throw;
            }

            /* clean up */
            for (int i = 0; i < (int)info_ptr->height; i++)
                if (row_pointers[i] != NULL) {
                    PDF_DEBUG_PRINT(("++ [%x] row_pointers(%d) delete.\n", 
                                (int)row_pointers[i], i));
                    delete[] row_pointers[i];
                }
            PDF_DEBUG_PRINT(("++ [%x] row_pointers delete.\n", 
                        (int)row_pointers));
            delete[] row_pointers;

        } else {
            png_bytep buf_ptr = new png_byte[len];

            /*
             * loading image with png_read_rows method and write data to
             * stream.
             */
            try {
                for (int i = 0; i < (int)info_ptr->height; i++) {
                    png_read_rows(png_ptr, (png_byte**)&buf_ptr, NULL, 1);
                    GetStream()->Write(buf_ptr, len);
                };
            } catch (PdfException& e) {
                delete[] buf_ptr;
                throw;
            }
            delete[] buf_ptr;
        }

        /* setting the info of the image. */
        fWidth = (unsigned int)info_ptr->width;
        fHeight = (unsigned int)info_ptr->height;
        fBitsPerComponent = (unsigned int)info_ptr->bit_depth;
    } catch (PdfException& e) {
        png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
        if (infile != NULL)
             fclose(infile);
        throw;
    }

    /* clean up */
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    if (infile != NULL)
        fclose(infile);

    fHasImage = true;
    return;
}

void
PdfPngImage::FreeImage()
{
    if (!fHasImage)
        return;
    GetStream()->FreeData();
    fWidth = 0;
    fHeight = 0;
    fBitsPerComponent = 0;
    fHasImage = false;
    fNumPallet = 0;

	PDF_DEBUG_PRINT(("++ [%x] fPallet delete.\n", (int)fPallet));	
    delete[] fPallet;
    fPallet = NULL;
}

FILE*
PdfPngImage::OpenImageFile(const char* filename)
{
    PDF_DEBUG_PRINT(("PdfPdfImage::OpenImageFile[%s]\n", filename));

    png_byte header[PNG_BYTES_TO_CHECK];
    memset(header, 0x00, PNG_BYTES_TO_CHECK);

    FILE *f = fopen(filename, "rb");
    if (f == NULL)
        throw PdfException(PDF_ERR_FILE_OPEN,
                "ERROR: cannot open %s.", filename);

    fread(header, PNG_BYTES_TO_CHECK, 1, f);
    if (png_sig_cmp(header, (png_size_t)0, PNG_BYTES_TO_CHECK) != 0) {
        fclose(f);

        throw PdfException(PDF_RUNTIME_ERROR,
            "ERROR: PdfPdfImage::OpenImageFile -- "
                "%s is not a valid PNG file.\n", filename);
    }

    return f;
}

/*----------------------------------------------------------------------------*/

#endif /* NOPNG */

