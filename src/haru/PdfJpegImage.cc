/*
 * << H aru -- Free PDF Library >> -- PdfJpegImage.cpp
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

#ifndef NOJPEG

#include <assert.h>
#include <stdio.h>
#include <new>
#include "libharu.h"

/*----- PdfJpegImage class ---------------------------------------------------*/

PdfJpegImage::PdfJpegImage(PdfDoc* doc)
        : PdfImage(doc)
{
    fWidth = 0;
    fHeight = 0;
    fBitsPerComponent = 0;
    fHasImage = false;
}

PdfJpegImage::~PdfJpegImage()
{
    FreeImage();
}

void
PdfJpegImage::LoadFromFile(const char* filename)
{
    FreeImage();

    /* read the jpeg-file to get the image parameters */
    struct jpeg_error_mgr errmgr;
    memset(&errmgr, 0x00, sizeof(jpeg_error_mgr));

    FILE* infile = fopen(filename, "rb");
    if (infile == NULL)
        throw PdfException(PDF_ERR_FILE_OPEN,
                "ERROR: cannot open(%s).", filename);

    size_t isize = sizeof(jpeg_decompress_struct);
    j_decompress_ptr pjpeg = new jpeg_decompress_struct;
    PDF_DEBUG_PRINT(("++ [%x] jpeg_decompress_struct new.\n", (int)pjpeg));
    if (pjpeg == NULL) {
        fclose(infile);
        throw PdfException(PDF_ERR_MALLOC,
            "ERROR: cannot allocate jpeg decompress sutruct.");
    }
    memset(pjpeg, 0x00, isize);

    try {
        pjpeg->err = jpeg_std_error(&errmgr);
        errmgr.error_exit = PdfJpegErrorExit;

        PDF_DEBUG_PRINT(("DEBUG: PdfJpegImage --jpeg_create_decompress\n"));
        jpeg_create_decompress(pjpeg);

        try {
            PDF_DEBUG_PRINT(("DEBUG: PdfJpegImage --jpeg_stdio_src\n"));
            jpeg_stdio_src(pjpeg, infile);

            PDF_DEBUG_PRINT(("DEBUG: PdfJpegImage --jpeg_read_header\n"));
            jpeg_read_header(pjpeg, FALSE);

            PDF_DEBUG_PRINT(("DEBUG: PdfJpegImage --jpeg_start_decompress\n"));
            jpeg_start_decompress(pjpeg);

            PDF_DEBUG_PRINT(("DEBUG: width=%d, height=%d, "
                "output_component=%d\n", pjpeg->output_width, 
                pjpeg->output_height, pjpeg->out_color_components));
            fWidth = (int)pjpeg->image_width;
            fHeight = (int)pjpeg->image_height;
            fBitsPerComponent = 8;
            AddFilter(PDF_FILTER_DCT_DECODE);
			 
			switch ((int)pjpeg->jpeg_color_space) { 
				case 1: 
					fColorSpace = PDF_CS_DEVICE_GRAY; 
					break; 
				case 3: 
					fColorSpace = PDF_CS_DEVICE_RGB; 
					break; 
				case 5: { 
					fColorSpace = PDF_CS_DEVICE_CMYK; 
					PdfArray* darray = new PdfArray(GetXref()); 
					AddElement("Decode", darray); 
					darray->Add(new PdfNumber(1)); 
					darray->Add(new PdfNumber(0)); 
					darray->Add(new PdfNumber(1)); 
					darray->Add(new PdfNumber(0)); 
					darray->Add(new PdfNumber(1)); 
					darray->Add(new PdfNumber(0)); 
					darray->Add(new PdfNumber(1)); 
					darray->Add(new PdfNumber(0)); 
				} 
					break; 
				default: 
					throw PdfException(PDF_RUNTIME_ERROR, 
								"error: unsupported format."); 
			} 
        } catch (PdfException& e) {
            PDF_DEBUG_PRINT(("DEBUG: PdfJpegImage --an error was "
                "occured(%s). \n", e.what()));
            jpeg_destroy_decompress(pjpeg);
            delete(pjpeg);
            throw;
        }
    } catch (PdfException& e) {
        fclose(infile);
        throw;
    }

    jpeg_destroy_decompress(pjpeg);
    PDF_DEBUG_PRINT(("++ [%x] jpeg_decompress_struct delete.\n", (int)pjpeg));
    delete pjpeg;

    /* load the image data to a Stream object */
    fseek(infile, 0, SEEK_SET);

    char buf[PDF_DEF_BUF_SIZE];
    for (;;) {
        int len = fread(buf, 1, PDF_DEF_BUF_SIZE, infile);
        if (len == 0)
            break;
        GetStream()->Write(buf, len);

        if (len < PDF_DEF_BUF_SIZE)
            break;
    }
    fclose(infile);
    fHasImage = true;
}

void
PdfJpegImage::FreeImage()
{
    if (!fHasImage)
        return;
    GetStream()->FreeData();
    fWidth = 0;
    fHeight = 0;
    fBitsPerComponent = 0;
    fHasImage = false;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

void
PdfJpegErrorExit(j_common_ptr cinfo) {
    char msg[JMSG_LENGTH_MAX + 1];

    /* Format message and throw PdfException with the message. */
    (*cinfo->err->format_message)(cinfo, msg);
    throw PdfException(PDF_RUNTIME_ERROR, "ERROR: PdfJpegImage[%s].", msg);
}

/*----------------------------------------------------------------------------*/
#endif /* NOJPEG */

