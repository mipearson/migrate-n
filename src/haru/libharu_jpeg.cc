/*
 * << H a r u -- Free PDF Library >> -- libharu_jpeg.cc
 *
 * Interface routines for "C".
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

#include <new>
#include "libharuc.h"

/*----- pdf_jpeg_image -------------------------------------------------------*/

pdf_jpeg_image
pdf_jpeg_image_new(pdf_doc doc)
{
    if (doc == NULL)
        return NULL;
    
    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        PdfJpegImage* image = new PdfJpegImage(doc_obj);
        return (pdf_jpeg_image)image;
    } catch (ALLOC_ERROR& e) {
        doc_obj->SetError(PDF_ERR_MALLOC);
    }
    return NULL;    
}

int
pdf_jpeg_image_load_from_file(pdf_jpeg_image image, char *filename)
{
    if (image == NULL || filename == NULL)
        return PDF_INVALID_PARAMETER;
    
    PdfJpegImage* image_obj = (PdfJpegImage*)image;
    try {
        image_obj->LoadFromFile(filename);
        return PDF_SUCCESS;
    } catch (PdfException& e) {
        image_obj->SetError(e.GetCode());
        return e.GetCode();
    } catch (ALLOC_ERROR& e) {
        image_obj->SetError(PDF_ERR_MALLOC);
        return PDF_ERR_MALLOC;
    }
}

void
pdf_jpeg_image_free_image(pdf_jpeg_image image)
{
    if (image == NULL)
        return;

    PdfJpegImage* image_obj = (PdfJpegImage*)image;
    image_obj->FreeImage();
}

/*---------------------------------------------------------------------------*/

#endif /* NOJPEG */

