/*
 * << H a r u -- Free PDF Library >> -- PdfUtils.cc
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

#include <new>
#include "libharu.h"

/*----- common routines -----------------------------------------------------*/

pdf_rgb_color 
PdfRGBColor(double r, double g, double b)
{
    pdf_rgb_color rgb;

    rgb.red = r;
    rgb.green = g;
    rgb.blue = b;

    return rgb;
}

pdf_box
PdfBox(int left, int bottom, int right, int top)
{
    pdf_box box;

    box.left = left;
    box.bottom = bottom;
    box.right = right;
    box.top = top;

    return box;
}

pdf_rect
PdfRect(double left, double bottom, double right, double top)
{
    pdf_rect rect;

    rect.left = left;
    rect.bottom = bottom;
    rect.right = right;
    rect.top = top;

    return rect;
}

pdf_point
PdfPoint(double x, double y)
{
    pdf_point point;

    point.x = x;
    point.y = y;

    return point;
}

pdf_text_matrix
PdfTextMatrix(double a, double b, double c, double d,double x, double y)
{
    pdf_text_matrix matrix;

    matrix.a = a;
    matrix.b = b;
    matrix.c = c;
    matrix.d = d;
    matrix.x = x;
    matrix.y = y;

    return matrix;
}

void
pdf_print_binary(const unsigned char* buf, int len, char* caption)
{
    printf("DEBUG: %s ", caption);
    for (int i = 0; i < len; i++)
        printf("%02X ", buf[i]);
    printf("\n");
}

#ifdef _NOT_SUPPORT_STD

int
throw_new_handler(size_t size)
{
    throw PDF_STD_EXCEPTION("bad alloc exception thrown");
    return 0;
}

void
throw_new_handler()
{
    throw PDF_STD_EXCEPTION("bad alloc exception thrown");
}

#endif

/*---------------------------------------------------------------------------*/

