/*
 * << H a r u --free pdf library >> -- PdfCIDFontDef.cpp
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

#include <assert.h>
#include <errno.h>
#include "libharu.h"

/*----- PdfCIDWItem class ----------------------------------------------------*/

PdfCIDWItem::~PdfCidWItem()
{
    PDF_DEBUG_PRINT(("++ [%x] PdfCidWItem delete. \n", (int)this));
}

/*----------------------------------------------------------------------------*/
/*----- PdfCidWItem1 class ---------------------------------------------------*/

PdfCidWItem1::PdfCidWItem1(unsigned int from, unsigned int to, 
        unsigned int value)
{
    PDF_DEBUG_PRINT(("++ [%x] PdfCidWItem1 new. \n", (int)this));
    fFrom = from;
    fTo = to;
    fValue = value;
}

/*----------------------------------------------------------------------------*/
/*----- PdfCidWItem2 class ---------------------------------------------------*/

PdfCidWItem2::PdfCidWItem2(unsigned int from, unsigned int count, 
        unsigned int* values)
{
    PDF_DEBUG_PRINT(("++ [%x] PdfCidWItem2 new. \n", (int)this));
    fFrom = from;
    fCount = count;

    PDF_DEBUG_PRINT(("++ [%x] PdfCidWItem::fWidths new. \n", (int)this));
    fWidths = new unsigned int[count];
    memcpy(fWidths, values, sizeof(unsigned int) * count);
}

PdfCidWItem2::~PdfCidWItem()
{
    PDF_DEBUG_PRINT(("++ [%x] PdfCidWItem::fWidths delete. \n", (int)this));
    delete fWidths();
}

