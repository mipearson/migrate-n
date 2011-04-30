/*
 * << H a r u free pdf library >> -- PdfMbFontDef_Ming.cc
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

#include "libharu_cntfonts.h"

/*---------------------------------------------------------------------------*/
/*----- Ming Font -----------------------------------------------------------*/

void
PdfMingFontDef::Init()
{
    PDF_DEBUG_PRINT(("PdfMingFontDef %X::Init().\n", (int)this));

    SetBaseFont("Ming");
    fAscent = 859;
    fDescent = -141;
    fCapHeight = 859;
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH;
    fFontBBox.left = 0;
    fFontBBox.top = 859;
    fFontBBox.right = 1000;
    fFontBBox.bottom = -141;
    fItalicAngle = 0;
    fStemV = 78;
    fMissingWidth = 500;
    AddWidths1(13648, 13742, 500);
}

/*---------------------------------------------------------------------------*/
/*----- Ming Font (Bold) ----------------------------------------------------*/

void
PdfMingBoldFontDef::Init()
{
    PdfMingFontDef::Init();

    SetBaseFont("Ming,Bold");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_FOURCE_BOLD;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/
/*----- Ming Font (Italic) --------------------------------------------------*/

void
PdfMingItalicFontDef::Init()
{
    PdfMingFontDef::Init();

    SetBaseFont("Ming,Italic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH;
    fItalicAngle = -11;
}

/*---------------------------------------------------------------------------*/
/*----- Ming Font (Bold-Italic) ---------------------------------------------*/

void
PdfMingBoldItalicFontDef::Init()
{
    PdfMingFontDef::Init();

    SetBaseFont("Ming,BoldItalic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_FOURCE_BOLD;
    fItalicAngle = -11;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/

