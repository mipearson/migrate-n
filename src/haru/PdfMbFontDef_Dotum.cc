/*
 * << H a r u free pdf library >> -- PdfMbFontDef_DotumChe.cc
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

#include "libharu_krfonts.h"

/*---------------------------------------------------------------------------*/
/*----- DotumChe Font -------------------------------------------------------*/

void
PdfDotumCheFontDef::Init()
{
    PDF_DEBUG_PRINT(("PdfDotumCheFontDef %X::Init().\n", (int)this));

    SetBaseFont("DotumChe");
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
    AddWidths1(8094, 8189, 500);
}

/*---------------------------------------------------------------------------*/
/*----- DotumChe Font (Bold) ------------------------------------------------*/

void
PdfDotumCheBoldFontDef::Init()
{
    PdfDotumCheFontDef::Init();

    SetBaseFont("DotumChe,Bold");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_FOURCE_BOLD;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/
/*----- DotumChe Font (Italic) ----------------------------------------------*/

void
PdfDotumCheItalicFontDef::Init()
{
    PdfDotumCheFontDef::Init();

    SetBaseFont("DotumChe,Italic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH;
    fItalicAngle = -11;
}

/*---------------------------------------------------------------------------*/
/*----- DotumChe Font (Bold-Italic) -----------------------------------------*/

void
PdfDotumCheBoldItalicFontDef::Init()
{
    PdfDotumCheFontDef::Init();

    SetBaseFont("DotumChe,BoldItalic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_FOURCE_BOLD;
    fItalicAngle = -11;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/

