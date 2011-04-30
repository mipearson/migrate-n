/*
 * << H a r u free pdf library >> -- PdfMbFontDef_BatangChe.cc
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
/*----- BatangChe Font ------------------------------------------------------*/

void
PdfBatangCheFontDef::Init()
{
    PDF_DEBUG_PRINT(("PdfBatangCheFontDef %X::Init().\n", (int)this));

    SetBaseFont("BatangChe");
    fAscent = 859;
    fDescent = -141;
    fCapHeight = 859;
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_SERIF;
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
/*----- BatangChe Font (Bold) -----------------------------------------------*/

void
PdfBatangCheBoldFontDef::Init()
{
    PdfBatangCheFontDef::Init();

    SetBaseFont("BatangChe,Bold");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + 
        PDF_FONT_FOURCE_BOLD + PDF_FONT_SERIF;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/
/*----- BatangChe Font (Italic) ---------------------------------------------*/

void
PdfBatangCheItalicFontDef::Init()
{
    PdfBatangCheFontDef::Init();

    SetBaseFont("BatangChe,Italic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_SERIF;
    fItalicAngle = -11;
}

/*---------------------------------------------------------------------------*/
/*----- BatangChe Font (Bold-Italic) ----------------------------------------*/

void
PdfBatangCheBoldItalicFontDef::Init()
{
    PdfBatangCheFontDef::Init();

    SetBaseFont("BatangChe,BoldItalic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + 
        PDF_FONT_FOURCE_BOLD + PDF_FONT_SERIF;
    fItalicAngle = -11;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/

