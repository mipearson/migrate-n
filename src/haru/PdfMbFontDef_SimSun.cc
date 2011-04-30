/*
 * << H a r u free pdf library >> -- PdfMbFontDef_SimSun.cc
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

#include "libharu_cnsfonts.h"

/*---------------------------------------------------------------------------*/
/*----- SimSun Font -----------------------------------------------------------*/

void
PdfSimSunFontDef::Init()
{
    PDF_DEBUG_PRINT(("PdfSimSunFontDef %X::Init().\n", (int)this));

    SetBaseFont("SimSun");
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
    AddWidths1(666, 699, 500);
    AddWidths1(814, 939, 500);
    AddWidths1(7716, 7716, 500);
}

/*---------------------------------------------------------------------------*/
/*----- SimSun Font (Bold) ----------------------------------------------------*/

void
PdfSimSunBoldFontDef::Init()
{
    PdfSimSunFontDef::Init();

    SetBaseFont("SimSun,Bold");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + 
		PDF_FONT_FOURCE_BOLD + PDF_FONT_SERIF;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/
/*----- SimSun Font (Italic) --------------------------------------------------*/

void
PdfSimSunItalicFontDef::Init()
{
    PdfSimSunFontDef::Init();

    SetBaseFont("SimSun,Italic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_SERIF;
    fItalicAngle = -11;
}

/*---------------------------------------------------------------------------*/
/*----- SimSun Font (Bold-Italic) ---------------------------------------------*/

void
PdfSimSunBoldItalicFontDef::Init()
{
    PdfSimSunFontDef::Init();

    SetBaseFont("SimSun,BoldItalic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + 
		PDF_FONT_FOURCE_BOLD + PDF_FONT_SERIF;
    fItalicAngle = -11;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/

