/*
 * << H a r u -- Free PDF Library >> -- libharu_cnsfonts.h
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
 *  Chinese (Traditional) fonts definition.
 *
 */

#ifndef _LIB_HARU_CNTFONTS_H 
#define _LIB_HARU_CNTFONTS_H 

#include "libharu.h"

/*---------------------------------------------------------------------------*/
/*----- PdfCMapETen_B5_H ----------------------------------------------------*/

class PdfCMap_ETen_B5_H : public PdfCMap
{
public:
        const char*         GetCMapName()   { return "ETen-B5-H"; };
        void                ParseText(const char* text, pdf_byte_type* btype);
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfMingFontDef ------------------------------------------------------*/

class PdfMingFontDef : public PdfCIDType2FontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfMingBoldFontDef --------------------------------------------------*/

class PdfMingBoldFontDef : public PdfMingFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfMingItalicFontDef ------------------------------------------------*/

class PdfMingItalicFontDef : public PdfMingFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfMingBoldItalicFontDef --------------------------------------------*/

class PdfMingBoldItalicFontDef : public PdfMingFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/

#endif /* _LIB_HARU_CNTFONTS_H */

