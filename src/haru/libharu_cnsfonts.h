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
 *  Chinese (Simplified) fonts definition.
 *
 */

#ifndef _LIB_HARU_CNSFONTS_H 
#define _LIB_HARU_CNSFONTS_H 

#include "libharu.h"

/*----- PdfCMap_GB_EUC_H ----------------------------------------------------*/

class PdfCMap_GB_EUC_H : public PdfCMap
{
public:
        const char*         GetCMapName()   { return "GB-EUC-H"; };
        void                ParseText(const char* text, pdf_byte_type* btype);
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfSimSunFontDef ----------------------------------------------------*/

class PdfSimSunFontDef : public PdfCIDType2FontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfSimSunBoldFontDef ------------------------------------------------*/

class PdfSimSunBoldFontDef : public PdfSimSunFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfSimSunItalicFontDef ----------------------------------------------*/

class PdfSimSunItalicFontDef : public PdfSimSunFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfSimSunBoldItalicFontDef ------------------------------------------*/

class PdfSimSunBoldItalicFontDef : public PdfSimSunFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfSimHeiFontDef ----------------------------------------------------*/

class PdfSimHeiFontDef : public PdfCIDType2FontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfSimHeiBoldFontDef ------------------------------------------------*/

class PdfSimHeiBoldFontDef : public PdfSimHeiFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfSimHeiItalicFontDef ----------------------------------------------*/

class PdfSimHeiItalicFontDef : public PdfSimHeiFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfSimHeiBoldItalicFontDef ------------------------------------------*/

class PdfSimHeiBoldItalicFontDef : public PdfSimHeiFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/

#endif /* _LIB_HARU_CNSFONTS_H */

