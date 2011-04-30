/*
 * << H a r u -- Free PDF Library >> -- libharu_krfonts.h
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
 * Korean fonts definition. 
 *
 */

#ifndef _LIB_HARU_KRFONTS_H 
#define _LIB_HARU_KRFONTS_H 

#include "libharu.h"

/*---------------------------------------------------------------------------*/
/*----- PdfCMap_KSCms_UHC_H -------------------------------------------------*/

class PdfCMap_KSCms_UHC_H : public PdfCMap
{
public:
        const char*         GetCMapName()   { return "KSCms-UHC-H"; };
        void                ParseText(const char* text, pdf_byte_type* btype);
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfCMap_KSC_EUC_H ---------------------------------------------------*/

class PdfCMap_KSC_EUC_H : public PdfCMap
{
public:
        const char*         GetCMapName()   { return "KSC-EUC-H"; };
        void                ParseText(const char* text, pdf_byte_type* btype);
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfBatangCheFontDef -------------------------------------------------*/

class PdfBatangCheFontDef : public PdfCIDType2FontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfBatangCheBoldFontDef ---------------------------------------------*/

class PdfBatangCheBoldFontDef : public PdfBatangCheFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfBatangCheItalicFontDef -------------------------------------------*/

class PdfBatangCheItalicFontDef : public PdfBatangCheFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfBatangCheBoldItalicFontDef ---------------------------------------*/

class PdfBatangCheBoldItalicFontDef : public PdfBatangCheFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfDotumCheFontDef --------------------------------------------------*/

class PdfDotumCheFontDef : public PdfCIDType2FontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfDotumCheBoldFontDef ----------------------------------------------*/

class PdfDotumCheBoldFontDef : public PdfDotumCheFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfDotumCheItalicFontDef --------------------------------------------*/

class PdfDotumCheItalicFontDef : public PdfDotumCheFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfDotumCheBoldItalicFontDef ----------------------------------------*/

class PdfDotumCheBoldItalicFontDef : public PdfDotumCheFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/

#endif /* _LIB_HARU_KRFONTS_H */

