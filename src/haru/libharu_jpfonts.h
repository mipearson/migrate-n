/*
 * << H a r u -- Free PDF Library >> -- libharu_jpfonts.h
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
 * Header file for Japanese fonts 
 *
 */

#ifndef _LIB_HARU_JPFONTS_H 
#define _LIB_HARU_JPFONTS_H 

#include "libharu.h"

/*---------------------------------------------------------------------------*/
/*----- PdfCMap_90msp_RKSJ_H ------------------------------------------------*/

class PdfCMap_90msp_RKSJ_H : public PdfCMap
{
public:
        const char*         GetCMapName()   { return "90msp-RKSJ-H"; };
        void                ParseText(const char* text, pdf_byte_type* btype);
protected:
        void                Init();
};  

/*---------------------------------------------------------------------------*/
/*----- PdfCMap_90ms_RKSJ_H -------------------------------------------------*/

class PdfCMap_90ms_RKSJ_H : public PdfCMap
{
public:
        const char*         GetCMapName()   { return "90ms-RKSJ-H"; };
        void                ParseText(const char* text, pdf_byte_type* btype);
protected:
        void                Init();
};  

/*---------------------------------------------------------------------------*/
/*----- PdfCMap_90ms_RKSJ_V -------------------------------------------------*/

class PdfCMap_90ms_RKSJ_V : public PdfCMap
{
public:
        const char*         GetCMapName()   { return "90ms-RKSJ-V"; };
        void                ParseText(const char* text, pdf_byte_type* btype);
        pdf_writing_mode    GetWritingMode()  { return PDF_WMODE_VERTICAL; };
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfCMap_EUC_H -------------------------------------------------------*/

class PdfCMap_EUC_H : public PdfCMap
{
public:
        const char*         GetCMapName()   { return "EUC-H"; };
        void                ParseText(const char* text, pdf_byte_type* btype);
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfCMap_EUC_V -------------------------------------------------------*/

class PdfCMap_EUC_V : public PdfCMap
{
public:
        const char*         GetCMapName()   { return "EUC-V"; };
        void                ParseText(const char* text, pdf_byte_type* btype);
        pdf_writing_mode    GetWritingMode()  { return PDF_WMODE_VERTICAL; }
protected:
        void                Init();
};
        
/*---------------------------------------------------------------------------*/
/*----- PdfGothicFontDef ----------------------------------------------------*/

class PdfGothicFontDef : public PdfCIDType2FontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfGothicBoldFontDef ------------------------------------------------*/

class PdfGothicBoldFontDef : public PdfGothicFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfGothicItalicFontDef ----------------------------------------------*/

class PdfGothicItalicFontDef : public PdfGothicFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfGothicBoldItalicFontDef ------------------------------------------*/

class PdfGothicBoldItalicFontDef : public PdfGothicFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfPGothicFontDef ---------------------------------------------------*/

class PdfPGothicFontDef : public PdfCIDType2FontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfPGothicBoldFontDef -----------------------------------------------*/

class PdfPGothicBoldFontDef : public PdfPGothicFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfPGothicItalicFontDef ---------------------------------------------*/

class PdfPGothicItalicFontDef : public PdfPGothicFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfPGothicBoldItalicFontDef -----------------------------------------*/

class PdfPGothicBoldItalicFontDef : public PdfPGothicFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfMincyoFontDef ----------------------------------------------------*/

class PdfMincyoFontDef : public PdfCIDType2FontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfMincyoBoldFontDef ------------------------------------------------*/

class PdfMincyoBoldFontDef : public PdfMincyoFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfMincyoItalicFontDef ----------------------------------------------*/

class PdfMincyoItalicFontDef : public PdfMincyoFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfMincyoBoldItalicFontDef ------------------------------------------*/

class PdfMincyoBoldItalicFontDef : public PdfMincyoFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfPMincyoFontDef ---------------------------------------------------*/

class PdfPMincyoFontDef : public PdfCIDType2FontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfPMincyoBoldFontDef -----------------------------------------------*/

class PdfPMincyoBoldFontDef : public PdfPMincyoFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfPMincyoItalicFontDef ---------------------------------------------*/

class PdfPMincyoItalicFontDef : public PdfPMincyoFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/
/*----- PdfPMincyoBoldItalicFontDef -----------------------------------------*/

class PdfPMincyoBoldItalicFontDef : public PdfPMincyoFontDef
{
public:
protected:
        void                Init();
};

/*---------------------------------------------------------------------------*/

#endif /* _LIB_HARU_JPFONTS_H */

