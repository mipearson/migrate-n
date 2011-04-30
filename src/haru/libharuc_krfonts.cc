/*
 * << H a r u --free pdf library >> -- libharuc_krfonts.cc
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

#include "libharuc_krfonts.h"
#include "libharu_krfonts.h"

/*-- Create CMap. --*/
    
pdf_cmap pdf_create_kr_cmap(pdf_predefined_kr_cmap map)
{
    try {
        switch (map) {
            case PDF_CMAP_KSC_MS_UHC_H:
                return new PdfCMap_KSCms_UHC_H();
            case PDF_CMAP_KSC_MS_UHC_V:
                return NULL;
            case PDF_CMAP_KSC_EUC_H:
                return new PdfCMap_KSC_EUC_H();
            case PDF_CMAP_KSC_EUC_V:
                return NULL;
        }
    } catch (...) {
        return NULL;
    }
    return NULL;
}

/*-- Create Korean font definition objects. --*/

pdf_cid_type2_fontdef pdf_create_kr_fontdef(pdf_kr_font font) {
    try {
        switch (font) {
            case PDF_FONT_BATANG_CHE:
                return new PdfBatangCheFontDef();
            case PDF_FONT_BATANG_CHE_BOLD:
                return new PdfBatangCheBoldFontDef();
            case PDF_FONT_BATANG_CHE_ITALIC:
                return new PdfBatangCheItalicFontDef();
            case PDF_FONT_BATANG_CHE_BOLD_ITALIC:
                return new PdfBatangCheBoldItalicFontDef();
            case PDF_FONT_DOTUM_CHE:
                return new PdfDotumCheFontDef();
            case PDF_FONT_DOTUM_CHE_BOLD:
                return new PdfDotumCheBoldFontDef();
            case PDF_FONT_DOTUM_CHE_ITALIC:
                return new PdfDotumCheItalicFontDef();
            case PDF_FONT_DOTUM_CHE_BOLD_ITALIC:
                return new PdfDotumCheBoldItalicFontDef();
        }
    } catch (...) {
       return NULL;
    }
    return NULL;
}

