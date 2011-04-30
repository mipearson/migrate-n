/*
 * << H a r u --free pdf library >> -- libharuc_jpfonts.cc
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

#include "libharuc_jpfonts.h"
#include "libharu_jpfonts.h"

/*-- Create CMap. --*/
    
pdf_cmap pdf_create_jp_cmap(pdf_predefined_jp_cmap map)
{
    try {
        switch (map) {
            case PDF_CMAP_90MSP_RKSJ_H:
                return new PdfCMap_90msp_RKSJ_H();
            case PDF_CMAP_90MS_RKSJ_H:
                return new PdfCMap_90ms_RKSJ_H();
            case PDF_CMAP_90MS_RKSJ_V:
                return new PdfCMap_90ms_RKSJ_V();
            case PDF_CMAP_EUC_H:
                return new PdfCMap_EUC_H();
            case PDF_CMAP_EUC_V:
                return new PdfCMap_EUC_V();
        }
    } catch (...) {
        return NULL;
    }
    return NULL;
}

/*-- Create Japanese font definition objects. --*/

pdf_cid_type2_fontdef pdf_create_jp_fontdef(pdf_jp_font font) {
    try {
        switch (font) {
            case PDF_FONT_GOTHIC:
                return new PdfGothicFontDef();
            case PDF_FONT_GOTHIC_BOLD:
                return new PdfGothicBoldFontDef();
            case PDF_FONT_GOTHIC_ITALIC:
                return new PdfGothicItalicFontDef();
            case PDF_FONT_GOTHIC_BOLD_ITALIC:
                return new PdfGothicBoldItalicFontDef();
            case PDF_FONT_PGOTHIC:
                return new PdfPGothicFontDef;
            case PDF_FONT_PGOTHIC_BOLD:
                return new PdfPGothicBoldFontDef;
            case PDF_FONT_PGOTHIC_ITALIC:
                return new PdfPGothicItalicFontDef;
            case PDF_FONT_PGOTHIC_BOLD_ITALIC:
                return new PdfPGothicBoldItalicFontDef;
            case PDF_FONT_MINCYO:
                return new PdfMincyoFontDef();
            case PDF_FONT_MINCYO_BOLD:
                return new PdfMincyoBoldFontDef();
            case PDF_FONT_MINCYO_ITALIC:
                return new PdfMincyoItalicFontDef();
            case PDF_FONT_MINCYO_BOLD_ITALIC:
                return new PdfMincyoBoldItalicFontDef();
            case PDF_FONT_PMINCYO:
                return new PdfPMincyoFontDef;
            case PDF_FONT_PMINCYO_BOLD:
                return new PdfPMincyoBoldFontDef;
            case PDF_FONT_PMINCYO_ITALIC:
                return new PdfPMincyoItalicFontDef;
            case PDF_FONT_PMINCYO_BOLD_ITALIC:
                return new PdfPMincyoBoldItalicFontDef;
        }
    } catch (...) {
       return NULL;
    }
    return NULL;
}

