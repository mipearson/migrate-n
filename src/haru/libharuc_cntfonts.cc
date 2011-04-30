/*
 * << H a r u --free pdf library >> -- libharuc_cntfonts.cc
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

#include "libharuc_cntfonts.h"
#include "libharu_cntfonts.h"

/*-- Create CMap. --*/
    
pdf_cmap pdf_create_cnt_cmap(pdf_predefined_cnt_cmap map)
{
    try {
        switch (map) {
            case PDF_CMAP_ETen_B5_H:
                return new PdfCMap_ETen_B5_H();
            case PDF_CMAP_ETen_B5_V:
                return NULL;
        }
    } catch (...) {
        return NULL;
    }
    return NULL;
}

/*-- Create Chinese font definition objects. --*/

pdf_cid_type2_fontdef pdf_create_cnt_fontdef(pdf_cnt_font font) {
    try {
        switch (font) {
            case PDF_FONT_MING:
                return new PdfMingFontDef();
            case PDF_FONT_MING_BOLD:
                return new PdfMingBoldFontDef();
            case PDF_FONT_MING_ITALIC:
                return new PdfMingItalicFontDef();
            case PDF_FONT_MING_BOLD_ITALIC:
                return new PdfMingBoldItalicFontDef();
        }
    } catch (...) {
       return NULL;
    }
    return NULL;
}

