/*
 * << H a r u --free pdf library >> -- libharuc_cnsfonts.cc
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

#include "libharuc_cnsfonts.h"
#include "libharu_cnsfonts.h"

/*-- Create CMap. --*/
    
pdf_cmap pdf_create_cns_cmap(pdf_predefined_cns_cmap map)
{
    try {
        switch (map) {
            case PDF_CMAP_GB_EUC_H:
                return new PdfCMap_GB_EUC_H();
            case PDF_CMAP_GB_EUC_V:
                return NULL;
        }
    } catch (...) {
        return NULL;
    }
    return NULL;
}

/*-- Create Chinese font definition objects. --*/

pdf_cid_type2_fontdef pdf_create_cns_fontdef(pdf_cns_font font) {
    try {
        switch (font) {
            case PDF_FONT_SIM_SUN:
                return new PdfSimSunFontDef();
            case PDF_FONT_SIM_SUN_BOLD:
                return new PdfSimSunBoldFontDef();
            case PDF_FONT_SIM_SUN_ITALIC:
                return new PdfSimSunItalicFontDef();
            case PDF_FONT_SIM_SUN_BOLD_ITALIC:
                return new PdfSimSunBoldItalicFontDef();
            case PDF_FONT_SIM_HEI:
                return new PdfSimHeiFontDef();
            case PDF_FONT_SIM_HEI_BOLD:
                return new PdfSimHeiBoldFontDef();
            case PDF_FONT_SIM_HEI_ITALIC:
                return new PdfSimHeiItalicFontDef();
            case PDF_FONT_SIM_HEI_BOLD_ITALIC:
                return new PdfSimHeiBoldItalicFontDef();
        }
    } catch (...) {
       return NULL;
    }
    return NULL;
}

