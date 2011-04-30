/*
 * << H a r u --free pdf library >> -- libharuc_krfonts.h
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

/*----------------------------------------------------------------------------*/

#ifndef _LIBHARUC_KRFONTS_H
#define _LIBHARUC_KRFONTS_H

#include <stdlib.h>
#include "libharuc.h"

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*----- "C" interfaces -------------------------------------------------------*/

enum pdf_predefined_kr_cmap_enum {
    PDF_CMAP_KSC_MS_UHC_H,
    PDF_CMAP_KSC_MS_UHC_V,  /* not supported. */
    PDF_CMAP_KSC_EUC_H,
    PDF_CMAP_KSC_EUC_V      /* not supported. */
};
typedef enum pdf_predefined_kr_cmap_enum pdf_predefined_kr_cmap;

enum pdf_kr_font_enum {
    PDF_FONT_BATANG_CHE = 0,
    PDF_FONT_BATANG_CHE_BOLD,
    PDF_FONT_BATANG_CHE_ITALIC,
    PDF_FONT_BATANG_CHE_BOLD_ITALIC,
    PDF_FONT_DOTUM_CHE,
    PDF_FONT_DOTUM_CHE_BOLD,
    PDF_FONT_DOTUM_CHE_ITALIC,
    PDF_FONT_DOTUM_CHE_BOLD_ITALIC
};
typedef enum pdf_kr_font_enum pdf_kr_font;

extern "C" {
#endif /* __cplusplus */

/*-- Create CMap. --*/
    
pdf_cmap pdf_create_kr_cmap(pdf_predefined_kr_cmap map);

/*-- Create Japanese font definition objects. --*/

pdf_cid_type2_fontdef pdf_create_kr_fontdef(pdf_kr_font font);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LIBHARUC_KRFONTS */

