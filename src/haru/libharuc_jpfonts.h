/*
 * << H a r u --free pdf library >> -- libharuc_jpfonts.h
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

#ifndef _LIBHARUC_JPFONTS_H
#define _LIBHARUC_JPFONTS_H

#include <stdlib.h>
#include "libharuc.h"

/*----------------------------------------------------------------------------*/
/*----- "C" interfaces -------------------------------------------------------*/

enum pdf_predefined_jp_cmap_enum {
    PDF_CMAP_90MSP_RKSJ_H,
    PDF_CMAP_90MS_RKSJ_H,
    PDF_CMAP_90MS_RKSJ_V,
    PDF_CMAP_EUC_H,
    PDF_CMAP_EUC_V
};
typedef enum pdf_predefined_jp_cmap_enum pdf_predefined_jp_cmap;

enum pdf_jp_font_enum {
    PDF_FONT_GOTHIC = 0,
    PDF_FONT_GOTHIC_BOLD,
    PDF_FONT_GOTHIC_ITALIC,
    PDF_FONT_GOTHIC_BOLD_ITALIC,
    PDF_FONT_PGOTHIC,
    PDF_FONT_PGOTHIC_BOLD,
    PDF_FONT_PGOTHIC_ITALIC,
    PDF_FONT_PGOTHIC_BOLD_ITALIC,
    PDF_FONT_MINCYO,
    PDF_FONT_MINCYO_BOLD,
    PDF_FONT_MINCYO_ITALIC,
    PDF_FONT_MINCYO_BOLD_ITALIC,
    PDF_FONT_PMINCYO,
    PDF_FONT_PMINCYO_BOLD,
    PDF_FONT_PMINCYO_ITALIC,
    PDF_FONT_PMINCYO_BOLD_ITALIC
};
typedef enum pdf_jp_font_enum pdf_jp_font;

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*-- Create CMap. --*/
    
pdf_cmap pdf_create_jp_cmap(pdf_predefined_jp_cmap map);

/*-- Create Japanese font definition objects. --*/

pdf_cid_type2_fontdef pdf_create_jp_fontdef(pdf_jp_font font);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LIBHARUC_JPFONTS */

