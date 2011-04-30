/*
 * << H a r u --free pdf library >> -- libharuc_cntfonts.h
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

#ifndef _LIBHARUC_CNTFONTS_H
#define _LIBHARUC_CNTFONTS_H

#include <stdlib.h>
#include "libharuc.h"

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*----- "C" interfaces -------------------------------------------------------*/

enum pdf_predefined_cnt_cmap_enum {
    PDF_CMAP_ETen_B5_H,
    PDF_CMAP_ETen_B5_V
};
typedef enum pdf_predefined_cnt_cmap_enum pdf_predefined_cnt_cmap;

enum pdf_cnt_font_enum {
    PDF_FONT_MING = 0,
    PDF_FONT_MING_BOLD,
    PDF_FONT_MING_ITALIC,
    PDF_FONT_MING_BOLD_ITALIC
};
typedef enum pdf_cnt_font_enum pdf_cnt_font;

extern "C" {
#endif /* __cplusplus */

/*-- Create CMap. --*/
    
pdf_cmap pdf_create_cnt_cmap(pdf_predefined_cnt_cmap map);

/*-- Create Chinese font definition objects. --*/

pdf_cid_type2_fontdef pdf_create_cnt_fontdef(pdf_cnt_font font);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LIBHARUC_CNTFONTS */

