/*
 * << H a r u --free pdf library >> -- libharuc_cnsfonts.h
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

#ifndef _LIBHARUC_CNSFONTS_H
#define _LIBHARUC_CNSFONTS_H

#include <stdlib.h>
#include "libharuc.h"

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*----- "C" interfaces -------------------------------------------------------*/

enum pdf_predefined_cns_cmap_enum {
    PDF_CMAP_GB_EUC_H,
    PDF_CMAP_GB_EUC_V
};
typedef enum pdf_predefined_cns_cmap_enum pdf_predefined_cns_cmap;

enum pdf_cns_font_enum {
    PDF_FONT_SIM_SUN = 0,
    PDF_FONT_SIM_SUN_BOLD,
    PDF_FONT_SIM_SUN_ITALIC,
    PDF_FONT_SIM_SUN_BOLD_ITALIC,
    PDF_FONT_SIM_HEI,
    PDF_FONT_SIM_HEI_BOLD,
    PDF_FONT_SIM_HEI_ITALIC,
    PDF_FONT_SIM_HEI_BOLD_ITALIC
};
typedef enum pdf_cns_font_enum pdf_cns_font;

extern "C" {
#endif /* __cplusplus */

/*-- Create CMap. --*/
    
pdf_cmap pdf_create_cns_cmap(pdf_predefined_cns_cmap map);

/*-- Create Chinese font definition objects. --*/

pdf_cid_type2_fontdef pdf_create_cns_fontdef(pdf_cns_font font);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LIBHARUC_CNSFONTS */

