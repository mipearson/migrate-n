/*
 * << H a r u --free pdf library >> -- libharu_MSCP.h
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

#ifndef _LIBHARUC_MSCP_H
#define _LIBHARUC_MSCP_H

#include <stdlib.h>
#include "libharuc.h"

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*----- "C" interfaces -------------------------------------------------------*/

extern "C" {
#endif /* __cplusplus */

pdf_encodingdef pdf_create_CP1251_encoding();

pdf_encodingdef pdf_create_CP1252_encoding();

pdf_encodingdef pdf_create_CP1253_encoding();

pdf_encodingdef pdf_create_CP1254_encoding();

pdf_encodingdef pdf_create_CP1255_encoding();

pdf_encodingdef pdf_create_CP1256_encoding();

pdf_encodingdef pdf_create_CP1257_encoding();

pdf_encodingdef pdf_create_CP1258_encoding();

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LIBHARUC_KOI8_H */

