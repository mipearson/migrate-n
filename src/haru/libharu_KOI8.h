/*
 * << H a r u --free pdf library >> -- libharu_KOI8.h
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

#ifndef _HARU_PDF_KOI8_H
#define _HARU_PDF_KOI8_H

#include <stdlib.h>
#include "libharu.h"

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*------ PdfEncoding_KOI8_R encoding --------------------------------------*/

class PdfEncoding_KOI8_R : public PdfWinAnsiEncoding
{
public:
        void            Init();
};

/*----------------------------------------------------------------------------*/

#endif /* __cplusplus */

#endif /* _HARU_PDF_KOI8_H */

