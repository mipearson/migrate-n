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

#ifndef _HARU_PDF_MSCP_H
#define _HARU_PDF_MSCP_H

#include <stdlib.h>
#include "libharu.h"

#ifdef __cplusplus

/*----------------------------------------------------------------------------*/
/*------ PdfEncoding_CP1251 encoding --------------------------------------*/

class PdfEncoding_CP1251 : public PdfWinAnsiEncoding
{
public:
        void            Init();
};

/*----------------------------------------------------------------------------*/
/*------ PdfEncoding_CP1252 encoding --------------------------------------*/

class PdfEncoding_CP1252 : public PdfWinAnsiEncoding
{
public:
        void            Init();
};

/*----------------------------------------------------------------------------*/
/*------ PdfEncoding_CP1253 encoding --------------------------------------*/

class PdfEncoding_CP1253 : public PdfWinAnsiEncoding
{
public:
        void            Init();
};

/*----------------------------------------------------------------------------*/
/*------ PdfEncoding_CP1254 encoding --------------------------------------*/

class PdfEncoding_CP1254 : public PdfWinAnsiEncoding
{
public:
        void            Init();
};

/*----------------------------------------------------------------------------*/
/*------ PdfEncoding_CP1255 encoding --------------------------------------*/

class PdfEncoding_CP1255 : public PdfWinAnsiEncoding
{
public:
        void            Init();
};

/*----------------------------------------------------------------------------*/
/*------ PdfEncoding_CP1256 encoding --------------------------------------*/

class PdfEncoding_CP1256 : public PdfWinAnsiEncoding
{
public:
        void            Init();
};

/*----------------------------------------------------------------------------*/
/*------ PdfEncoding_CP1257 encoding --------------------------------------*/

class PdfEncoding_CP1257 : public PdfWinAnsiEncoding
{
public:
        void            Init();
};

/*----------------------------------------------------------------------------*/
/*------ PdfEncoding_CP1258 encoding --------------------------------------*/

class PdfEncoding_CP1258 : public PdfWinAnsiEncoding
{
public:
        void            Init();
};

/*----------------------------------------------------------------------------*/

#endif /* __cplusplus */

#endif /* _HARU_PDF_MSCP_H */

