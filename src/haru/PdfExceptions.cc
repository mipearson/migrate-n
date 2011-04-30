/*
 * << H a r u --free pdf library >> -- PdfStreams.cpp
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

#include <stdarg.h>
#include <stdio.h>
#include "libharu.h"

PdfException::PdfException(int code, const char* fmt, ...)
{
    va_list args;
    
    va_start(args, fmt);
#ifdef __WIN32__
    _vsnprintf(fErrorBuf, 512, fmt, args);
#else
    vsnprintf(fErrorBuf, 512, fmt, args);
#endif
    va_end(args);
    fCode = code;
}

const char*
#ifdef _NO_EXCEPT_LIST
PdfException::what() const
#else
PdfException::what() const throw()
#endif
{
    return fErrorBuf;
}

