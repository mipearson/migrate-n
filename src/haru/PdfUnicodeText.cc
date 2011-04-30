/*
 * << H a r u --free pdf library >> -- PdfUnicodeText.cc
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

#include "libharu.h"

/*------ PdfUnicodeText ------------------------------------------------------*/

PdfUnicodeText::PdfUnicodeText()
    : PdfBinary()
{
    fCMap = NULL;
    fDef = NULL;
    fValue = NULL;
}

PdfUnicodeText::PdfUnicodeText(PdfCMap* cmap)
    : PdfBinary()
{
    fCMap = cmap;
    fDef = NULL;
    fValue = NULL;
}

PdfUnicodeText::PdfUnicodeText(PdfEncodingDef* def)
    : PdfBinary()
{
    fCMap = NULL;
    fDef = def;
    fValue = NULL;
}

PdfUnicodeText::~PdfUnicodeText()
{
    PDF_DEBUG_PRINT(("++ [%x] fValue delete[].\n", (int)fValue));
    delete[] fValue;
}

void
PdfUnicodeText::SetText(const char* text)
{
    if (fCMap == NULL && fDef == NULL)
        throw PdfException(PDF_ERR_INVALID_CMAP_OBJECT, "ERROR: "
                "Encoding object is not set or is invalid.");
    
    int len = strlen(text) * 2 + PDF_UNICODE_HEADER_LEN;
    if (fValue != NULL) {
        delete[] fValue;
        fValue = NULL;
    }
    fValue = new unsigned char[len];

    if (fCMap != NULL)
        fCMap->ToUnicode(text, fValue, &len);
    else {
        fDef->ToUnicode(text, fValue, &len);
    }
    SetData(fValue, len, true);
}

const char*
PdfUnicodeText::GetText()
{
    // not implemented yet.
    return "";
}

/*----------------------------------------------------------------------------*/

