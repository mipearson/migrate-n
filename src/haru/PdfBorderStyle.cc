/*
 * << H a r u --free pdf library >> -- PdfBorderStyle.cpp
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

const static pdf_bs_subtype BS_DEF_SUBTYPE = PDF_BS_SOLID;
const static double BS_DEF_WIDTH = 1;
const unsigned int BS_DEF_DASH_ON = 3;

/*----- PdfBorderStype -------------------------------------------------------*/

PdfBorderStyle::PdfBorderStyle(PdfXref* xref)
            : PdfDictionary(xref)
{
    AddElement("Type", new PdfName("Border"));
    fSubtype = BS_DEF_SUBTYPE;
    fWidth = BS_DEF_WIDTH;
    fDashOn = BS_DEF_DASH_ON;
    fDashOff = 0;
    fPhase = 0;
}

void
PdfBorderStyle::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    if (fSubtype == BS_DEF_SUBTYPE)
        RemoveElement("S");
    else
        switch (fSubtype) {
            case PDF_BS_SOLID: 
                AddElement("S", new PdfName("S"));
                break;
            case PDF_BS_DASHED: 
                AddElement("S", new PdfName("D"));
                break;
            case PDF_BS_BEVELED: 
                AddElement("S", new PdfName("B"));
                break;
            case PDF_BS_INSET: 
                AddElement("S", new PdfName("I"));
                break;
            case PDF_BS_UNDERLINED: 
                AddElement("S", new PdfName("U"));
                break;
            default: 
                RemoveElement("S");
        }

    if (fWidth == BS_DEF_WIDTH)
        RemoveElement("W");
    else
        AddElement("W", new PdfReal(fWidth));

    if (fDashOn == BS_DEF_DASH_ON && fDashOff == 0 && fPhase == 0)
        RemoveElement("D");
    else {
        PdfArray* array = new PdfArray(GetXref());
        try {
            array->Add(new PdfNumber(fDashOn));
            if (fDashOff != 0 && fPhase != 0)
                array->Add(new PdfNumber(fDashOff));
            if (fPhase != 0)
                array->Add(new PdfNumber(fPhase));  
        } catch (PdfException& e) {
            delete array;
            throw;
        }
        AddElement("D", array);
    }
        
    PdfDictionary::InternalWriteStream(out, e);
}

void
PdfBorderStyle::SetDash(unsigned int on, unsigned int off, unsigned int phase)
{
    fDashOn = on;
    fDashOff = off;
    fPhase = phase;
}

/*----------------------------------------------------------------------------*/

