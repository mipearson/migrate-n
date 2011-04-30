/*
 * << H a r u --free pdf library >> -- PdfDestination.cpp
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

#include <assert.h>
#include "libharu.h"

/*----------------------------------------------------------------------------*/

const char* PDF_DESTINATION_TYPE_NAMES[] = {"XYZ",
                                            "Fit",
                                            "FitH",
                                            "FitV",
                                            "FitR",
                                            "FitB",
                                            "FitBH",
                                            "FitBV"};

/*----------------------------------------------------------------------------*/

pdf_destination_type
StrToDestinationType(const char* type)
{
    for (int i = 0; i < PDF_DST_EOF; i++) {
        if (strcmp(type, PDF_DESTINATION_TYPE_NAMES[i]) == 0)
            return (pdf_destination_type)i;
    }
    return PDF_DST_EOF;
}

const char*
DestinationTypeToSTr(pdf_destination_type type)
{
    if (PDF_XYZ <= type && type < PDF_DST_EOF)
        return PDF_DESTINATION_TYPE_NAMES[(int)type];
    else
        return NULL;
}

/*----------------------------------------------------------------------------*/
/*----- PdfDestination class -------------------------------------------------*/


PdfDestination::PdfDestination(PdfPage* page, bool create_as_shared)
    : PdfArray(page->GetXref())
{
    assert(page);
    fPage = page;
    fShared = create_as_shared;
}

void
PdfDestination::CheckShared()
{
    if (fShared && GetObjectType() == PDF_OBJ_TYPE_DIRECT)
        GetXref()->AddObject(this);
}

bool
PdfDestination::CheckValid()
{
    if (fPage == NULL)
        return false;
    if (GetCount() == 0)
        SetFit();
    CheckShared();
    return true;
}

pdf_destination_type
PdfDestination::Type()
{
    PdfObject* type = GetItem(1);

    if (type == NULL || type->GetClass() == ocName) {
        const char* type_name = ((PdfName*)type)->GetValue();
        return StrToDestinationType(type_name);
    } else
        return PDF_DST_EOF;
}

double
PdfDestination::GetParam(int index)
{
    PdfObject* param = GetItem(index + 2);
    if (param != NULL)
        return ((PdfReal*)param)->GetValue();
    else
        return PDF_PARAM_NODEF;
}

void
PdfDestination::SetXYZ(double left, double top, double zoom)
{
    if (zoom < 0 || zoom > PDF_MAX_ZOOMSIZE)
        zoom = 1;

    CheckShared();
    Clear();
    Add(fPage);
    Add(new PdfName(DestinationTypeToSTr(PDF_XYZ)));
    Add(new PdfReal(left));
    Add(new PdfReal(top));
    Add(new PdfReal(zoom));
}

void
PdfDestination::SetFit()
{
    CheckShared();
    Clear();
    Add(fPage);
    Add(new PdfName(DestinationTypeToSTr(PDF_FIT)));
}

void
PdfDestination::SetFitH(double top)
{
    CheckShared();
    Clear();
    Add(fPage);
    Add(new PdfName(DestinationTypeToSTr(PDF_FIT_H)));
    Add(new PdfReal(top));
}
void
PdfDestination::SetFitV(double left)
{
    CheckShared();
    Clear();
    Add(fPage);
    Add(new PdfName(DestinationTypeToSTr(PDF_FIT_V)));
    Add(new PdfReal(left));
}

void
PdfDestination::SetFitR(double left, double bottom, double right, double top)
{
    CheckShared();
    Clear();
    Add(fPage);
    Add(new PdfName(DestinationTypeToSTr(PDF_FIT_R)));
    Add(new PdfReal(left));
    Add(new PdfReal(bottom));
    Add(new PdfReal(right));
    Add(new PdfReal(top));
}

void
PdfDestination::SetFitB()
{
    CheckShared();
    Clear();
    Add(fPage);
    Add(new PdfName(DestinationTypeToSTr(PDF_FIT_B)));
}

void
PdfDestination::SetFitBH(double top)
{
    CheckShared();
    Clear();
    Add(fPage);
    Add(new PdfName(DestinationTypeToSTr(PDF_FIT_BH)));
    Add(new PdfReal(top));
}
void
PdfDestination::SetFitBV(double left)
{
    CheckShared();
    Clear();
    Add(fPage);
    Add(new PdfName(DestinationTypeToSTr(PDF_FIT_BV)));
    Add(new PdfReal(left));
}

/*----------------------------------------------------------------------------*/

