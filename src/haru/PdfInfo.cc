/*
 * << H a r u --free pdf library >> -- PdfInfo.cpp
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

const static pdf_date INIT_DATE = {0, 1, 1, 0, 0, 0, '-', 0, 0};

/*------ PdfInfo -------------------------------------------------------------*/

PdfInfo::PdfInfo(PdfXref* xref)
        : PdfDictionary(xref)
{
    PDF_DEBUG_PRINT(("PdfInfo::PdfInfo: this:[%x] \n", (int)this));
}

PdfInfo::PdfInfo(int objectID, PdfXref* xref)
        : PdfDictionary(objectID, xref)
{
    PDF_DEBUG_PRINT(("PdfInfo::PdfInfo: this:[%x] objectID:[%d]\n",
        (int)this, objectID));
}

void
PdfInfo::Init()
{
    GetXref()->AddObject(this);
}

pdf_date
PdfInfo::CreationDate()
{
    const char* tmpStr = GetTextValue("CreationDate");

    if (tmpStr != NULL)
        return StrToPdfDate(tmpStr);
    else
        return INIT_DATE;
}

pdf_date
PdfInfo::ModDate()
{
    const char* tmpStr = GetTextValue("ModDate");

    if (tmpStr != NULL)
        return StrToPdfDate(tmpStr);
    else
        return INIT_DATE;
}

void
PdfInfo::SetAuthor(const char* value, PdfEncodingDef* encoding)
{
    if (encoding == NULL) 
        AddElement("Author", new PdfText(value));
    else
        SetTextAsUnicode("Author", value, encoding);
}

void
PdfInfo::SetCreator(const char* value, PdfEncodingDef* encoding)
{
    if (encoding == NULL)
        AddElement("Creator", new PdfText(value));
    else
        SetTextAsUnicode("Creator", value, encoding);
}

void
PdfInfo::SetProducer(const char* value, PdfEncodingDef* encoding)
{
    if (encoding == NULL)
        AddElement("Producer", new PdfText(value));
    else
        SetTextAsUnicode("Producer", value, encoding);
}

void
PdfInfo::SetTitle(const char* value, PdfEncodingDef* encoding)
{
    if (encoding == NULL)
        AddElement("Title", new PdfText(value));
    else
        SetTextAsUnicode("Title", value, encoding);
}

void
PdfInfo::SetSubject(const char* value, PdfEncodingDef* encoding)
{
    if (encoding == NULL)
        AddElement("Subject", new PdfText(value));
    else
        SetTextAsUnicode("Subject", value, encoding);
}

void
PdfInfo::SetKeywords(const char* value, PdfEncodingDef* encoding)
{
    if (encoding == NULL)
        AddElement("Keywords", new PdfText(value));
    else
        SetTextAsUnicode("Keywords", value, encoding);
}

void
PdfInfo::SetCreationDate(pdf_date value)
{
    char tmpStr[24];

    if (PdfDateToStr(value, tmpStr, 24)) {
        PdfText *date = new PdfText(tmpStr);
        AddElement("CreationDate", date);
    }
}

void
PdfInfo::SetModDate(pdf_date value)
{
    char tmpStr[24];

    if (PdfDateToStr(value, tmpStr, 24)) {
        PdfText *date = new PdfText(tmpStr);
        AddElement("ModDate", date);
    }
}

pdf_date
PdfInfo::StrToPdfDate(const char* value)
{
    pdf_date tmpDate = INIT_DATE;

    if (strlen(value) < 16) 
        throw PdfException(PDF_INVALID_PARAMETER, 
                "error: PdfInfo::StrToPdfDate --invalid value[%s]", value);

    int ret = sscanf(value, "D:%4d%2d%2d%2d%2d%2d%c%2d'%2d'", &tmpDate.year,
        &tmpDate.month, &tmpDate.day, &tmpDate.hour, &tmpDate.minutes,
        &tmpDate.seconds, &tmpDate.ind, &tmpDate.off_hour,
        &tmpDate.off_minutes);
    if (ret < 6)
        throw PdfException(PDF_INVALID_PARAMETER, 
                "error: PdfInfo::StrToPdfDate --invalid value[%s]", value);

    if (tmpDate.year < 0 ||
        tmpDate.month < 1 || 12 < tmpDate.month ||
        tmpDate.day < 1 || 31 < tmpDate.day ||
        tmpDate.hour < 0 || 23 < tmpDate.hour ||
        tmpDate.minutes < 0 || 59 < tmpDate.minutes ||
        tmpDate.seconds < 0 || 59 < tmpDate.seconds ||
        (tmpDate.ind != '+' && tmpDate.ind != '-' && tmpDate.ind != 'z') ||
        tmpDate.off_hour < 0 || 23 < tmpDate.off_hour ||
        tmpDate.off_minutes < 0 || 59 < tmpDate.off_minutes) {
        throw PdfException(PDF_INVALID_PARAMETER, 
                "error: PdfInfo::StrToPdfDate --invalid value[%s]", value);
    }

    return tmpDate;
}

bool
PdfInfo::PdfDateToStr(const pdf_date date, char* value, int length)
{
    char tmpStr[24];

    memset(tmpStr, 0x00, 24);
#ifdef __WIN32__
    int ret = _snprintf(tmpStr, 24, "D:%04d%02d%02d%02d%02d%02d%c%02d'%02d'",
        date.year, date.month, date.day, date.hour, date.minutes,
        date.seconds, date.ind, date.off_hour, date.off_minutes);
#else 
    int ret = snprintf(tmpStr, 24, "D:%04d%02d%02d%02d%02d%02d%c%02d'%02d'",
        date.year, date.month, date.day, date.hour, date.minutes,
        date.seconds, date.ind, date.off_hour, date.off_minutes);
#endif
    memset(value, 0x00, length);
    strncpy(value, tmpStr, length - 1);

    if (ret != 23)
        throw PdfException(PDF_INVALID_PARAMETER, 
                "error: PdfInfo::PdfDateToStr --invalid date[%s:%d].", 
                tmpStr, ret);

    return true;
}

void
PdfInfo::SetTextAsUnicode(const char* key, const char* value, 
        PdfEncodingDef* encoding)
{
    if (encoding == NULL)
        throw PdfException(PDF_INVALID_PARAMETER, "ERROR: "
                "PdfOutlineItem::SetTextAsUnicode -- encoding object is null.");

    if (!encoding->IsValid())
        GetXref()->GetDoc()->RegisterObject(encoding);

    PdfUnicodeText* ut = new PdfUnicodeText(encoding);
    AddElement(key, ut);
    ut->SetText(value);
}
    
void
PdfInfo::SetMbText(const char* key, const char* value, PdfCMap* cmap)
{
    if (cmap == NULL)
        throw PdfException(PDF_INVALID_PARAMETER,
                "ERROR: PdfOutlineItem::SetMbText --cmap object is null.");

    if (!cmap->IsValid())
        GetXref()->GetDoc()->RegisterObject(cmap);
    
    PdfUnicodeText* ut = new PdfUnicodeText(cmap);
    AddElement(key, ut);
    ut->SetText(value);
} 

/*----------------------------------------------------------------------------*/

