/*
 * << H a r u -- Free PDF Library >> -- PdfXObjects.cpp
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

/*----- constants ------------------------------------------------------------*/

const char* PDF_COLOR_SPACE_NAMES[] = {"DeviceGray",
                                       "DeviceRGB",
                                       "DeviceCMYK",
                                       "CalGray",
                                       "CalRGB",
                                       "Lab",
                                       "ICCBased",
                                       "Separation",
                                       "DeviceN",
                                       "Indexed",
                                       "Pattern"};

/*----------------------------------------------------------------------------*/
/*----- PdfXObjectBase class -------------------------------------------------*/

PdfXObjectBase::PdfXObjectBase(PdfDoc* doc)
        : PdfStream(doc->Xref())
{
	memset(fName, 0x00, PDF_LIMIT_MAX_NAME + 1);

    try {
        AddElement("Type", new PdfName("XObject"));
    } catch (...) {}
}

void
PdfXObjectBase::SetName(const char* name)
{
	strncpy(fName, name, PDF_LIMIT_MAX_NAME);
}

PdfXObjectBase::~PdfXObjectBase()
{
}

/*----------------------------------------------------------------------------*/
/*----- PdfXObjectMgr class --------------------------------------------------*/

PdfXObjectMgr::PdfXObjectMgr(PdfXref* xref)
{
    assert(xref);

    fXref = xref;
    fList = NULL;
}

PdfXObjectMgr::~PdfXObjectMgr()
{
    delete fList;
}

void
PdfXObjectMgr::RegisterXObject(PdfXObject* xobject, const char* name)
{
    assert(xobject);

    fXref->AddObject(xobject);

    if (xobject->Name()[0] != 0x00)
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: PdfXObjectMgr::RegisterXObject "
                "--cannot register this object[%s].", xobject->Name());
        
    if (!xobject->IsValidObject())
        throw PdfException(PDF_RUNTIME_ERROR,
            "ERROR: PdfXObjectMgr::RegisterXObject --Invalid XObject.");

    if (GetXObject(name) != NULL)
        throw PdfException(PDF_RUNTIME_ERROR,   
            "ERROR: PdfXObjectMgr::RegisterXObject "
                "duplicate xobject registration[%s].", xobject->Name());

    PdfNumber *length = new PdfNumber(0);
    fXref->AddObject(length);
    xobject->AddElement("Length", length);

    if (fList == NULL)
        fList = new PdfList();
    fList->AddItem(xobject);
    xobject->SetName(name);
}

PdfXObject*
PdfXObjectMgr::GetXObject(const char* name)
{
    if (fList == NULL)
        return NULL;
    
    for (int i = 0; i < fList->CountItems(); i++) {
        PdfXObject* f = GetXObject(i);

        const char* fname = f->Name();
        if (strcmp(name, fname) == 0) {
            PDF_DEBUG_PRINT(("PdfXObjectMgr found xobject %s[%d]\n", fname, i));
            return f;
        }
    }
    return NULL;
}

/*----------------------------------------------------------------------------*/
/*----- PdfImage class -------------------------------------------------------*/

PdfImage::PdfImage(PdfDoc* doc)
            : PdfXObject(doc)
{
    AddElement("Subtype", new PdfName("Image"));
    fColorSpace = PDF_CS_EOF;
}

void
PdfImage::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    AddElement("Width", new PdfNumber(fWidth));
    AddElement("Height", new PdfNumber(fHeight));
    AddElement("BitsPerComponent", new PdfNumber(fBitsPerComponent));
    if (fColorSpace != PDF_CS_INDEXED) {
        AddElement("ColorSpace", 
                new PdfName(PdfGetColorSpaceName(fColorSpace)));
	}

    PDF_DEBUG_PRINT(("PdfImage::InternalWriteStream\n"));
    PdfStream::InternalWriteStream(out, e);
}

/*----------------------------------------------------------------------------*/
/*------ utilities -----------------------------------------------------------*/

const char*
PdfGetColorSpaceName(pdf_color_space cs)
{
    return PDF_COLOR_SPACE_NAMES[(int)cs];
}

/*----------------------------------------------------------------------------*/

