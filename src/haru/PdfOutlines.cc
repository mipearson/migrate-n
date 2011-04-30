/*
 * << H a r u --free pdf library >> -- PdfOutlines.cpp
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
/*----- PdfOutlineBase class -------------------------------------------------*/

PdfOutlineBase::PdfOutlineBase(PdfXref* xref)
        : PdfDictionary(xref)
{
    assert(xref != NULL);
    PDF_DEBUG_PRINT(("PdfOutlineBase::PdfOutlineBase\n"));

    fOpened = false;
    xref->AddObject(this);
}

void
PdfOutlineBase::AddChild(PdfOutlineItem* item)
{
    assert(item != NULL);
    assert(item->Parent() == NULL);
    assert(item != this);

    if (First() == NULL)
        AddElement("First", item);

    PdfOutlineItem* last = Last();
    if (last != NULL) {
        last->AddElement("Next", item);
        item->AddElement("Prev", last);
    }
    AddElement("Last", item);
    item->AddElement("Parent", this);
}

int
PdfOutlineBase::ChildCount()
{
    PdfOutlineItem* item = First();
    int count = 0;

    while (item != NULL) {
        count++;
        if (item->Opened())
            count += item->ChildCount();
        item = item->Next();
    }

    return count;
}

void
PdfOutlineBase::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    int num_children = (int)ChildCount();
    if (num_children == 0)
        RemoveElement("Count");
    else
    if (fOpened)
        AddElement("Count", new PdfNumber(num_children));
    else
        AddElement("Count", new PdfNumber(-num_children));

    PDF_DEBUG_PRINT(("PdfOutlineBase::InternalWriteStream\n"));
    PdfDictionary::InternalWriteStream(out, e);
}

/*----------------------------------------------------------------------------*/
/*----- PdfOutlineRoot class -------------------------------------------------*/

PdfOutlineRoot::PdfOutlineRoot(PdfXref *xref)
        : PdfOutlineBase(xref)
{
    AddElement("Type", new PdfName("Outline"));
    fOpened = true;
}

/*----------------------------------------------------------------------------*/
/*----- PdfOutlineItem class -------------------------------------------------*/

PdfOutlineItem::PdfOutlineItem(PdfOutline *parent)
        : PdfOutlineBase(parent->GetXref())
{
    AddElement("Type", new PdfName("Outline"));
    PDF_DEBUG_PRINT(("PdfOutlineItem::PdfOutlineItem.1\n"));
    parent->AddChild(this);
    PDF_DEBUG_PRINT(("PdfOutlineItem::PdfOutlineItem.2\n"));
}

const char*
PdfOutlineItem::Title()
{
    PdfText* title = (PdfText*)GetValue("Title");
    if (title == NULL)
        return NULL;
    else 
        return title->GetValue();
}

void
PdfOutlineItem::SetTitle(const char* title, PdfEncodingDef* encoding)
{
    if (encoding == NULL)
        AddElement("Title", new PdfText(title));
    else {
        if (!encoding->IsValid())
            GetXref()->GetDoc()->RegisterObject(encoding);

        PdfUnicodeText* ut = new PdfUnicodeText(encoding);
        AddElement("Title", ut);
        ut->SetText(title);
    }
}

void
PdfOutlineItem::SetTitleMb(const char* title, PdfCMap* cmap)
{
    if (cmap == NULL)
        throw PdfException(PDF_INVALID_PARAMETER, 
                "ERROR: PdfOutlineItem::SetTitleMb --cmap object is null.");

    if (!cmap->IsValid())
        GetXref()->GetDoc()->RegisterObject(cmap);
    
    PdfUnicodeText* ut = new PdfUnicodeText(cmap);
    AddElement("Title", ut);
    ut->SetText(title);
}

void
PdfOutlineItem::SetDestination(PdfDestination* dst)
{
    PDF_DEBUG_PRINT(("PdfOutlineItem::SetDestination\n"));
    if (dst == NULL) {
        RemoveElement("Dest");
        return;
    }

    if (dst->CheckValid() == false) 
        throw PdfException(PDF_RUNTIME_ERROR, "ERROR: PdfOutlineItem::"
                "SetDestination --Invalid destination.");

    AddElement("Dest", dst);
}

/*----------------------------------------------------------------------------*/

