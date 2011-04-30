/*
 * << H a r u --free pdf library >> -- PdfPages.cpp
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
/*----- PdfPageBase class ----------------------------------------------------*/

PdfPageBase::PdfPageBase(PdfXref *xref)
        : PdfDictionary(xref)
{
    fParent = NULL;
    fProcSet = 0;
}

PdfPageBase::~PdfPageBase()
{
}

void
PdfPageBase::Init()
{
    GetXref()->AddObject(this);
}

PdfObject*
PdfPageBase::FindElement(const char* name)
{
    PdfObject *ret = NULL;
    PdfPageBase *page = this;

    /* find element recursivvely */
    while (ret == NULL) {
        ret = page->GetValue(name);
        if (ret == NULL) {
            page = page->Parent();
            if (page == NULL) {
                PDF_DEBUG_PRINT(("WARNING: PdfPageBase::FindElement "
                    "%s not found.\n", name));
                return(NULL);
            }
        }
    }
    return(ret);
}

pdf_box
PdfPageBase::GetElementRect(const char* name)
{
    pdf_box ret = {0, 0, 0, 0};

    PdfObject *rect = FindElement(name);
    if (rect == NULL || rect->GetClass() != ocArray) {
        PDF_DEBUG_PRINT(("WARNING: PdfPage::GetElementRect %s invalid.\n",
            name));
        return ret;
    }
    ret.left = ((PdfArray*)rect)->GetAsInteger(0);
    ret.bottom = ((PdfArray*)rect)->GetAsInteger(1);
    ret.right = ((PdfArray*)rect)->GetAsInteger(2);
    ret.top = ((PdfArray*)rect)->GetAsInteger(3);
    return ret;
}

void
PdfPageBase::SetElementRect(const char* name, pdf_box rect)
{
    PdfArray *newRect = new PdfArray(GetXref());

    newRect->Add(new PdfNumber(rect.left));
    newRect->Add(new PdfNumber(rect.bottom));
    newRect->Add(new PdfNumber(rect.right));
    newRect->Add(new PdfNumber(rect.top));
    AddElement(name, newRect);
}

pdf_box
PdfPageBase::MediaBox()
{
    return GetElementRect("MediaBox");
}

void
PdfPageBase::SetMediaBox(pdf_box rect)
{
    SetElementRect("MediaBox", rect);
}

pdf_box
PdfPageBase::CropBox()
{
    return GetElementRect("CropBox");
}

void
PdfPageBase::SetCropBox(pdf_box rect)
{
    SetElementRect("CropBox", rect);
}

PdfDictionary*
PdfPageBase::Resources()
{
    PdfObject *obj = FindElement("Resources");

    if (obj == NULL) {
        PDF_DEBUG_PRINT(("DEBUG: Adding resource at[%d].\n", GetObjectID()));
        obj = new PdfDictionary(this->GetXref());
        SetResources((PdfDictionary*)obj);
    }

    assert(obj->GetClass() == ocDictionary);
    return (PdfDictionary*)obj;
}

PdfDictionary*
PdfPageBase::GetResource(const char* element_name)
{
    /* find element of resource dictionary from page(s) objects */

    PdfDictionary* res = Resources();
    assert(res);

    PdfObject* element = res->GetValue(element_name);
    if (element == NULL) {

        if (strcmp("Font", element_name) != 0 &&
            strcmp("ColorSpace", element_name) != 0 &&
            strcmp("XObject", element_name) != 0 &&
            strcmp("ExtGState", element_name) != 0 &&
            strcmp("Pattern", element_name) != 0 &&
            strcmp("Properties", element_name) != 0 &&
            strcmp("Shading", element_name) != 0) {
            throw PdfException(PDF_RUNTIME_ERROR,
                    "ERROR: Invalid resource element name "
                    "[%s].", element_name);
        }

        PDF_DEBUG_PRINT(("DEBUG: Adding %s element at[%d].\n",
            element_name, GetObjectID()));
        element = new PdfDictionary(this->GetXref());
        res->AddElement(element_name, element);
    }
    return (PdfDictionary*)element;
}

void
PdfPageBase::SetResources(PdfDictionary *res)
{
    if (res == NULL) {
        RemoveElement("Resources");
        return;
    }

    if (res->GetClass() != ocDictionary)
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: PdfPage::SetResources invarid resources");

    AddElement("Resources", res);
}

int
PdfPageBase::Rotate()
{
    PdfObject *obj = FindElement("Rotate");
    if (obj == NULL || obj->GetClass() != ocNumber) {
        PDF_DEBUG_PRINT(("WARNING: PdfPageBase::GetRotate not found/n"));
        return 0;
    }
    return(((PdfNumber*)obj)->GetValue());
}

void
PdfPageBase::SetRotate(int rotate)
{
    if (rotate > 360 || rotate < 0 || rotate % 90 != 0)
        throw PdfException(PDF_RUNTIME_ERROR, "ERROR: PdfPageBase::SetRotate "
                "invarid value[%d].", rotate);

    PdfNumber *obj = new PdfNumber(rotate);
    AddElement("Rotate", obj);
}

int
PdfPageBase::CountFonts()
{
    PdfDictionary* font_res = GetResource("Font");
    assert(font_res);

    return font_res->GetCount();
}

bool
PdfPageBase::AddFont(PdfFont *font, const char* fontname)
{
    assert(font);
    assert(fontname);

    PdfDictionary* font_res = GetResource("Font");
    assert(font_res);

    if (font_res->GetValue(fontname) != NULL)
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: duplicate font-name[%s].", fontname);

    font_res->AddElement(fontname, font);
    return true;
}

const char*
PdfPageBase::GetFontName(PdfFont *font)
{
    PdfDictionary* font_res = GetResource("Font");
    assert(font_res);

    for (int i = 0; i < font_res->GetCount(); i++) {
        const char* key_name = font_res->GetKeyValue(i);
        PdfObject* obj = font_res->GetValue(key_name);
        if (obj == font)
            return key_name;
    }

    return NULL;
}

bool
PdfPageBase::AddXObject(PdfXObject* xobject, const char* name)
{
    assert(xobject);
    assert(name);

    PdfDictionary* xobj_res = GetResource("XObject");
    assert(xobj_res);

    if (xobj_res->GetValue(name) != NULL)
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: duplicate XObject-name[%s].", name);

    const char* subtype = xobject->GetNameValue("Subtype");
    
    if (name != NULL && strcmp(subtype, "Image") == 0) {
        switch (((PdfImage*)xobject)->ColorSpace()) {
            case PDF_CS_DEVICE_GRAY:
            case PDF_CS_CAL_GRAY: 
                AddProcSet(PDF_PROCSET_IMAGEB);
                break;
            case PDF_CS_INDEXED: 
                AddProcSet(PDF_PROCSET_IMAGEI);
                break;
            default: 
                AddProcSet(PDF_PROCSET_IMAGEC);
        }
    }
    
    xobj_res->AddElement(name, xobject);
    return true;
}

const char*
PdfPageBase::GetXObjectName(PdfXObject* xobject)
{
    PdfDictionary* xobj_res = GetResource("XObject");
    assert(xobj_res);

    for (int i = 0; i < xobj_res->GetCount(); i++) {
        const char* key_name = xobj_res->GetKeyValue(i);
        PdfObject* obj = xobj_res->GetValue(key_name);
        if (xobject == obj)
            return key_name;
    }

    return NULL;
}

int
PdfPageBase::GetIndex()
{
    if (Parent() != NULL)
        return Parent()->IndexOf(this);

    return -1;
}

void
PdfPageBase::SetSize(int width, int height)
{
    if (width <= 0 || height <= 0) {
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: invalid page size(%d,%d).", width, height);
    }

    SetMediaBox(PdfBox(0, 0, width, height));
}

void
PdfPageBase::InternalAddProcSet(int procset) 
{
    PdfArray *array = (PdfArray*)Resources()->GetValue("ProcSet");
    if (array == NULL) {
        array = new PdfArray(this->GetXref());
        Resources()->AddElement("ProcSet", array);
    }
    PDF_DEBUG_PRINT(("PdfPageBase::InternalAddProcSet (%d)\n", procset));
    
    if ((procset & PDF_PROCSET_PDF) == PDF_PROCSET_PDF)
        array->Add(new PdfName("PDF"));
    if ((procset & PDF_PROCSET_TEXT) == PDF_PROCSET_TEXT)
        array->Add(new PdfName("Text"));
    if ((procset & PDF_PROCSET_IMAGEB) == PDF_PROCSET_IMAGEB)
        array->Add(new PdfName("ImageB"));
    if ((procset & PDF_PROCSET_IMAGEC) == PDF_PROCSET_IMAGEC)
        array->Add(new PdfName("ImageC"));
    if ((procset & PDF_PROCSET_IMAGEI) == PDF_PROCSET_IMAGEI)
        array->Add(new PdfName("ImageI"));

    fProcSet |= procset;
}

/*----------------------------------------------------------------------------*/
/*----- PdfPages class -------------------------------------------------------*/

PdfPages::PdfPages(PdfXref* xref, PdfFontMgr* fmgr, PdfXObjectMgr* omgr)
        : PdfPageBase(xref)
{
    fFontMgr = fmgr;
    fXObjectMgr = omgr;
    fPagesType = ptHasNoKids;
}

void
PdfPages::SetParent(PdfPages* parent)
{
    if (fParent != NULL) 
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: this page[%d] has parent object yet.", (long)this);

    AddElement("Parent", parent);
    fParent = parent;
    fFontMgr = fParent->FontMgr();
    fXObjectMgr = fParent->XObjectMgr();
}

bool
PdfPages::AddKids(PdfPageBase* page)
{
    PdfArray *fKids;

    assert(page != NULL);

    fKids = Kids();
    fKids->Add(page);
    page->AddElement("Parent", this);

    return true;
}

bool
PdfPages::InsertKids(PdfPageBase* page, int index)
{
    PdfArray *fKids;

    assert(page != NULL);

    fKids = Kids();
    fKids->Insert(page, index);

    page->AddElement("Parent", this);

    return true;
}

void
PdfPages::Init()
{
    SetType("Pages");
    AddElement("Kids", new PdfArray(GetXref()));
    PdfPageBase::Init();
}

PdfPage*
PdfPages::AddPage(int index)
{
    /* check the type of the page is collect. */
    if (index < -1 ||
        index > Kids()->GetCount() ||
        fPagesType == ptHasPages) {
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: Invalid operation cannot add Page.");
    } else
       fPagesType = ptHasPage;

    PdfPage* page = new PdfPage(GetXref());
    page->Init();

    if (index < 0)
        AddKids(page);
    else
        InsertKids(page, index);
    page->SetParent(this);
    GetXref()->GetDoc()->SetCurrentPage(page);
    GetXref()->GetDoc()->SetCurrentPages(this);

    return page;
}

PdfPages*
PdfPages::AddPages(int index)
{
    /* check the type of the page is collect. */
    if (index < -1 ||
        index > Kids()->GetCount() || fPagesType == ptHasPage) {
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: Invalid operation cannot add Pages.");
    } else
        fPagesType = ptHasPages;

    PdfPages* pages = new PdfPages(GetXref());
    pages->Init();

    if (index < 0)
        AddKids(pages);
    else
        InsertKids(pages, index);
    pages->SetParent(this);
    GetXref()->GetDoc()->SetCurrentPages(pages);

    return pages;
}

unsigned int
PdfPages::GetKidsCount()
{
    PdfArray* kids = Kids();
    if (kids == NULL)
        return 0;
    else
        return kids->GetCount();
}

int
PdfPages::IndexOf(PdfPageBase* page)
{
    PdfArray* kids = Kids();
    if (kids == NULL)
        return -1;
    else
        return kids->IndexOf(page);
}

unsigned int
PdfPages::GetPageCount()
{
     PdfArray* kids = Kids();
     unsigned int count = 0;

     if (kids == NULL)
        return 0;
     else {
        for (int i = 0; i < kids->GetCount(); i++) {
            PdfPageBase* page = (PdfPageBase*)kids->GetItem(i);
            count += page->GetPageCount();
        }
     }
     return count;
}

void
PdfPages::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    AddElement("Count", new PdfNumber(GetPageCount()));
    PdfDictionary::InternalWriteStream(out, e);
}

/*----------------------------------------------------------------------------*/
/*----- PdfPage class --------------------------------------------------------*/

void
PdfPage::SetParent(PdfPages *parent)
{
    assert(fParent == NULL);

    SetType("Page");
    AddElement("Parent", parent);
    if (Resources() == NULL)
        AddElement("Resources", new PdfDictionary(GetXref()));

    fParent = parent;

    /* if media-box is not inherited, set the default value */
    if (FindElement("MediaBox") == NULL) {
        pdf_box box = {0, 0, PDF_DEFAULT_PAGE_WIDTH, PDF_DEFAULT_PAGE_HEIGHT};

        SetMediaBox(box);
    }
}

PdfContents*
PdfPage::Canvas()
{
    // if the page does not have a canvas. create a canvas.
    if (fCanvas == NULL) {
        fCanvas = new PdfContents(this);
        try {
            fCanvas->Init();
        } catch (...) {
            delete fCanvas;
            throw;
        }
    }

    AddProcSet(PDF_PROCSET_PDF);
    return fCanvas;
}

PdfLinkAnnot*
PdfPage::AddLink(pdf_rect rect, PdfDestination* dest, pdf_annot_hl_mode mode)
{
    PdfLinkAnnot* link = new PdfLinkAnnot(GetXref());

    AddAnnotation(link);
    link->SetRect(rect);
    link->SetDest(dest);
    if (mode != PDF_ANNOT_HL_EOF)
        link->SetHightlightMode(mode);

    return link;
}

PdfTextAnnot*
PdfPage::AddTextAnnot(pdf_rect rect)
{
    PdfTextAnnot* annot = new PdfTextAnnot(GetXref());

    AddAnnotation(annot);
    annot->SetRect(rect);

    return annot;
}

void
PdfPage::AddAnnotation(PdfAnnotation* annot)
{
    /* make annotation to indirect object. */
    if (annot->GetObjectType() != PDF_OBJ_TYPE_INDIRECT)
        GetXref()->AddObject(annot);

    PdfArray* array = (PdfArray*)GetValue("Annots");
    if (array == NULL) {
        array = new PdfArray(GetXref());
        AddElement("Annots", array);
    }
    array->Add(annot);
}

/*----------------------------------------------------------------------------*/

