/*
 * << H a r u --free pdf library >> -- PdfCatalog.cpp
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

const static char* PAGE_LAYOUT_NAMES[] = {
    "SinglePage",
    "OneColumn",
    "TwoColumnLeft",
    "TwoColumnRight",
    NULL
};
const static char* PAGE_MODE_NAMES[] = {
    "UseNone",
    "UseOutlines",
    "UseThumbs",
    "FullScreen",
    NULL
};

/*----------------------------------------------------------------------------*/
/*----- PdfCatalog class -----------------------------------------------------*/

PdfCatalog::PdfCatalog(PdfXref *xref)
        : PdfDictionary(xref)
{
    fPageLayout = PDF_DEF_PAGE_LAYOUT;
    fPageMode = PDF_DEF_PAGE_MODE;
    fNonFullScreenPageMode = PDF_DEF_PAGE_MODE;
    fOpenAction = NULL;
    fViewerPreferences = 0;
}

PdfCatalog::PdfCatalog(int objectID, PdfXref *xref)
        : PdfDictionary(objectID, xref)
{
    fPageLayout = PDF_DEF_PAGE_LAYOUT;
    fPageMode = PDF_DEF_PAGE_MODE;
    fNonFullScreenPageMode = PDF_DEF_PAGE_MODE;
    fOpenAction = NULL;
}

PdfOutlineRoot*
PdfCatalog::Outlines()
{
    PdfOutlineRoot* root = (PdfOutlineRoot*)GetValue("Outlines");
    if (root == NULL) {
        root = new PdfOutlineRoot(GetXref());
        AddElement("Outlines", root);
    }
    return root;
}

void
PdfCatalog::Init()
{
    GetXref()->AddObject(this);
    SetType("Catalog");
}

PdfPages*
PdfCatalog::Pages()
{
    PdfPages *pages = (PdfPages*)GetValue("Pages");

    if (pages == NULL)
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: Pages element not found.");
    return pages;
}

void
PdfCatalog::SetOpenAction(PdfDestination* action)
{
    if (action != NULL && action->CheckValid() == false)
        throw PdfException(PDF_RUNTIME_ERROR, "ERROR: PdfCatalog::"
                "SetOpenAction --Invalid destination.");    
    if (fOpenAction != NULL && fOpenAction != action && 
            fOpenAction->GetObjectType() == PDF_OBJ_TYPE_DIRECT)
        delete fOpenAction;
    fOpenAction = action;
}

void
PdfCatalog::AddPageLabel(unsigned int pagenum, pdf_page_num_style style,
        unsigned int firstpage, const char* prefix)
{
    PdfPageLabel* label = new PdfPageLabel(GetXref());

    label->SetNumberingStyle(style);
    if (prefix != NULL)
        label->SetLabelPrefix(prefix);
    if (firstpage != 1)
        label->SetFirstPage(firstpage);

    InternalAddPageLabel(pagenum, label);
}

void
PdfCatalog::InternalAddPageLabel(unsigned int pagenum, PdfPageLabel* label)
{
    if (label == NULL)
        return;

    PdfDictionary* pagelabels = (PdfDictionary*)GetValue("PageLabels");

    if (pagelabels == NULL) {
        pagelabels = new PdfDictionary(GetXref());
        AddElement("PageLabels", pagelabels);
    }

    PdfArray* nums = (PdfArray*)pagelabels->GetValue("Nums");

    if (nums == NULL) {
        nums = new PdfArray(GetXref());
        pagelabels->AddElement("Nums", nums);
    }

    nums->Add(new PdfNumber(pagenum));
    nums->Add(label);
}

void
PdfCatalog::ClearPageLabel()
{
    RemoveElement("PageLabels");
}

void
PdfCatalog::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    if (fPageLayout == PDF_DEF_PAGE_LAYOUT)
        RemoveElement("PageLayout");
    else
        AddElement("PageLayout", new PdfName(PAGE_LAYOUT_NAMES[fPageLayout]));

    if (fPageMode == PDF_DEF_PAGE_MODE)
        RemoveElement("PageMode");
    else
        AddElement("PageMode", new PdfName(PAGE_MODE_NAMES[fPageMode]));

    if (fNonFullScreenPageMode == PDF_DEF_PAGE_MODE)
        RemoveElement("NonFullScreenPageMode");
    else
        AddElement("NonFullScreenPageMode", 
            new PdfName(PAGE_MODE_NAMES[fNonFullScreenPageMode]));

    if (fOpenAction == NULL)
        RemoveElement("OpenAction");
    else
        AddElement("OpenAction", fOpenAction);

    if (fViewerPreferences == 0) {
        RemoveElement("ViewerPreferences");
    } else {
        PdfDictionary *dict_preferences =
                (PdfDictionary*)GetValue("ViewerPreferences");
        if (dict_preferences == NULL) {
            dict_preferences = new PdfDictionary(GetXref());
            AddElement("ViewerPreferences", dict_preferences);
        }

        if (IsSetViewerPreference(PDF_HIDE_TOOLBAR))
            dict_preferences->AddElement("HideToolbar", new PdfBoolean(true));
        else
            dict_preferences->RemoveElement("HideToolbar");

        if (IsSetViewerPreference(PDF_HIDE_MENUBAR))
            dict_preferences->AddElement("HideMenubar", new PdfBoolean(true));
        else
            dict_preferences->RemoveElement("HideMenubar");

        if (IsSetViewerPreference(PDF_HIDE_WINDOW_UI))
            dict_preferences->AddElement("HideWindowUI", new PdfBoolean(true));
        else
            dict_preferences->RemoveElement("HideWindowUI");

        if (IsSetViewerPreference(PDF_FIT_WINDOW))
            dict_preferences->AddElement("FitWindow", new PdfBoolean(true));
        else
            dict_preferences->RemoveElement("FitWindow");

        if (IsSetViewerPreference(PDF_CENTER_WINDOW))
            dict_preferences->AddElement("CenterWindow", new PdfBoolean(true));
        else
            dict_preferences->RemoveElement("CenterWindow");
    }

    PdfDictionary::InternalWriteStream(out, e);
}

/*----------------------------------------------------------------------------*/
/*----- PdfPageLabel class ---------------------------------------------------*/

PdfPageLabel::PdfPageLabel(PdfXref *xref)
        : PdfDictionary(xref)
{
    fNumberingStyle = PDF_PAGE_NUM_DECIMAL;
    fLabelPrefix = NULL;
    fFirstPage = PDF_DEF_PAGE_NUM;
}

PdfPageLabel::~PdfPageLabel()
{
    delete[] fLabelPrefix;
}

void
PdfPageLabel::SetLabelPrefix(const char* value)
{
    if (fLabelPrefix != NULL) 
        delete[] fLabelPrefix;

    if (value == NULL) {
        fLabelPrefix = NULL;
        return;
    }   
    
    int len = strlen(value);
    fLabelPrefix = new char[len + 1];
    strncpy(fLabelPrefix, value, len + 1);
}

void
PdfPageLabel::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    switch (fNumberingStyle) {
        case PDF_PAGE_NUM_DECIMAL: 
            AddElement("S", new PdfName("D"));
            break;
        case PDF_PAGE_NUM_UPPER_ROMAN:
            AddElement("S", new PdfName("R"));
            break;
        case PDF_PAGE_NUM_LOWER_ROMAN: 
            AddElement("S", new PdfName("r"));
            break;
        case PDF_PAGE_NUM_UPPER_LETTERS: 
            AddElement("S", new PdfName("A"));
            break;
        case PDF_PAGE_NUM_LOWER_LETTERS:
            AddElement("S", new PdfName("a"));
            break;
    }
    
    if (fLabelPrefix == NULL)
        RemoveElement("P");
    else
        AddElement("P", new PdfText(fLabelPrefix));

    if (fFirstPage == PDF_DEF_PAGE_NUM)
        RemoveElement("St");
    else
        AddElement("St", new PdfNumber(fFirstPage));

    PdfDictionary::InternalWriteStream(out, e);
}

/*----------------------------------------------------------------------------*/

