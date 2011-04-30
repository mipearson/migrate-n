/*
 * << H a r u --free pdf library >> -- PdfAnnotations.cpp
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

const static char* PDF_ANNOT_TYPE_NAMES[] = {
                                        "Text",
                                        "Link",
                                        "Sound",
                                        "FreeText",
                                        "Stamp",
                                        "Square",
                                        "Circle",
                                        "StrikeOut",
                                        "Highlight",
                                        "Underline",
                                        "Ink",
                                        "FileAttachment",
                                        "Popup"
                                        };

const static char* PDF_ANNOT_ICON_NAMES_NAMES[] = {
                                        "Comment",
                                        "Key",
                                        "Note",
                                        "Help",
                                        "NewParagraph",
                                        "Paragraph",
                                        "Insert"
                                        };


/*----- PdfAnnotation --------------------------------------------------------*/

PdfAnnotation::PdfAnnotation(PdfXref* xref, pdf_annot_type subtype)
            : PdfDictionary(xref)
{
    fSubtype = subtype;
    AddElement("Type", new PdfName("Annot"));
    AddElement("Subtype", new PdfName(SubtypeText(fSubtype)));
}

void
PdfAnnotation::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    if (GetValue("Rect") == NULL)
        throw PdfException(PDF_RUNTIME_ERROR, "ERROR: PdfAnnotation::"
                "WriteToStream Rect element must be set.");
    PdfDictionary::InternalWriteStream(out, e);
}

void
PdfAnnotation::SetRect(pdf_rect rect)
{
    if (!PDF_CHECK_RECT(rect))
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: Invalid rect values (%f,%f,%f,%f).",
                rect.left, rect.bottom, rect.right, rect.top);

    fRect = rect;
    PdfArray *array = new PdfArray(GetXref());
    try {
        array->Add(rect);
    } catch (PdfException& e) {
        delete array;
    }
    AddElement("Rect", array);
}

const char*
PdfAnnotation::SubtypeText(pdf_annot_type subtype)
{
    return PDF_ANNOT_TYPE_NAMES[(int)subtype];
}

void
PdfAnnotation::SetBorder(double horiz, double vert, double width)
{
    if (horiz < 0 || vert < 0 || width < 0)
        throw PdfException(PDF_RUNTIME_ERROR, 
            "PdfLinkAnnot::SetBorder --invalid value [%f,%f,%f]", 
                horiz, vert, width);
    
    if (horiz == 0 && vert == 0 && width == 1) {
        RemoveElement("Border");
        return;
    }
    
    PdfArray* array = new PdfArray(GetXref());

    array->Add(new PdfReal(horiz)); 
    array->Add(new PdfReal(vert));  
    array->Add(new PdfReal(width));

    fBorder[0] = horiz;
    fBorder[1] = vert;
    fBorder[2] = width;

    AddElement("Border", array);
}

void
PdfAnnotation::GetBorder(double* horiz, double* vert, double* width)
{
    *horiz = fBorder[0];
    *vert = fBorder[1];
    *width = fBorder[2];
}

PdfBorderStyle*
PdfAnnotation::GetBorderStyle()
{
    PdfBorderStyle* bs = (PdfBorderStyle*)GetValue("BS");

    if (bs == NULL) {
        bs = new PdfBorderStyle(GetXref());
        AddElement("BS", bs);
    }
    
    return bs;
}

/*----------------------------------------------------------------------------*/
/*------ PdfLinkAnnot --------------------------------------------------------*/

void
PdfLinkAnnot::SetDest(PdfDestination* dest)
{
    if (dest == NULL)
        RemoveElement("Dest");
    else {
        if (dest->CheckValid() == false)
            throw PdfException(PDF_INVALID_PARAMETER, "PdfLinkAnnot::SetDest");
        AddElement("Dest", dest);
    }
}

PdfDestination*
PdfLinkAnnot::Dest()
{
    return (PdfDestination*)GetValue("Dest");
}

void
PdfLinkAnnot::SetHightlightMode(pdf_annot_hl_mode mode)
{
    switch (mode) {
        case PDF_ANNOT_NO_HIGHTLIGHT:
            AddElement("H", new PdfName("N"));
            break;
        case PDF_ANNOT_INVERT_BOX:
            /* default value */
            RemoveElement("H");
            /* AddElement("H", new PdfName("I")); */
            break;
        case PDF_ANNOT_INVERT_BORDER:
            AddElement("H", new PdfName("O"));
            break;
        case PDF_ANNOT_DOWN_APPEARANCE:
            AddElement("H", new PdfName("P"));
            break;
        default:
            throw PdfException(PDF_RUNTIME_ERROR, "Invalid Annotation"
                    "HighlightMode[%d].", (int)mode);
    }
}
/*----------------------------------------------------------------------------*/
/*------ PdfTextAnnot --------------------------------------------------------*/

PdfTextAnnot::PdfTextAnnot(PdfXref* xref)
    : PdfAnnotation(xref, PDF_ANNOT_TEXT_NOTES)
{
    fContents = NULL;
    fOpened = false;
    fIcon = PDF_ANNOT_ICON_EOF;
}

void
PdfTextAnnot::SetContents(const char* text, PdfEncodingDef* encoding)
{
    if (encoding == NULL) {
        PdfText *t = new PdfText();
        t->SetValue(text);
        AddElement("Contents", t);
    } else {
        if (!encoding->IsValid())
            GetXref()->GetDoc()->RegisterObject(encoding);
        
        PdfUnicodeText* ut = new PdfUnicodeText(encoding);
        AddElement("Contents", ut);
        ut->SetText(text);
    }
}

void
PdfTextAnnot::SetContentsMb(const char* text, PdfCMap* cmap)
{
    if (text == NULL || strlen(text) == 0)
        RemoveElement("Contents");

    PdfUnicodeText *ut = new PdfUnicodeText(cmap);
    AddElement("Contents", ut);
    ut->SetText(text);   
}

void
PdfTextAnnot::InternalWriteStream(PdfStreamBase* out,
        PdfEncryptor* e)
{
    if (fContents == NULL)
        fContents = new PdfText("");

    PdfBoolean* opened = new PdfBoolean(fOpened);
    AddElement("Open", opened);

    if (fIcon >= 0 && fIcon < PDF_ANNOT_ICON_EOF) {
        PdfName* icon = new PdfName(PDF_ANNOT_ICON_NAMES_NAMES[(int)fIcon]);
        AddElement("Name", icon);
    }

    PdfAnnotation::InternalWriteStream(out, e);
}

