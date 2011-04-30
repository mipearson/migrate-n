/*
 * << H a r u --free pdf library >> -- PdfMbFonts.cc
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
#include <errno.h>
#include "libharu.h"

/*----- PdfType0Font class ---------------------------------------------------*/

PdfType0Font::PdfType0Font(PdfXref* xref)
        : PdfFontBase(xref) 
{
    fCMap = NULL;
    fDescendantFont = NULL;
    fValid = false;
}

PdfType0Font::~PdfType0Font()
{
}

void
PdfType0Font::SetAttributes(const char* name, PdfCIDFont* font, PdfCMap* cmap)
{
    if (font == NULL || cmap == NULL)
        throw PdfException(PDF_INVALID_PARAMETER, 
                "error: PdfType0Font::SetAttributes().");

    if (font->GetObjectType() == PDF_OBJ_TYPE_DIRECT)
        throw PdfException(PDF_ERR_INVALID_OBJECT,
                "error: PdfType0Font::SetAttributes() "
                "--font must be indirect object.");
    
    if (fCMap != NULL) 
        throw PdfException(PDF_ERR_INVALID_OPERATION,
                "error: PdfType0Font::SetAttributes() "
                "-- cannot set cmap.");

    fName = new char[strlen(name) + 1];
    PDF_DEBUG_PRINT(("++ [%x] fName new.\n", (int)fName));
    strcpy(fName, name);
    
    fDescendantFont = font;
    fCMap = cmap;
    SetWritingMode(cmap->GetWritingMode());
    
    AddElement("Type", new PdfName("Font"));
    AddElement("Subtype", new PdfName("Type0"));
    AddElement("Encoding", new PdfName(cmap->GetCMapName()));
    AddElement("BaseFont", new PdfName(font->BaseFont()));

    PdfArray* array = new PdfArray(GetXref());
    AddElement("DescendantFonts", array);
    array->Add(font);
    fValid = true;
}

pdf_text_width
PdfType0Font::TextWidths(const char* text, unsigned int* widths)
{
    pdf_text_width ret;
    
    ret.width = 0;
    ret.numchars = 0;
    ret.numwords = 0;
    
    if (!fValid || text == NULL || (ret.numchars = strlen(text)) == 0)
        return ret;

    pdf_byte_type* btype = new pdf_byte_type[ret.numchars];
    PDF_DEBUG_PRINT(("++ [%x] pdf_byte_type new[].\n", (int)btype));
    try { 
        fCMap->ParseText(text, btype);
        const unsigned char* ptext = (unsigned char*)text;
        for (int i = 0; i < ret.numchars; i++, ptext++) {
            pdf_cid cid;
            unsigned int w = 0;
        
            if (btype[i] == PDF_BYTE_TYPE_SINGLE) {
                if (fWritingMode == PDF_WMODE_HORIZONTAL) {
                    cid = fCMap->GetCID((unsigned int)*ptext);
                    w = fDescendantFont->CIDWidth(cid);
                } else
                    w = -fDescendantFont->FontDef()->DW2()[1];
                ret.numwords++;
            } else if (btype[i] == PDF_BYTE_TYPE_LEAD) {
                if (fWritingMode == PDF_WMODE_HORIZONTAL) {
                    unsigned int code = (unsigned int)*ptext * 256 + 
                        *(ptext + 1);
                    cid = fCMap->GetCID(code);
                    w = fDescendantFont->CIDWidth(cid);
                } else 
                    w = -fDescendantFont->FontDef()->DW2()[1];
                ret.numwords++;
            }  
            ret.width += w;
            if (widths != NULL)
                widths[i] = w;
        }
        PDF_DEBUG_PRINT(("++ [%x] pdf_byte_type delete[].\n", (int)btype));
        delete[] btype;
    } catch (...) {
        PDF_DEBUG_PRINT(("++ [%x] pdf_byte_type delete[].\n", (int)btype));
        delete[] btype;
    }

    return ret;
}

unsigned int
PdfType0Font::MeasureText(const char* text, double width, double fontsize, 
        double charspace, double wordspace, double* realwidth)
{
    if (!fValid) {
        pdf_text_width ret;

        ret.width = 0;
        ret.numchars = 0;
        ret.numwords = 0;

        return ret.width;
    }
    return 0;
}

int
PdfType0Font::Ascent()
{
    return (GetValid() ? fDescendantFont->Ascent() : 0);
}

int
PdfType0Font::Descent()
{
    return (GetValid() ? fDescendantFont->Descent() : 0);
}

/*----------------------------------------------------------------------------*/
/*----- PdfCIDFontDef class --------------------------------------------------*/

PdfCIDType2FontDef::PdfCIDType2FontDef()
    : PdfAutoPtrObject()
{
    PDF_DEBUG_PRINT(("++ [%x] PdfCIDType2FontDef new.\n", (int)this));
    fBaseFont = NULL;
    fAscent = 0;
    fDescent = 0;
    fCapHeight = 0;
    fFlags = 0;
    fFontBBox = PdfBox(0, 0, 0, 0);
    fItalicAngle = 0;
    fStemV = 0;
    fAvgWidth = 0;
    fLeading = 0;
    fMaxWidth = 0;
    fMissingWidth = 0;
    fStemH = 0;
    fDW = 1000;
    fDW2[0] = 880;
    fDW2[1] = -1000;
    fWidths = NULL;
}

void
PdfCIDType2FontDef::SetParam(char** dst, const char* src)
{
    PDF_DEBUG_PRINT(("++ [%x] dst delete[].\n", (int)*dst));
    delete[] *dst;
    *dst = NULL;
    
    if (src == NULL) 
        return;
    
    size_t len = strlen(src) + 1;
    *dst = new char[len];
    PDF_DEBUG_PRINT(("++ [%x] dst new[]. \n", (int)*dst));
    strncpy(*dst, src, len);
}

PdfCIDType2FontDef::~PdfCIDType2FontDef()
{
    if (fWidths != NULL) {
        for (int i = 0; i < fWidths->CountItems(); i++) {
            pdf_cid_width* w = (pdf_cid_width*)fWidths->ItemAt(i);
            PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width->widths delete. \n", 
                        (int)w->widths));
            delete[] w->widths;
            PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width delete. \n", (int)w)); 
            delete w;
        }
    }
    PDF_DEBUG_PRINT(("++ [%x] fWidths delete. \n", (int)fWidths));
    delete fWidths;
    PDF_DEBUG_PRINT(("++ [%x] fBaseFont delete[]. \n", (int)fBaseFont));
    delete[] fBaseFont;
    PDF_DEBUG_PRINT(("++ [%x] PdfCIDType2FontDef delete. \n", (int)this));
}

unsigned int
PdfCIDType2FontDef::CIDWidth(pdf_cid cid)
{
    if (fWidths == NULL)
        return fDW;

    for (int i = 0; i < fWidths->CountItems(); i++) {
        pdf_cid_width* w = (pdf_cid_width*)fWidths->ItemAt(i);

        if (w->type == PDF_CID_W_TYPE_FROM_TO) {
            if (w->from_cid <= cid && w->to_cid >= cid) 
                return w->widths[0];
        } else {
            if (w->from_cid <= cid && w->to_cid >= cid)
                return w->widths[cid - w->from_cid];
        }
    }

    /* Not found in pdf_cid_width array. */
    return fDW;
}

int
PdfCIDType2FontDef::NumWidths()
{
    return (fWidths == NULL) ? 0 : fWidths->CountItems();
}

void
PdfCIDType2FontDef::AddWidths1(pdf_cid fromcid, pdf_cid tocid,
        unsigned int width)
{
    if (fWidths == NULL) {
        fWidths = new PdfList();
        PDF_DEBUG_PRINT(("++ [%x] PdfCIDType2FontDef::fWidths new. \n", 
                    (int)fWidths));
    }
    pdf_cid_width* w = new pdf_cid_width;
    PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width new. \n", (int)w));
    
    try {
        w->type = PDF_CID_W_TYPE_FROM_TO;
        w->from_cid = fromcid;
        w->to_cid = tocid;
        w->widths = new unsigned int[1];
        PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width->widths new. \n", 
                    (int)w->widths));
        w->widths[0] = width;
    } catch (...) {
        PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width delete. \n", (int)w));
        delete[] w;
        throw;
    }
        
    try {
        fWidths->AddItem(w);
    } catch (...) {
        PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width delete. \n", (int)w->widths));
        delete w->widths;
        PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width delete. \n", (int)w));
        delete[] w;
        throw;
    }
}

void
PdfCIDType2FontDef::AddWidths2(pdf_cid fromcid, pdf_cid tocid,
        const unsigned int* widths)
{
    if (fWidths == NULL) {
        fWidths = new PdfList();
        PDF_DEBUG_PRINT(("++ [%x] PdfCIDType2FontDef::fWidths new. \n", 
                    (int)fWidths));
    }
    pdf_cid_width* w = new pdf_cid_width;
    PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width new. \n", (int)w));

    try {
        w->type = PDF_CID_W_TYPE_FROM_ARRAY;
        w->from_cid = fromcid;
        w->to_cid = tocid;
        w->widths = new unsigned int[tocid - fromcid + 1];
        PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width->widths new. \n", 
                    (int)w->widths));
        for (int i = 0; i < tocid - fromcid + 1; i++) {
            w->widths[i] = *widths;
            widths++;
        }
    } catch (...) {
        PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width delete. \n", (int)w));
        delete[] w;
        throw;
    }

    try {
        fWidths->AddItem(w);
    } catch (...) {
        PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width->widths delete. \n", 
                    (int)w->widths));
        delete w->widths;
        PDF_DEBUG_PRINT(("++ [%x] pdf_cid_width delete. \n", (int)w));
        delete[] w;
        throw;
    }
}

/*----------------------------------------------------------------------------*/
/*----- PdfCIDFont class -----------------------------------------------------*/

PdfCIDType2Font::PdfCIDType2Font(PdfXref* xref)
            : PdfFontBase(xref)
{
    fValid = false;
    fFontDef = NULL;
}

PdfCIDType2Font::~PdfCIDType2Font()
{
}

void
PdfCIDType2Font::SetAttributes(PdfCIDFontDef* fontdef)
{
    if (fontdef == NULL) 
        throw PdfException(PDF_INVALID_PARAMETER, 
            "error: PdfCIDType2Font::SetAttributes().");

    fFontDef = fontdef;
    
    AddElement("Type", new PdfName("Font"));
    AddElement("Subtype", new PdfName("CIDFontType0"));
    AddElement("BaseFont", new PdfName(fontdef->BaseFont()));
    AddElement("DW", new PdfNumber(fontdef->DW()));

    PdfArray* dw2 = new PdfArray(GetXref());
    AddElement("DW2", dw2);
    dw2->Add(new PdfNumber(fontdef->DW2()[0]));
    dw2->Add(new PdfNumber(fontdef->DW2()[1]));
    
    PdfArray* w_array = new PdfArray(GetXref());
    AddElement("W", w_array);
    
    /* Create W array. 
     * "W"attribute constst of forrowing data format.
     * type1. <CID of first char> <CID if last char> <char width>
     * type2. <CID od first char> "["<array of char widths>"]"
     */
    for (int i = 0; i< fontdef->NumWidths(); i++) {
        pdf_cid_width* w = fontdef->GetWidths(i);

        if (w->type == PDF_CID_W_TYPE_FROM_ARRAY) {
            /* type2 */
            w_array->Add(new PdfNumber((int)w->from_cid));
            PdfArray* sub_array = new PdfArray(GetXref());
            w_array->Add(sub_array);
            for (int j = 0; j <= w->to_cid - w->from_cid; j++) 
                sub_array->Add(new PdfNumber((int)w->widths[j]));
        } else {
            /* type1 */
            w_array->Add(new PdfNumber((int)w->from_cid));
            w_array->Add(new PdfNumber((int)w->to_cid));
            w_array->Add(new PdfNumber((int)w->widths[0]));
        }
    }

    CreateDescriptor();
    fValid = true;
}

void
PdfCIDType2Font::CreateDescriptor()
{
    /* create font descriptor as a indirect object. */
    fDescriptor = new PdfDictionary(GetXref());
    GetXref()->AddObject(fDescriptor);
    AddElement("FontDescriptor", fDescriptor);

    fDescriptor->AddElement("Type", new PdfName("FontDescriptor"));
    fDescriptor->AddElement("FontName", new PdfName(fFontDef->BaseFont()));
    fDescriptor->AddElement("Ascent", new PdfNumber(fFontDef->Ascent()));
    fDescriptor->AddElement("Descent", new PdfNumber(fFontDef->Descent()));
    fDescriptor->AddElement("CapHeight", new PdfNumber(fFontDef->CapHeight()));
    fDescriptor->AddElement("MissingWidth", 
            new PdfNumber(fFontDef->MissingWidth()));
    fDescriptor->AddElement("Flags", new PdfNumber(fFontDef->Flags()));

    PdfArray *bbox = new PdfArray(GetXref());
    fDescriptor->AddElement("FontBBox", bbox);
    bbox->Add(fFontDef->FontBBox());
            
    fDescriptor->AddElement("ItalicAngle", 
            new PdfNumber(fFontDef->ItalicAngle()));
    fDescriptor->AddElement("StemV", new PdfNumber(fFontDef->StemV()));
}   

int
PdfCIDType2Font::Ascent()
{
    return (GetValid() ? fFontDef->Ascent() : 0);
}

int
PdfCIDType2Font::Descent()
{
    return (GetValid() ? fFontDef->Descent() : 0);
}

