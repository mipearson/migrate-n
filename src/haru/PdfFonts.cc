/*
 * << H a r u --free pdf library >> -- PdfFonts.cpp
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
#include <new>
#include "libharu.h"

/*----------------------------------------------------------------------------*/
/*----- PdfFontBase class ----------------------------------------------------*/

PdfFontBase::PdfFontBase(PdfXref *xref)
        : PdfStream(xref)
{
    fName = NULL;
    fWritingMode = PDF_WMODE_HORIZONTAL;
}

PdfFontBase::~PdfFontBase()
{
    delete[] fName;
    PDF_DEBUG_PRINT(("++ [%x] fName delete.\n", (int)fName));
}

pdf_text_width
PdfFontBase::TextWidth(const char* text)
{
    return  TextWidths(text, NULL);
}

pdf_text_width
PdfFontBase::TextWidths(const char* text, unsigned int* widths)
{
    pdf_text_width ret;

    ret.width = 0;
    ret.numchars = 0;
    ret.numwords = 0;

    return ret;
}

unsigned int
PdfFontBase::MeasureText(const char* text, double width,
    double fontsize, double charspace, double wordspace, double* realwdth)
{
    return 0;
}

/*----------------------------------------------------------------------------*/
/*----- PdfFontMgr class -----------------------------------------------------*/

PdfFontMgr::PdfFontMgr(PdfXref *xref)
{
    assert(xref);

    fXref = xref;
    fList = NULL;
}

PdfFontMgr::~PdfFontMgr()
{
    delete fList;
}

int
PdfFontMgr::RegisterFont(PdfFont* font)
{
    assert(font);

    if (fList == NULL)
        fList = new PdfList();

    if (font->Name() == NULL) 
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: PdfFontMgr::RegisterFont: "
                "invalid font --name is null.");
    
    if (GetFont(font->Name()) != NULL) 
        throw PdfException(PDF_RUNTIME_ERROR, 
                "ERROR: PdfFontMgr::RegisterFont: "
                "duplicate font registration[%s]", font->Name());

    if (font->GetObjectType() != PDF_OBJ_TYPE_INDIRECT) 
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: PdfFontMgr::RegisterFont: "
                "font must be indirect object[%s]", font->Name());
    
    fList->AddItem(font);
    return fList->CountItems() - 1;
}

PdfFont*
PdfFontMgr::GetFont(const char* name)
{
    if (fList == NULL || name == NULL)
        return NULL;
    
    for (int i = 0; i < fList->CountItems(); i++) {
        PdfFont* f = GetFont(i);
        const char* fname = f->Name();
        if (strcmp(name, fname) == 0) {
            PDF_DEBUG_PRINT(("PdfFontMgr found font %s[%d]\n", fname, i));
            return f;
        }
    }
    return NULL;
}

/*----------------------------------------------------------------------------*/
/*----- PdfType1Font class ---------------------------------------------------*/

PdfType1Font::PdfType1Font(PdfXref *xref)
       : PdfFontBase(xref)
{
    fValid = false;
    fWidths = NULL;
    fFontDef = NULL;
    fEncoding = NULL;
    fDescriptor = NULL;
    fFirstChar = 0;
    fLastChar = 0;
    fMissingWidth = 0;
	fFlags = 0;
}

void
PdfType1Font::SetAttributes(const char* name, PdfType1FontDef* fontdef, 
        PdfEncodingDef* encoding)
{
    assert(fValid == false);
    assert(fontdef != false);
    
    fFontDef = fontdef;
    fEncoding = encoding;
    PdfArray* width_array = NULL;
	fFlags = fFontDef->Flags();
    
    fName = new char[strlen(name) + 1];
    PDF_DEBUG_PRINT(("++ [%x] fName new.\n", (int)fName));
        
    strcpy(fName, name);
    AddElement("Type", new PdfName("Font"));
    
    if (!IsBase14Font()) {
        width_array = new PdfArray(GetXref());
        AddElement("Widths", width_array);
    }

    /* FONT_SPECIFIC encoding */
    if (fEncoding == NULL || fEncoding->BaseEncoding() == PDF_FONT_SPECIFIC) {
        fFirstChar = fFontDef->FirstChar();
        fLastChar = fFontDef->LastChar();

        PDF_DEBUG_PRINT(("PdfType1Font --FirstChar=%d, LastChar=%d \n",
                    fFirstChar, fLastChar));

		if (fFontDef->DefaultEncoding() == PDF_FONT_SPECIFIC) {
			fFlags ^= PDF_FONT_STD_CHARSET;
			fFlags |= PDF_FONT_SYMBOLIC;
		} 
    } else {
        fFirstChar = fEncoding->FirstChar();
        fLastChar = fEncoding->LastChar();
    }

    fWidths = new int[fLastChar - fFirstChar + 1];
    PDF_DEBUG_PRINT(("++ [%x] PdfType1Font::fWidths new.\n", (int)fWidths));

    if (fEncoding == NULL || fEncoding->BaseEncoding() == PDF_FONT_SPECIFIC) {
        for (int i = fFirstChar, j = 0; i <= fLastChar; i++, j++) {

            fWidths[j] = fFontDef->Widths(i);
            if (!IsBase14Font())
                width_array->Add(new PdfNumber(fWidths[j]));
        }
    } else {
        for (unsigned int i = fEncoding->FirstChar(), j = 0;
            i <= fEncoding->LastChar(); i++, j++) {
            const char* char_name = fEncoding->GetCharName(i);

            fWidths[j] = fFontDef->Widths(char_name);
            PDF_DEBUG_PRINT(("PdfType1Font --char_name=%s, code=%d, "
                        "width=%d\n", char_name, i, fWidths[j]));
            if (!IsBase14Font())
                width_array->Add(new PdfNumber(fWidths[j]));
        }
    }
    fMissingWidth = fFontDef->MissingWidth();

    AddElement("Subtype", new PdfName("Type1"));
    AddElement("BaseFont", new PdfName(fFontDef->BaseFont()));
    AddElement("FirstChar", new PdfNumber(fFirstChar));
    AddElement("LastChar", new PdfNumber(fLastChar));

    if (fEncoding != NULL) {
        PdfObject* obj = fEncoding->GetEncoding(this->GetXref());
        if (obj != NULL)
            AddElement("Encoding", obj);
    }

    if (!IsBase14Font())
        CreateDescriptor();
    else
        fDescriptor = NULL;

    fValid = true;
}

PdfType1Font::~PdfType1Font()
{
    PDF_DEBUG_PRINT(("++ [%x] PdfType1Font::fWidths delete.\n", (int)fWidths));
    delete[] fWidths;
}

pdf_text_width
PdfType1Font::TextWidth(const char* text)
{
    return TextWidths(text, NULL);
}

pdf_text_width
PdfType1Font::TextWidths(const char* text, unsigned int* widths)
{
    int len = strlen(text);
    pdf_text_width ret;
    unsigned int* tmp_widths = widths;

    ret.width = 0;
    ret.numchars = 0;
    ret.numwords = 0;
    const unsigned char* p = (const unsigned char*)text;

    for (int i = 0; i < len; i++, p++) {
        ret.width += CharWidth(*p);

        if (tmp_widths != NULL) {
            *tmp_widths = CharWidth(*p);
            tmp_widths++;
        }
        ret.numchars++;
        if (*p == ' ')
            ret.numwords++;
    }

    return ret;
}

unsigned int
PdfType1Font::MeasureText(const char* text, double width,
        double fontsize, double charspace, double wordspace, double* realwidth)
{
    double w = 0;
    int ln = strlen(text);
    int tmp_ln = 0;
    const char* tmp_char = text;

    for (int i = 0; i < ln; i++) {
        w += CharWidth(*tmp_char) * fontsize / 1000;
        if (w > width)
            if (tmp_ln > 0)
                return  (unsigned int)tmp_ln;
            else {
                if (realwidth != NULL)
                    *realwidth = w;
                return  (unsigned int)i;
            }
        if (realwidth != NULL)
            *realwidth = w;
        if (i > 0)
            w += charspace;
        if (*tmp_char == ' ') {
            w += wordspace;
            tmp_ln = i + 1;
        }
        tmp_char++;
    }

    return (int)ln;
}

int
PdfType1Font::Ascent()
{
    return (GetValid() ? fFontDef->Ascent() : 0);
}

int
PdfType1Font::Descent()
{
    return (GetValid() ? fFontDef->Descent() : 0);
}

bool
PdfType1Font::CreateDescriptor()
{
    /* if the font data is shared by other font objects, also FontDescripter
     * data will be shared. 
     */
    if (fFontDef->Descriptor() != NULL) {
        fDescriptor = fFontDef->Descriptor();
        AddElement("FontDescriptor", fDescriptor);
        return false;
    }
    
    /* create font descriptor as a indirect object. */
    fDescriptor = new PdfDictionary(GetXref());
    GetXref()->AddObject(fDescriptor);
    AddElement("FontDescriptor", fDescriptor);

    fDescriptor->AddElement("Type", new PdfName("FontDescriptor"));
    fDescriptor->AddElement("Ascent", new PdfNumber(fFontDef->Ascent()));
    fDescriptor->AddElement("CapHeight", new PdfNumber(fFontDef->CapHeight()));
    fDescriptor->AddElement("Descent", new PdfNumber(fFontDef->Descent()));

    PdfArray *bbox = new PdfArray(GetXref());
    fDescriptor->AddElement("FontBBox", bbox);
    bbox->Add(fFontDef->FontBBox());
    fDescriptor->AddElement("FontName", new PdfName(fFontDef->FontName()));

    fDescriptor->AddElement("ItalicAngle", 
            new PdfNumber(fFontDef->ItalicAngle()));
    fDescriptor->AddElement("StemV", new PdfNumber(fFontDef->StemV()));
    fDescriptor->AddElement("Flags", new PdfNumber(fFlags));

    if (strlen(fFontDef->CharSet()) > 0)
        fDescriptor->AddElement("CharSet", new PdfName(fFontDef->CharSet()));
    if (fFontDef->XHeight() != 0)
        fDescriptor->AddElement("XHeight", new PdfNumber(fFontDef->XHeight()));

    if (fFontDef->FontData() != NULL) {
        /* create font file stream as a indirect object. */
        PdfStream *font_file = new PdfStream(GetXref());
        GetXref()->AddObject(font_file);
        font_file->AddFilter(PDF_FILTER_DEFLATE);
        fDescriptor->AddElement("FontFile", font_file);

        fFontDef->FontData()->WriteToStream(font_file->GetStream());
        font_file->AddElement("Length1", new PdfNumber(fFontDef->Length1()));
        font_file->AddElement("Length2", new PdfNumber(fFontDef->Length2()));
        font_file->AddElement("Length3", new PdfNumber(fFontDef->Length3()));
    }

    fFontDef->SetDescriptor(fDescriptor);

    return true;
}

/*----------------------------------------------------------------------------*/

