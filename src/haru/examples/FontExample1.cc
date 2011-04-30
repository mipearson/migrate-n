/*
 * << H a r u -- Free PDF Library >> -- FontExample1.cc
 *
 * Copyright (c) 1999-2002 Takeshi Kanno <takeshi_kanno@est.hi-ho.ne.jp>
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

int main ()
{
    const char* page_title = "Font Example1";
    PdfDoc* doc = new PdfDoc();

    try {
        /* Start creating PDF document. */
        doc->NewDoc();
        
        /* Create font definition objects and register to the PdfDoc object. */

        /* Helvetica */
        doc->AddType1Font(new PdfHelveticaFontDef());
        doc->AddType1Font(new PdfHelveticaBoldFontDef());
        doc->AddType1Font(new PdfHelveticaObliqueFontDef());
        doc->AddType1Font(new PdfHelveticaBoldObliqueFontDef());

        /* Courier */
        doc->AddType1Font(new PdfTimesRomanFontDef());
        doc->AddType1Font(new PdfTimesBoldFontDef());
        doc->AddType1Font(new PdfTimesItalicFontDef()); 
        doc->AddType1Font(new PdfTimesBoldItalicFontDef());
    
        /* Courier */
        doc->AddType1Font(new PdfCourierFontDef());
        doc->AddType1Font(new PdfCourierBoldFontDef());
        doc->AddType1Font(new PdfCourierObliqueFontDef()); 
        doc->AddType1Font(new PdfCourierBoldObliqueFontDef());

        /* Symbol */
        doc->AddType1Font(new PdfSymbolFontDef());

        /* ZapfDingbats */
        doc->AddType1Font(new PdfZapfDingbatsFontDef()); 

        /* Add a new page object. */
        PdfPage* page = doc->AddPage();
        PdfContents* canvas = page->Canvas();
        
        /* Print the lines of the page. */
        canvas->SetLineWidth(1);
        canvas->Rectangle(50, 50, canvas->Width() - 100, 
                canvas->Height() - 110);
        canvas->Stroke();

        /* Print the title of the page (with positioning center). */
        canvas->SetFontAndSize("Helvetica", 24);
        double w = canvas->TextWidth(page_title);
        canvas->BeginText();
        canvas->MoveTextPos((canvas->Width() - w) / 2, canvas->Height() - 50);
        canvas->ShowText(page_title);
        canvas->EndText();

        /* Output subtitle. */
        canvas->BeginText();
        canvas->MoveTextPos(60, canvas->Height() - 80);
        canvas->SetFontAndSize("Helvetica", 16);
        canvas->ShowText("<Standerd Type1 fonts samples>");
        canvas->EndText();
        
        /* Output texts in order using all fonts. */
        canvas->BeginText();
        canvas->MoveTextPos(60, canvas->Height() - 105);
        for (unsigned int i = 0; i < doc->FontMgr()->CountFonts(); i++) {
            const char* samp_text = "abcdefgABCDEFG12345!#$%&+-@?";
            
            PdfFont* f = doc->FontMgr()->GetFont(i);
        
            /* print a label of text */ 
            canvas->SetFontAndSize("Helvetica", 9);
            canvas->ShowText(f->Name());
            canvas->MoveTextPos(0, -18);
            
            /* print a sample text. */
            canvas->SetFontAndSize(f, 20);
            canvas->ShowText(samp_text);
            canvas->MoveTextPos(0, -20);
        }
        canvas->EndText();
        
        /* Save the document to a file */
        doc->WriteToFile("FontExample1.pdf");
        doc->FreeDoc();
    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
    }
    delete doc;
    return 0;
}
