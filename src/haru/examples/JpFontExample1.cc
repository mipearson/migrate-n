/*
 * << H a r u -- Free PDF Library >> -- JpFontExample1.cc
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
 *   Japanese Font example for SJIS encoding.
 */

#include "libharu.h"
#include "libharu_jpfonts.h"
#include "errno.h"

int main ()
{
    char samp_text[2048];
    const char* page_title = "JPFont Example";
    PdfDoc* doc = new PdfDoc();

    try {
        /* Start creating PDF document. */
        doc->NewDoc();
        doc->AddType1Font(new PdfHelveticaFontDef());

        /* Create CMap objects. */

        /* CMap for sjis encording (fix width font). */
        PdfCMap* map_sjis1 = new PdfCMap_90ms_RKSJ_H();

        /* Once the CMap object register to PdfDoc object by AddType0Font
         * method, the CMap object will be disposed automatically. But 
         * if AddType0Font failed, the CMap object may never be disposed.
         * So It is recommended to call RegisterCMap method to register
         * CMap object to PdfDoc object. */
        doc->RegisterObject(map_sjis1);

        /* Gothic */
        doc->AddType0Font(new PdfGothicFontDef(), "gothic-sjis", map_sjis1);
        
        doc->AddType0Font(new PdfGothicBoldFontDef(), 
                "gothic-bold-sjis", map_sjis1);
        doc->AddType0Font(new PdfGothicItalicFontDef(), 
                "gothic-italic-sjis", map_sjis1);
        doc->AddType0Font(new PdfGothicBoldItalicFontDef(), 
                "gothic-bolditalic", map_sjis1);

        /* Mincyo */
        doc->AddType0Font(new PdfMincyoFontDef(), "mincyo-sjis", map_sjis1);
        doc->AddType0Font(new PdfMincyoBoldFontDef(), 
                "mincyo-bold-sjis", map_sjis1);
        doc->AddType0Font(new PdfMincyoItalicFontDef(), 
                "mincyo-italic-sjis", map_sjis1);
        doc->AddType0Font(new PdfMincyoBoldItalicFontDef(), 
                "mincyo-bolditalic-sjis", map_sjis1);
    
        /* CMap for sjis encoding (proportional font) */
        PdfCMap* map_sjis2 = new PdfCMap_90msp_RKSJ_H();
        doc->RegisterObject(map_sjis2);

        /* Proportional Gothic */
        doc->AddType0Font(new PdfPGothicFontDef(), "pgothic-sjis", map_sjis2);
        doc->AddType0Font(new PdfPGothicBoldFontDef(), 
                "pgothic-bold-sjis", map_sjis2);
        doc->AddType0Font(new PdfPGothicItalicFontDef(), 
                "pgothic-italic-sjis", map_sjis2);
        doc->AddType0Font(new PdfPGothicBoldItalicFontDef(), 
                "pgothic-bolditalic", map_sjis2);

        /* Proportional Mincyo */
        doc->AddType0Font(new PdfPMincyoFontDef(), "pmincyo-sjis", map_sjis2);
        doc->AddType0Font(new PdfPMincyoBoldFontDef(), 
                "pmincyo-bold-sjis", map_sjis2);
        doc->AddType0Font(new PdfPMincyoItalicFontDef(), 
                "pmincyo-italic-sjis", map_sjis2);
        doc->AddType0Font(new PdfPMincyoBoldItalicFontDef(), 
                "pmincyo-bolditalic-sjis", map_sjis2);

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
        canvas->ShowText("<Japanese fonts samples>");
        canvas->EndText();
       
        /* Load sample text. */
        FILE* f = fopen("mbtext/sjis.txt", "rb");
        if (f == NULL) {
            perror("Cannot open [mbtext/sjis.txt].");
            return errno;
        }
        if (fgets(samp_text, 2048, f) == NULL) {
            perror("Cannot read from [mbtext/sjis.txt].");
            return errno;
        }

        /* Print the text in order using all fonts. */
        canvas->BeginText();
        canvas->MoveTextPos(60, canvas->Height() - 105);
        for (unsigned int i = 1; i < doc->FontMgr()->CountFonts(); i++) {
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
        doc->WriteToFile("JpFontExample1.pdf");
        doc->FreeDoc();
    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
    }
    delete doc;
    return 0;
}

