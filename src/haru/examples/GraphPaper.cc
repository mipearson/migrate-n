/*
 * << H a r u --free pdf library >> -- GraphPaper.cpp
 *
 * Copyright (c) 1999-2002 Takeshi Kanno <takeshi_kanno@est.hi-ho.ne.jp>
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 * It is provided "as is" without express or implied warranty.
 *
 */

#include "libharu.h"

int main()
{

    /* Create a PdfDoc object and start making a PDF document. */
    PdfDoc* doc = new PdfDoc();
    try {
        doc->NewDoc();

        /* Make helvetica font object. */
        PdfType1FontDef* fd1 = new PdfHelveticaFontDef();

        /* Add font with PdfStandardEncoding and name it "helvetica" */
        doc->AddType1Font(fd1, "helvetica");
    
        for (int i = 0; i < 2; i++) {
        
        /* Add new page and set the size of the page */
        PdfPage* page = doc->AddPage();
        
        /* Draw graphs on the canvas. */
        PdfContents* canvas = page->Canvas();
        canvas->AddFilter(PDF_FILTER_DEFLATE);

        canvas->SetLineWidth(0.1);
        canvas->SetFontAndSize("helvetica", 8);

        /* Draw horizontal line */  
        int y = 60;
        canvas->SetDash(2, 0, 1);
        while (y <= page->Height() - 50) {
            canvas->MoveTo(50, y);
            canvas->LineTo(page->Width() - 50, y);
            canvas->Stroke();
            
            y += 20;
        }
        
        y = 50;
        canvas->SetDash(0, 0, 0);
        while (y <= page->Height() - 50) {
            char label[11];
            
            canvas->MoveTo(50, y);
            canvas->LineTo(page->Width() - 50, y);
            canvas->Stroke();

#ifdef __WIN32__
            _snprintf(label, 11, "%d", y);
#else
            snprintf(label, 11, "%d", y);
#endif
            canvas->BeginText();
            canvas->MoveTextPos(50 - canvas->TextWidth(label) - 2, y - 3);
            canvas->ShowText(label);
            canvas->EndText();
            
            y += 20;
        }

        /* Draw virtical line */
        int x = 60;
        canvas->SetDash(2, 0, 1);
        while (x <= page->Width() - 50) {
            canvas->MoveTo(x, 50);
            canvas->LineTo(x, canvas->Height() - 50);
            canvas->Stroke();
    
            x += 20;
        }

        x = 50;
        canvas->SetDash(0, 0, 0);
        while (x <= page->Width() - 50) {
            char label[11];
            
            canvas->MoveTo(x, 50);
            canvas->LineTo(x, canvas->Height() - 50);
            canvas->Stroke();

            if (x > 50) {
#ifdef __WIN32__
                _snprintf(label, 11, "%d", x);
#else
                snprintf(label, 11, "%d", x);
#endif
                canvas->BeginText();
                canvas->MoveTextPos(x - canvas->TextWidth(label) / 2, 40);
                canvas->ShowText(label);
                canvas->EndText();
            }
            x += 20;
        }

        }

        /* Save the document to the file. */
        doc->WriteToFile("GraphPaper.pdf");

    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
        delete doc;
        return 1;
    }
    delete doc;
    return 0;
}


