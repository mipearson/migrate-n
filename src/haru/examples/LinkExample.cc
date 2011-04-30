/*
 * << H a r u --free pdf library >> -- LinkExample.cpp
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

/*------ Simple example for Link-Annotation ----------------------------------*/

void PrintPage(PdfContents* canvas, int page);

void PrintPage(PdfContents* canvas, int page)
{
    char buf[50];
    
    canvas->BeginText();
    canvas->SetFontAndSize("Helvetica", 20);
    canvas->MoveTextPos(50, 150);
#ifdef __WIN32__
    _snprintf(buf, 50, "Page:%d", page);
#else
    snprintf(buf, 50, "Page:%d", page);
#endif
    canvas->ShowText(buf);
    canvas->EndText();
}

int main()
{
    /* Create a PdfDoc object and start to make new PDF file. */
    PdfDoc* doc = new PdfDoc();
    try {
        doc->NewDoc();

        /* Create font object and register it to the PdfDoc object.*/
        doc->AddType1Font(new PdfHelveticaFontDef());

        /* Add new page(index page) and set the size of the page. */
        PdfPage* index_page = doc->AddPage();
        index_page->SetSize(200, 200);

        /* Add new page (page 2). */
        PdfPage* page2 = doc->RootPages()->AddPage(); 
        page2->SetSize(200, 200);
        PrintPage(page2->Canvas(), 2);  

        /* Add new page (page 3). */
        PdfPage* page3 = doc->RootPages()->AddPage(); 
        page3->SetSize(200, 200);  
        PrintPage(page3->Canvas(), 3);  

        /* Add new page (page 4). */
        PdfPage* page4 = doc->RootPages()->AddPage(); 
        page4->SetSize(200, 200);  
        PrintPage(page4->Canvas(), 4);

        /* make Destination objects for Link-Annotation. */
        PdfDestination* dest1 = new PdfDestination(page2, true);
        PdfDestination* dest2 = new PdfDestination(page3);
        PdfDestination* dest3 = new PdfDestination(page4);

        /*
         * Create Link-Annotation object with these Destination objects
         * on index page 
         */
        PdfContents* canvas = index_page->Canvas();
        canvas->BeginText();
        canvas->SetFontAndSize("Helvetica", 15);
        canvas->MoveTextPos(50, 150);
        canvas->SetTextLeading(25);
    
        int text_width = (int)canvas->TextWidth("Jump to PageX");   
        pdf_rect rect = PdfRect(canvas->TextPos().x, 
                                canvas->TextPos().y,
                                canvas->TextPos().x + text_width,
                                canvas->TextPos().y + 20);
        index_page->AddLink(rect, dest1);
        canvas->ShowText("Jump to Page2");
        canvas->MoveToNextLine();
        
        rect.bottom = canvas->TextPos().y;
        rect.top = rect.bottom + 20;
        index_page->AddLink(rect, dest2);
        
        canvas->ShowText("Jump to Page3");
        canvas->MoveToNextLine();
        
        rect.bottom = canvas->TextPos().y;
        rect.top = rect.bottom + 20;
        index_page->AddLink(rect, dest3);
        
        canvas->ShowText("Jump to Page4");
        canvas->MoveToNextLine();

        canvas->EndText();

        /* Save the document to the file. */
        doc->WriteToFile("LinkExample.pdf"); 
    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());                 
        delete doc;
        return 1;
    }
    delete doc;
    return 0;
}

