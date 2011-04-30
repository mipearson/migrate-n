/*
 * << H a r u --free pdf library >> -- OutlineExample.cpp
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

void PrintPage(PdfContents* canvas, int page);

void PrintPage(PdfContents* canvas, int page)
{
    char buf[50];
    
    canvas->BeginText();
    canvas->SetFontAndSize("Helvetica", 20);
    canvas->MoveTextPos(50, 250);
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
    PdfDoc *doc = new PdfDoc();

    try {
        doc->NewDoc();

        /* Create font object and register it to the PdfDoc object.*/
        doc->AddType1Font(new PdfHelveticaFontDef());

        /* Set page mode to use outlines. */
        doc->Catalog()->SetPageMode(PDF_USE_OUTLINES);

        /* Add 3 pages to the document. */
        PdfPage *page1 = doc->AddPage();
        page1->SetSize(200, 300);
        PrintPage(page1->Canvas(), 1);
    
        PdfPage *page2 = doc->AddPage();
        page2->SetSize(200, 300);
        PrintPage(page2->Canvas(), 2);
    
        PdfPage *page3 = doc->AddPage();
        page3->SetSize(200, 300);
        PrintPage(page3->Canvas(), 3);
    
        /* create outline root. */
        PdfOutlineRoot *outline_root = doc->Outlines();
        PdfOutlineItem *outline_item = new PdfOutlineItem(outline_root);
        outline_item->SetTitle("OutlineRoot");
        outline_item->SetOpened(true);
    
        /* create outline items. */
        PdfOutlineItem *outline_item1 = new PdfOutlineItem(outline_item);
        outline_item1->SetTitle("page1");

        PdfOutlineItem *outline_item2 = new PdfOutlineItem(outline_item);
        outline_item2->SetTitle("page2");

        PdfOutlineItem *outline_item3 = new PdfOutlineItem(outline_item);
        outline_item3->SetTitle("page3");
    
        /* create destination objects on each pages 
         * and link it to outline items. 
        */  
        PdfDestination *dst = new PdfDestination(page1);
        dst->SetXYZ(0, page1->Height(), 1);
        outline_item1->SetDestination(dst);
        doc->Catalog()->SetOpenAction(dst);

        dst = new PdfDestination(page2);
        dst->SetXYZ(0, page2->Height(), 1);
        outline_item2->SetDestination(dst);

        dst = new PdfDestination(page3);
        dst->SetXYZ(0, page3->Height(), 1);
        outline_item3->SetDestination(dst);

        /* save the document to a file. */
        doc->WriteToFile("OutlineExample.pdf");
    } catch (PDF_STD_EXCEPTION& error) {
        fprintf(stderr, "%s\n", error.what());
    }

    /* clean up. */
    delete doc;
    
    return 0;
}
