/*
 * << H a r u --free pdf library >> -- ViewStyle.cc
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

PdfPage* AddPage(PdfDoc *doc, int page_num);

int main(int argc, char** argv)
{
    char buf[1024];
    pdf_page_mode page_mode = PDF_USE_NONE;
    pdf_page_layout layout = PDF_SINGLE_PAGE;
    
    /* select page-mode */  
    printf("Select PageMode \n 0: Use None\n 1: Use Outlines\n"
            " 2: Use Thumbs\n 3: Full Screen\n x: Exit demo\n : ");
    for (;;) {
        if (fgets(buf, 1024, stdin) != NULL) {
            bool selected = true;
            
            switch (buf[0]) {
                case '0': page_mode = PDF_USE_NONE;
                          break;
                case '1': page_mode = PDF_USE_OUTLINES;
                          break;
                case '2': page_mode = PDF_USE_THUMBS;
                          break;
                case '3': page_mode = PDF_FULL_SCREEN;
                          break;
                case 'x': return 0;
                default:  printf("Select PageMode \n 0: Use None\n"
                                  " 1: Use Outlines\n 2: Use Thumbs\n"
                                  " 3: Full Screen\n x: Exit demo\n : ");
                          selected = false;
            }
            if (selected)
                break;
        }
    }
    
    /* select PageLayout */
    printf("Select PageLayout \n 0: Single Page\n 1: One Column\n"
            " 2: Two Column Left\n 3: Two Column Right\n x: Exit demo\n : ");
    for (;;) {
        if (fgets(buf, 1024, stdin) != NULL) {
            bool selected = true;
            
            switch (buf[0]) {
                case '0': layout = PDF_SINGLE_PAGE;
                          break;
                case '1': layout = PDF_ONE_COLUMN;
                          break;
                case '2': layout = PDF_TWO_COLUMN_LEFT;
                          break;
                case '3': layout = PDF_TWO_COLUMN_RIGHT;
                          break;
                case 'x': return 0;
                default: printf("Select PageLayout \n 0: Single Page\n"
                                 " 1: One Column\n 2: Tow Column Left\n"
                                 " 3: Two Column Light\n x: Exit demo\n : ");
                         selected = false;
            }
            if (selected)
                break;
        }
    }
    
    PdfDoc *doc = new PdfDoc();

    try {
        /* Create a new PDF document. */
        doc->NewDoc();

        /* Set page-mode to selected value. */
        doc->Catalog()->SetPageMode(page_mode);

        /* Set page-layout to selected value. */
        doc->Catalog()->SetPageLayout(layout);
        
        /* Add Times-Roman font to the document. */
        doc->AddType1Font(new PdfTimesRomanFontDef());

        /* Add 4 pages to the document */ 
        PdfPage* page1 = AddPage(doc, 1);
        PdfPage* page2 = AddPage(doc, 2);
        PdfPage* page3 = AddPage(doc, 3);
        PdfPage* page4 = AddPage(doc, 4);
        
        /* Create outlines. */
        PdfOutlineRoot *outline_root = doc->Outlines();
        PdfOutlineItem *outline_item1 = new PdfOutlineItem(outline_root);
        outline_item1->SetTitle("Page1");
        PdfOutlineItem *outline_item2 = new PdfOutlineItem(outline_root);
        outline_item2->SetTitle("Page2");
        PdfOutlineItem *outline_item3 = new PdfOutlineItem(outline_root);
        outline_item3->SetTitle("Page3");
        PdfOutlineItem *outline_item4 = new PdfOutlineItem(outline_root);
        outline_item4->SetTitle("Page4");
    
        /* Atatch the page object to outline-item. */
        PdfDestination *dst1 = new PdfDestination(page1);
        outline_item1->SetDestination(dst1);
        PdfDestination *dst2 = new PdfDestination(page2);
        outline_item2->SetDestination(dst2);
        PdfDestination *dst3 = new PdfDestination(page3);
        outline_item3->SetDestination(dst3);
        PdfDestination *dst4 = new PdfDestination(page4);
        outline_item4->SetDestination(dst4);
        
        /* Save the document as "Hello.pdf" */
        doc->WriteToFile("ViewStyle.pdf");

    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
    }

    /* Clean up */
    delete doc;

    return 0;
}

PdfPage* AddPage(PdfDoc *doc, int page_num) {

    char buf[10];

    sprintf(buf, "Page %d", page_num);
    
    /* Add a page to the document. */
    PdfPage *page = doc->AddPage();

    /* Get the canvas object of the page. */
    PdfContents *canvas = page->Canvas();

    /*Set current font to "Times-Roman" and set the size to 80. */
    canvas->SetFontAndSize("Times-Roman", 80);

    /* Get textWidth and calcurate the position of the test. */
    double w = canvas->TextWidth(buf);
    double xpos = (canvas->Width() - w) / 2;
    double ypos = (canvas->Height() - canvas->FontSize()) / 2;
    
    /* Start to print text. */
    canvas->BeginText();

    /* Move the position of the text to 150,400. */
    canvas->MoveTextPos(xpos, ypos);

    /* Print "Hello". */
    canvas->ShowText(buf);

    /* Finish to print text. */
    canvas->EndText();

    return page;
}

