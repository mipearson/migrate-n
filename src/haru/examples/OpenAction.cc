/*
 * << H a r u --free pdf library >> -- OpenAction.cpp
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

/*------ Simple example for OpenAction --------------------------------------*/

int main()
{
    const char* text1 = "When this file is opened, the size of page is ";
    const char* text2 = "adjusted to the size of viewer window.";
    
    /* Create a PdfDoc object and start to make new PDF file. */
    PdfDoc* doc = new PdfDoc();
    try {
        doc->NewDoc();

        /* Create font object and register it to the PdfDoc object.*/
        doc->AddType1Font(new PdfHelveticaFontDef());

        /* Add new page. */
        PdfPage* page = doc->RootPages()->AddPage(); 
        page->SetSize(600, 800);  

        /* Make Destination object for OpenAction. */
        PdfDestination* dest = new PdfDestination(page);
        dest->SetFit();

        /* Set the Destination object to the catalog of the document. */
        doc->Catalog()->SetOpenAction(dest);

        /* Get contents of the page. */
        PdfContents* canvas = page->Canvas();
        
        /* Print a message. */
        canvas->BeginText();
        canvas->SetFontAndSize("Helvetica", 20);
        double w = canvas->TextWidth(text1);
        canvas->MoveTextPos((canvas->Width() - w) / 2, 400);
        canvas->SetTextLeading(25);
        canvas->ShowText(text1);
        canvas->ShowTextNextLine(text2);
        canvas->EndText();
        
        /* Save the document to the file. */
        doc->WriteToFile("OpenAction.pdf"); 
    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());                 
        delete doc;
        return 1;
    }
    delete doc;
    return 0;
}

