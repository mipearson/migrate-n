/*
 * << H a r u --free pdf library >> -- Hello.cc
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

int main(int argc, char** argv)
{
    PdfDoc *doc = new PdfDoc();

    try {
        /* Create a new PDF document. */
        doc->NewDoc();

        /* Add Times-Roman font to the document. */
        doc->AddType1Font(new PdfTimesRomanFontDef());
    
        /* Add a page to the document. */
        PdfPage *page = doc->AddPage();

        /* Get the canvas object of the page. */
        PdfContents *canvas = page->Canvas();

        /*Set current font to "Times-Roman" and set the size to 80. */
        canvas->SetFontAndSize("Times-Roman", 80);

        /* Start to print text. */
        canvas->BeginText();

        /* Move the position of the text to 150,400. */
        canvas->MoveTextPos(150, 400);

        /* Print "Hello". */
        canvas->ShowText("Hello");

        /* Finish to print text. */
        canvas->EndText();

        /* Save the document as "Hello.pdf" */
        doc->WriteToFile("Hello.pdf");

    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
    }

    /* Clean up */
    delete doc;

    return 0;
}

