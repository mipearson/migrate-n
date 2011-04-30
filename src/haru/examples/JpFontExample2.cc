/*
 * << H a r u -- Free PDF Library >> -- JpFontExample2.cc
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
 * << Example for virtical writing of Japanese text. >>
 *
 */

#include "libharu.h"
#include "libharu_jpfonts.h"
#include "errno.h"

int main ()
{
    char samp_text[2048];

    PdfDoc* doc = new PdfDoc();

    try {
        /* Start creating PDF document. */
        doc->NewDoc();

        /* Create CMap for EUC encoding(vertical writing). */
        PdfCMap* map2 = new PdfCMap_EUC_V();

        /* Create font definition objects and register to the PdfDoc object. */
        PdfCIDFontDef* def = new PdfGothicFontDef();
        doc->AddType0Font(def, "gothic-euc-v", map2);

        /* Add a new page object. */
        PdfPage* page = doc->AddPage();
        page->SetSize(100, 430);
        PdfContents* canvas = page->Canvas();

        /* Load sample text. */
        FILE* f = fopen("mbtext/euc.txt", "rb");
        if (f == NULL) {
                perror("Cannot open [mbtext/euc.txt].");
                    return errno;
        }
        if (fgets(samp_text, 2048, f) == NULL) {
                perror("Cannot read from [mbtext/euc.txt].");
                    return errno;
        }

        /* Print the text. */
        canvas->BeginText();
        canvas->SetFontAndSize("gothic-euc-v", 20);
        canvas->MoveTextPos(35, 410);
        canvas->ShowText(samp_text);
        canvas->MoveTextPos(30, 0);
        canvas->ShowText(samp_text);
        canvas->EndText();

        /* Save the document to a file */
        doc->WriteToFile("JpFontExample2.pdf");
        doc->FreeDoc();
    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
    }
    delete doc;
    return 0;
}

