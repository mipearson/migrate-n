/*
 * << H a r u --free pdf library >> -- Permisson.cpp
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

/*------ Simple example for restricting operation. ---------------------------*/

const static char* text = "User cannot print and copy this document.";
const static char* owner_passwd = "owner";
const static char* user_passwd = "";

int main()
{
    /* Create a PdfDoc object and start to make new PDF file. */
    PdfDoc* doc = new PdfDoc();
    try {
        doc->NewDoc();

        /* Create font object and register it to the PdfDoc object.*/
        doc->AddType1Font(new PdfHelveticaFontDef());

        /* Add new page. */
        PdfPage* page = doc->RootPages()->AddPage(); 
        page->SetSize(600, 800);  

        /* Get contents of the page. */
        PdfContents* canvas = page->Canvas();

        /* Compless the contents stream with deflator */
        canvas->AddFilter(PDF_FILTER_DEFLATE);
        
        /* Print a message. */
        canvas->BeginText();
        canvas->SetFontAndSize("Helvetica", 20);
        double w = canvas->TextWidth(text);
        canvas->MoveTextPos((canvas->Width() - w) / 2, 400);
        canvas->ShowText(text);
        canvas->EndText();

        /* Set owner an user password. */
        doc->SetPassword(owner_passwd, user_passwd);

        /* Set user access permission. Only viewing is allowed. */
        doc->SetPermission(PDF_ENABLE_READ);
        
        /* Save the document to the file. */
        doc->WriteToFile("Permission.pdf"); 
    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());                 
        delete doc;
        return 1;
    }
    delete doc;
    return 0;
}

