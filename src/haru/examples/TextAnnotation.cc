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
#include "libharu_ISO8859.h"

int main(int argc, char** argv)
{
    PdfDoc* doc = new PdfDoc();
    pdf_rect rect1 = {50, 350, 150, 400};
    pdf_rect rect2 = {210, 350, 350, 400};
    pdf_rect rect3 = {50, 250, 150, 300};
    pdf_rect rect4 = {210, 250, 350, 300};
    pdf_rect rect5 = {50, 150, 150, 200};
    pdf_rect rect6 = {210, 150, 350, 200};
    pdf_rect rect7 = {50, 50, 150, 100};
    pdf_rect rect8 = {210, 50, 350, 100};

    try {
        /* Create a new PDF document. */
        doc->NewDoc();

        /* Add Times-Roman font to the document. */
        doc->AddType1Font(new PdfTimesRomanFontDef());
    
        /* Add a page to the document. */
        PdfPage* page = doc->AddPage();
        page->SetSize(400, 500);
        
        /* Get the canvas object of the page. */
        PdfContents* canvas = page->Canvas();
        canvas->BeginText();
        canvas->SetFontAndSize("Times-Roman", 16);
        canvas->MoveTextPos(130, 450);
        canvas->ShowText("Annotation Example");
        canvas->EndText();

        PdfTextAnnot* annot = page->AddTextAnnot(rect1);
        annot->SetContents("Annotation with Comment Icon. \n"
                "This annotation set to be opened initially.");
        annot->SetIcon(PDF_ANNOT_ICON_COMMENT);
        annot->SetOpened(true);
        
        annot = page->AddTextAnnot(rect2);
        annot->SetContents("Annotation with Key Icon");
        annot->SetIcon(PDF_ANNOT_ICON_PARAGRAPH);
        
        annot = page->AddTextAnnot(rect3);
        annot->SetContents("Annotation with Note Icon");
        annot->SetIcon(PDF_ANNOT_ICON_NOTE);
        
        annot = page->AddTextAnnot(rect4);
        annot->SetContents("Annotation with Help Icon");
        annot->SetIcon(PDF_ANNOT_ICON_HELP);
        
        annot = page->AddTextAnnot(rect5);
        annot->SetContents("Annotation with NewParagraph Icon");
        annot->SetIcon(PDF_ANNOT_ICON_NEW_PARAGRAPH);
        
        annot = page->AddTextAnnot(rect6);
        annot->SetContents("Annotation with Paragraph Icon");
        annot->SetIcon(PDF_ANNOT_ICON_PARAGRAPH);
        
        annot = page->AddTextAnnot(rect7);
        annot->SetContents("Annotation with Insert Icon");
        annot->SetIcon(PDF_ANNOT_ICON_INSERT);
       
        PdfEncodingDef* encoding = new PdfEncoding_ISO8859_2();
        annot = page->AddTextAnnot(rect8);
        annot->SetContents("Annotation with ISO8859 text гдежзий", encoding);
      
        canvas->SetFontAndSize("Times-Roman", 11);
        
        canvas->BeginText();
        canvas->MoveTextPos(rect1.left + 35, rect1.top - 20);
        canvas->ShowText("Comment Icon.");
        canvas->EndText();
        
        canvas->BeginText();
        canvas->MoveTextPos(rect2.left + 35, rect2.top - 20);
        canvas->ShowText("Key Icon");
        canvas->EndText();
        
        canvas->BeginText();
        canvas->MoveTextPos(rect3.left + 35, rect3.top - 20);
        canvas->ShowText("Note Icon.");
        canvas->EndText();
       
        canvas->BeginText();
        canvas->MoveTextPos(rect4.left + 35, rect4.top - 20);
        canvas->ShowText("Help Icon");
        canvas->EndText();
       
        canvas->BeginText();
        canvas->MoveTextPos(rect5.left + 35, rect5.top - 20);
        canvas->ShowText("NewParagraph Icon");
        canvas->EndText();
       
        canvas->BeginText();
        canvas->MoveTextPos(rect6.left + 35, rect6.top - 20);
        canvas->ShowText("Paragraph Icon");
        canvas->EndText();
       
        canvas->BeginText();
        canvas->MoveTextPos(rect7.left + 35, rect7.top - 20);
        canvas->ShowText("Insert Icon");
        canvas->EndText();
       
        canvas->BeginText();
        canvas->MoveTextPos(rect8.left + 35, rect8.top - 20);
        canvas->ShowText("Text Icon(ISO8859-2 text)");
        canvas->EndText();
        
        /* Save the document */
        doc->WriteToFile("TextAnnotation.pdf");

    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
    }

    /* Clean up */
    delete doc;

    return 0;
}

