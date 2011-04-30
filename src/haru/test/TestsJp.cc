/*
 * << H a r u --free pdf library >> -- Tests.cc
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
 * << Test program for libharu >>
 */

#include "libharu.h"
#include "libharu_ISO8859.h"
#include "libharu_KOI8.h"
#include "libharu_jpfonts.h"

int main(int argc, char** argv)
{
    const char SAMP_TEXT[] = "ABCDEFGabcdefg\"#$%&'()123*+--=";
    const char SAMP_TEXT2[] = {
        0x82, 0xB1, 0x82, 0xEA, 0x82, 0xCD, 0x53, 0x4A,
        0x49, 0x53, 0x83, 0x52, 0x81, 0x5B, 0x83, 0x68, 
        0x82, 0xCC, 0x93, 0xFA, 0x96, 0x7B, 0x8C, 0xEA, 
        0x83, 0x65, 0x83, 0x4C, 0x83, 0x58, 0x83, 0x67, 
        0x82, 0xC5, 0x82, 0xB7, 0x81, 0x42, 0x00
    };

    const char SAMP_TEXT3[] = {
        0xd3, 0xd4, 0xd5, 0xd6, 0xd9, 0xda, 0xdb, 0xdc,
        0xf3, 0xf4, 0xf5, 0xf6, 0xf9, 0xfa, 0xfb, 0xfc,
        0x00
    };
    
    PdfDoc *doc = new PdfDoc();

    try {
        PDF_DEBUG_PRINT(("doc->NewDoc()\n"));
        doc->NewDoc();

        PDF_DEBUG_PRINT(("doc->AddType1Font()\n"));
        doc->AddType1Font(new PdfTimesRomanFontDef());
        doc->AddType1Font(new PdfHelveticaFontDef());
    
        /* Add a page to the document. */
        PDF_DEBUG_PRINT(("doc->AddPage()\n"));
        PdfPage *page = doc->AddPage();

        /* Set page size to 400 x 600. */
        PDF_DEBUG_PRINT(("page->SetSize()\n"));
        page->SetSize(400, 600);

        /* Get the canvas object of the page. */
        PDF_DEBUG_PRINT(("page->Canvas()\n"));
        PdfContents *canvas = page->Canvas();

        /*Set current font to "Times-Roman" and set the size to 80. */
        PDF_DEBUG_PRINT(("canvas->SetFontAndSize()\n"));
        canvas->SetFontAndSize("Times-Roman", 20);

        /* Start to print text. */
        PDF_DEBUG_PRINT(("canvas->BeginText()\n"));
        canvas->BeginText();

        /* Move the position of the text to 150,400. */
        PDF_DEBUG_PRINT(("canvas->MoveTextPos()\n"));
        canvas->MoveTextPos(40, 500);

        /* Print text. */
        PDF_DEBUG_PRINT(("canvas->ShowText()\n"));
        canvas->ShowText(SAMP_TEXT);
        canvas->SetTextLeading(25);
        canvas->SetFontAndSize("Helvetica", 20);
        canvas->ShowTextNextLine(SAMP_TEXT);

        /* TextWidths */
        PDF_DEBUG_PRINT(("canvas->TextWidths(SAMP_TEXT)\n"));
        float w = canvas->TextWidth(SAMP_TEXT);
        
        /* Finish to print text. */
        PDF_DEBUG_PRINT(("canvas->EndText()\n"));
        canvas->EndText();

        /* Draw under line for text. */
        PDF_DEBUG_PRINT(("Drawing line\n"));
        canvas->SetLineWidth(0.5);
        canvas->MoveTo(40, 500 - 30);
        canvas->LineTo(40 + w, 500 - 30);
        canvas->Stroke();

        /* CMap for sjis encording (fix width font). */
        PDF_DEBUG_PRINT(("Creating cmap object\n"));
        PdfCMap* map_sjis1 = new PdfCMap_90msp_RKSJ_H();
        doc->RegisterObject(map_sjis1);
        
        /* Create PGothic font and register it to PdfDoc object */
        PDF_DEBUG_PRINT(("Creating Type0 Font.\n"));
        doc->AddType0Font(new PdfPGothicFontDef(), "pgothic-sjis", map_sjis1);

        /* Start to print text. */
        PDF_DEBUG_PRINT(("canvas->BeginText()\n"));
        canvas->BeginText();
        canvas->SetFontAndSize("pgothic-sjis", 14);

        /* Move the position of the text to 150,400. */
        PDF_DEBUG_PRINT(("canvas->MoveTextPos()\n"));
        canvas->MoveTextPos(40, 440);

        /* Print text. */
        PDF_DEBUG_PRINT(("canvas->ShowText()\n"));
        canvas->ShowText(SAMP_TEXT2);
        
        /* TextWidths */
        PDF_DEBUG_PRINT(("canvas->TextWidths(SAMP_TEXT2)\n"));
        w = canvas->TextWidth(SAMP_TEXT2);

        canvas->EndText();

        /* Draw under line for text. */
        PDF_DEBUG_PRINT(("Drawing line\n"));
        canvas->SetLineWidth(0.5);
        canvas->SetRGBStroke(255, 0, 0);
        canvas->MoveTo(40, 440 - 5);
        canvas->LineTo(40 + w, 440 - 5);
        canvas->Stroke();

        /* Load type1 font file */
        PDF_DEBUG_PRINT(("Loading Type1 font.\n"));
        PdfType1FontDef* fd2 = new PdfType1FontDef();
        
        try {
#ifdef _WIN32
            fd2->LoadFromFile("..\\examples\\type1\\a010013l.afm", 
                    "..\\examples\\type1\\a010013l.pfb");
#else
            fd2->LoadFromFile("../examples/type1/a010013l.afm", 
                    "../examples/type1/a010013l.pfb");
#endif
        } catch (PDF_STD_EXCEPTION& e) {
            fprintf(stderr, "ERROR: failed to load a010013l font.\n");
            throw;
        }

        /* create encodingdef */
        PdfEncodingDef* encoding = new PdfEncoding_ISO8859_2();
        PDF_DEBUG_PRINT(("doc->AddType1Font\n"));
        doc->AddType1Font(fd2, "iso8859-2 Encoding", encoding);

        PdfEncodingDef* encoding2 = new PdfEncoding_KOI8_R();

        /* Draw iso8859-2 Encoding text. */
        PDF_DEBUG_PRINT(("doc->SetFontAndSize\n"));
        canvas->SetFontAndSize("iso8859-2 Encoding", 30);
        canvas->BeginText();
        canvas->MoveTextPos(40, 400);
        canvas->ShowText(SAMP_TEXT3);
        canvas->EndText();

        /* make outline with unicode */
        PDF_DEBUG_PRINT(("doc->Outlines\n"));
        PdfOutlineRoot *outline_root = doc->Outlines();
        PdfOutlineItem *outline_item = new PdfOutlineItem(outline_root);
        PDF_DEBUG_PRINT(("doc->SetTitleMb\n"));
        outline_item->SetTitleMb(SAMP_TEXT2, map_sjis1);
        outline_item->SetOpened(true);
        outline_item = new PdfOutlineItem(outline_root);
        outline_item->SetTitle(SAMP_TEXT3, encoding);
        outline_item = new PdfOutlineItem(outline_root);
        outline_item->SetTitle(SAMP_TEXT3);
        outline_item = new PdfOutlineItem(outline_root);
        outline_item->SetTitle(SAMP_TEXT3, encoding2);

        /* Save the document as "Test.pdf" */
        PDF_DEBUG_PRINT(("canvas->WriteToFile()\n"));
        doc->WriteToFile("Test.pdf");

        /* Set RGBFill */
        canvas->SetRGBFill(96, 192, 255);
        canvas->SetRGBStroke(255, 192, 96);
        canvas->Rectangle(40, 360, 140, 20);
        canvas->FillStroke();

        canvas->SetCMYKFill(0, 1, 0, 0);
        canvas->SetCMYKStroke(0, 0, 1, 0);
        canvas->Rectangle(200, 360, 140, 20);
        canvas->FillStroke();

        /* Set Infomation of the documant */
        PdfInfo* info = doc->Info();
        info->SetAuthor(SAMP_TEXT3, encoding);
        info->SetTitleMb(SAMP_TEXT2, map_sjis1);
    
        /* Save the document with using GetBuf method. */
        PDF_DEBUG_PRINT(("canvas->GetBuf()\n"));
        PdfMemStream* strm = new PdfMemStream();
        doc->WriteToStream(strm);
        unsigned int siz = strm->BufSize();
        
        FILE* f = fopen("Test2.pdf", "wb");
        for (int i = 0; i < strm->GetBufCount(); i++) {
            unsigned int siz = 0;
            const void* buffer = strm->GetBuf(i, &siz);
            fwrite(buffer, siz, 1, f);
        }
        fclose(f);

    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
    }

    /* Clean up */
    delete doc;

    return 0;
}

