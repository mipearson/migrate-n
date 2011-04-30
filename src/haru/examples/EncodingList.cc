/*
 * << H a r u --free pdf library >> -- EncodingList.cpp
 *
 * Copyright (c) 1999-2004 Takeshi Kanno <takeshi_kanno@est.hi-ho.ne.jp>
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
#include "libharu_MSCP.h"
#include "libharu_KOI8.h"

/*
 * Character lists for various encodings.  
 */

void DrawGraph(PdfContents* canvas);
void DrawFonts(PdfContents* canvas);

static const int PAGE_WIDTH = 420;
static const int PAGE_HEIGHT = 400;
static const int CELL_WIDTH = 20;
static const int CELL_HEIGHT = 20;
static const int CELL_HEADER = 10;

int main()
{

    /* Create a PdfDoc object and start making a PDF document. */
    PdfDoc* doc = new PdfDoc();
    try {
        doc->NewDoc();
        doc->Catalog()->SetPageMode(PDF_USE_OUTLINES);

        /* Make helvetica font object. */
        PdfType1FontDef* fd1 = new PdfHelveticaFontDef();

        /* Add font with PdfStandardEncoding and name it "helvetica" */
        doc->AddType1Font(fd1, "helvetica");

        /* Load courier font objects. */
        PdfType1FontDef* fd2 = new PdfType1FontDef();
        try {
        #ifdef __WIN32__
            fd2->LoadFromFile("type1\\a010013l.afm", "type1\\a010013l.pfb");
        #else
            fd2->LoadFromFile("type1/a010013l.afm", "type1/a010013l.pfb");
        #endif
        } catch (PdfException& e) {
            fprintf(stderr, "ERROR: failed to load URWGothicL font.\n");
            throw;
        }

        /* Add font with PdfStandardEncoding and name it "StandardEncoding" */
        doc->AddType1Font(fd2, "StandardEncoding", new PdfStandardEncoding());

        /* Add font with PdfWinAnsiEncoding and name it "WinAnsiEncoding" */
        doc->AddType1Font(fd2, "WinAnsiEncoding", new PdfWinAnsiEncoding());

        /* Add font with PdfMacRomanEncoding and name it "MacRomanEncoding" */
        doc->AddType1Font(fd2, "MacRomanEncoding", new PdfMacRomanEncoding());

        /*
         * add font with various encodings.
         */
        doc->AddType1Font(fd2, "iso8859-2 Encoding",
                new PdfEncoding_ISO8859_2());
        
        doc->AddType1Font(fd2, "iso8859-3 Encoding",
                new PdfEncoding_ISO8859_3());
        
        doc->AddType1Font(fd2, "iso8859-4 Encoding",
                new PdfEncoding_ISO8859_4());
        
        doc->AddType1Font(fd2, "iso8859-5 Encoding",
                new PdfEncoding_ISO8859_5());
        
        doc->AddType1Font(fd2, "iso8859-9 Encoding",
                new PdfEncoding_ISO8859_9());
        
        doc->AddType1Font(fd2, "iso8859-10 Encoding",
                new PdfEncoding_ISO8859_10());
        
        doc->AddType1Font(fd2, "iso8859-13 Encoding",
                new PdfEncoding_ISO8859_13());
        
        doc->AddType1Font(fd2, "iso8859-15 Encoding",
                new PdfEncoding_ISO8859_15());
        
        doc->AddType1Font(fd2, "iso8859-16 Encoding",
                new PdfEncoding_ISO8859_16());
       
        doc->AddType1Font(fd2, "CP1251 Encoding", new PdfEncoding_CP1251());

        doc->AddType1Font(fd2, "CP1252 Encoding", new PdfEncoding_CP1252());

        doc->AddType1Font(fd2, "CP1254 Encoding", new PdfEncoding_CP1254());

        doc->AddType1Font(fd2, "CP1257 Encoding", new PdfEncoding_CP1257());

        doc->AddType1Font(fd2, "KOI8-R Encoding", new PdfEncoding_KOI8_R());

        /* Add Symbolic font */
        doc->AddType1Font(new PdfSymbolFontDef(), "Symbol Set");

        /* Add ZapfDingbats font */
        doc->AddType1Font(new PdfZapfDingbatsFontDef(), "ZapfDingbats Set");

        /* Create outline root */
        PdfOutlineRoot *outline_root = doc->Outlines();
        PdfOutlineItem *outline_item = new PdfOutlineItem(outline_root);
        outline_item->SetTitle("Encoding list");
        outline_item->SetOpened(true);
   
        for (unsigned int i = 1; i < doc->FontMgr()->CountFonts(); i++) {
            /* Add new page and set the size of the page */
            PdfPage* page = doc->AddPage();
            page->SetSize(PAGE_WIDTH, PAGE_HEIGHT);

            /* Add the page to outline. */
            PdfOutlineItem *outline_item1 = new PdfOutlineItem(outline_item);
            PdfDestination *dst = new PdfDestination(page);
            dst->SetXYZ(0, page->Height(), 1);
            outline_item1->SetDestination(dst);

            /* Get font name from FontManager object. */
            PdfFont* f = doc->FontMgr()->GetFont(i);
            const char* fontname = f->Name();
            outline_item1->SetTitle(fontname);

            /* Draw graphs on the canvas. */
            PdfContents* canvas = page->Canvas();
            canvas->AddFilter(PDF_FILTER_DEFLATE);
            DrawGraph(canvas);

            canvas->BeginText();
            canvas->SetFontAndSize("helvetica", 20);
            canvas->MoveTextPos(40, PAGE_HEIGHT - 50);
            canvas->ShowText(f->Name());
            canvas->EndText();

            /* Draw the gryphs of the font. */
            canvas->SetFontAndSize(fontname, 14);
            DrawFonts(canvas);
        }

        /* Save the document to the file. */
        doc->WriteToFile("EncodingList.pdf");

    } catch (PdfException& e) {
        fprintf(stderr, "%s\n", e.what());
        delete doc;
        return 1;
    }
    delete doc;
    return 0;
}

void DrawGraph(PdfContents* canvas)
{
    /* Draw 16 X 15 cells */
    char buf[50];

    /* Draw vertical lines. */
    canvas->SetLineWidth(0.5);
    canvas->SetFontAndSize("helvetica", 15);
    for (int i = 0; i <= 17; i++) {
        int x = i * CELL_WIDTH + 40;

        canvas->MoveTo(x, PAGE_HEIGHT - 60);
        canvas->LineTo(x, 40);
        canvas->Stroke();

        if (i > 0 && i <= 16) {
            canvas->BeginText();
            canvas->MoveTextPos(x + 5, PAGE_HEIGHT - 75);
#ifdef __WIN32__
            _snprintf(buf, 5, "%X", i - 1);
#else
            snprintf(buf, 5, "%X", i - 1);
#endif
            canvas->ShowText(buf);
            canvas->EndText();
        }
    }

    /* Draw horizontal lines. */
    for (int j = 0; j <= 15; j++) {
       int y = j * CELL_HEIGHT + 40;

        canvas->MoveTo(40, y);
        canvas->LineTo(PAGE_WIDTH - 40, y);
        canvas->Stroke();

        if (j < 14) {
            canvas->BeginText();
            canvas->MoveTextPos(45, y + 5);
#ifdef __WIN32__
            _snprintf(buf, 5, "%X", 15 - j);
#else
            snprintf(buf, 5, "%X", 15 - j);
#endif
            canvas->ShowText(buf);
            canvas->EndText();
        }
    }
}

void DrawFonts(PdfContents* canvas)
{
    /* Draw all character from 0x20 to 0xFF to the canvas. */
    for (int i = 1; i < 17; i++) {
        for (int j = 1; j < 17; j++) {
            unsigned char buf[2];

            buf[1] = 0x00;
            int y = PAGE_HEIGHT - 55 - ((i - 1) * CELL_HEIGHT);
            int x = j * CELL_WIDTH + 50;

            buf[0] = (i - 1) * 16 + (j - 1);
            if (buf[0] >= 32) {
                canvas->BeginText();
                double d = x - canvas->TextWidth((char*)buf) / 2;
                canvas->MoveTextPos(d, y);
                canvas->ShowText((char*)buf);
                canvas->EndText();
            }
        }
    }
}

