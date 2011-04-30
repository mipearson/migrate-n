/*
 * << H a r u -- Free PDF Library >> -- PngExample.pdf
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
#ifndef NOPNG

void DrawImage(PdfDoc* doc, PdfContents* canvas, char* filename, 
        double x, double y, char* text);

void DrawImage(PdfDoc* doc, PdfContents* canvas, char* filename, 
        double x, double y, char* text) 
{
#ifdef __WIN32__
    const char* FILE_SEPARATOR = "\\";  
#else
    const char* FILE_SEPARATOR = "/";
#endif
    char filename1[255];
    
    strcpy(filename1, "pngsuite");
    strcat(filename1, FILE_SEPARATOR);
    strcat(filename1, filename);

    PdfPngImage* image = new PdfPngImage(doc);
    try {
        /* Load the image. */
        image->LoadFromFile(filename1);

        /* Add the deflater filter. */
        image->AddFilter(PDF_FILTER_DEFLATE);

    } catch (...) {
        delete image;
        throw;
    }   
    /* Register the image to the PdfDoc object. */
    doc->AddXObject(image);
        
    /* Draw image to the canvas. */
    canvas->GSave();
    canvas->Concat(image->Width(), 0, 0, image->Height(), x, y);
    canvas->ExecuteXObject(image);
    canvas->GRestore();

    /* Print the text. */
    canvas->BeginText();
    canvas->SetTextLeading(16);
    canvas->MoveTextPos(x, y);
    canvas->ShowTextNextLine(filename);
    canvas->ShowTextNextLine(text);
    canvas->EndText();
}

int main()
{
    PdfDoc *doc = new PdfDoc();
    try {
        doc->NewDoc();
        PdfPage *page = doc->AddPage();
        page->SetSize(550, 650);

        PdfDestination* dst = new PdfDestination(page);
        dst->SetXYZ(0, page->Height(), 1);
        doc->Catalog()->SetOpenAction(dst);
        
        PdfContents *canvas = page->Canvas();       

        doc->AddType1Font(new PdfHelveticaFontDef());
        
        canvas->BeginText();
        canvas->SetFontAndSize("Helvetica", 20);
        canvas->MoveTextPos(220, canvas->Height() - 70);
        canvas->ShowText("PngExample");
        canvas->EndText();
        
        canvas->SetFontAndSize("Helvetica", 12);
        
        DrawImage(doc, canvas, "basn0g01.png", 100, canvas->Height() - 150,
                "1bit grayscale.");
        DrawImage(doc, canvas, "basn0g02.png", 200, canvas->Height() - 150,
                "2bit grayscale.");
        DrawImage(doc, canvas, "basn0g04.png", 300, canvas->Height() - 150,
                "4bit grayscale.");
        DrawImage(doc, canvas, "basn0g08.png", 400, canvas->Height() - 150,
                "8bit grayscale.");
    
        DrawImage(doc, canvas, "basn2c08.png", 100, canvas->Height() - 250,
                "8bit color.");
        DrawImage(doc, canvas, "basn2c16.png", 200, canvas->Height() - 250,
                "16bit color.");

        DrawImage(doc, canvas, "basn3p01.png", 100, canvas->Height() - 350,
                "01bit pallet.");
        DrawImage(doc, canvas, "basn3p02.png", 200, canvas->Height() - 350,
                "2bit pallet.");
        DrawImage(doc, canvas, "basn3p04.png", 300, canvas->Height() - 350,
                "4bit pallet.");
        DrawImage(doc, canvas, "basn3p08.png", 400, canvas->Height() - 350,
                "8bit pallet.");
        
        DrawImage(doc, canvas, "basn4a08.png", 100, canvas->Height() - 450,
                "8bit alpha.");
        DrawImage(doc, canvas, "basn4a16.png", 200, canvas->Height() - 450,
                "16bit alpha.");

        DrawImage(doc, canvas, "basn6a08.png", 100, canvas->Height() - 550,
                "8bit alpha.");
        DrawImage(doc, canvas, "basn6a16.png", 200, canvas->Height() - 550,
                "16bit alpha.");
        
        doc->WriteToFile("PngExample.pdf");
    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
        exit(1);
    }
    delete doc;

    return 0;
}

#else

int main()
{
    printf("WARNING: if you want to run this example, \n"
           "make Haru without NOPNG option.\n");
    return 0;
}

#endif

