/*
 * << H a r u -- Free PDF Library >> -- JpegExample.pdf
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
#ifndef NOJPEG

double XPos = 20;
double YPos = 500;

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
    
    strcpy(filename1, "images");
    strcat(filename1, FILE_SEPARATOR);
    strcat(filename1, filename);

    PdfJpegImage* image = new PdfJpegImage(doc);
    try {
        /* Load the image. */
        image->LoadFromFile(filename1);

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
        
    } catch (PDF_STD_EXCEPTION& e) {
        delete image;
        throw;
    }   
}

int main()
{
    PdfDoc *doc = new PdfDoc();
    try {
        doc->NewDoc();
        PdfPage *page = doc->AddPage();
        page->SetSize(650, 500);

        PdfDestination* dst = new PdfDestination(page);
        dst->SetXYZ(0, page->Height(), 1);
        doc->Catalog()->SetOpenAction(dst);
        
        PdfContents *canvas = page->Canvas();       

        doc->AddType1Font(new PdfHelveticaFontDef());
        
        canvas->BeginText();
        canvas->SetFontAndSize("Helvetica", 20);
        canvas->MoveTextPos(270, canvas->Height() - 70);
        canvas->ShowText("JpegExample");
        canvas->EndText();
        
        canvas->SetFontAndSize("Helvetica", 12);
        
        DrawImage(doc, canvas, "rgb.jpg", 70, canvas->Height() - 410,
                "24bit color image");
        DrawImage(doc, canvas, "gray.jpg", 340, canvas->Height() - 410,
                "8bit grayscale image");
        
        doc->WriteToFile("JpegExample.pdf");
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
           "make Haru without NOJPEG option.\n");
    return 0;
}

#endif

