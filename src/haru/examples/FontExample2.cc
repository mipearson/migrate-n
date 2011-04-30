/*
 * << H a r u -- Free PDF Library >> -- FontExample2.cc
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
#include <math.h>

void ShowStripesPattern(PdfContents* canvas, double x, double y);
void ShowDescription(PdfContents* canvas, double x, double y, char* text);

int main ()
{
    const char* page_title = "Font Example2";
    const char* samp_text = "abcdefgABCDEFG123!#$%&@?";
    unsigned int cnt = (unsigned int)strlen(samp_text);

    PdfDoc* doc = new PdfDoc();

    try {
        /* Start creating PDF document. */
        doc->NewDoc();
        
        /* Helvetica */
        doc->AddType1Font(new PdfHelveticaFontDef());

        /* Add a new page object. */
        PdfPage* page = doc->AddPage();
        PdfContents* canvas = page->Canvas();
        
        /* Print the lines of the page. */
        canvas->SetLineWidth(1);
        canvas->Rectangle(50, 50, canvas->Width() - 100, 
                canvas->Height() - 110);
        canvas->Stroke();

        /* Print the title of the page (with positioning center). */
        canvas->SetFontAndSize("Helvetica", 24);
        double w = canvas->TextWidth(page_title);
        canvas->BeginText();
        canvas->MoveTextPos((canvas->Width() - w) / 2, canvas->Height() - 50);
        canvas->ShowText(page_title);
        canvas->EndText();

        /* Output subtitle. */
        canvas->BeginText();
        canvas->MoveTextPos(60, canvas->Height() - 80);
        canvas->SetFontAndSize("Helvetica", 16);
        canvas->ShowText("<Various properties of fonts>");

        /* Font size */
        double fsize = 8;
        while (fsize < 60) {
            char buf[50];
            
            /* Set style and size of font. */ 
            canvas->SetFontAndSize("Helvetica", fsize);

            /* Set the position of the text. */
            canvas->MoveTextPos(0, -10 - fsize);

            /* Measure the number of characters which included in the page. */
            strcpy(buf, samp_text);
            unsigned int len = canvas->MeasureText(samp_text, 
                    canvas->Width() - 120);

            /* Truncate the text. */
            buf[len] = 0x00;
            
            canvas->ShowText(buf);

            /* print the description. */
            canvas->MoveTextPos(0, -10);
            canvas->SetFontAndSize("Helvetica", 8);
#ifdef __WIN32__
            _snprintf(buf, 50, "Fontsize=%.0f", fsize);
#else
            snprintf(buf, 50, "Fontsize=%.0f", fsize);
#endif
            canvas->ShowText(buf);
            
            fsize *= 1.5;
        }

        /* Font color */
        canvas->SetFontAndSize("Helvetica", 18);
        canvas->MoveTextPos(0, -50);
        for (unsigned int i = 0; i < cnt; i++) {
            char buf[2];
            double r = (double)i / (double)cnt;
            double g = 1 - ((double)i / (double)cnt);
            buf[0] = samp_text[i];
            buf[1] = 0x00;
            
            canvas->SetRGBFill(r, g, 0.0);
            canvas->ShowText(buf);
        }
        canvas->MoveTextPos(0, -25);

        for (unsigned int j = 0; j < cnt; j++) {
            char buf[2];
            double r = (double)j / (double)cnt;
            double b = 1 - ((double)j / (double)cnt);
            buf[0] = samp_text[j];
            buf[1] = 0x00;
            
            canvas->SetRGBFill(r, 0.0, b);
            canvas->ShowText(buf);
        }
        canvas->MoveTextPos(0, -25);
                
        for (unsigned int k = 0; k < cnt; k++) {
            char buf[2];
            double b = (double)k / (double)cnt;
            double g = 1 - ((double)k / (double)cnt);
            buf[0] = samp_text[k];
            buf[1] = 0x00;

            canvas->SetRGBFill(0.0, g, b);
            canvas->ShowText(buf);
        }

        double save_ypos = canvas->TextPos().y - 20;
        canvas->EndText();
        
        /* Font rendering mode */
        canvas->SetFontAndSize("Helvetica", 40);
        canvas->SetRGBFill(0.5, 0.5, 0.0);
        canvas->SetLineWidth(1.5);
        
        /* PDF_FILL */
        ShowDescription(canvas, 60, save_ypos - 50, "RenderingMode=PDF_FILL");
        canvas->SetTextRenderingMode(PDF_FILL);
        canvas->BeginText();
        canvas->MoveTextPos(60, save_ypos - 50);
        canvas->ShowText("ABCabc123");
        canvas->EndText();

        /* PDF_STROKE */
        ShowDescription(canvas, 60, save_ypos - 110,
                "RenderingMode=PDF_STROKE");
        canvas->SetTextRenderingMode(PDF_STROKE);
        canvas->BeginText();
        canvas->MoveTextPos(60, save_ypos - 110);
        canvas->ShowText("ABCabc123");
        canvas->EndText();

        /* PDF_FILL_THEN_STROKE */
        ShowDescription(canvas, 60, save_ypos - 170,
                "RenderingMode=PDF_FILL_THEN_STROKE");
        canvas->SetTextRenderingMode(PDF_FILL_THEN_STROKE);
        canvas->BeginText();
        canvas->MoveTextPos(60, save_ypos - 170);
        canvas->ShowText("ABCabc123");
        canvas->EndText();

        /* PDF_FILL_CLIPPING */
        ShowDescription(canvas, 60, save_ypos - 230, 
                "RenderingMode=PDF_FILL_CLIPPING");
        canvas->GSave();
        canvas->SetTextRenderingMode(PDF_FILL_CLIPPING);
        canvas->BeginText();
        canvas->MoveTextPos(60, save_ypos - 230);
        canvas->ShowText("ABCabc123");
        canvas->EndText();
        ShowStripesPattern(canvas, 60, save_ypos - 230);
        canvas->GRestore();
        
        /* PDF_STROKE_CLIPPING */
        ShowDescription(canvas, 60, save_ypos - 290, 
                "RenderingMode=PDF_STROKE_CLIPPING");
        canvas->GSave();
        canvas->SetTextRenderingMode(PDF_STROKE_CLIPPING);
        canvas->BeginText();
        canvas->MoveTextPos(60, save_ypos - 290);
        canvas->ShowText("ABCabc123");
        canvas->EndText();
        ShowStripesPattern(canvas, 60, save_ypos - 290);
        canvas->GRestore();

        /* PDF_FILL_STROKE_CLIPPING */
        ShowDescription(canvas, 60, save_ypos - 350,
                "RenderingMode=PDF_FILL_STROKE_CLIPPING");
        canvas->GSave();
        canvas->SetTextRenderingMode(PDF_FILL_STROKE_CLIPPING);
        canvas->BeginText();
        canvas->MoveTextPos(60, save_ypos - 350);
        canvas->ShowText("ABCabc123");
        canvas->EndText();
        ShowStripesPattern(canvas, 60, save_ypos - 350);
        canvas->GRestore();

        /* Reset text attributes */ 
        canvas->SetTextRenderingMode(PDF_FILL);
        canvas->SetRGBFill(0, 0, 0);
        canvas->SetFontAndSize("Helvetica", 30);
        
        /* Rotating text */
        double angle = 30;                   /* A rotation of 30 degrees. */
        double rad = angle / 180 * 3.141592; /* Calcurate the radian value. */
    
        ShowDescription(canvas, 300, save_ypos - 100, "Rotating text");
        canvas->BeginText();
        canvas->SetTextMatrix(cos(rad), sin(rad), -sin(rad), cos(rad), 
                310, save_ypos - 100);
        canvas->ShowText("ABCabc123");
        canvas->EndText();

        /* Skewing text. */
        ShowDescription(canvas, 300, save_ypos - 160, "Skewing text");
        canvas->BeginText();
        
        double angle1 = 10;
        double angle2 = 20;
        double rad1 = angle1 / 180 * 3.141592;
        double rad2 = angle2 / 180 * 3.141592;

        canvas->SetTextMatrix(1, tan(rad1), tan(rad2), 1, 
                300, save_ypos - 160);
        canvas->ShowText("ABCabc123");
        canvas->EndText();

        /* Scaling text (X direction) */
        ShowDescription(canvas, 300, save_ypos - 220, 
                "Scaling text (X direction)");
        canvas->BeginText();
        canvas->SetTextMatrix(1.5, 0, 0, 1, 300, save_ypos - 220);
        canvas->ShowText("ABCabc12");
        canvas->EndText();

        /* Scaling text (Y direction) */
        ShowDescription(canvas, 300, save_ypos - 300, 
                "Scaling text (Y direction)");
        canvas->BeginText();
        canvas->SetTextMatrix(1, 0, 0, 2, 300, save_ypos - 300);
        canvas->ShowText("ABCabc123");
        canvas->EndText();
        
        /* Save the document to a file */
        doc->WriteToFile("FontExample2.pdf");
        doc->FreeDoc();
    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
    }
    delete doc;
    return 0;
}

void ShowStripesPattern(PdfContents* canvas, double x, double y)
{
    int iy = 0;
    while (iy < 50) {
        canvas->SetRGBStroke(0.0, 0.0, 0.5);
        canvas->SetLineWidth(1);
        canvas->MoveTo(x, y + iy);
        canvas->LineTo(x + canvas->TextWidth("ABCabc123"), y + iy);
        canvas->Stroke();
        iy += 3;
    }
    canvas->SetLineWidth(2.5);
}

void ShowDescription(PdfContents* canvas, double x, double y, char* text)
{
    double fsize = canvas->FontSize();
    const char* fname = canvas->FontName();
    pdf_rgb_color c = canvas->RGBFill();
    
    canvas->BeginText();
    canvas->MoveTextPos(x, y - 12);
    canvas->SetRGBFill(0, 0, 0);
    canvas->SetTextRenderingMode(PDF_FILL);
    canvas->SetFontAndSize("Helvetica", 10);
    canvas->ShowText(text);
    canvas->EndText();
    
    canvas->SetFontAndSize(fname, fsize);
    canvas->SetRGBFill(c);
}

