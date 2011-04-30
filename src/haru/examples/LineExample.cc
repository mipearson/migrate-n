/*
 * << H a r u -- Free PDF Library >> -- LineExample.cc
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

void DrawLine(PdfContents* canvas, double x, double y, const char* label);
void DrawRect(PdfContents* canvas, double x, double y, const char* label);

void DrawLine(PdfContents* canvas, double x, double y, const char* label)
{
    canvas->SetFontAndSize("Helvetica", 10);
    canvas->BeginText();
    canvas->MoveTextPos(x, y - 10);
    canvas->ShowText(label);
    canvas->EndText();
    
    canvas->MoveTo(x, y - 15);
    canvas->LineTo(x + 220, y - 15);
    canvas->Stroke();
}

void DrawLine2(PdfContents* canvas, double x, double y, const char* label)
{
    canvas->SetFontAndSize("Helvetica", 10);
    canvas->BeginText();
    canvas->MoveTextPos(x, y);
    canvas->ShowText(label);
    canvas->EndText();

    canvas->MoveTo(x + 30, y - 25);
    canvas->LineTo(x + 160, y - 25);
    canvas->Stroke();
}   

void DrawRect(PdfContents* canvas, double x, double y, const char* label)
{
    canvas->SetFontAndSize("Helvetica", 10);
    canvas->BeginText();
    canvas->MoveTextPos(x, y - 10);
    canvas->ShowText(label);
    canvas->EndText();
    
    canvas->Rectangle(x, y - 40, 220, 25);
}

int main ()
{
    const char* page_title = "Line Example";
    PdfDoc* doc = new PdfDoc();

    try {
        /* Start creating PDF document. */
        doc->NewDoc();
        
        /* Add Helvetica Font. */
        doc->AddType1Font(new PdfHelveticaFontDef());
        doc->AddType1Font(new PdfHelveticaObliqueFontDef());

        /* Add a new page object. */
        PdfPage* page = doc->AddPage();
        PdfContents* canvas = page->Canvas();
        
        /* Print the lines of the page. */
        canvas->SetLineWidth(1);
        canvas->Rectangle(50, 50, canvas->Width() - 100, canvas->Height() - 110);
        canvas->Stroke();

        /* Print the title of the page (with positioning center). */
        canvas->SetFontAndSize("Helvetica", 24);
        double w = canvas->TextWidth(page_title);
        canvas->BeginText();
        canvas->MoveTextPos((canvas->Width() - w) / 2, canvas->Height() - 50);
        canvas->ShowText(page_title);
        canvas->EndText();

        /* Draw verious widths of lines. */
        canvas->SetLineWidth(0);
        DrawLine(canvas, 60, 770, "line width = 0");
        
        canvas->SetLineWidth(1.0);
        DrawLine(canvas, 60, 740, "line width = 1.0");

        canvas->SetLineWidth(2.0);
        DrawLine(canvas, 60, 710, "line width = 2.0");

        /* Line dash pattern */ 
        canvas->SetLineWidth(1.0);
        
        canvas->SetDash(3, 0, 0);
        DrawLine(canvas, 60, 680, "SetDash(3, 0, 0) -- "
                "3 units on, 3 units off,...");

        canvas->SetDash(4, 0, 2);
        DrawLine(canvas, 60, 650, "SetDash(4, 2, 0) -- "
                "2 on, 4 off, 4 on, 4 off,...");

        canvas->SetDash(3, 7, 2);
        DrawLine(canvas, 60, 620, "SetDash(3, 7, 2) -- "
                "2 off, 3 on, 7 off, 3 on, 7 off,..");

        canvas->SetLineWidth(30);
        canvas->SetDash(0, 0, 0);
        canvas->SetRGBStroke(0.0, 0.5, 0.0);

        /* Line Cap Style */
        canvas->SetLineCap(PDF_BUTT_END);
        DrawLine2(canvas, 60, 570, "PDF_BUTT_END");
        
        canvas->SetLineCap(PDF_ROUND_END);
        DrawLine2(canvas, 60, 505, "PDF_ROUND_END");
        
        canvas->SetLineCap(PDF_PROJECTING_SCUARE_END);
        DrawLine2(canvas, 60, 440, "PDF_PROJECTING_SCUARE_END");

        /* Line Join Style */
        canvas->SetLineWidth(30);
        canvas->SetRGBStroke(0.0, 0.0, 0.5);
        
        canvas->SetLineJoin(PDF_MITER_JOIN);
        canvas->MoveTo(120, 300);
        canvas->LineTo(160, 340);
        canvas->LineTo(200, 300);
        canvas->Stroke();
        
        canvas->BeginText();
        canvas->MoveTextPos(60, 360);
        canvas->ShowText("PDF_MITER_JOIN");
        canvas->EndText();

        canvas->SetLineJoin(PDF_ROUND_JOIN);
        canvas->MoveTo(120, 195);
        canvas->LineTo(160, 235);
        canvas->LineTo(200, 195);
        canvas->Stroke();
        
        canvas->BeginText();
        canvas->MoveTextPos(60, 255);
        canvas->ShowText("PDF_ROUND_JOIN");
        canvas->EndText();

        canvas->SetLineJoin(PDF_BEVEL_JOIN);
        canvas->MoveTo(120, 90);
        canvas->LineTo(160, 130);
        canvas->LineTo(200, 90);
        canvas->Stroke();
        
        canvas->BeginText();
        canvas->MoveTextPos(60, 150);
        canvas->ShowText("PDF_BEVEL_JOIN");
        canvas->EndText();

        /* Draw Rectangle */
        canvas->SetLineWidth(2);
        canvas->SetRGBStroke(0, 0, 0);
        canvas->SetRGBFill(0.75, 0.0, 0.0);
        
        DrawRect(canvas, 300, 770, "Stroke");
        canvas->Stroke();
        
        DrawRect(canvas, 300, 720, "Fill");
        canvas->Fill();

        DrawRect(canvas, 300, 670, "Fill then Stroke");
        canvas->FillStroke();

        /* Clip Rect */
        canvas->GSave();  /* Save the current graphic state */
        DrawRect(canvas, 300, 620, "Clip Rectangle");
        canvas->Clip();
        canvas->Stroke();
        canvas->SetFontAndSize("Helvetica-Oblique", 13);
        
        canvas->BeginText();
        canvas->MoveTextPos(290, 600);
        canvas->SetTextLeading(12);
        canvas->ShowText("Clip Clip Clip Clip Clip Clipi Clip Clip Clip");
        canvas->ShowTextNextLine("Clip Clip Clip Clip Clip Clip Clip Clip Clip");
        canvas->ShowTextNextLine("Clip Clip Clip Clip Clip Clip Clip Clip Clip");
        canvas->EndText();
        canvas->GRestore();

        /* Curve Example(CurveTo2) */
        canvas->SetRGBFill(0, 0, 0);

        double x = 330;
        double y = 440;
        double x1 = 430;
        double y1 = 530;
        double x2 = 480;
        double y2 = 470;

        canvas->BeginText();
        canvas->MoveTextPos(300, 540);
        canvas->ShowText("CurveTo2(x1, y1, x2. y2)");
        canvas->EndText();
        
        canvas->BeginText();
        canvas->MoveTextPos(x + 5, y - 5);
        canvas->ShowText("Current point");
        canvas->MoveTextPos(x1 - x, y1 - y);
        canvas->ShowText("(x1, y1)");
        canvas->MoveTextPos(x2 - x1, y2 - y1);
        canvas->ShowText("(x2, y2)");
        canvas->EndText();

        canvas->SetDash(3, 0, 0);
        canvas->SetLineWidth(0.5);
        canvas->MoveTo(x1, y1);
        canvas->LineTo(x2, y2);
        canvas->Stroke();
        
        canvas->SetDash(0, 0, 0);
        canvas->SetLineWidth(1.5);
        
        canvas->MoveTo(x, y);
        canvas->CurveTo2(x1, y1, x2, y2);
        canvas->Stroke();
    
        /* Curve Example(CurveTo3) */
        y -= 150;
        y1 -= 150;
        y2 -= 150;
        
        canvas->BeginText();
        canvas->MoveTextPos(300, 390);
        canvas->ShowText("CurveTo3(x1, y1, x2. y2)");
        canvas->EndText();
        
        canvas->BeginText();
        canvas->MoveTextPos(x + 5, y - 5);
        canvas->ShowText("Current point");
        canvas->MoveTextPos(x1 - x, y1 - y);
        canvas->ShowText("(x1, y1)");
        canvas->MoveTextPos(x2 - x1, y2 - y1);
        canvas->ShowText("(x2, y2)");
        canvas->EndText();

        canvas->SetDash(3, 0, 0);
        canvas->SetLineWidth(0.5);
        canvas->MoveTo(x, y);
        canvas->LineTo(x1, y1);
        canvas->Stroke();
        
        canvas->SetDash(0, 0, 0);
        canvas->SetLineWidth(1.5);
        canvas->MoveTo(x, y);
        canvas->CurveTo3(x1, y1, x2, y2);
        canvas->Stroke();

        /* Curve Example(CurveTo) */
        y -= 150;
        y1 -= 160;
        y2 -= 130;
        x2 += 10;
        double x3 = 480;
        double y3 = 90;

        canvas->BeginText();
        canvas->MoveTextPos(300, 240);
        canvas->ShowText("CurveTo(x1, y1, x2. y2, x3, y3)");
        canvas->EndText();
        
        canvas->BeginText();
        canvas->MoveTextPos(x + 5, y - 5);
        canvas->ShowText("Current point");
        canvas->MoveTextPos(x1 - x, y1 - y);
        canvas->ShowText("(x1, y1)");
        canvas->MoveTextPos(x2 - x1, y2 - y1);
        canvas->ShowText("(x2, y2)");
        canvas->MoveTextPos(x3 - x2, y3 - y2);
        canvas->ShowText("(x3, y3)");
        canvas->EndText();

        canvas->SetDash(3, 0, 0);
        canvas->SetLineWidth(0.5);
        canvas->MoveTo(x, y);
        canvas->LineTo(x1, y1);
        canvas->Stroke();
        canvas->MoveTo(x2, y2);
        canvas->LineTo(x3, y3);
        canvas->Stroke();

        canvas->SetDash(0, 0, 0);
        canvas->SetLineWidth(1.5);
        canvas->MoveTo(x, y);
        canvas->CurveTo(x1, y1, x2, y2, x3, y3);
        canvas->Stroke();
        
        /* Save the document to a file */
        doc->WriteToFile("LineExample.pdf");
        doc->FreeDoc();
    } catch (PDF_STD_EXCEPTION& e) {
        fprintf(stderr, "%s\n", e.what());
    }
    delete doc;
    return 0;
}
