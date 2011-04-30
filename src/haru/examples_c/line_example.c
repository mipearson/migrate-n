/*
 * << H a r u -- Free PDF Library >> -- line_example.c
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

#include "libharuc.h"

#define CHECK_ERROR(ret) \
{ \
    int err = ret; \
    if (err != PDF_SUCCESS) { \
        fprintf(stderr, "ERROR (%d) line:%d.\n", err, __LINE__); \
        pdf_doc_free(doc); \
        return err; \
    } \
}

#define CHECK_ERROR2(obj) \
{ \
    void* ret = obj; \
    if (ret == NULL) { \
        fprintf(stderr, "ERROR line:%d.\n", __LINE__); \
        pdf_doc_free(doc); \
        return -1; \
    } \
}

void draw_line(pdf_contents canvas, float x, float y, const char* label);
void draw_rect(pdf_contents canvas, float x, float y, const char* label);

void draw_line(pdf_contents canvas, float x, float y, const char* label)
{
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, x, y - 10);
    pdf_contents_show_text(canvas, label);
    pdf_contents_end_text(canvas);
    
    pdf_contents_move_to(canvas, x, y - 15);
    pdf_contents_line_to(canvas, x + 220, y - 15);
    pdf_contents_stroke(canvas);
}

void draw_line2(pdf_contents canvas, float x, float y, const char* label)
{
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, x, y);
    pdf_contents_show_text(canvas, label);
    pdf_contents_end_text(canvas);

    pdf_contents_move_to(canvas, x + 30, y - 25);
    pdf_contents_line_to(canvas, x + 160, y - 25);
    pdf_contents_stroke(canvas);
}   

void draw_rect(pdf_contents canvas, float x, float y, const char* label)
{
    pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, x, y - 10);
    pdf_contents_show_text(canvas, label);
    pdf_contents_end_text(canvas);
    
    pdf_contents_rectangle(canvas, x, y - 40, 220, 25);
}

int main ()
{
    int ret = 0;
	pdf_type1_fontdef font1_def;
	pdf_type1_fontdef font2_def;
	pdf_page page;
	float w;
	pdf_contents canvas;
	float x = 330;
	float y = 440;
	float x1 = 430;
	float y1 = 530;
	float x2 = 480;
	float y2 = 470;
	float x3;
	float y3;
	
    const char* page_title = "Line Example";
    pdf_doc doc = pdf_doc_new();

    /* Start creating PDF document. */
    ret = pdf_doc_new_doc(doc);
    CHECK_ERROR(ret);
        
    /* Add Helvetica Font. */
    font1_def = pdf_create_type1_fontdef(PDF_FONT_HELVETICA);
    CHECK_ERROR2(font1_def);

    ret = pdf_doc_add_type1font(doc, font1_def, NULL, NULL);
    CHECK_ERROR(ret);

    font2_def = 
        pdf_create_type1_fontdef(PDF_FONT_HELVETICA_OBLIQUE);
    CHECK_ERROR2(font2_def);

    ret = pdf_doc_add_type1font(doc, font2_def, NULL, NULL);
    CHECK_ERROR(ret);

    /* Add a new page object. */
    page = pdf_doc_add_page(doc);
    CHECK_ERROR2(page);

    canvas = pdf_page_get_canvas(page);
    CHECK_ERROR2(canvas);
        
    /* Print the lines of the page. */
    ret = pdf_contents_set_line_width(canvas, 1);
    CHECK_ERROR(ret);

    ret = pdf_contents_rectangle(canvas, 50, 50, 
            pdf_contents_get_width(canvas) - 100, 
            pdf_contents_get_height(canvas) - 110);
    CHECK_ERROR(ret);
    
    ret = pdf_contents_stroke(canvas);
    CHECK_ERROR(ret);

    /* Print the title of the page (with positioning center). */
    ret = pdf_contents_set_font_and_size(canvas, "Helvetica", 24);
    CHECK_ERROR(ret);

    w = pdf_contents_get_text_width(canvas, page_title, NULL, NULL);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_text_pos(canvas, 
            (pdf_contents_get_width(canvas) - w) / 2, 
            pdf_contents_get_height(canvas) - 50);
    CHECK_ERROR(ret);

    ret = pdf_contents_show_text(canvas, page_title);
    CHECK_ERROR(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR(ret);

    /* Draw verious widths of lines. */
    ret = pdf_contents_set_line_width(canvas, 0);
    CHECK_ERROR(ret);

    draw_line(canvas, 60, 770, "line width = 0");
        
    ret = pdf_contents_set_line_width(canvas, 1.0);
    CHECK_ERROR(ret);

    draw_line(canvas, 60, 740, "line width = 1.0");

    ret = pdf_contents_set_line_width(canvas, 2.0);
    CHECK_ERROR(ret);

    draw_line(canvas, 60, 710, "line width = 2.0");
    CHECK_ERROR(ret);


    /* Line dash pattern */ 
    ret = pdf_contents_set_line_width(canvas, 1.0);
    CHECK_ERROR(ret);

    ret = pdf_contents_set_dash(canvas, 3, 0, 0);
    CHECK_ERROR(ret);

    draw_line(canvas, 60, 680, "set_dash(3, 0, 0) -- "
            "3 units on, 3 units off,...");

    ret = pdf_contents_set_dash(canvas, 4, 0, 2);
    CHECK_ERROR(ret);

    draw_line(canvas, 60, 650, "set_dash(4, 2, 0) -- "
            "2 on, 4 off, 4 on, 4 off,...");

    ret = pdf_contents_set_dash(canvas, 3, 7, 2);
    CHECK_ERROR(ret);

    draw_line(canvas, 60, 620, "set_dash(3, 7, 2) -- "
            "2 off, 3 on, 7 off, 3 on, 7 off,..");

    ret = pdf_contents_set_line_width(canvas, 30);
    CHECK_ERROR(ret);

    ret = pdf_contents_set_dash(canvas, 0, 0, 0);
    CHECK_ERROR(ret);

    ret = pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0.5, 0));
    CHECK_ERROR(ret);

    /* Line Cap Style */
    ret = pdf_contents_set_line_cap(canvas, PDF_BUTT_END);
    CHECK_ERROR(ret);

    draw_line2(canvas, 60, 570, "PDF_BUTT_END");
        
    ret =  pdf_contents_set_line_cap(canvas, PDF_ROUND_END);
    CHECK_ERROR(ret);

    draw_line2(canvas, 60, 505, "PDF_ROUND_END");
        
    ret = pdf_contents_set_line_cap(canvas, PDF_PROJECTING_SCUARE_END);
    CHECK_ERROR(ret);

    draw_line2(canvas, 60, 440, "PDF_PROJECTING_SCUARE_END");

    /* Line Join Style */
    ret = pdf_contents_set_line_width(canvas, 30);
    CHECK_ERROR(ret);

    ret = pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0.5));
    CHECK_ERROR(ret);
        
    ret = pdf_contents_set_line_join(canvas, PDF_MITER_JOIN);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_to(canvas, 120, 300);
    CHECK_ERROR(ret);

    ret = pdf_contents_line_to(canvas, 160, 340);
    CHECK_ERROR(ret);

    ret = pdf_contents_line_to(canvas, 200, 300);
    CHECK_ERROR(ret);

    ret = pdf_contents_stroke(canvas);
    CHECK_ERROR(ret);
        
    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, 360);
    CHECK_ERROR(ret);

    ret = pdf_contents_show_text(canvas, "PDF_MITER_JOIN");
    CHECK_ERROR(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_set_line_join(canvas, PDF_ROUND_JOIN);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_to(canvas, 120, 195);
    CHECK_ERROR(ret);

    ret = pdf_contents_line_to(canvas, 160, 235);
    CHECK_ERROR(ret);

    ret = pdf_contents_line_to(canvas, 200, 195);
    CHECK_ERROR(ret);

    ret = pdf_contents_stroke(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, 255);
    CHECK_ERROR(ret);

    ret = pdf_contents_show_text(canvas, "PDF_ROUND_JOIN");
    CHECK_ERROR(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_set_line_join(canvas, PDF_BEVEL_JOIN);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_to(canvas, 120, 90);
    CHECK_ERROR(ret);

    ret = pdf_contents_line_to(canvas, 160, 130);
    CHECK_ERROR(ret);

    ret = pdf_contents_line_to(canvas, 200, 90);
    CHECK_ERROR(ret);

    ret = pdf_contents_stroke(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, 150);
    CHECK_ERROR(ret);

    ret = pdf_contents_show_text(canvas, "PDF_BEVEL_JOIN");
    CHECK_ERROR(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR(ret);

    /* Draw rectangle */
    pdf_contents_set_line_width(canvas, 2);
    pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0));
    pdf_contents_set_rgb_fill(canvas, PdfRGBColor(0.75, 0, 0));
        
    draw_rect(canvas, 300, 770, "stroke");
    pdf_contents_stroke(canvas);
        
    draw_rect(canvas, 300, 720, "Fill");
    pdf_contents_fill(canvas);

    draw_rect(canvas, 300, 670, "Fill then stroke");
    pdf_contents_fill_stroke(canvas);

    /* clip Rect */
    pdf_contents_gsave(canvas);  /* Save the current graphic state */
    draw_rect(canvas, 300, 620, "clip rectangle");
    pdf_contents_clip(canvas);
    pdf_contents_stroke(canvas);
    pdf_contents_set_font_and_size(canvas, "Helvetica-Oblique", 13);
        
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, 290, 600);
    pdf_contents_set_text_leading(canvas, 12);
    pdf_contents_show_text(canvas, 
            "clip clip clip clip clip clipi clip clip clip");
    pdf_contents_show_text_next_line(canvas,
            "clip clip clip clip clip clip clip clip clip");
    pdf_contents_show_text_next_line(canvas,
            "clip clip clip clip clip clip clip clip clip");
    pdf_contents_end_text(canvas);
    pdf_contents_grestore(canvas);

    /* Curve Example(curve_to2) */
    pdf_contents_set_rgb_fill(canvas, PdfRGBColor(0, 0, 0));

    x = 330;
    y = 440;
    x1 = 430;
    y1 = 530;
    x2 = 480;
    y2 = 470;

    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, 300, 540);
    pdf_contents_show_text(canvas, "curve_to2(x1, y1, x2. y2)");
    pdf_contents_end_text(canvas);
        
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, x + 5, y - 5);
    pdf_contents_show_text(canvas, "Current point");
    pdf_contents_move_text_pos(canvas, x1 - x, y1 - y);
    pdf_contents_show_text(canvas, "(x1, y1)");
    pdf_contents_move_text_pos(canvas, x2 - x1, y2 - y1);
    pdf_contents_show_text(canvas, "(x2, y2)");
    pdf_contents_end_text(canvas);

    pdf_contents_set_dash(canvas, 3, 0, 0);
    pdf_contents_set_line_width(canvas, 0.5);
    pdf_contents_move_to(canvas, x1, y1);
    pdf_contents_line_to(canvas, x2, y2);
    pdf_contents_stroke(canvas);
        
    pdf_contents_set_dash(canvas, 0, 0, 0);
    pdf_contents_set_line_width(canvas, 1.5);
        
    pdf_contents_move_to(canvas, x, y);
    pdf_contents_curve_to2(canvas, x1, y1, x2, y2);
    pdf_contents_stroke(canvas);
    
    /* Curve Example(curve_to3) */
    y -= 150;
    y1 -= 150;
    y2 -= 150;
        
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, 300, 390);
    pdf_contents_show_text(canvas, "curve_to3(x1, y1, x2. y2)");
    pdf_contents_end_text(canvas);
        
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, x + 5, y - 5);
    pdf_contents_show_text(canvas, "Current point");
    pdf_contents_move_text_pos(canvas, x1 - x, y1 - y);
    pdf_contents_show_text(canvas, "(x1, y1)");
    pdf_contents_move_text_pos(canvas, x2 - x1, y2 - y1);
    pdf_contents_show_text(canvas, "(x2, y2)");
    pdf_contents_end_text(canvas);

    pdf_contents_set_dash(canvas, 3, 0, 0);
    pdf_contents_set_line_width(canvas, 0.5);
    pdf_contents_move_to(canvas, x, y);
    pdf_contents_line_to(canvas, x1, y1);
    pdf_contents_stroke(canvas);
        
    pdf_contents_set_dash(canvas, 0, 0, 0);
    pdf_contents_set_line_width(canvas, 1.5);
    pdf_contents_move_to(canvas, x, y);
    pdf_contents_curve_to3(canvas, x1, y1, x2, y2);
    pdf_contents_stroke(canvas);

    /* Curve Example(curve_to) */
    y -= 150;
    y1 -= 160;
    y2 -= 130;
    x2 += 10;
    x3 = 480;
    y3 = 90;

    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, 300, 240);
    pdf_contents_show_text(canvas, "curve_to(x1, y1, x2. y2, x3, y3)");
    pdf_contents_end_text(canvas);
        
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, x + 5, y - 5);
    pdf_contents_show_text(canvas, "Current point");
    pdf_contents_move_text_pos(canvas, x1 - x, y1 - y);
    pdf_contents_show_text(canvas, "(x1, y1)");
    pdf_contents_move_text_pos(canvas, x2 - x1, y2 - y1);
    pdf_contents_show_text(canvas, "(x2, y2)");
    pdf_contents_move_text_pos(canvas, x3 - x2, y3 - y2);
    pdf_contents_show_text(canvas, "(x3, y3)");
    pdf_contents_end_text(canvas);

    pdf_contents_set_dash(canvas, 3, 0, 0);
    pdf_contents_set_line_width(canvas, 0.5);
    pdf_contents_move_to(canvas, x, y);
    pdf_contents_line_to(canvas, x1, y1);
    pdf_contents_stroke(canvas);
    pdf_contents_move_to(canvas, x2, y2);
    pdf_contents_line_to(canvas, x3, y3);
    pdf_contents_stroke(canvas);

    pdf_contents_set_dash(canvas, 0, 0, 0);
    pdf_contents_set_line_width(canvas, 1.5);
    pdf_contents_move_to(canvas, x, y);
    pdf_contents_curve_to(canvas, x1, y1, x2, y2, x3, y3);
    pdf_contents_stroke(canvas);
        
    /* Save the document to a file */
    pdf_doc_write_to_file(doc, "line_example.pdf");
    pdf_doc_free(doc);
    return 0;
}
