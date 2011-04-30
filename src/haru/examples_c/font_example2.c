/*
 * << H a r u --free pdf library >> -- font_example2.c
 *
 * Copyright (c) 1999-2003 Takeshi Kanno <takeshi_kanno@est.hi-ho.ne.jp>
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 * It is provided "as is" without express or implied warranty.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "libharuc.h"

#define CHECK_ERROR(ret) \
{ \
    if (ret != PDF_SUCCESS) { \
        fprintf(stderr, "ERROR (%d) line:%d.\n", ret, __LINE__); \
        pdf_doc_free(doc); \
        return ret; \
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

#define CHECK_ERROR3(ret) \
{ \
    if (ret != PDF_SUCCESS) { \
        fprintf(stderr, "ERROR (%d) line:%d.\n", ret, __LINE__); \
        return ret; \
    } \
}

static const char* samp_text = "abcdefgABCDEFG123!#$%&@?";

int font_size_demo(pdf_contents canvas);
int font_color_demo(pdf_contents canvas);
int font_rendering_mode_demo(pdf_contents canvas, float ypos);
int text_matrix_demo(pdf_contents canvas, float ypos);
int show_description(pdf_contents canvas, float x, float y, const char *text);
int show_stripe_pattern(pdf_contents canvas, float x, float y);

int main()
{
    const char *page_title = "Font Example2";
    int ret;
    pdf_doc doc = pdf_doc_new();
    pdf_page page;
    pdf_contents canvas;
    float w;
    float page_height;
    float page_width;
    float save_ypos;
    pdf_type1_fontdef fd = NULL;

    /* Start creating PDF document. */
    ret = pdf_doc_new_doc(doc);
    CHECK_ERROR(ret);

    /* Helvetica */
    fd = pdf_create_type1_fontdef(PDF_FONT_HELVETICA);
    CHECK_ERROR2(fd);

    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);

    /* Add a new page object. */
    page = pdf_doc_add_page(doc);
    CHECK_ERROR2(page);

    canvas = pdf_page_get_canvas(page);
    CHECK_ERROR2(page);

    ret = pdf_stream_add_filter(canvas, PDF_FILTER_DEFLATE);
    CHECK_ERROR(ret);
    
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);

    /* Print the lines of the page. */
    ret = pdf_contents_set_line_width(canvas, 1);
    CHECK_ERROR(ret);
    
    ret = pdf_contents_rectangle(canvas, 50, 50, page_width - 100, 
            page_height - 110);
    CHECK_ERROR(ret);

    ret = pdf_contents_stroke(canvas);
    CHECK_ERROR(ret);

    /* Print the title of the page (with positioning center). */       
    ret = pdf_contents_set_font_and_size(canvas, "Helvetica", 24);
    CHECK_ERROR(ret);

    w = pdf_contents_get_text_width(canvas, page_title, NULL, NULL);
    
    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_text_pos(canvas, (page_width - w) / 2,
            page_height - 50);
    CHECK_ERROR(ret);

    ret = pdf_contents_show_text(canvas, page_title);
    CHECK_ERROR(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR(ret);

    /* Output subtitle. */
    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, page_height - 80);
    CHECK_ERROR(ret);

    ret = pdf_contents_set_font_and_size(canvas, "Helvetica", 16);
    CHECK_ERROR(ret);

    ret = pdf_contents_show_text(canvas, "<Various properties of fonts>");
    CHECK_ERROR(ret);

    ret = font_size_demo(canvas);
    CHECK_ERROR(ret);

    ret = font_color_demo(canvas);
    CHECK_ERROR(ret);

    save_ypos = pdf_contents_get_text_pos(canvas).y;
    
    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR(ret);

    ret = font_rendering_mode_demo(canvas, save_ypos);
    CHECK_ERROR(ret);

    ret = text_matrix_demo(canvas, save_ypos);
    CHECK_ERROR(ret);

    ret = pdf_doc_write_to_file(doc, "font_example2.pdf");
    CHECK_ERROR(ret);

    pdf_doc_free(doc);

    return 0;
}

int font_size_demo(pdf_contents canvas)
{
    float fsize = 8;
    int ret;

    float page_width = pdf_contents_get_width(canvas);

    while (fsize < 60) {
        char buf[50];
        unsigned int len;

        /* Set style and size of font. */
        ret = pdf_contents_set_font_and_size(canvas, "Helvetica", fsize);
        CHECK_ERROR3(ret);
        
        /* Set the position of the text. */
        ret = pdf_contents_move_text_pos(canvas, 0, -10 - fsize);

        /* Measure the number of characters which included in the page. */
        strcpy(buf, samp_text);
        len = pdf_contents_measure_text(canvas, buf, page_width - 120, NULL);

        /* Truncate the text. */
        buf[len] = 0x00;

        ret = pdf_contents_show_text(canvas, buf);
        CHECK_ERROR3(ret);
        if (ret != PDF_SUCCESS) return ret;
        
        /* print the description. */
        ret = pdf_contents_move_text_pos(canvas, 0, -10);
        CHECK_ERROR3(ret);
        
        ret = pdf_contents_set_font_and_size(canvas, "Helvetica", 8);
        CHECK_ERROR3(ret);
        
        snprintf(buf, 50, "Fontsize=%.0f", fsize);

        ret = pdf_contents_show_text(canvas, buf);
        CHECK_ERROR3(ret);

        fsize *= 1.5;
    }

    return PDF_SUCCESS;
}

int font_color_demo(pdf_contents canvas)
{
    int ret;
    unsigned int i;
    unsigned int cnt = (unsigned int)strlen(samp_text);
    
    ret = pdf_contents_set_font_and_size(canvas, "Helvetica", 18);
    CHECK_ERROR3(ret);

    ret = pdf_contents_move_text_pos(canvas, 0, -50);
    CHECK_ERROR3(ret);

    for (i = 0; i < cnt; i++) {
        char buf[2];
        float r = (float)i / (float)cnt;
        float g = 1 - ((float)i / (float)cnt);
        buf[0] = samp_text[i];
        buf[1] = 0x00;

        ret = pdf_contents_set_rgb_fill(canvas, PdfRGBColor(r, g, 0));
        CHECK_ERROR3(ret);

        ret = pdf_contents_show_text(canvas, buf);
        CHECK_ERROR3(ret);
    }
    ret = pdf_contents_move_text_pos(canvas, 0, -25);
    CHECK_ERROR3(ret);
    
    for (i = 0; i < cnt; i++) {
        char buf[2];
        float r = (float)i / (float)cnt;
        float b = 1 - ((float)i / (float)cnt);
        buf[0] = samp_text[i];
        buf[1] = 0x00;

        ret = pdf_contents_set_rgb_fill(canvas, PdfRGBColor(r, 0, b));
        CHECK_ERROR3(ret);

        ret = pdf_contents_show_text(canvas, buf);
        CHECK_ERROR3(ret);
    }
    ret =pdf_contents_move_text_pos(canvas, 0, -25);
    CHECK_ERROR3(ret);

    for (i = 0; i < cnt; i++) {
        char buf[2];
        float b = (float)i / (float)cnt;
        float g = 1 - ((float)i / (float)cnt);
        buf[0] = samp_text[i];
        buf[1] = 0x00;

        ret = pdf_contents_set_rgb_fill(canvas, PdfRGBColor(0, g, b));
        CHECK_ERROR3(ret);

        ret = pdf_contents_show_text(canvas, buf);
        CHECK_ERROR3(ret);
    }
    ret = pdf_contents_move_text_pos(canvas, 0, -25);
    CHECK_ERROR3(ret);
    
    return PDF_SUCCESS;
}       

int show_stripe_pattern(pdf_contents canvas, float x, float y)
{
    int ret;
    int iy = 0;
    float tw = pdf_contents_get_text_width(canvas, "ABCabc123", NULL, NULL);
    
    while (iy < 50) {
        ret = pdf_contents_set_rgb_stroke(canvas, PdfRGBColor(0, 0, 0.5));
        CHECK_ERROR3(ret);

        ret = pdf_contents_set_line_width(canvas, 1);
        CHECK_ERROR3(ret);

        ret = pdf_contents_move_to(canvas, x, y + iy);
        CHECK_ERROR3(ret);

        ret = pdf_contents_line_to(canvas, x + tw, y + iy);
        CHECK_ERROR3(ret);

        ret = pdf_contents_stroke(canvas);
        CHECK_ERROR3(ret);

        iy += 3;
    }
    ret = pdf_contents_set_line_width(canvas, 2.5);
    CHECK_ERROR3(ret);

    return PDF_SUCCESS;
}

int show_description(pdf_contents canvas, float x, float y, const char *text)
{
    float fsize = pdf_contents_get_font_size(canvas);
    const char *fname = pdf_contents_get_font_name(canvas);
    pdf_rgb_color c = pdf_contents_get_rgb_fill(canvas);
    int ret;

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_move_text_pos(canvas, x, y - 12);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_rgb_fill(canvas, PdfRGBColor(0, 0, 0));
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_rendering_mode(canvas, PDF_FILL);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_font_and_size(canvas, "Helvetica", 10);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, text);
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_font_and_size(canvas, fname, fsize);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_rgb_fill(canvas, c);
    CHECK_ERROR3(ret);

    return PDF_SUCCESS;
}

int font_rendering_mode_demo(pdf_contents canvas, float ypos)
{
    int ret;

    ret = pdf_contents_set_font_and_size(canvas, "Helvetica", 40);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_rgb_fill(canvas, PdfRGBColor(0.5, 0.5, 0));
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_line_width(canvas, 1.5);
    CHECK_ERROR3(ret);

    /* PDF_FILL */
    ret = show_description(canvas, 60, ypos - 50, "RenderingMode=PDF_FILL");
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_rendering_mode(canvas, PDF_FILL);
    CHECK_ERROR3(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, ypos - 50);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, "ABCabc123");
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);
 
    /* PDF_STROKE */
    ret = show_description(canvas, 60, ypos - 110, "RenderingMode=PDF_STROKE");
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_rendering_mode(canvas, PDF_STROKE);
    CHECK_ERROR3(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, ypos - 110);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, "ABCabc123");
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);

    /* PDF_FILL_THEN_STROKE */
    ret = show_description(canvas, 60, ypos - 170, 
            "RenderingMode=PDF_FILL_THEN_STROKE");
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_rendering_mode(canvas, PDF_FILL_THEN_STROKE);
    CHECK_ERROR3(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, ypos - 170);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, "ABCabc123");
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);

    /* PDF_FILL_CLIPPING */
    ret = show_description(canvas, 60, ypos - 230, 
            "RenderingMode=PDF_FILL_CLIPPING");
    CHECK_ERROR3(ret);

    ret = pdf_contents_gsave(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_rendering_mode(canvas, PDF_FILL_CLIPPING);
    CHECK_ERROR3(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, ypos - 230);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, "ABCabc123");
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);

    ret = show_stripe_pattern(canvas, 60, ypos - 230);
    CHECK_ERROR3(ret);

    ret = pdf_contents_grestore(canvas);
    CHECK_ERROR3(ret);

    /* PDF_STROKE_CLIPPING */
    ret = show_description(canvas, 60, ypos - 290, 
            "RenderingMode=PDF_STROKE_CLIPPING");
    CHECK_ERROR3(ret);

    ret = pdf_contents_gsave(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_rendering_mode(canvas, PDF_STROKE_CLIPPING);
    CHECK_ERROR3(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, ypos - 290);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, "ABCabc123");
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);

    ret = show_stripe_pattern(canvas, 60, ypos - 290);
    CHECK_ERROR3(ret);

    ret = pdf_contents_grestore(canvas);
    CHECK_ERROR3(ret);

    /* PDF_FILL_STROKE_CLIPPING */
    ret = show_description(canvas, 60, ypos - 350,
            "RenderingMode=PDF_FILL_STROKE_CLIPPING");
    CHECK_ERROR3(ret);

    ret = pdf_contents_gsave(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_rendering_mode(canvas, 
            PDF_FILL_STROKE_CLIPPING);
    CHECK_ERROR3(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, ypos - 350);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, "ABCabc123");
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);

    ret = show_stripe_pattern(canvas, 60, ypos - 350);
    CHECK_ERROR3(ret);

    ret = pdf_contents_grestore(canvas);
    CHECK_ERROR3(ret);

    return PDF_SUCCESS;
}

int text_matrix_demo(pdf_contents canvas, float ypos)
{
    int ret;
    float angle = 30;                   /* A rotation of 30 degrees. */
    float rad = angle / 180 * 3.141592; /* Calcurate the radian value. */
    float angle1 = 10;
    float angle2 = 20;
    float rad1 = angle1 / 180 * 3.141592;
    float rad2 = angle2 / 180 * 3.141592;

    /* Reset text attributes */
    ret = pdf_contents_set_text_rendering_mode(canvas, PDF_FILL);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_rgb_fill(canvas, PdfRGBColor(0, 0, 0));
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_font_and_size(canvas, "Helvetica", 30);
    CHECK_ERROR3(ret);

    /* Rotating text */
    ret = show_description(canvas, 300, ypos - 100, "Rotating text");
    CHECK_ERROR3(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_matrix(canvas, cos(rad), sin(rad), -sin(rad), 
            cos(rad), 310, ypos - 100);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, "ABCabc123");
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);
    
    /* Skewing text. */     
    ret = show_description(canvas, 300, ypos - 160, "Skewing text");
    CHECK_ERROR3(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_matrix(canvas, 1, tan(rad1), tan(rad2), 1, 
            300, ypos - 160);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, "ABCabc123");
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);
    
    /* Scaling text (X direction) */
    ret = show_description(canvas, 300, ypos - 220, 
            "Scaling text (X direction)");
    CHECK_ERROR3(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_matrix(canvas, 1.5, 0, 0, 1, 300, ypos - 220);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, "ABCabc12");
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);

    /* Scaling text (Y direction) */
    ret = show_description(canvas, 300, ypos - 300,
            "Scaling text (Y direction)");
    CHECK_ERROR3(ret);

    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR3(ret);

    ret = pdf_contents_set_text_matrix(canvas, 1, 0, 0, 2, 300, ypos - 300);
    CHECK_ERROR3(ret);

    ret = pdf_contents_show_text(canvas, "ABCabc123");
    CHECK_ERROR3(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR3(ret);

    return PDF_SUCCESS;
}


