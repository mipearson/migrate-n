/*
 * << H a r u --free pdf library >> -- encoding_list.c
 *
 * Copyright (c) 1999-2003 Takeshi Kanno <takeshi_kanno@est.hi-ho.ne.jp>
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.
 * It is provided "as is" without express or implied warranty.
 *
 * 2003.04.20 created.
 *
 */

#include <stdio.h>
#include "libharuc.h"
#include "libharuc_ISO8859.h"

static const int PAGE_WIDTH = 420;
static const int PAGE_HEIGHT = 400;
static const int CELL_WIDTH = 20;
static const int CELL_HEIGHT = 20;
static const int CELL_HEADER = 10;

void draw_fonts(pdf_contents canvas, const char *font_name)
{
    int i, j;

    for (i = 1; i < 17; i++) {
        for (j = 1; j < 17; j++) {
            unsigned char buf[2];

            int y = PAGE_HEIGHT - 55 - ((i - 1) * CELL_HEIGHT);
            int x = j * CELL_WIDTH + 50;

            buf[1] = 0x00;
            buf[0] = (i - 1) * 16 + (j - 1);
            if (buf[0] >= 32) {
                double d;
                
                pdf_contents_begin_text(canvas);
                d = x - pdf_contents_get_text_width(canvas, (char*)buf, 
                        NULL, NULL) / 2;
                pdf_contents_move_text_pos(canvas, d, y);
                pdf_contents_show_text(canvas, (char*)buf);
                pdf_contents_end_text(canvas);
            }
        }
    }
}

void draw_graph(pdf_contents canvas)
{
    char buf[50];
    int i, j;
    
    /* Draw vertical lines. */
    pdf_contents_set_line_width(canvas, 0.5);
    pdf_contents_set_font_and_size(canvas, "helvetica", 14);
    
    for (i = 0; i <= 17; i++) {
        int x = i * CELL_WIDTH + 40;

        pdf_contents_move_to(canvas, x, PAGE_HEIGHT - 60);
        pdf_contents_line_to(canvas, x, 40);
        pdf_contents_stroke(canvas);

        if (i > 0 && i <= 16) {
            pdf_contents_begin_text(canvas);
            pdf_contents_move_text_pos(canvas, x + 5, PAGE_HEIGHT - 75);
            snprintf(buf, 5, "%X", i - 1);
            pdf_contents_show_text(canvas, buf);
            pdf_contents_end_text(canvas);
        }
    }

    /* Draw horizontal lines. */
    for (j = 0; j <= 15; j++) {
        int y = j * CELL_HEIGHT + 40;
        
        pdf_contents_move_to(canvas, 40, y);
        pdf_contents_line_to(canvas, PAGE_WIDTH - 40, y);
        pdf_contents_stroke(canvas);

        if (j < 14) {
            pdf_contents_begin_text(canvas);
            pdf_contents_move_text_pos(canvas, 45, y + 5);
            snprintf(buf, 5, "%X", 15 - j);
            pdf_contents_show_text(canvas, buf);
            pdf_contents_end_text(canvas);
        }
    }
}

void draw_page(pdf_page page, const char *font_name)
{
    pdf_contents canvas = pdf_page_get_canvas(page);
    
    pdf_stream_add_filter(canvas, PDF_FILTER_DEFLATE);
    
    pdf_page_set_size(page, PAGE_WIDTH, PAGE_HEIGHT);
    draw_graph(canvas);

    pdf_contents_begin_text(canvas);
    pdf_contents_set_font_and_size(canvas, "helvetica", 20);
    pdf_contents_move_text_pos(canvas, 40, PAGE_HEIGHT - 50);
    pdf_contents_show_text(canvas, font_name);
    pdf_contents_end_text(canvas);

    pdf_contents_set_font_and_size(canvas, font_name, 14);
    draw_fonts(canvas, font_name);
}

int main(int argc, char** argv)
{
    #ifdef _WIN32
    const char* AFM_PATH = "..\\examples\\type1\\a010013l.afm";
    const char* PFA_PATH = "..\\examples\\type1\\a010013l.pfb";
    #else
    const char* AFM_PATH = "../examples/type1/a010013l.afm";
    const char* PFA_PATH = "../examples/type1/a010013l.pfb";
    #endif
    int errcd = 0;
    pdf_doc doc;
    pdf_type1_fontdef fd;
    int i;
    pdf_page page;
    
    /* Create a PdfDoc object and start making a PDF document. */
    doc = pdf_doc_new();
    if (doc == NULL) {
        fprintf(stderr, "Error: Cannot create PdfDoc Object.\n");
        return -1;
    }   
    
    pdf_doc_new_doc(doc);

    fd = pdf_create_type1_fontdef(PDF_FONT_HELVETICA);
    pdf_doc_add_type1font(doc, fd, "helvetica", NULL);

    fd = pdf_type1_fontdef_new(NULL);
    errcd = pdf_type1_fontdef_load_from_file(fd, AFM_PATH, PFA_PATH);
    if (errcd != PDF_SUCCESS) {
        fprintf(stderr, "ERROR: failed to load courier font(%s,%s).\n",
                AFM_PATH, PFA_PATH);
        return -1;
    }
    
    pdf_doc_add_type1font(doc, fd, "StandardEncoding", NULL);
    pdf_doc_add_type1font(doc, fd, "WinAnsiEncoding", 
            pdf_create_encodingdef(PDF_EN_WIN_ANSI_ENCODING)); 
    pdf_doc_add_type1font(doc, fd, "MacRomanEncoding", 
            pdf_create_encodingdef(PDF_EN_MAC_ROMAN_ENCODING));

    pdf_doc_add_type1font(doc, fd, "iso8859-2 Encoding",
            pdf_create_iso8859_2_encoding());

    /* Add Symbolic font */
    fd = pdf_create_type1_fontdef(PDF_FONT_SYMBOL);
    pdf_doc_add_type1font(doc, fd, "Symbol Set", NULL);

    /* Add ZapfDingbats font */
    fd = pdf_create_type1_fontdef(PDF_FONT_ZAP_DINGBATS);
    pdf_doc_add_type1font(doc, fd, "ZapfDingbats Set", NULL);

    for (i = 1; i < pdf_doc_count_fonts(doc); i++)
    {
        page = pdf_doc_add_page(doc);
        draw_page(page, pdf_doc_get_font_name(doc, i));
    }   

    errcd = pdf_doc_write_to_file(doc, "encoding_list.pdf");
    if (errcd != 0)
        fprintf(stderr, "ERROR %d.\n", errcd);
    pdf_doc_free(doc);

    return 0;
}

