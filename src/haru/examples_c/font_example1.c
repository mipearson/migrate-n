/*
 * << H a r u --free pdf library >> -- font_example.c
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
#include <stdlib.h>
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

int add_fonts(pdf_doc doc);

int main()
{
    const char *page_title = "Font Example1";
    int ret;
    pdf_doc doc = pdf_doc_new();
    pdf_page page;
    pdf_contents canvas;
    float w;
    float page_height;
    float page_width;
    int i;

    /* Create a new document object. */
    ret = pdf_doc_new_doc(doc);
    CHECK_ERROR(ret);

    /* Add predefined fonts */
    ret = add_fonts(doc);
    CHECK_ERROR(ret);

    /* Add new page. */
    page = pdf_doc_add_page(doc);
    CHECK_ERROR2(page);

    /* Get the contents object. */
    canvas = pdf_page_get_canvas(page);
    CHECK_ERROR2(canvas);

    ret = pdf_stream_add_filter(canvas, PDF_FILTER_DEFLATE);
    CHECK_ERROR(ret);
    
    page_height = pdf_contents_get_height(canvas);
    page_width = pdf_contents_get_width(canvas);

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

    /* Print the subtitle. */
    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, page_height - 80);
    CHECK_ERROR(ret);

    ret = pdf_contents_set_font_and_size(canvas, "Helvetica", 16);
    CHECK_ERROR(ret);

    ret = pdf_contents_show_text(canvas, "<Standerd Type1 fonts samples>");
    CHECK_ERROR(ret);

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR(ret);
    
    ret = pdf_contents_begin_text(canvas);
    CHECK_ERROR(ret);

    ret = pdf_contents_move_text_pos(canvas, 60, page_height - 105);
    CHECK_ERROR(ret);

    /* Print the sample text with using each predefined fonts. */
    for (i = 0; i < 14; i++) {
        const char* samp_text = "abcdefgABCDEFG12345!#$%&+-@?";
        const char *font_name = pdf_doc_get_font_name(doc, i);
    
        if (font_name == NULL) 
            break;

        /* Print the description label. */
        ret = pdf_contents_set_font_and_size(canvas, "Helvetica", 9);
        CHECK_ERROR(ret);

        ret = pdf_contents_show_text(canvas, font_name);
        CHECK_ERROR(ret);

        ret = pdf_contents_move_text_pos(canvas, 0, -18);
        CHECK_ERROR(ret);
        
        /* Prit the sample text. */
        ret = pdf_contents_set_font_and_size(canvas, font_name, 20);
        CHECK_ERROR(ret);

        ret = pdf_contents_show_text(canvas, samp_text);
        CHECK_ERROR(ret);

        ret = pdf_contents_move_text_pos(canvas, 0, -20);
        CHECK_ERROR(ret);
    }

    ret = pdf_contents_end_text(canvas);
    CHECK_ERROR(ret);

    /* Save the document as PDF file. */
    ret = pdf_doc_write_to_file(doc, "font_example1.pdf");
    CHECK_ERROR(ret);

    /* Cleanup. */
    pdf_doc_free(doc);
    return 0;
}

int add_fonts(pdf_doc doc)
{
    int ret;
    pdf_type1_fontdef fd = NULL;

    fd = pdf_create_type1_fontdef(PDF_FONT_HELVETICA);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);

    fd = pdf_create_type1_fontdef(PDF_FONT_HELVETICA_BOLD);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);

    fd = pdf_create_type1_fontdef(PDF_FONT_HELVETICA_OBLIQUE);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_HELVETICA_BOLD_OBLIQUE);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_TIMES_ROMAN);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_TIMES_BOLD);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_TIMES_ITALIC);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_TIMES_BOLD_ITALIC);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_COURIRE);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_COURIRE_BOLD);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_COURIRE_OBLIQUE);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_COURIRE_BOLD_OBLIQUE);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_SYMBOL);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, 
            pdf_create_encodingdef(PDF_EN_SYMBOL_FONT_ENCODING));
    CHECK_ERROR(ret);
    
    fd = pdf_create_type1_fontdef(PDF_FONT_ZAP_DINGBATS);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, 
            pdf_create_encodingdef(PDF_EN_ZAP_DINGBATS_FONT_ENCODING));
    CHECK_ERROR(ret);
    
    return PDF_SUCCESS;
}

