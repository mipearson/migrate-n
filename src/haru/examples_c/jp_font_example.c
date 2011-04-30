/*
 * << H a r u -- Free PDF Library >> -- jp_font_example1.cc
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
 *   Japanese Font example for SJIS encoding.
 */

#include "libharuc.h"
#include "libharuc_jpfonts.h"
#include "errno.h"

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

int main ()
{
    char samp_text[2048];
    const char* page_title = "JPFont Example";
    pdf_doc doc = pdf_doc_new();
    pdf_type1_fontdef fd = NULL;
    int ret;
    double w;
    unsigned int i;
    pdf_cid_type2_fontdef jpfont[16];
    pdf_cmap map_sjis[2];
    pdf_page page;
    pdf_contents canvas;
	FILE* f;
    
    /* Start creating PDF document. */
    ret = pdf_doc_new_doc(doc);
    CHECK_ERROR(ret);

    fd = pdf_create_type1_fontdef(PDF_FONT_HELVETICA);
    CHECK_ERROR2(fd);
    
    ret = pdf_doc_add_type1font(doc, fd, NULL, NULL);
    CHECK_ERROR(ret);

    /* Create CMap objects. */

    /* CMap for sjis encording (fix width font). */
    map_sjis[0] = pdf_create_jp_cmap(PDF_CMAP_90MS_RKSJ_H);
    CHECK_ERROR2(map_sjis[0]);

    /* Once the CMap object register to PdfDocObject by AddType0Font
     * method, the CMap object will be disposed automatically. But 
     * if AddType0Font failed, the CMap object may never be disposed.
     * So It is recommended to call pdf_doc_register_object method to register
     * CMap object to PdfDoc object. */
    pdf_doc_register_object(doc, map_sjis[0]);

    jpfont[0] = pdf_create_jp_fontdef(PDF_FONT_GOTHIC);
    jpfont[1] = pdf_create_jp_fontdef(PDF_FONT_GOTHIC_BOLD);
    jpfont[2] = pdf_create_jp_fontdef(PDF_FONT_GOTHIC_ITALIC);
    jpfont[3] = pdf_create_jp_fontdef(PDF_FONT_GOTHIC_BOLD_ITALIC);
    jpfont[4] = pdf_create_jp_fontdef(PDF_FONT_MINCYO);
    jpfont[5] = pdf_create_jp_fontdef(PDF_FONT_MINCYO_BOLD);
    jpfont[6] = pdf_create_jp_fontdef(PDF_FONT_MINCYO_ITALIC);
    jpfont[7] = pdf_create_jp_fontdef(PDF_FONT_MINCYO_BOLD_ITALIC);

    /* Gothic */
    pdf_doc_add_type0font(doc, jpfont[0], "gothic", map_sjis[0]);
    pdf_doc_add_type0font(doc, jpfont[1], "gothic-bold", map_sjis[0]);
    pdf_doc_add_type0font(doc, jpfont[2], "gothic-italic", map_sjis[0]);
    pdf_doc_add_type0font(doc, jpfont[3], "gothic-bolditalic", map_sjis[0]);

    /* Mincyo */
    pdf_doc_add_type0font(doc, jpfont[4], "mincyo", map_sjis[0]);
    pdf_doc_add_type0font(doc, jpfont[5], "mincyo-bold", map_sjis[0]);
    pdf_doc_add_type0font(doc, jpfont[6], "mincyo-italic", map_sjis[0]);
    pdf_doc_add_type0font(doc, jpfont[7], "mincyo-bolditalic", map_sjis[0]);
    
    /* CMap for sjis encoding (proportional font) */
    map_sjis[1] = pdf_create_jp_cmap(PDF_CMAP_90MSP_RKSJ_H);
    pdf_doc_register_object(doc, map_sjis[1]);
    
    jpfont[8] = pdf_create_jp_fontdef(PDF_FONT_PGOTHIC);
    jpfont[9] = pdf_create_jp_fontdef(PDF_FONT_PGOTHIC_BOLD);
    jpfont[10] = pdf_create_jp_fontdef(PDF_FONT_PGOTHIC_ITALIC);
    jpfont[11] = pdf_create_jp_fontdef(PDF_FONT_PGOTHIC_BOLD_ITALIC);
    jpfont[12] = pdf_create_jp_fontdef(PDF_FONT_PMINCYO);
    jpfont[13] = pdf_create_jp_fontdef(PDF_FONT_PMINCYO_BOLD);
    jpfont[14] = pdf_create_jp_fontdef(PDF_FONT_PMINCYO_ITALIC);
    jpfont[15] = pdf_create_jp_fontdef(PDF_FONT_PMINCYO_BOLD_ITALIC);

    /* Proportional Gothic */
    pdf_doc_add_type0font(doc, jpfont[8], "pgothic", map_sjis[1]);
    pdf_doc_add_type0font(doc, jpfont[9], "pgothic-bold", map_sjis[1]);
    pdf_doc_add_type0font(doc, jpfont[10], "pgothic-italic", map_sjis[1]);
    pdf_doc_add_type0font(doc, jpfont[11], "pgothic-bolditalic", map_sjis[1]);

    /* Proportional Mincyo */
    pdf_doc_add_type0font(doc, jpfont[12], "pmincyo", map_sjis[1]);
    pdf_doc_add_type0font(doc, jpfont[13], "pmincyo-bold", map_sjis[1]);
    pdf_doc_add_type0font(doc, jpfont[14], "pmincyo-italic", map_sjis[1]);
    pdf_doc_add_type0font(doc, jpfont[15], "pmincyo-bolditalic", map_sjis[1]);
    
    /* Add a new page object. */
    page = pdf_doc_add_page(doc);
    canvas = pdf_page_get_canvas(page);
        
    /* Print the lines of the page. */
    pdf_contents_set_line_width(canvas, 1);
    pdf_contents_rectangle(canvas, 50, 50, 
            pdf_contents_get_width(canvas) - 100, 
            pdf_contents_get_height(canvas) - 110);

    pdf_contents_stroke(canvas);

    /* Print the title of the page (with positioning center). */
    pdf_contents_set_font_and_size(canvas, "Helvetica", 24);
    w = pdf_contents_get_text_width(canvas, page_title, NULL, NULL);
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, (pdf_contents_get_width(canvas) - w) / 2,
            pdf_contents_get_height(canvas) - 50);
    pdf_contents_show_text(canvas, page_title);
    pdf_contents_end_text(canvas);

    /* Output subtitle. */
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, 60, 
            pdf_contents_get_height(canvas) - 80);
    pdf_contents_set_font_and_size(canvas, "Helvetica", 16);
    pdf_contents_show_text(canvas, "<Japanese fonts samples>");
    pdf_contents_end_text(canvas);
       
    /* Load sample text. */
    f = fopen("../examples/mbtext/sjis.txt", "rb");
    if (f == NULL) {
        perror("Cannot open [../examples/mbtext/sjis.txt].");
        return errno;
    }
    if (fgets(samp_text, 2048, f) == NULL) {
        perror("Cannot read from [../examples/mbtext/sjis.txt].");
        return errno;
    }

    /* Print the text in order using all fonts. */
    pdf_contents_begin_text(canvas);
    pdf_contents_move_text_pos(canvas, 60, 
            pdf_contents_get_height(canvas) - 105);
    for (i = 1; i < pdf_doc_count_fonts(doc); i++) {
        const char* fname = pdf_doc_get_font_name(doc, i);
        
        /* print a label of text */ 
        pdf_contents_set_font_and_size(canvas, "Helvetica", 9);
        pdf_contents_show_text(canvas, fname);
        pdf_contents_move_text_pos(canvas, 0, -18);
            
        /* print a sample text. */
        pdf_contents_set_font_and_size(canvas, fname, 20);
        pdf_contents_show_text(canvas, samp_text);
        pdf_contents_move_text_pos(canvas, 0, -20);
    }
    pdf_contents_end_text(canvas);
        
    /* Save the document to a file */
    pdf_doc_write_to_file(doc, "jp_font_example1.pdf");
    
    pdf_doc_free_doc(doc);
    pdf_doc_free(doc);
    return 0;
}

