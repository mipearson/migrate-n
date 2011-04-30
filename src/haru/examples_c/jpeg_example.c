/*
 * << H a r u -- Free PDF Library >> -- jpeg_example.c
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

#define CHECK_ERROR(func) \
{ \
    int ret_cd = func; \
    if (ret_cd != PDF_SUCCESS) { \
        fprintf(stderr, "ERROR (%d) line:%d.\n", ret_cd, __LINE__); \
        pdf_doc_free(doc); \
        return ret_cd; \
    } \
}

#define CHECK_ERROR2(func) \
{ \
    void* ret = func; \
    if (ret == NULL) { \
        fprintf(stderr, "ERROR line:%d.\n", __LINE__); \
        pdf_doc_free(doc); \
        return -1; \
    } \
}

int draw_image(pdf_doc doc, pdf_contents canvas, char *filename,
        float x, float y, char *text)
{
#ifdef _WIN32
    const char* FILE_SEPARATOR = "\\";
#else
    const char* FILE_SEPARATOR = "/";
#endif
    char filename1[255];
    pdf_jpeg_image image;

    strcpy(filename1, "../examples/images");
    strcat(filename1, FILE_SEPARATOR);
    strcat(filename1, filename);

    image = pdf_jpeg_image_new(doc);
    CHECK_ERROR2(image);

    CHECK_ERROR(pdf_jpeg_image_load_from_file(image, filename1));
    
    CHECK_ERROR(pdf_doc_add_xobject(doc, image, NULL));

    CHECK_ERROR(pdf_contents_gsave(canvas));

    CHECK_ERROR(pdf_contents_concat(canvas, pdf_image_get_width(image),
                0, 0, pdf_image_get_height(image), x, y));

    CHECK_ERROR(pdf_contents_execute_xobject2(canvas, image));

    CHECK_ERROR(pdf_contents_grestore(canvas));

    CHECK_ERROR(pdf_contents_begin_text(canvas));
    
    CHECK_ERROR(pdf_contents_set_text_leading(canvas, 16));
    
    CHECK_ERROR(pdf_contents_move_text_pos(canvas, x, y));
    
    CHECK_ERROR(pdf_contents_show_text_next_line(canvas, filename));
    
    CHECK_ERROR(pdf_contents_show_text_next_line(canvas, text));
    
    CHECK_ERROR(pdf_contents_end_text(canvas));
    
    return PDF_SUCCESS;
}

int main()
{
    pdf_page page;
    pdf_contents canvas;
    pdf_type1_fontdef fd;
    int ret;
    
    pdf_doc doc = pdf_doc_new();

    if (doc == NULL) {
        fprintf(stderr, "ERROR cannot create pdf_doc object.\n");
        return -1;
    }

    CHECK_ERROR(pdf_doc_new_doc(doc));
    
    page = pdf_doc_add_page(doc);
    CHECK_ERROR2(page);

    canvas = pdf_page_get_canvas(page);
    CHECK_ERROR2(canvas);

    fd = pdf_create_type1_fontdef(PDF_FONT_HELVETICA);
    CHECK_ERROR2(fd);

    CHECK_ERROR(pdf_doc_add_type1font(doc, fd, NULL, NULL));

    CHECK_ERROR(pdf_contents_begin_text(canvas));

    CHECK_ERROR(pdf_contents_set_font_and_size(canvas, "Helvetica", 20));

    CHECK_ERROR(pdf_contents_move_text_pos(canvas, 220, 
                pdf_contents_get_height(canvas) - 70));

    CHECK_ERROR(pdf_contents_show_text(canvas, "PngExample"));

    CHECK_ERROR(pdf_contents_end_text(canvas));

    CHECK_ERROR(pdf_contents_set_font_and_size(canvas, "Helvetica", 12));

    ret = draw_image(doc, canvas, "basn0g01.jpeg", 100, 
            pdf_contents_get_height(canvas) - 150, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn0g02.jpeg", 200, 
            pdf_contents_get_height(canvas) - 150, "2bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn0g04.jpeg", 300, 
            pdf_contents_get_height(canvas) - 150, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn0g08.jpeg", 400, 
            pdf_contents_get_height(canvas) - 150, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn2c08.jpeg", 100, 
            pdf_contents_get_height(canvas) - 250, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn2c16.jpeg", 200, 
            pdf_contents_get_height(canvas) - 250, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn3p01.jpeg", 100, 
            pdf_contents_get_height(canvas) - 350, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn3p02.jpeg", 200, 
            pdf_contents_get_height(canvas) - 350, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn3p04.jpeg", 300, 
            pdf_contents_get_height(canvas) - 350, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn3p08.jpeg", 400, 
            pdf_contents_get_height(canvas) - 350, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn4a08.jpeg", 100, 
            pdf_contents_get_height(canvas) - 450, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn4a16.jpeg", 200, 
            pdf_contents_get_height(canvas) - 450, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn6a08.jpeg", 100, 
            pdf_contents_get_height(canvas) - 550, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    ret = draw_image(doc, canvas, "basn6a16.jpeg", 200, 
            pdf_contents_get_height(canvas) - 550, "1bit grayscale.");
    if (ret != PDF_SUCCESS)
        return -1;

    CHECK_ERROR(pdf_doc_write_to_file(doc, "jpeg_example.pdf"));

    pdf_doc_free(doc);

    return 0;
}

    
