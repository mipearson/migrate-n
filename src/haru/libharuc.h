/*
 * << H a r u -- Free PDF Library >> -- libharuc.h
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

#ifndef _LIB_HARUC_H 
#define _LIB_HARUC_H

#include "libharu.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
 
/*----------------------------------------------------------------------------*/
/*----- "C" interfaces -------------------------------------------------------*/

typedef void *pdf_object;
typedef void *pdf_auto_ptr_object;
typedef void *pdf_page_label;
typedef pdf_object pdf_name;
typedef pdf_object pdf_dictionary;
typedef pdf_object pdf_stream;
typedef pdf_dictionary pdf_font;
typedef pdf_font pdf_type1font;
typedef pdf_font pdf_type0font;
typedef pdf_font pdf_cid_type2_font;
typedef pdf_dictionary pdf_xobject;
typedef pdf_dictionary pdf_page_base;
typedef pdf_dictionary pdf_link_anot;
typedef pdf_page_base pdf_pages;
typedef pdf_page_base pdf_page;
typedef pdf_xobject pdf_image;
typedef pdf_image pdf_png_image;
typedef pdf_image pdf_jpeg_image;
typedef pdf_auto_ptr_object pdf_encodingdef;
typedef pdf_auto_ptr_object pdf_type1_fontdef;
typedef pdf_auto_ptr_object pdf_cmap;
typedef pdf_auto_ptr_object pdf_cid_fontdef;
typedef pdf_cid_fontdef pdf_cid_type2_fontdef;
typedef void *pdf_doc;

#define PDF_SUCCESS 0
#define PDF_FAILURE 1

#define PDF_TRUE    1
#define PDF_FALSE   0

/*---------------------------------------------------------------------------*/
/*----- PdfStream class -----------------------------------------------------*/

int pdf_stream_add_filter(pdf_stream stream, pdf_filter filter);

/*---------------------------------------------------------------------------*/
/*----- PdfType1FontDef class -----------------------------------------------*/

pdf_type1_fontdef pdf_type1_fontdef_new(const char* basefont);

int pdf_type1_fontdef_load_from_file(pdf_type1_fontdef fontdef, 
		const char* afmfile, const char* fontfile);

int pdf_type1_fontdef_clear(pdf_type1_fontdef fontdef);

enum pdf_base_14_font_enum {
    PDF_FONT_HELVETICA = 0,
    PDF_FONT_HELVETICA_BOLD,
    PDF_FONT_HELVETICA_OBLIQUE,
    PDF_FONT_HELVETICA_BOLD_OBLIQUE,
    PDF_FONT_TIMES_ROMAN,
    PDF_FONT_TIMES_BOLD,
    PDF_FONT_TIMES_ITALIC,
    PDF_FONT_TIMES_BOLD_ITALIC,
    PDF_FONT_COURIRE,
    PDF_FONT_COURIRE_BOLD,
    PDF_FONT_COURIRE_OBLIQUE,
    PDF_FONT_COURIRE_BOLD_OBLIQUE,
    PDF_FONT_SYMBOL,
    PDF_FONT_ZAP_DINGBATS
};
typedef enum pdf_base_14_font_enum pdf_base_14_font;
	
pdf_type1_fontdef pdf_create_type1_fontdef(pdf_base_14_font font);

/*---------------------------------------------------------------------------*/
/*----- PdfEncodingDef class ------------------------------------------------*/

enum pdf_predefined_encoding_enum {
    PDF_EN_STANDARD_ENCODING,
    PDF_EN_MAC_ROMAN_ENCODING,
    PDF_EN_WIN_ANSI_ENCODING,
    PDF_EN_SYMBOL_FONT_ENCODING,
    PDF_EN_ZAP_DINGBATS_FONT_ENCODING
};
typedef enum pdf_predefined_encoding_enum pdf_predefined_encoding;

pdf_encodingdef pdf_create_encodingdef(pdf_predefined_encoding encoding);

/*---------------------------------------------------------------------------*/
/*----- PdfCMap class -------------------------------------------------------*/

unsigned int pdf_cmap_get_cid(pdf_cmap map, unsigned int code);

/*---------------------------------------------------------------------------*/
/*----- PdfCIDType2FontDef class --------------------------------------------*/

unsigned int pdf_cid_type2_fontdef_cid_width(pdf_cid_type2_font font, 
        unsigned int cid);

const char* pdf_cid_type2_fontdef_basefont(pdf_cid_type2_font font);

int pdf_cid_type2_fontdef_ascent(pdf_cid_type2_font font);

int pdf_cid_type2_fontdef_descent(pdf_cid_type2_font font);

int pdf_cid_type2_fontdef_capheight(pdf_cid_type2_font font);

unsigned int pdf_cid_type2_fontdef_flags(pdf_cid_type2_font font);

pdf_box pdf_cid_type2_fontdef_font_bbox(pdf_cid_type2_font font);

int pdf_cid_type2_fontdef_italic_angle(pdf_cid_type2_font font);

unsigned int pdf_cid_type2_fontdef_stemv(pdf_cid_type2_font font);

int pdf_cid_type2_fontdef_dw(pdf_cid_type2_font font);

int pdf_cid_type2_fontdef_missing_width(pdf_cid_type2_font font);

/*---------------------------------------------------------------------------*/
/*----- PdfCatalog class ----------------------------------------------------*/

typedef pdf_dictionary pdf_catalog;
typedef pdf_dictionary pdf_outline_root;
typedef pdf_dictionary pdf_destination;

pdf_destination pdf_catalog_get_openaction(pdf_catalog catalog);

pdf_outline_root pdf_catalog_get_outlines(pdf_catalog catalog);

pdf_page_layout pdf_catalog_get_page_layout(pdf_catalog catalog);

pdf_page_mode pdf_catalog_get_page_mode(pdf_catalog catalog);

pdf_page_mode pdf_catalog_get_nonfullscreen_page_mode(pdf_catalog catalog);

int pdf_catalog_get_viewer_preferences(pdf_catalog catalog);

int pdf_catalog_set_page_layout(pdf_catalog catalog, pdf_page_layout layout);

int pdf_catalog_set_page_mode(pdf_catalog catalog, pdf_page_mode mode);

int pdf_catalog_set_nonfullscreen_page_mode(pdf_catalog catalog, 
		pdf_page_mode mode);

int pdf_catalog_set_viewer_preferences(pdf_catalog catalog, 
		int viewer_preferences);

int pdf_catalog_set_openaction(pdf_catalog catalog, pdf_destination action);

int pdf_catalog_add_page_label(pdf_catalog catalog, unsigned int pagenum, 
		pdf_page_num_style style, unsigned int firstpage, const char* prefix);

int pdf_catalog_clear_page_label(pdf_catalog catalog);

/*---------------------------------------------------------------------------*/
/*----- PdfPages class ------------------------------------------------------*/

pdf_box pdf_pages_get_media_box(pdf_pages pages);

pdf_box pdf_pages_get_crop_box(pdf_pages pages);

pdf_dictionary pdf_pages_get_resources(pdf_pages pages);

int pdf_pages_get_rotate(pdf_pages pages);

  /* BEERLI addition for missing set_rotate function*/
  void pdf_pages_set_rotate(pdf_pages pages, int degree);
  void pdf_page_set_rotate(pdf_page pages, int degree);


pdf_pages pdf_pages_get_parent(pdf_pages pages);

pdf_dictionary pdf_pages_get_resource(pdf_pages pages, 
		const char *element_name);

int pdf_pages_set_media_box(pdf_pages pages, pdf_box rect);

int pdf_pages_set_crop_box(pdf_pages pages, pdf_box rect);

int pdf_pages_set_resources(pdf_page page, pdf_dictionary resources);

int pdf_pages_set_size(pdf_pages pages, int width, int height);

pdf_pages pdf_pages_add_pages(pdf_pages pages, int index);

pdf_page pdf_pages_add_page(pdf_pages pages, int index);

/*---------------------------------------------------------------------------*/
/*----- PdfPage class -------------------------------------------------------*/

typedef pdf_stream pdf_contents;

pdf_box pdf_page_get_media_box(pdf_page page);

pdf_box pdf_page_get_crop_box(pdf_page page);

pdf_dictionary pdf_page_get_resources(pdf_page page);

int pdf_page_get_rotate(pdf_page page);

pdf_pages pdf_page_get_parent(pdf_page page);

pdf_dictionary pdf_page_get_resource(pdf_page page, const char *element_name);

int pdf_page_set_media_box(pdf_page page, pdf_box rect);

int pdf_page_set_crop_box(pdf_page page, pdf_box rect);

int pdf_page_set_resources(pdf_page page, pdf_dictionary resources);

int pdf_page_set_size(pdf_page page, int width, int height);

pdf_contents pdf_page_get_canvas(pdf_page page);

int pdf_page_get_width(pdf_page page);

int pdf_page_get_height(pdf_page page);

pdf_link_anot pdf_page_add_link(pdf_page page, pdf_rect rect, 
		pdf_destination dest, pdf_annot_hl_mode mode);

/*---------------------------------------------------------------------------*/
/*----- PdfContents class ---------------------------------------------------*/

pdf_page pdf_contents_get_page(pdf_contents contents);

int pdf_contents_get_width(pdf_contents contents);

int pdf_contents_get_height(pdf_contents contents);

pdf_point pdf_contents_get_cur_point(pdf_contents contents);

pdf_point pdf_contents_get_text_pos(pdf_contents contents);

double pdf_contents_get_test_leading(pdf_contents contents);

double pdf_contents_get_font_size(pdf_contents contents);

const char *pdf_contents_get_font_name(pdf_contents contents);

double pdf_contents_get_char_space(pdf_contents contents);

double pdf_contents_get_word_space(pdf_contents contents);

double pdf_contents_get_text_width(pdf_contents contents, const char *text,
		int *numchars, int *numwords);

int pdf_contents_get_char_widths(pdf_contents contents, const char*text,
		double *widths);

unsigned int pdf_contents_measure_text(pdf_contents contents, 
		const char *text, double width, double *realwidth);

pdf_rgb_color pdf_contents_get_rgb_fill(pdf_contents contents);

pdf_rgb_color pdf_contents_get_rgb_stroke(pdf_contents contents);

double pdf_contents_get_gray_fill(pdf_contents contents);

double pdf_contents_get_gray_stroke(pdf_contents contents);

/*--- General graphics state ---*/

int pdf_contents_set_line_width(pdf_contents contents, double linewidth);

int pdf_contents_set_line_cap(pdf_contents contents, 
		pdf_line_cap_style linecap);

int pdf_contents_set_line_join(pdf_contents contents,
		pdf_line_join_style linejoin);

int pdf_sontents_set_miter_limit(pdf_contents contents, double miterlimit);

int pdf_contents_set_dash(pdf_contents contents, unsigned int on, 
		unsigned int off, unsigned int phase);

int pdf_contents_set_flat(pdf_contents contents, unsigned int flatness);

/*--- Special graphic state ---*/

int pdf_contents_gsave(pdf_contents contents);

int pdf_contents_grestore(pdf_contents contents);

int pdf_contents_concat(pdf_contents contents, double a, double b, double c,
				double d, double e, double f);

/*--- Path construction ---*/

int pdf_contents_move_to(pdf_contents contents, double x, double y);

int pdf_contents_line_to(pdf_contents contents, double x, double y);

int pdf_contents_curve_to(pdf_contents contents, double x1, double y1,
				double x2, double y2, double x3, double y3);

int pdf_contents_curve_to2(pdf_contents contents, double x2, double y2,
				double x3, double y3);

int pdf_contents_curve_to3(pdf_contents contents, double x1, double y1,
				double x3, double y3);

int pdf_contents_close_path(pdf_contents contents);

int pdf_contents_rectangle(pdf_contents contents, double x, double y,
				double width, double height);

/*--- Path painting ---*/

int pdf_contents_stroke(pdf_contents contents);

int pdf_contents_close_path_stroke(pdf_contents contents);

int pdf_contents_fill(pdf_contents contents);

int pdf_contents_eofill(pdf_contents contents);

int pdf_contents_fill_stroke(pdf_contents contents);

int pdf_contents_eofill_stroke(pdf_contents contents);

int pdf_contents_close_path_fill_stroke(pdf_contents contents);

int pdf_contents_close_path_eofill_stroke(pdf_contents contents);

int pdf_contents_end_path(pdf_contents contents);

/*--- Clipping path ---*/

int pdf_contents_clip(pdf_contents contents);

int pdf_contents_eoclip(pdf_contents contents);

/*--- Text object ---*/

int pdf_contents_begin_text(pdf_contents contents);

int pdf_contents_end_text(pdf_contents contents);

/*--- Text state ---*/

int pdf_contents_set_char_space(pdf_contents contents, double value);

int pdf_contents_set_word_space(pdf_contents contents, double value);

int pdf_contents_set_horizontal_scalling(pdf_contents contents, double value);

int pdf_contents_set_text_leading(pdf_contents contents, double value);

int pdf_contents_set_font_and_size(pdf_contents contents,
				const char *fontname, double size);

int pdf_contents_set_text_rendering_mode(pdf_contents contents,
				pdf_text_rendering_mode mode);

int pdf_contents_set_text_raise(pdf_contents contents, double value);

/*--- Text positioning ---*/

int pdf_contents_move_text_pos(pdf_contents contents, double x, double y);

int pdf_contents_move_text_pos2(pdf_contents contents, double x, double y);

int pdf_contents_set_text_matrix(pdf_contents contents, double a, double b, 
		double c, double d, double x, double y);

int pdf_contents_move_to_next_line(pdf_contents contents);

/*--- Text showing ---*/

int pdf_contents_show_text(pdf_contents contents, const char *text);

int pdf_contents_show_text_next_line(pdf_contents contents, const char *text);

/*--- color showing ---*/

int pdf_contents_set_gray_fill(pdf_contents contents, double gray);

int pdf_contents_set_gray_stroke(pdf_contents contents, double gray);

int pdf_contents_set_rgb_fill(pdf_contents contents, pdf_rgb_color c);

int pdf_contents_set_rgb_stroke(pdf_contents contents, pdf_rgb_color c);

int pdf_contents_set_cmyk_fill(pdf_contents contents, double c, double m, 
				double y, double k);

int pdf_contents_set_cmyk_stroke(pdf_contents contents, double c, double m, 
				double y, double k);

/*--- XObjects ---*/

int pdf_contents_execute_xobject(pdf_contents contents, const char *name);

int pdf_contents_execute_xobject2(pdf_contents contents,
				pdf_xobject xobject);

/*---------------------------------------------------------------------------*/
/*----- PdfInfo class -------------------------------------------------------*/

typedef pdf_dictionary pdf_info;

const char* pdf_info_get_author(pdf_info info);

const char* pdf_info_get_creator(pdf_info info);

const char* pdf_info_get_producer(pdf_info info);

const char* pdf_info_get_title(pdf_info info);

const char* pdf_info_get_subject(pdf_info info);

const char* pdf_info_get_keywords(pdf_info info);

int pdf_info_set_author(pdf_info info, const char* value);

int pdf_info_set_creator(pdf_info info, const char* value);

int pdf_info_set_producer(pdf_info info, const char* value);

int pdf_info_set_title(pdf_info info, const char* value);

int pdf_info_set_subject(pdf_info info, const char* value);

int pdf_info_set_keywords(pdf_info info, const char* value);

int pdf_info_set_author_mb(pdf_info info, const char* value, pdf_cmap map);

int pdf_info_set_creator_mb(pdf_info info, const char* value, pdf_cmap map);

int pdf_info_set_producer_mb(pdf_info info, const char* value, pdf_cmap map);

int pdf_info_set_title_mb(pdf_info info, const char* value, pdf_cmap map);

int pdf_info_set_subject_mb(pdf_info info, const char* value, pdf_cmap map);

int pdf_info_set_keywords_mb(pdf_info info, const char* value, pdf_cmap map);

/*----------------------------------------------------------------------------*/
/*----- PdfOutline class -----------------------------------------------------*/

typedef pdf_dictionary pdf_outline;
typedef pdf_outline pdf_outline_item;

int pdf_outline_has_child(pdf_outline outline);

pdf_outline pdf_outline_first(pdf_outline outline);

pdf_outline pdf_outline_last(pdf_outline outline);

int pdf_outline_has_child(pdf_outline outline);

const char *pdf_outline_get_title(pdf_outline_item item);

int pdf_outline_set_title(pdf_outline_item item, const char *title);

int pdf_outline_set_title_mb(pdf_outline_item item, const char *title,
        pdf_cmap map);

/*----------------------------------------------------------------------------*/
/*----- PdfImage class -------------------------------------------------------*/

double pdf_image_get_height(pdf_image image);

double pdf_image_get_width(pdf_image image);

/*----------------------------------------------------------------------------*/
/*----- PdfPngImage class ----------------------------------------------------*/

pdf_png_image pdf_png_image_new(pdf_doc doc);

int pdf_png_image_load_from_file(pdf_png_image image, char *filename);

void pdf_png_image_free_image(pdf_png_image image);

/*----------------------------------------------------------------------------*/
/*----- PdfJpegPngImage class ------------------------------------------------*/

pdf_jpeg_image pdf_jpeg_image_new(pdf_doc doc);

int pdf_jpeg_image_load_from_file(pdf_jpeg_image image, char *filename);

void pdf_jpeg_image_free_image(pdf_jpeg_image image);

/*----------------------------------------------------------------------------*/
/*----- PdfDoc class ---------------------------------------------------------*/

pdf_doc pdf_doc_new(void);

void pdf_doc_free(pdf_doc doc);

int pdf_doc_new_doc(pdf_doc doc);

int pdf_doc_write_to_file(pdf_doc doc, const char *filename);

int pdf_doc_free_doc(pdf_doc doc);

int pdf_doc_add_type1font(pdf_doc doc, pdf_type1_fontdef font, 
		const char *name, pdf_encodingdef encoding);

int pdf_doc_add_type0font(pdf_doc doc, pdf_cid_type2_fontdef font,
        const char *name, pdf_cmap map);

int pdf_doc_add_xobject(pdf_doc doc, pdf_xobject obj, const char *name);

int pdf_doc_register_object(pdf_doc doc, pdf_auto_ptr_object obj);

pdf_info pdf_doc_get_info(pdf_doc doc);

pdf_page pdf_doc_add_page(pdf_doc doc);

pdf_pages pdf_doc_root_pages(pdf_doc doc);

pdf_page pdf_doc_current_page(pdf_doc doc);

int pdf_doc_has_doc(pdf_doc doc);

int pdf_doc_get_error(pdf_doc doc);

int pdf_doc_count_fonts(pdf_doc doc);

const char *pdf_doc_get_font_name(pdf_doc doc, int index);

int pdf_doc_set_password(pdf_doc doc, const char* owner_passwd, 
        const char* user_passwd);

int pdf_doc_set_permission(pdf_doc doc, int permission);

/*----------------------------------------------------------------------------*/
	
#ifdef __cplusplus
}
#endif

#endif /* _LIB_HARU_H */

