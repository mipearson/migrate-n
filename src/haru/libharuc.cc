/*
 * << H a r u -- Free PDF Library >> -- libharu.cc
 *
 * Interface routines for "C".
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

#include <new>
#include "libharuc.h"

/*----- PdfStream -----------------------------------------------------------*/

int pdf_stream_add_filter(pdf_stream stream, pdf_filter filter)
{
    if (stream == NULL)
        return PDF_INVALID_PARAMETER;
    else
        ((PdfStream*)stream)->AddFilter(filter);
    return PDF_SUCCESS;
}

/*---------------------------------------------------------------------------*/
/*----- PdfType1FontDef -----------------------------------------------------*/

pdf_type1_fontdef pdf_type1_fontdef_new(const char* basefont)
{
    try {
        if (basefont == NULL)
            return new PdfType1FontDef();
        else
            return new PdfType1FontDef(basefont);
    } catch (...) {
        return NULL;
    }
}

int 
pdf_type1_fontdef_load_from_file(pdf_type1_fontdef fontdef,
        const char* afmfile, const char* fontfile)
{
    if (fontdef == NULL)
        return PDF_INVALID_PARAMETER;
    try {
        ((PdfType1FontDef*)fontdef)->LoadFromFile(afmfile, fontfile);
        return PDF_SUCCESS;
    } catch (PdfException &e) {
        return e.GetCode();
    } catch (...) {
        return PDF_UNKNOWN_ERROR;
    }
}

int 
pdf_type1_fontdef_clear(pdf_type1_fontdef fontdef)
{
    if (fontdef == NULL)
        return PDF_INVALID_PARAMETER;
    try {
        ((PdfType1FontDef*)fontdef)->Clear();
        return PDF_SUCCESS;
    } catch (...) {
        return PDF_UNKNOWN_ERROR;
    }
}

pdf_type1_fontdef
pdf_create_type1_fontdef(pdf_base_14_font font)
{
    try {
        switch (font) {
        case PDF_FONT_HELVETICA:
            return new PdfHelveticaFontDef();
        case PDF_FONT_HELVETICA_BOLD:
            return new PdfHelveticaBoldFontDef();
        case PDF_FONT_HELVETICA_OBLIQUE:
            return new PdfHelveticaObliqueFontDef();
        case PDF_FONT_HELVETICA_BOLD_OBLIQUE:
            return new PdfHelveticaBoldObliqueFontDef();
        case PDF_FONT_TIMES_ROMAN:
            return new PdfTimesRomanFontDef();
        case PDF_FONT_TIMES_BOLD:
            return new PdfTimesBoldFontDef();
        case PDF_FONT_TIMES_ITALIC:
            return new PdfTimesItalicFontDef();
        case PDF_FONT_TIMES_BOLD_ITALIC:
            return new PdfTimesBoldItalicFontDef();
        case PDF_FONT_COURIRE:
            return new PdfCourierFontDef();
        case PDF_FONT_COURIRE_BOLD:
            return new PdfCourierBoldFontDef();
        case PDF_FONT_COURIRE_OBLIQUE:
            return new PdfCourierObliqueFontDef();
        case PDF_FONT_COURIRE_BOLD_OBLIQUE:
            return new PdfCourierBoldObliqueFontDef();
        case PDF_FONT_SYMBOL:
            return new PdfSymbolFontDef();
        case PDF_FONT_ZAP_DINGBATS:
            return new PdfZapfDingbatsFontDef();
        default: return NULL;
        }
    } catch (...) {
        return NULL;
    }
}
        
/*---------------------------------------------------------------------------*/
/*----- PdfEncodingDef ------------------------------------------------------*/

pdf_encodingdef 
pdf_create_encodingdef(pdf_predefined_encoding encoding)
{
    try {
        switch (encoding) {
        case PDF_EN_STANDARD_ENCODING:
            return new PdfStandardEncoding();
        case PDF_EN_MAC_ROMAN_ENCODING:
            return new PdfMacRomanEncoding();
        case PDF_EN_WIN_ANSI_ENCODING:
            return new PdfWinAnsiEncoding();
        default:
            return NULL;
        }
    } catch (...) {
        return NULL;
    }
}

/*---------------------------------------------------------------------------*/
/*----- PdfCMap class -------------------------------------------------------*/

unsigned int 
pdf_cmap_get_cid(pdf_cmap map, unsigned int code)
{
    return (map != NULL) ? ((PdfCMap*)map)->GetCID(code) : 0;
}

/*---------------------------------------------------------------------------*/
/*----- PdfCIDTypeFontDef class ---------------------------------------------*/

unsigned int 
pdf_cid_type2_fontdef_cid_width(pdf_cid_type2_fontdef font, unsigned int cid)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->CIDWidth(cid) : 0;
}

const char* 
pdf_cid_type2_fontdef_basefont(pdf_cid_type2_fontdef font)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->BaseFont() : NULL;
}

int 
pdf_cid_type2_fontdef_ascent(pdf_cid_type2_fontdef font)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->Ascent() : NULL;
}

int 
pdf_cid_type2_fontdef_descent(pdf_cid_type2_fontdef font)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->Descent() : 0;
}

int 
pdf_cid_type2_fontdef_capheight(pdf_cid_type2_fontdef font)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->CapHeight() : 0;
}

unsigned int 
pdf_cid_type2_fontdef_flags(pdf_cid_type2_fontdef font)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->Flags() : 0;
}

pdf_box 
pdf_cid_type2_fontdef_font_bbox(pdf_cid_type2_fontdef font)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->FontBBox() : 
        PdfBox(0, 0, 0, 0);
}

int 
pdf_cid_type2_fontdef_italic_angle(pdf_cid_type2_fontdef font)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->ItalicAngle() : 0;
}

unsigned int 
pdf_cid_type2_fontdef_stemv(pdf_cid_type2_fontdef font)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->StemV() : 0;
}

int 
pdf_cid_type2_fontdef_dw(pdf_cid_type2_fontdef font)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->DW() : 0;
}

int 
pdf_cid_type2_fontdef_missing_width(pdf_cid_type2_fontdef font)
{
    return (font != NULL) ? ((PdfCIDType2FontDef*)font)->MissingWidth() : 0;
}
        
/*---------------------------------------------------------------------------*/
/*----- pdf_catalog ---------------------------------------------------------*/

pdf_destination
pdf_catalog_get_openaction(pdf_catalog catalog)
{
    return (catalog != NULL) ? 
        (pdf_destination)((PdfCatalog*)catalog)->OpenAction() : NULL;
}
                           
pdf_outline_root
pdf_catalog_get_outlines(pdf_catalog catalog)
{
    return (catalog != NULL) ?
        (pdf_outline_root)((PdfCatalog*)catalog)->Outlines() : NULL;
}

pdf_page_layout
pdf_catalog_get_page_layout(pdf_catalog catalog)
{
    return (catalog != NULL) ? 
        ((PdfCatalog*)catalog)->PageLayout() : PDF_SINGLE_PAGE;
}

pdf_page_mode
pdf_catalog_get_page_mode(pdf_catalog catalog)
{
    return (catalog != NULL) ? 
        ((PdfCatalog*)catalog)->PageMode() : PDF_USE_NONE;
}

pdf_page_mode
pdf_catalog_get_nonfullscreen_page_mode(pdf_catalog catalog)
{
    return (catalog != NULL) ? 
        ((PdfCatalog*)catalog)->NonFullScreenPageMode() : PDF_USE_NONE;
}

int
pdf_catalog_get_viewer_preferences(pdf_catalog catalog)
{
    return (catalog != NULL) ?
        ((PdfCatalog*)catalog)->ViewerPreferences() : 0;
}

int
pdf_catalog_set_page_layout(pdf_catalog catalog, pdf_page_layout layout)
{
    if (catalog == NULL)
        return PDF_INVALID_PARAMETER;
    else 
        ((PdfCatalog*)catalog)->SetPageLayout(layout);
    return PDF_SUCCESS;
}   

int
pdf_catalog_set_page_mode(pdf_catalog catalog, pdf_page_mode mode)
{
    if (catalog == NULL)
        return PDF_INVALID_PARAMETER;
    else
        ((PdfCatalog*)catalog)->SetPageMode(mode);
    return PDF_SUCCESS;
}   

int
pdf_catalog_set_nonfullscreen_page_mode(pdf_catalog catalog, 
        pdf_page_mode mode)
{
    if (catalog == NULL)
        return PDF_INVALID_PARAMETER;
    else
        ((PdfCatalog*)catalog)->SetNonFullScreenPageMode(mode);
    return PDF_SUCCESS;
}   

int
pdf_catalog_set_viewer_preferences(pdf_catalog catalog, 
        int viewer_preferences)
{
    if (catalog == NULL)
        return PDF_INVALID_PARAMETER;
    else 
        ((PdfCatalog*)catalog)->SetViewerPreferences(viewer_preferences);
    return PDF_SUCCESS;
}

int
pdf_catalog_set_openaction(pdf_catalog catalog, pdf_destination action)
{
    if (catalog == NULL)
        return PDF_INVALID_PARAMETER;
    else
        ((PdfCatalog*)catalog)->SetOpenAction((PdfDestination*)action);
    return PDF_SUCCESS;
}

int
pdf_catalog_add_page_label(pdf_catalog catalog, unsigned int pagenum, 
        pdf_page_num_style style, unsigned int firstpage, const char* prefix)
{
    if (catalog == NULL)
        return PDF_INVALID_PARAMETER;
    try {
        ((PdfCatalog*)catalog)->AddPageLabel(pagenum, style, firstpage, prefix);
        return PDF_SUCCESS;
    } catch (PdfException &e) {
        ((PdfCatalog*)catalog)->SetError(e.GetCode());
        return e.GetCode();
    } catch (...) {
        ((PdfCatalog*)catalog)->SetError(PDF_UNKNOWN_ERROR);
        return PDF_UNKNOWN_ERROR;
    }
}

int
pdf_catalog_clear_page_label(pdf_catalog catalog)
{
    if (catalog == NULL)
        return PDF_INVALID_PARAMETER;
    try {
        ((PdfCatalog*)catalog)->ClearPageLabel();
        return PDF_SUCCESS;
    } catch (PdfException &e) {
        ((PdfCatalog*)catalog)->SetError(e.GetCode());
        return e.GetCode();
    } catch (...) {
        ((PdfCatalog*)catalog)->SetError(PDF_UNKNOWN_ERROR);
        return PDF_UNKNOWN_ERROR;
    }
}

/*---------------------------------------------------------------------------*/
/*----- pdf_pages -----------------------------------------------------------*/

pdf_box 
pdf_pages_get_media_box(pdf_pages pages)
{
    if (pages == NULL)
        return PdfBox(0, 0, 0, 0);
    return ((PdfPageBase*)pages)->MediaBox();
}

pdf_box 
pdf_pages_get_crop_box(pdf_pages pages)
{
    if (pages == NULL)
        return PdfBox(0, 0, 0, 0);
    return ((PdfPageBase*)pages)->CropBox();
}

pdf_dictionary
pdf_pages_get_resources(pdf_pages pages)
{
    if (pages == NULL)
        return NULL;
    return ((PdfPageBase*)pages)->Resources();
}

int
pdf_pages_get_rotate(pdf_pages pages)
{
    if (pages == NULL)
        return -1;
    return ((PdfPageBase*)pages)->Rotate();
}


// BEERLI addition for missing set_rotate interface

void
pdf_pages_set_rotate(pdf_pages pages, int degree)
{
    if (pages != NULL)
      ((PdfPageBase*)pages)->SetRotate(degree);
}

void
pdf_page_set_rotate(pdf_page page, int degree)
{
  pdf_pages_set_rotate(page, degree);
}


pdf_pages pdf_pages_get_parent(pdf_pages pages)
{
    return (pages != NULL) ? ((PdfPageBase*)pages)->Parent() : NULL;
}

pdf_pages pdf_pages_get_resource(pdf_pages pages, 
        const char *element_name)
{
    if (pages == NULL)
        return NULL;
    return ((PdfPageBase*)pages)->GetResource(element_name);
}

int
pdf_pages_set_media_box(pdf_pages pages, pdf_box rect)
{
    if (pages == NULL)
        return PDF_INVALID_PARAMETER;
    
    try {
        ((PdfPages*)pages)->SetMediaBox(rect);
        return PDF_SUCCESS;
    } catch (PdfException &e) {
        return e.GetCode();
    } catch (...) {
        return PDF_UNKNOWN_ERROR;
    }
}

int
pdf_pages_set_crop_box(pdf_pages pages, pdf_box rect)
{
    if (pages == NULL)
        return PDF_INVALID_PARAMETER;

    try {
        ((PdfPages*)pages)->SetCropBox(rect);
        return PDF_SUCCESS;
    } catch (PdfException &e) {
        return e.GetCode();
    } catch (...) {
        return PDF_UNKNOWN_ERROR;
    }
}

int
pdf_pages_set_resources(pdf_pages pages, pdf_dictionary resources)
{
    if (pages == NULL)
        return PDF_INVALID_PARAMETER;
    
    try {
        ((PdfPages*)pages)->SetResources((PdfDictionary*)resources);
        return PDF_SUCCESS;
    } catch (PdfException &e) {
        return e.GetCode();
    } catch (...) {
        return PDF_UNKNOWN_ERROR;
    }
}

int
pdf_pages_set_size(pdf_pages pages, int width, int height)
{
    if (pages == NULL)
        return PDF_INVALID_PARAMETER;
    
    try {
        ((PdfPageBase*)pages)->SetSize(width, height);
        return PDF_SUCCESS;
    } catch (PdfException &e) {
        return e.GetCode();
    } catch (...) {
        return PDF_UNKNOWN_ERROR;
    }
}

pdf_pages
pdf_pages_add_pages(pdf_pages pages, int index)
{
    if (pages == NULL)
        return NULL;
    
    try {
        return ((PdfPages*)pages)->AddPages(index);
    } catch (PdfException &e) {
        ((PdfPages*)pages)->SetError(e.GetCode());
    } catch (...) {
        ((PdfPages*)pages)->SetError(PDF_UNKNOWN_ERROR);
    }
    return NULL;
}

pdf_page
pdf_pages_add_page(pdf_pages pages, int index)
{
    if (pages == NULL)
        return NULL;
    
    try {
        return ((PdfPages*)pages)->AddPage(index);
    } catch (PdfException &e) {
        ((PdfPages*)pages)->SetError(e.GetCode());
    } catch (...) {
        ((PdfPages*)pages)->SetError(PDF_UNKNOWN_ERROR);
    }
    return NULL;
}

/*---------------------------------------------------------------------------*/
/*----- pdf_page ------------------------------------------------------------*/

pdf_dictionary 
pdf_page_get_resource(pdf_page page, const char *element_name)
{
    return pdf_pages_get_resource(page, element_name);
}

int
pdf_page_set_media_box(pdf_page page, pdf_box rect)
{
    return pdf_pages_set_media_box(page, rect);
}

int
pdf_page_set_crop_box(pdf_page page, pdf_box rect)
{
    return pdf_pages_set_crop_box(page, rect);
}

int
pdf_page_set_resources(pdf_page page, pdf_dictionary resources)
{
    return pdf_pages_set_resources(page, resources);
}
    
int
pdf_page_set_size(pdf_page page, int width, int height)
{
    return pdf_pages_set_size(page, width, height);
}

pdf_contents
pdf_page_get_canvas(pdf_page page)
{
    if (page == NULL)
        return NULL;
    try {
        return ((PdfPage*)page)->Canvas();
    } catch (...) {
        return NULL;
    }
}

int
pdf_page_get_width(pdf_page page)
{
    return (page != NULL) ? ((PdfPage*)page)->Width() : -1;
}

int
pdf_page_get_height(pdf_page page)
{
    return (page == NULL) ? ((PdfPage*)page)->Height() : -1;
}

pdf_link_anot pdf_page_add_link(pdf_page page, pdf_rect rect,
        pdf_destination dest, pdf_annot_hl_mode mode)
{
    if (page == NULL)
        return NULL;
    try {
        return ((PdfPage*)page)->AddLink(rect, (PdfDestination*)dest, mode);
    } catch (...) {
        return NULL;
    }
}

pdf_box
pdf_page_get_media_box(pdf_page page)
{
    return pdf_pages_get_media_box(page);
}

pdf_box
pdf_page_get_crop_box(pdf_page page)
{
    return pdf_pages_get_crop_box(page);
}

pdf_dictionary 
pdf_page_get_resources(pdf_page page)
{
    return pdf_pages_get_resources(page);
}

int
pdf_page_get_rotate(pdf_page page)
{
    return pdf_pages_get_rotate(page);
}

pdf_pages
pdf_page_get_parent(pdf_page page)
{
    return pdf_pages_get_parent(page);
}

/*---------------------------------------------------------------------------*/
/*----- pdf_contents --------------------------------------------------------*/

/*
 * This macro protect C program from aborting program by exception.
 */

#define PDF_CONTENTS_SAFE_CALL(CONTENTS, FUNC, ARGS) \
    if (CONTENTS == NULL) \
        return PDF_INVALID_PARAMETER; \
    try { \
        ((PdfContents*)CONTENTS)->FUNC ARGS; \
        return PDF_SUCCESS; \
    } catch (PdfException &e) { \
        ((PdfDictionary*)CONTENTS)->SetError(e.GetCode()); \
        return e.GetCode(); \
    } catch (...) { \
        ((PdfDictionary*)CONTENTS)->SetError(PDF_UNKNOWN_ERROR); \
        return PDF_UNKNOWN_ERROR; \
    } \

pdf_page
pdf_contents_get_page(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->Page() : NULL;
}

int 
pdf_contents_get_width(pdf_contents contents)
{
    return (contents != NULL && ((PdfContents*)contents)->Page() != NULL) ?
        ((PdfContents*)contents)->Width() : -1;
}

int 
pdf_contents_get_height(pdf_contents contents)
{
    return (contents != NULL && ((PdfContents*)contents)->Page() != NULL) ?
        ((PdfContents*)contents)->Height() : -1;
}

pdf_point
pdf_contents_get_cur_point(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->CurPoint() : 
        PdfPoint(0, 0);
}

pdf_point
pdf_contents_get_text_pos(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->TextPos() : 
        PdfPoint(0, 0);
}

double
pdf_contents_get_test_leading(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->TextLeading() : 0;
}

double
pdf_contents_get_font_size(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->FontSize() : 0;
}

const char 
*pdf_contents_get_font_name(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->FontName() : NULL;
}

double
pdf_contents_get_char_space(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->CharSpace() : 0;
}

double
pdf_contents_get_word_space(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->WordSpace() : 0;
}

double
pdf_contents_get_text_width(pdf_contents contents, const char *text,
        int *numchars, int *numwords)
{
    if (contents == NULL) {
        *numchars = 0;
        *numwords = 0;
        return 0;
    } else 
        return ((PdfContents*)contents)->TextWidth(text, numchars, numwords);
}

int 
pdf_contents_get_char_widths(pdf_contents contents, const char*text, 
        double *widths)
{
    PDF_CONTENTS_SAFE_CALL(contents, CharWidths, (text, widths));
}

unsigned int
pdf_contents_measure_text(pdf_contents contents, const char *text, 
        double width, double *realwidth)
{
    return (contents != NULL) ? 
        ((PdfContents*)contents)->MeasureText(text, width, realwidth) : 0;
}

pdf_rgb_color
pdf_contents_get_rgb_fill(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->RGBFill() : 
        PdfRGBColor(0, 0, 0);
}

pdf_rgb_color
pdf_contents_get_rgb_stroke(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->RGBStroke() :
        PdfRGBColor(0, 0, 0);
}

double
pdf_contents_get_gray_fill(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->GrayFill() : 0;
}

double
pdf_contents_get_gray_stroke(pdf_contents contents)
{
    return (contents != NULL) ? ((PdfContents*)contents)->GrayStroke() : 0;
}

int 
pdf_contents_set_line_width(pdf_contents contents, double linewidth)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetLineWidth, (linewidth));
}

int
pdf_contents_set_line_cap(pdf_contents contents, pdf_line_cap_style linecap)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetLineCap, (linecap));
}

int
pdf_contents_set_line_join(pdf_contents contents, pdf_line_join_style linejoin)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetLineJoin, (linejoin));
}

int
pdf_contents_set_miter_limit(pdf_contents contents, double miterlimit)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetMiterLimit, (miterlimit));
}

int
pdf_contents_set_dash(pdf_contents contents, unsigned int on, 
        unsigned int off, unsigned int phase)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetDash, (on, off, phase));
}

int
pdf_contents_set_flat(pdf_contents contents, unsigned int flatness)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetFlat, (flatness));
}

int
pdf_contents_gsave(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, GSave, ());
}

int
pdf_contents_grestore(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, GRestore, ());
}

int
pdf_contents_concat(pdf_contents contents, double a, double b, double c,
        double d, double e, double f)
{
    PDF_CONTENTS_SAFE_CALL(contents, Concat, (a, b, c, d, e, f));
}

int
pdf_contents_move_to(pdf_contents contents, double x, double y)
{
    PDF_CONTENTS_SAFE_CALL(contents, MoveTo, (x, y));
}

int
pdf_contents_line_to(pdf_contents contents, double x, double y)
{
    PDF_CONTENTS_SAFE_CALL(contents, LineTo, (x, y));
}

int
pdf_contents_curve_to(pdf_contents contents, double x1, double y1, double x2, 
        double y2, double x3, double y3)
{
    PDF_CONTENTS_SAFE_CALL(contents, CurveTo, (x1, y1, x2, y2, x3, y3));
}

int
pdf_contents_curve_to2(pdf_contents contents, double x2, double y2, double x3, 
        double y3)
{
    PDF_CONTENTS_SAFE_CALL(contents, CurveTo2, (x2, y2, x3, y3));
}

int
pdf_contents_curve_to3(pdf_contents contents, double x1, double y1, double x3,
        double y3)
{
    PDF_CONTENTS_SAFE_CALL(contents, CurveTo3, (x1, y1, x3, y3));
}

int
pdf_contents_close_path(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, ClosePath, ());
}

int
pdf_contents_rectangle(pdf_contents contents, double x, double y, double width, 
        double height)
{
    PDF_CONTENTS_SAFE_CALL(contents, Rectangle, (x, y, width, height));
}

int
pdf_contents_stroke(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, Stroke, ());
}

int
pdf_contents_fill(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, Fill, ());
}

int
pdf_contents_eofill(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, Eofill, ());
}

int
pdf_contents_fill_stroke(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, FillStroke, ());
}

int
pdf_contents_eofill_stroke(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, EofillStroke, ());
}

int
pdf_contents_close_path_fill_stroke(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, ClosePathFillStroke, ());
}

int
pdf_contents_close_path_eofill_stroke(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, ClosePathEofillStroke, ());
}

int
pdf_contents_end_path(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, EndPath, ());
}

int
pdf_contents_clip(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, Clip, ());
}

int
pdf_contents_eoclip(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, EoClip, ());
}

int
pdf_contents_begin_text(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, BeginText, ());
}

int
pdf_contents_end_text(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, EndText, ());
}

int
pdf_contents_set_char_space(pdf_contents contents, double value)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetCharSpace, (value));
}

int
pdf_contents_set_word_space(pdf_contents contents, double value)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetWordSpace, (value));
}

int
pdf_contents_set_horizontal_scalling(pdf_contents contents, double value)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetHorizontalScalling, (value));
}

int
pdf_contents_set_text_leading(pdf_contents contents, double value)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetTextLeading, (value));
}

int
pdf_contents_set_font_and_size(pdf_contents contents, const char *fontname, 
        double size)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetFontAndSize, (fontname, size));
}

int
pdf_contents_set_text_rendering_mode(pdf_contents contents, 
        pdf_text_rendering_mode mode)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetTextRenderingMode, (mode));
}

int
pdf_contents_set_text_raise(pdf_contents contents, double value)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetTextRaise, (value));
}

int
pdf_contents_move_text_pos(pdf_contents contents, double x, double y)
{
    PDF_CONTENTS_SAFE_CALL(contents, MoveTextPos, (x, y));
}

int
pdf_contents_move_text_pos2(pdf_contents contents, double x, double y)
{
    PDF_CONTENTS_SAFE_CALL(contents, MoveTextPos2, (x, y));
}

int pdf_contents_set_text_matrix(pdf_contents contents, double a, double b,
        double c, double d, double x, double y)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetTextMatrix, (a, b, c, d, x, y));
}

int
pdf_contents_move_to_next_line(pdf_contents contents)
{
    PDF_CONTENTS_SAFE_CALL(contents, MoveToNextLine, ());
}

int
pdf_contents_show_text(pdf_contents contents, const char *text)
{
    PDF_CONTENTS_SAFE_CALL(contents, ShowText, (text));
}

int
pdf_contents_show_text_next_line(pdf_contents contents, const char *text)
{
    PDF_CONTENTS_SAFE_CALL(contents, ShowTextNextLine, (text));
}

int
pdf_contents_set_gray_fill(pdf_contents contents, double gray)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetGrayFill, (gray));
}

int
pdf_contents_set_gray_stroke(pdf_contents contents, double gray)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetGrayStroke, (gray));
}

int
pdf_contents_set_rgb_fill(pdf_contents contents, pdf_rgb_color c)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetRGBFill, (c));
}

int
pdf_contents_set_rgb_stroke(pdf_contents contents, pdf_rgb_color c)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetRGBStroke, (c));
}

int
pdf_contents_set_cmyk_fill(pdf_contents contents, double c, double m, 
        double y, double k)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetCMYKFill, (c, y, m, k));
}

int
pdf_contents_set_cmyk_stroke(pdf_contents contents, double c, double m, 
        double y, double k)
{
    PDF_CONTENTS_SAFE_CALL(contents, SetCMYKStroke, (c, y, m, k));
}

int
pdf_contents_execute_xobject(pdf_contents contents, const char *name)
{
    PDF_CONTENTS_SAFE_CALL(contents, ExecuteXObject, (name));
}

int
pdf_contents_execute_xobject2(pdf_contents contents, pdf_xobject xobject)
{
    if (xobject == NULL)
        return PDF_INVALID_PARAMETER;
    PDF_CONTENTS_SAFE_CALL(contents, ExecuteXObject, ((PdfXObject*)xobject));
}

/*---------------------------------------------------------------------------*/
/*----- pdf_info ------------------------------------------------------------*/

#define PDF_INFO_SAFE_CALL(INFO, FUNC, ARGS) \
    if (INFO == NULL) \
        return PDF_INVALID_PARAMETER; \
    try { \
        ((PdfInfo*)INFO)->FUNC ARGS; \
        return PDF_SUCCESS; \
    } catch (PdfException &e) { \
        ((PdfDictionary*)INFO)->SetError(e.GetCode()); \
        return e.GetCode(); \
    } catch (...) { \
        ((PdfDictionary*)INFO)->SetError(PDF_UNKNOWN_ERROR); \
        return PDF_UNKNOWN_ERROR; \
    } \

const char*
pdf_info_get_author(pdf_info info)
{
    return (info != NULL) ? ((PdfInfo*)info)->Author() : NULL;
}

const char*
pdf_info_get_creator(pdf_info info)
{
    return (info != NULL) ? ((PdfInfo*)info)->Creator() : NULL;
}

const char*
pdf_info_get_producer(pdf_info info)
{
    return (info != NULL) ? ((PdfInfo*)info)->Producer() : NULL;
}

const char*
pdf_info_get_title(pdf_info info)
{
    return (info != NULL) ? ((PdfInfo*)info)->Title() : NULL;
}

const char*
pdf_info_get_subject(pdf_info info)
{
    return (info != NULL) ? ((PdfInfo*)info)->Subject() : NULL;
}

const char*
pdf_info_get_keywords(pdf_info info)
{
    return (info != NULL) ? ((PdfInfo*)info)->Keywords() : NULL;
}

int
pdf_info_set_author(pdf_info info, const char* value)
{
    PDF_INFO_SAFE_CALL(info, SetAuthor, (value));
}

int
pdf_info_set_creator(pdf_info info, const char* value)
{
    PDF_INFO_SAFE_CALL(info, SetCreator, (value));
}

int
pdf_info_set_producer(pdf_info info, const char* value)
{
    PDF_INFO_SAFE_CALL(info, SetProducer, (value));
}

int
pdf_info_set_title(pdf_info info, const char* value)
{
    PDF_INFO_SAFE_CALL(info, SetTitle, (value));
}

int
pdf_info_set_subject(pdf_info info, const char* value)
{
    PDF_INFO_SAFE_CALL(info, SetSubject, (value));
}

int
pdf_info_set_keywords(pdf_info info, const char* value)
{
    PDF_INFO_SAFE_CALL(info, SetKeywords, (value));
}

/*---------------------------------------------------------------------------*/
/*----- pdf_image -----------------------------------------------------------*/

double
pdf_image_get_width(pdf_image image)
{
    if (image == NULL)
        return 0;

    PdfImage* image_obj = (PdfImage*)image;
    return image_obj->Width();
}

double
pdf_image_get_height(pdf_image image)
{
    if (image == NULL)
        return 0;

    PdfImage* image_obj = (PdfImage*)image;
    return image_obj->Height();
}

/*---------------------------------------------------------------------------*/
/*----- pdf_doc -------------------------------------------------------------*/

pdf_doc pdf_doc_new()
{
    try {
        PdfDoc* doc = new PdfDoc();
        return doc;
    } catch (...) {
        return NULL;
    }
}

void pdf_doc_free(pdf_doc doc)
{
    delete (PdfDoc*)doc;
}

int pdf_doc_new_doc(pdf_doc doc)
{
    if (doc == NULL)
        return PDF_INVALID_PARAMETER;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->NewDoc();
        doc_obj->SetError(0);
    } catch (PdfException &e) {
        doc_obj->SetError(e.GetCode());
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
    }
    return doc_obj->LastError();
}

int pdf_doc_free_doc(pdf_doc doc)
{
    if (doc == NULL)
        return PDF_INVALID_PARAMETER;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->FreeDoc();
        doc_obj->SetError(0);
    } catch (PdfException &e) {
        doc_obj->SetError(e.GetCode());
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
    }
    return doc_obj->LastError();
}

int pdf_doc_write_to_file(pdf_doc doc, const char *filename)
{
    if (doc == NULL)
        return PDF_INVALID_PARAMETER;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->WriteToFile(filename);
        doc_obj->SetError(0);
    } catch (PdfException &e) {
        doc_obj->SetError(e.GetCode());
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
    }
    return doc_obj->LastError();
}

int
pdf_doc_add_type1font(pdf_doc doc, pdf_type1_fontdef font, 
        const char *name, pdf_encodingdef encoding)
{
    if (doc == NULL)
        return PDF_INVALID_PARAMETER;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->AddType1Font((PdfType1FontDef*)font, name, 
                (PdfEncodingDef*)encoding);
        doc_obj->SetError(0);
    } catch (PdfException &e) {
        doc_obj->SetError(e.GetCode());
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
    }
    return doc_obj->LastError();
}

int
pdf_doc_add_type0font(pdf_doc doc, pdf_cid_fontdef font, 
        const char *name, pdf_cmap map)
{
    if (doc == NULL)
        return PDF_INVALID_PARAMETER;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->AddType0Font((PdfCIDFontDef*)font, name, 
                (PdfCMap*)map);
        doc_obj->SetError(0);
    } catch (PdfException &e) {
        doc_obj->SetError(e.GetCode());
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
    }
    return doc_obj->LastError();
}

pdf_info
pdf_doc_get_info(pdf_doc doc)
{
    if (doc == NULL)
        return NULL;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->SetError(0);
        return doc_obj->Info(); 
    } catch (PdfException &e) {
        doc_obj->SetError(e.GetCode());
        return NULL;
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
        return NULL;
    }
}

int pdf_doc_add_xobject(pdf_doc doc, pdf_xobject obj, const char *name)
{
    if (doc == NULL)
        return PDF_INVALID_PARAMETER;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->AddXObject((PdfXObject*)obj, name);
        doc_obj->SetError(0);
    } catch (PdfException &e) {
        doc_obj->SetError(e.GetCode());
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
    }
    return doc_obj->LastError();
}

int pdf_doc_register_object(pdf_doc doc, pdf_auto_ptr_object obj)
{
    if (doc == NULL)
        return PDF_INVALID_PARAMETER;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->RegisterObject((PdfAutoPtrObject*)obj);
        doc_obj->SetError(0);
    } catch (PdfException &e) {
        doc_obj->SetError(e.GetCode());
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
    }
    return doc_obj->LastError();
}
    
pdf_page
pdf_doc_add_page(pdf_doc doc)
{
    if (doc == NULL)
        return NULL;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->SetError(0);
        return doc_obj->AddPage();
    } catch (PdfException &e) {
        doc_obj->SetError(e.GetCode());
        return NULL;
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
        return NULL;
    }
}

pdf_pages 
pdf_doc_root_pages(pdf_doc doc)
{
    if (doc == NULL)
        return NULL;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->SetError(0);
        return doc_obj->RootPages();
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
        return NULL;
    }
}

pdf_page
pdf_doc_current_page(pdf_doc doc)
{
    if (doc == NULL)
        return NULL;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->SetError(0);
        return doc_obj->CurrentPage();
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
        return NULL;
    }
}

int
pdf_doc_has_doc(pdf_doc doc)
{
    if (doc == NULL) return PDF_FALSE;
    return (((PdfDoc*)doc)->HasDoc() == true) ? PDF_TRUE : PDF_FALSE;
}

const char
*pdf_doc_get_font_name(pdf_doc doc, int index)
{
    if (doc == NULL) return NULL;
    PdfFont* font = ((PdfDoc*)doc)->FontMgr()->GetFont(index);
    return (font != NULL) ? font->Name() : NULL;
}

int
pdf_doc_count_fonts(pdf_doc doc)
{
    return (doc != NULL) ? 
        (int)((PdfDoc*)doc)->FontMgr()->CountFonts() : -1;
}

int
pdf_doc_get_error(pdf_doc doc)
{
    return (doc != NULL) ? ((PdfDoc*)doc)->LastError() : 0;
}

int
pdf_doc_set_password(pdf_doc doc, const char* owner_passwd, 
        const char* user_passwd)
{
    if (doc == NULL)
        return PDF_INVALID_PARAMETER;

    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->SetPassword(owner_passwd, user_passwd);
    } catch (PdfException& e) {
        doc_obj->SetError(e.GetCode());
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
    }

    return doc_obj->LastError();
}

int
pdf_doc_set_permission(pdf_doc doc, int permission)
{
    if (doc == NULL)
        return PDF_INVALID_PARAMETER;
    
    PdfDoc* doc_obj = (PdfDoc*)doc;
    try {
        doc_obj->SetPermission(permission);
    } catch (PdfException& e) {
        doc_obj->SetError(e.GetCode());
    } catch (...) {
        doc_obj->SetError(PDF_UNKNOWN_ERROR);
    }

    return doc_obj->LastError();
}

