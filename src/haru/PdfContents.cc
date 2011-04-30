/*
 * << H a r u --free pdf library >> -- PdfContents.cpp
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

#include <assert.h>
#include "libharu.h"

/*----- PdfContents class ----------------------------------------------------*/

PdfContents::PdfContents(PdfPage *page)
        : PdfStream(page->GetXref())
{
    /* initialize attributes to the default */
    fPage = page;
    fFontMgr = fPage->FontMgr();
    fXObjectMgr = fPage->XObjectMgr();
    
    fWordSpace = PDF_DEF_WORDSPACE;
    fCharSpace = PDF_DEF_CHARSPACE;
    fFontSize = PDF_DEF_FONTSIZE;
    fHScalling = PDF_DEF_HSCALING;
    fTextLeading = PDF_DEF_LEADING;
    fRenderingMode = PDF_DEF_RENDERING_MODE;
    fTextRaise = PDF_DEF_RAISE;
    fLineWidth = PDF_DEF_LINEWIDTH;
    fLineCap = PDF_DEF_LINECAP;
    fLineJoin = PDF_DEF_LINEJOIN;
    fMiterLimit = PDF_DEF_MITERLIMIT;
    fDashOn = 0;
    fDashOff = 0;
    fDashPhase = 0;
    fFlatness = 0;
    
    fCurPoint.x = 0;
    fCurPoint.y = 0;
    fTextPoint.x = 0;
    fTextPoint.y = 0;
    fMatrix = PdfTextMatrix(1, 0, 0, 1, 0, 0);
    fGMode = PDF_GMODE_PAGE_DESCRIPTION;

    fRGBFill.red = 0;
    fRGBFill.green = 0;
    fRGBFill.blue = 0;
    fRGBStroke.red = 0;
    fRGBStroke.green = 0;
    fRGBStroke.blue = 0;
    fGrayFill = 0;
    fGrayStroke = 0;
}

void
PdfContents::Init()
{
    GetXref()->AddObject(this);
    fPage->AddElement("Contents", this);

    PdfNumber *length = new PdfNumber(0);
    GetXref()->AddObject(length);
    AddElement("Length", length);
}

/* graphics mode checking macro */
#define GMODE_ERROR(method_name)  { \
    throw PdfException(PDF_RUNTIME_ERROR, "ERROR: %s invalid graphics " \
            "mode(%d).", method_name, (int)fGMode); \
}

/* font checking macro */



double
PdfContents::TextWidth(const char* text, int* numchars, int* numwords)
{
    if (fFont == NULL) 
        throw PdfException(PDF_RUNTIME_ERROR, "ERROR: PdfContents::TextWidth "
                "current font has not been set.");

    pdf_text_width pw = fFont->TextWidth(text);

    double w1 = pw.width * fFontSize / 1000;
    double w2 = (pw.numchars - 1) * fCharSpace;
    double w3 = (pw.numwords - 1) * fWordSpace;

    PDF_DEBUG_PRINT(("PdfContents::TextWidth pw.width=%d, pw.numchars=%d, "
               "pw.numwords=%d\n", pw.width, pw.numchars, pw.numwords));
    
    if (numchars != NULL)
        *numchars = pw.numchars;
    
    if (numwords != NULL)
        *numwords = pw.numwords;

    return w1 + w2 + w3;
}

void
PdfContents::CharWidths(const char* chars, double* widths)
{
    assert(chars);
    assert(widths);

    if (fFont == NULL) 
        throw PdfException(PDF_RUNTIME_ERROR, "ERROR: PdfContents::CharWidths "
                "current font has not been set.");

    int len = strlen(chars);
    unsigned int* w = new unsigned int[len];
    unsigned int* pw = w;
    double* pf = widths;

    fFont->TextWidths(chars, w);
    for (int i = 0; i < len; i++) {
        *pf = *pw * fFontSize / 1000;
        pf++;
        pw++;
    }

    delete[] w;
}

unsigned int
PdfContents::MeasureText(const char* text, double width, double* realwidth)
{
    assert(text);
    assert(width > 0);

    if (fFont == NULL) 
        throw PdfException(PDF_RUNTIME_ERROR, "ERROR: "
                "PdfContents::MeasureText current font has not been set.");

    return fFont->MeasureText(text, width, fFontSize, fCharSpace,
            fWordSpace, realwidth);
}

void
PdfContents::TextOut(double x, double y, const char* text)
{
    BeginText();
    MoveTextPos(x, y);
    ShowText(text);
    EndText();
}

void
PdfContents::TextRect(const char* text, pdf_rect rect, int max_len)
{
    char *tmp;
    double w = rect.right - rect.left;
    int len = strlen(text);
    if (len == 0)
        return;

    /* copy the text to temporary buffer */
    tmp = new char[len + 1];
    int j = 0;
    for (int i = 0; i < len; i++) {
        if (text[i] != 0x0a) {
            tmp[j] = text[i];
            j++;
        }
    }
    len = j;

    /* replace \n charactor to 0x00 */
    for (int k = 0; k < len; k++)
        if (tmp[k] == 0x0d)
            tmp[k] = 0x00;

    try {
        char* buf = new char[max_len + 1];

        try {
            BeginText();
            double ascent = Font()->Ascent() * FontSize() / 1000;
            MoveTextPos(rect.left, rect.top - ascent);

            char* pstr = tmp;
            char* end = tmp + len;

            while (pstr < end) {
                if (*pstr == 0x00)
                    MoveToNextLine();
                else {
                    /* calcurate how many charactor of the text can be 
                     * included in a line
                     */
                    int cnt = MeasureText(pstr, w);
                    if (cnt <= max_len) {
                        memcpy(buf, pstr, cnt);
                        buf[cnt] = 0x00;
                    } else {
                        memcpy(buf, pstr, max_len);
                        buf[max_len] = 0x00;
                    }
                    ShowText(buf);
                    MoveToNextLine();

                    if ((TextPos()).y < rect.bottom + TextLeading())
                        break;
                    pstr += cnt;
                }
            }

            EndText();
        } catch (...) {
            delete[] buf;
            throw;
        }
        delete[] buf;
    } catch (...) {
        delete[] tmp;
        throw;
    }
    delete[] tmp;
}

/*----- General graphics state -------------------------------------*/

void
PdfContents::SetLineWidth(double linewidth)
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION && fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("SetLineWidth")

    if (linewidth < 0)
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
                "ERROR: SetLineWidth invalid value(%d).", (int)linewidth);

    *GetStream() << linewidth
             << " w\012";
    fLineWidth = linewidth;
    return;
}

void
PdfContents::SetLineCap(pdf_line_cap_style linecap)
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION && fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("SetLineCap")

    *GetStream() << (unsigned int)linecap
             << " J\012";
    fLineCap = linecap;
    return;
}

void
PdfContents::SetLineJoin(pdf_line_join_style linejoin)
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION && fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("SetLineJoin")

    *GetStream() << (unsigned int)linejoin
             << " j\012";
    fLineJoin = linejoin;
    return;
}

void
PdfContents::SetMiterLimit(double miterlimit)
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION && fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("SetMiterLimit")

    if (miterlimit < 1)
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
            "ERROR: PdfContents::SetLineJoin (%d)"
            ")--linejoin must be more than or equal to 1.", (int)miterlimit);

    *GetStream() << miterlimit
             << " M\012";
    fMiterLimit = miterlimit;
    return;
}

void
PdfContents::SetDash(unsigned int on, unsigned int off, unsigned int phase)
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION && fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("SetDash")

    fDashOn = on;
    if (on == off || on == 0)
        fDashOff = 0;
    else
        fDashOff = off;
    fDashPhase = phase;

    *GetStream() << '[';
    if (fDashOn > 0) {
        *GetStream() << fDashOn;
        if (fDashOff > 0) {
            *GetStream() << ' ' << fDashOff;
        }
    }
    *GetStream() << ']'
             << fDashPhase
             << " d\012";
    return;
}

void
PdfContents::SetFlat(unsigned int flatness)
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION && fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("SetFlat")

    if (flatness > 100)
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
            "ERROR: PdfContents::SetFlat flatness must "
            "be between 0 to 100.(%d).\n", flatness);

    fFlatness = flatness;
    *GetStream() << flatness
             << " i\012";
    return;
}

/*----- Special graphics state -------------------------------------*/

void
PdfContents::GSave()
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("GSave")

    *GetStream() << "q\012";
    return;
}

void
PdfContents::GRestore()
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("GRestore")

    *GetStream() << "Q\012";
    return;
}

void
PdfContents::Concat(double a, double b, double c, double d, double e, double f)
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("Concat")

    *GetStream() << a << ' '
             << b << ' '
             << c << ' '
             << d << ' '
             << e << ' '
             << f << " cm\012";
    return;
}

/*----- Path construction ------------------------------------------*/

void
PdfContents::MoveTo(double x, double y)
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION && fGMode != PDF_GMODE_PATH_OBJECT)
        GMODE_ERROR("MoveTo")

    *GetStream() << x << ' '
             << y << " m\012";
    fCurPoint.x = x;
    fCurPoint.y = y;
    fStartPoint.x = x;
    fStartPoint.y = y;
    fGMode = PDF_GMODE_PATH_OBJECT;
    return;
}

void
PdfContents::LineTo(double x, double y)
{
    if (fGMode != PDF_GMODE_PATH_OBJECT)
        GMODE_ERROR("LineTo")

    *GetStream() << x << ' '
             << y << " l\012";
    fCurPoint.x = x;
    fCurPoint.y = y;
    return;
}

void
PdfContents::CurveTo(double x1, double y1, double x2,
        double y2, double x3, double y3)
{
    if (fGMode != PDF_GMODE_PATH_OBJECT)
        GMODE_ERROR("CurveTo")

    *GetStream() << x1 << ' '
             << y1 << ' '
             << x2 << ' '
             << y2 << ' '
             << x3 << ' '
             << y3 << " c\012";
    fCurPoint.x = x3;
    fCurPoint.y = y3;
    return;
}

void
PdfContents::CurveTo2(double x2, double y2, double x3, double y3)
{
    if (fGMode != PDF_GMODE_PATH_OBJECT)
        GMODE_ERROR("CurveTo2")

    *GetStream() << x2 << ' '
             << y2 << ' '
             << x3 << ' '
             << y3 << " v\012";
    fCurPoint.x = x3;
    fCurPoint.y = y3;
    return;
}

void
PdfContents::CurveTo3(double x1, double y1, double x3,
                double y3)
{
    if (fGMode != PDF_GMODE_PATH_OBJECT)
        GMODE_ERROR("CurveTo3")

    *GetStream() << x1 << ' '
             << y1 << ' '
             << x3 << ' '
             << y3 << " y\012";
    fCurPoint.x = x3;
    fCurPoint.y = y3;
    return;
}

void
PdfContents::ClosePath()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT)
        GMODE_ERROR("ClosePath")

    *GetStream() << "h\012";
    fCurPoint = fStartPoint;
    return;
}

void
PdfContents::Rectangle(double x, double y, double width,
                double height)
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION && fGMode != PDF_GMODE_PATH_OBJECT)
        GMODE_ERROR("Rectangle")

    *GetStream() << x << ' '
             << y << ' '
             << width << ' '
             << height << " re\012";
    fCurPoint.x = x;
    fCurPoint.y = y;
    fGMode = PDF_GMODE_PATH_OBJECT;
    return;
}

/*----- Path painting ----------------------------------------------*/

void
PdfContents::Stroke()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT && fGMode != PDF_GMODE_CLIPPING_PATH)
        GMODE_ERROR("Stroke")

    *GetStream() << "S\012";
    fGMode = PDF_GMODE_PAGE_DESCRIPTION;
    return;
}

void
PdfContents::ClosePathStroke()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT && fGMode != PDF_GMODE_CLIPPING_PATH)
        GMODE_ERROR("ClosePathStroke")

    *GetStream() << "s\012";
    fGMode = PDF_GMODE_PAGE_DESCRIPTION;
    return;
}

void
PdfContents::Fill()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT && fGMode != PDF_GMODE_CLIPPING_PATH)
        GMODE_ERROR("Fill")

    *GetStream() << "f\012";
    fGMode = PDF_GMODE_PAGE_DESCRIPTION;
    return;
}

void
PdfContents::Eofill()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT && fGMode != PDF_GMODE_CLIPPING_PATH)
        GMODE_ERROR("Eofill")

    *GetStream() << "f*\012";
    fGMode = PDF_GMODE_PAGE_DESCRIPTION;
    return;
}

void
PdfContents::FillStroke()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT && fGMode != PDF_GMODE_CLIPPING_PATH)
        GMODE_ERROR("FillStroke")

    *GetStream() << "B\012";
    fGMode = PDF_GMODE_PAGE_DESCRIPTION;
    return;
}

void
PdfContents::EofillStroke()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT && fGMode != PDF_GMODE_CLIPPING_PATH)
        GMODE_ERROR("EofillStroke")

    *GetStream() << "B*\012";
    fGMode = PDF_GMODE_PAGE_DESCRIPTION;
    return;
}

void
PdfContents::ClosePathFillStroke()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT && fGMode != PDF_GMODE_CLIPPING_PATH)
        GMODE_ERROR("ClosePathFillStroke")

    *GetStream() << "b\012";
    fGMode = PDF_GMODE_PAGE_DESCRIPTION;
    return;
}

void
PdfContents::ClosePathEofillStroke()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT && fGMode != PDF_GMODE_CLIPPING_PATH)
        GMODE_ERROR("ClosePathEofillStroke")

    *GetStream() << "b*\012";
    fGMode = PDF_GMODE_PAGE_DESCRIPTION;
    return;
}

void
PdfContents::EndPath()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT && fGMode != PDF_GMODE_CLIPPING_PATH)
        GMODE_ERROR("EndPath")

    *GetStream() << "n\012";
    fCurPoint.x = 0;  //
    fCurPoint.y = 0;  //
    fGMode = PDF_GMODE_PAGE_DESCRIPTION;
    return;
}

/*----- Clipping paths ---------------------------------------------*/

void
PdfContents::Clip()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT)
        GMODE_ERROR("Clip")

    *GetStream() << "W\012";
    fGMode = PDF_GMODE_CLIPPING_PATH;
    return;
}

void
PdfContents::EoClip()
{
    if (fGMode != PDF_GMODE_PATH_OBJECT)
        GMODE_ERROR("Eolip")

    *GetStream() << "W*\012";
    fGMode = PDF_GMODE_CLIPPING_PATH;
    return;
}

/*----- Text object ------------------------------------------------*/

void
PdfContents::BeginText()
{
    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("BeginText")

    fPage->AddProcSet(PDF_PROCSET_TEXT);

    *GetStream() << "BT\012";
    fTextPoint.x = 0;
    fTextPoint.y = 0;
    fMatrix = PdfTextMatrix(1, 0, 0, 1, 0, 0);

    fGMode = PDF_GMODE_TEXT_OBJECT;
    return;
}

void
PdfContents::EndText()
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("EndText")

    *GetStream() << "ET\012";

    fGMode = PDF_GMODE_PAGE_DESCRIPTION;
    return;
}

/*----- Text state -------------------------------------------------*/

void
PdfContents::SetCharSpace(double value)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetCharSpace")

    if (value < PDF_MIN_CHARSPACE || value > PDF_MAX_CHARSPACE)
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
            "ERROR: SetCharSpace out of range(%f)", value);

    *GetStream() << value << " "
             << "Tc\012";
    fCharSpace = value;
    return;
}

void
PdfContents::SetWordSpace(double value)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetWordSpace")

    if (value < PDF_MIN_WORDSPACE || value > PDF_MAX_WORDSPACE)
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
            "ERROR: SetWordSpace out of range(%f)", value);

    *GetStream() << value << " "
             << "Tw\012";
    fWordSpace = value;
    return;
}

void
PdfContents::SetHorizontalScalling(double value)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetHorizontalScalling")

    if (value < PDF_MIN_HORIZONTALSCALING ||
        value > PDF_MAX_HORIZONTALSCALING) {
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
                "ERROR: SetHorizontalScalling out of range(%f)", value);
    }

    *GetStream() << value << " "
             << "Tz\012";
    fHScalling = value;
    return;
}

void
PdfContents::SetTextLeading(double value)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetTextLeading")

    *GetStream() << value << " "
             << "TL\012";
    fTextLeading = value;
    return;
}

void
PdfContents::SetFontAndSize(const char* fontname, double size)
{
    assert(fFontMgr);

    PdfFont* font = fFontMgr->GetFont(fontname);
    if (font == NULL) {
        throw PdfException(PDF_RUNTIME_ERROR,
            "ERROR: font[%s] is not found.", fontname);
    }

    SetFontAndSize(font, size);
}

void
PdfContents::SetFontAndSize(PdfFont* font, double size)
{
    if (font == NULL) 
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: Invalid paramator[font].");
    if (size <= 0)
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: Invalid font size[%f].", size);
    
    fFont = font;
    const char* local_fontname = fPage->GetFontName(font);
    if (local_fontname != NULL)
        WriteEscapeName(GetStream(), local_fontname);
    else {
        char tmp_fontname[20];

#ifdef __WIN32__
        _snprintf(tmp_fontname, 20, "F%d", (int)fPage->CountFonts());
#else
        snprintf(tmp_fontname, 20, "F%d", (int)fPage->CountFonts());
#endif

        fPage->AddFont(font, tmp_fontname);
        WriteEscapeName(GetStream(), tmp_fontname);
    }

    *GetStream() << " "
             << size
             << " Tf\012";

    fFontSize = size;
    return;
}

void
PdfContents::SetTextRenderingMode(pdf_text_rendering_mode mode)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetTextRenderingMode")

    *GetStream() << (int)mode
             << " Tr\012";
    fRenderingMode = mode;
    return;
}

void
PdfContents::SetTextRaise(double value)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetTextRaise")

    *GetStream() << value
             << " Ts\012";
    fTextRaise = value;
    return;
}

/*--- Text positioning -----------------------------------------------------*/

void
PdfContents::MoveTextPos(double tx, double ty)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("MoveTextPos")

    *GetStream() << tx << " "
             << ty
             << " Td\012";
    fMatrix.x += tx;
    fMatrix.y += ty;
    fTextPoint.x = fMatrix.x;
    fTextPoint.y += ty;
    return;
}

void
PdfContents::MoveTextPos2(double tx, double ty)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("MoveTextPos2")

    *GetStream() << tx << " "
             << ty
             << " TD\012";
    fMatrix.x += tx;
    fMatrix.y += ty;
    fTextPoint.x = fMatrix.x;
    fTextPoint.y += ty;
    fTextLeading = -ty;
    return;
}

void
PdfContents::SetTextMatrix(double a, double b, double c,
                        double d, double x, double y)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("SetTextMatrix")

    *GetStream() << a << " "
             << b << " "
             << c << " "
             << d << " "
             << x << " "
             << y
             << " Tm\012";
    fMatrix.a = a;
    fMatrix.b = b;
    fMatrix.c = c;
    fMatrix.d = d;
    fMatrix.x = x;
    fMatrix.y = y;
    fTextPoint.x = fMatrix.x;
    fTextPoint.y = fMatrix.y;
    return;
}

void
PdfContents::MoveToNextLine()
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("MoveToNextLine")

    *GetStream() << "T*\012";

    /* calculate the reference point of text */
    fMatrix.x += fTextLeading * fMatrix.b;
    fMatrix.y -= fTextLeading * fMatrix.a;

    /* set text point to the start of new line. */
    fTextPoint.x = fMatrix.x;
    fTextPoint.y = fMatrix.y;

    PDF_DEBUG_PRINT(("PdfContents::MoveToNextLine -- x=%f, y=%f\n"
                , fTextPoint.x, fTextPoint.y));

    return;
}

/*--- Text showing ---------------------------------------------------------*/
void
PdfContents::ShowText(const char* text)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("ShowText")

    if (text == NULL || strlen(text) < 1)
        return;

    WriteEscapeText(GetStream(), text);
    *GetStream() << " Tj\012";

    /* set text point to the end of this line */
    double w = TextWidth(text);

    if (fFont->WritingMode() == PDF_WMODE_HORIZONTAL) {
        fTextPoint.x += w * fMatrix.a;
        fTextPoint.y += w * fMatrix.b;
    } else {
        fTextPoint.y -= w * fMatrix.a;
        fTextPoint.x -= w * fMatrix.b;
    }
    return;
}

void
PdfContents::ShowTextNextLine(const char* text)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("ShowTextNextLine")

    if (text == NULL || strlen(text) < 1)
        return;

    /* calculate the reference point of text. */
    fMatrix.x += fTextLeading * fMatrix.b;
    fMatrix.y -= fTextLeading * fMatrix.a;

    /* set text point to the start of new line. */
    fTextPoint.x = fMatrix.x;
    fTextPoint.y = fMatrix.y;

    WriteEscapeText(GetStream(), text);
    *GetStream() << " '\012";

    /* set text point to the end of this line */
    double w = TextWidth(text);

    if (fFont->WritingMode() == PDF_WMODE_HORIZONTAL) {
        fTextPoint.x += w * fMatrix.a;
        fTextPoint.y += w * fMatrix.b;
    } else {
        fTextPoint.y -= w * fMatrix.a;
        fTextPoint.x -= w * fMatrix.b;
    }

    return;
}

void
PdfContents::ShowTextNextLine(double aw, double ac,
        const char* text)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT)
        GMODE_ERROR("ShowTextNextLine")

    if (text == NULL || strlen(text) < 1)
        return;

    fWordSpace = aw;
    fCharSpace = ac;

    /* calculate the reference point of text */
    fMatrix.x += fTextLeading * fMatrix.b;
    fMatrix.y -= fTextLeading * fMatrix.a;

    /* set text point to the start of new line. */
    fTextPoint.x = fMatrix.x;
    fTextPoint.y = fMatrix.y;

    *GetStream() << aw << " "
             << ac << " ";
    WriteEscapeText(GetStream(), text);
    *GetStream() << " \"\012";

    /* set text point to the end of this line */
    double w = TextWidth(text);

    if (fFont->WritingMode() == PDF_WMODE_HORIZONTAL) {
        fTextPoint.x += w * fMatrix.a;
        fTextPoint.y += w * fMatrix.b;
    } else {
        fTextPoint.y -= w * fMatrix.a;
        fTextPoint.x -= w * fMatrix.b;
    }

    return;
}

/*--- Color showing --------------------------------------------------------*/

void
PdfContents::SetGrayFill(double gray)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetGrayFill")

    if (gray < 0 || gray > 1)
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
            "ERROR: SetGrayFill out of range(%f)", gray);

    *GetStream() << gray
             << " g\012";

    fGrayFill = gray;
    
    return;
}

void
PdfContents::SetGrayStroke(double gray)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetGrayStroke")

    if (gray < 0 || gray > 1)
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
            "ERROR: SetGrayStroke out of range(%f).", gray);

    *GetStream() << gray
             << " G\012";

    fGrayStroke = gray;
    
    return;
}

void
PdfContents::SetRGBFill(int r,  int g, int b)
{
    SetRGBFill((double)r / 255, (double)g / 255, (double)b / 255);
}

void
PdfContents::SetRGBFill(double r, double g, double b)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetRGBFill")

    if (r < 0 || r > 1 ||
        g < 0 || g > 1 ||
        b < 0 || b > 1) {
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
            "ERROR: SetRGBFill out of range(%f:%f:%f).", r, g, b);
    }

    fRGBFill.red = r;
    fRGBFill.green = g;
    fRGBFill.blue = b;

    *GetStream() << r
             << " "
             << g
             << " "
             << b
             << " rg\012";

    return;
}

void
PdfContents::SetRGBStroke(int r, int g, int b)
{
    SetRGBStroke((double)r / 255, (double)g / 255, (double)b / 255);
}

void
PdfContents::SetRGBStroke(double r, double g, double b)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetRGBStroke")

    if (r < 0 || r > 1 ||
        g < 0 || g > 1 ||
        b < 0 || b > 1) {
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
                "ERROR: SetRGBStroke out of range(%f:%f:%f).", r, g, b);
    }

    fRGBStroke.red = r;
    fRGBStroke.green = g;
    fRGBStroke.blue = b;

    *GetStream() << r
             << " "
             << g
             << " "
             << b
             << " RG\012";

    return;
}

void
PdfContents::SetCMYKFill(double c, double m, double y, double k)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetCMYKFill")

    if (c < 0 || c > 1 ||
        m < 0 || m > 1 ||
        y < 0 || y > 1 ||
        k < 0 || k > 1) {
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
            "ERROR: SetCYMKFill out of range(%f:%f:%f:%f).", c, m, y, k);
    }

    *GetStream() << c
             << " "
             << m
             << " "
             << y
             << " "
             << k
             << " k\012";

    return;
}

void
PdfContents::SetCMYKStroke(double c, double m, double y, double k)
{
    if (fGMode != PDF_GMODE_TEXT_OBJECT && fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("SetCMYKStroke")

    if (c < 0 || c > 1 ||
        m < 0 || m > 1 ||
        y < 0 || y > 1 ||
        k < 0 || k > 1) {
        throw PdfException(PDF_ERR_OUT_OF_RANGE,
            "ERROR: SetCMYKStroke out of range(%f:%f:%f:%f).", c, m, y, k);
    }

    *GetStream() << c
             << " "
             << m
             << " "
             << y
             << " "
             << k
             << " K\012";

    return;
}

/*----------------------------------------------------------------------------*/
/*------ XObjects ------------------------------------------------------------*/

void
PdfContents::ExecuteXObject(const char *name)
{
    assert(fXObjectMgr);

    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("ExecuteXObject")

    /* find x-object from XObjectManager by name */
    PdfXObject* obj = fXObjectMgr->GetXObject(name);
    if (obj == NULL)
        throw PdfException(PDF_RUNTIME_ERROR,
            "ERROR: XObject[%s] is not found.", name);
    
    ExecuteXObject(obj);
}

void
PdfContents::ExecuteXObject(PdfXObject* obj)
{
    assert(fXObjectMgr);

    if (fGMode != PDF_GMODE_PAGE_DESCRIPTION)
        GMODE_ERROR("ExecuteXObject");

    /* check the xobject is valid or not. */
    if (obj == NULL || obj->Name() == NULL) 
        throw PdfException(PDF_RUNTIME_ERROR,
                "ERROR: Invalid XObject.");

    /* get the local-name which is used in page-object */
    const char* local_name = fPage->GetXObjectName(obj);

    if (local_name != NULL)
        WriteEscapeName(GetStream(), local_name);
    else {
        fPage->AddXObject(obj, obj->Name());
        WriteEscapeName(GetStream(), obj->Name());
    }

    *GetStream() << " Do\012";
}

/*----------------------------------------------------------------------------*/

