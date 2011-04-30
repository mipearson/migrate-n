/*
 * << H a r u -- Free PDF Library >> -- PdfWMFLoader.cc
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
 * This is an experimental implementation.
 * Therefore, many features of a WMF format are not implemented yet.
 *
 */

#include <errno.h>
#include "libharu_wmf.h"
#include "libharu.h"

PdfWMFLoader::PdfWMFLoader()
{
    fFile = NULL;
    fXOrg = 0;
    fYOrg = 0;
    fCnt = 0;
    fPenStyle = WMF_PEN_STYLE_PS_SOLID;
    fPolyFillMode = WMF_POLY_FILL_MODE_ALTERNATE;
    fTextColor[0] = 0;
    fTextColor[1] = 0;
    fTextColor[2] = 0;
    fXScale = 0;
    fYScale = 0;
    fXOffset = 0;
    fYOffset = 0;
}

void
PdfWMFLoader::DrawWMFImage(const char* filename, PdfContents* canvas,
        pdf_rect rect)
{   
    if (canvas->GMode() != PDF_GMODE_PAGE_DESCRIPTION)
        throw PdfException(PDF_ERR_INVALID_GRAPHIC_MODE, 
                "ERROR: Invalid graphics mode %d.", canvas->GMode());
    
    fFile = fopen(filename, "rb");
    if (fFile == NULL)
        throw PdfException(errno, "ERROR: cannot open %s.", filename);

    try {
        LoadMetafileHeader();

        double width = rect.right - rect.left;
        double height = abs((int)(rect.bottom - rect.top));
        double image_width = fHeader.Right - fHeader.Left;
        double image_height = fHeader.Bottom - fHeader.Top;
       
        fprintf(stderr, "%f,%f,%f,%f\n",width,height,image_width,image_height);
        fXScale = width / image_width;
        fYScale = height / image_height;
        fXOffset = rect.left;
        if (rect.bottom > rect.top)
            fYOffset = canvas->Height() - rect.bottom;
        else
            fYOffset = canvas->Height() - rect.top;
        fprintf(stderr, "%f,%f,%f,%f\n",fXScale, fYScale, fXOffset, fYOffset);
        canvas->SetLineCap(PDF_PROJECTING_SCUARE_END);

        while (DrawMetafileRecord(canvas));
    } catch (...) {
        fclose(fFile);
        throw;
    }
    fclose(fFile);
}

PLACEABLEMETAHEADER 
PdfWMFLoader::LoadMetafileHeader(const char* filename)
{
    fFile = fopen(filename, "rb");
    if (fFile == NULL)
        throw PdfException(errno, "ERROR: cannot open %s.", filename);

    try {
        LoadMetafileHeader();
    } catch (...) {
        fclose(fFile);
        throw;
    }
    fclose(fFile);
    return fHeader;
}

void 
PdfWMFLoader::LoadMetafileHeader()
{
    const unsigned long WMF_HEADER_KEY = 0x9AC6CDD7;

    if (fread(&fHeader, 22, 1, fFile) != 1)
        throw PdfException(ferror(fFile), "ERROR: cannot read from file.");
    else if (fHeader.Key != WMF_HEADER_KEY || fHeader.Handle != 0 ||
            fHeader.Reserved != 0) {
            throw PdfException (PDF_RUNTIME_ERROR, 
                    "ERROR invalid meta-file fHeader %04X %02X.\n",
                    (unsigned int)fHeader.Key,
                    (unsigned int)fHeader.Handle);
    } else {
        PDF_DEBUG_PRINT(("Header Check OK.\n"));
        PDF_DEBUG_PRINT(("Left=%d\n", fHeader.Left));
        PDF_DEBUG_PRINT(("Top=%d\n", fHeader.Top));
        PDF_DEBUG_PRINT(("Right=%d\n", fHeader.Right));
        PDF_DEBUG_PRINT(("Bottom=%d\n", fHeader.Bottom));
        PDF_DEBUG_PRINT(("Inch=%u\n\n", fHeader.Inch));

        unsigned short chk = 0;
        unsigned short* ptr = (unsigned short*)&fHeader;
        for (int i = 0; i < 10; i++, ptr++)
            chk ^= *ptr;

        if (chk != fHeader.Checksum)
            throw PdfException(PDF_ERR_INVALID_WMF_FILE, 
                "Invalid MetaFile Header %X <> %X.\n",
                fHeader.Checksum, chk);
    }

    if (fread(&fWMFHeader.FileType, 2, 1, fFile) != 1 ||
            fread(&fWMFHeader.HeaderSize, 2, 1, fFile) != 1 ||
            fread(&fWMFHeader.Version, 2, 1, fFile) != 1 ||
            fread(&fWMFHeader.FileSize, 4, 1, fFile) != 1 ||
            fread(&fWMFHeader.NumOfObjects, 2, 1, fFile) != 1 ||
            fread(&fWMFHeader.MaxRecordSize, 4, 1, fFile) != 1 ||
            fread(&fWMFHeader.NumOfParams, 2, 1, fFile) != 1) {
        throw PdfException(ferror(fFile), "ERROR: cannot read from file.");
    } else {
        PDF_DEBUG_PRINT(("FileType=%u.\n", (unsigned int)fWMFHeader.FileType));
        PDF_DEBUG_PRINT(("Header Size=%u\n", 
                    (unsigned int)fWMFHeader.HeaderSize));
        PDF_DEBUG_PRINT(("Version=0x%04X\n", 
                    (unsigned int)fWMFHeader.Version));
        PDF_DEBUG_PRINT(("FileSize=%u\n", 
                    (unsigned int)fWMFHeader.FileSize));
        PDF_DEBUG_PRINT(("NumOfObjects=%u\n", 
                    (unsigned int)fWMFHeader.NumOfObjects));
        PDF_DEBUG_PRINT(("MaxRecordSize=%u\n", 
                    (unsigned int)fWMFHeader.MaxRecordSize));
        PDF_DEBUG_PRINT(("NumOfParams=0x%04X\n", fWMFHeader.NumOfParams));
        if (fWMFHeader.NumOfParams != 0) 
            throw PdfException (PDF_ERR_INVALID_WMF_FILE, 
                    "ERROR: Invalid NumOfParams.");
    }

    fCnt = 0;
}

bool
PdfWMFLoader::DrawMetafileRecord(PdfContents* canvas)
{
    unsigned long size;
    unsigned short func;
    short param;

    ReadFile(&size, sizeof(size));
    ReadFile(&func, sizeof(func));

    fCnt++;
    PDF_DEBUG_PRINT(("[param%d]\n", fCnt));
    PDF_DEBUG_PRINT(("\tsize=%d\n", size));
    PDF_DEBUG_PRINT(("\tfunc=0x%04X\n", func));
    
    if (size < 3)
        throw PdfException(PDF_ERR_INVALID_WMF_FILE,
                "ERROR: Invalid metafile record-size(%d)", size);
    size -= 3;

    if (func == 0x0000 && size == 0)
        return false;

    switch(func) {
        case 0x020B:
            size -= SetWindowOrg(canvas);
            break;
        case 0x02FC:
            size -= CreateBrushIndirect(canvas);
            break;
        case 0x02FA:
            size -= CreatePenIndirect(canvas);
            break;
        case 0x0324:
            size -= Polygon(canvas, size);
            break;
        case 0x0325:
            size -= PolyLine(canvas, size);
            break;
        case 0x0538:
            size -= PolyPolygon(canvas);
            break;
        case 0x041B:
            size -= Rectangle(canvas);
            break;
        case 0x0106:
            size -= SetPolyFillMode(canvas);
            break;
        case 0x0521:
            size -= TextOut(canvas, size);
            break;
        case 0x02FB:
            size -= CreateFontIndirect(canvas);
            break;
        case 0x0209:
            size -= SetTextColor(canvas);
            break;
        default:
            PDF_DEBUG_PRINT(("\tunsupported function=0x%04X\n", func));
    }
    
    while (size > 0) {
        ReadFile(&param, sizeof(param));
        PDF_DEBUG_PRINT(("\tvoid param%d\n", param)); 
        size--;
    }
    return true;
}

unsigned int
PdfWMFLoader::SetWindowOrg(PdfContents* canvas)
{
    ReadFile(&fYOrg, sizeof(fYOrg));
    ReadFile(&fXOrg, sizeof(fXOrg));
    PDF_DEBUG_PRINT(("\tSetWindowOrg(x,y) = (%d,%d)\n", fYOrg, fXOrg));

    return 2;
}

unsigned int 
PdfWMFLoader::CreateBrushIndirect(PdfContents* canvas)
{
    unsigned char rgb[4];
    unsigned short param;

    ReadFile(&param, sizeof(param));
    PDF_DEBUG_PRINT(("\tvoid param%d\n", param));
    ReadFile(rgb, 4);
    
    PDF_DEBUG_PRINT(("\tSetRGBFill(%d, %d, %d)\n", rgb[0], rgb[1], rgb[2]));
    canvas->SetRGBFill(rgb[0], rgb[1], rgb[2]);

    return 3;
}

unsigned int 
PdfWMFLoader::CreatePenIndirect(PdfContents* canvas)
{
    unsigned char rgb[4];
    unsigned short param;

    ReadFile(&fPenStyle, sizeof(fPenStyle));

    ReadFile(&param, sizeof(param));
    PDF_DEBUG_PRINT(("\tSetLineWidth(%d)\n", param));
    canvas->SetLineWidth(param);

    ReadFile(&param, sizeof(param));
    PDF_DEBUG_PRINT(("\tvoid param%d\n", param));

    ReadFile(&rgb, sizeof(rgb));
    PDF_DEBUG_PRINT(("\tSetRGBStroke(%d, %d, %d)\n", rgb[0], rgb[1], rgb[2]));
    canvas->SetRGBStroke(rgb[0], rgb[1], rgb[2]);

    return 5;
}

unsigned int 
PdfWMFLoader::Polygon(PdfContents* canvas, unsigned int size)
{
    unsigned short param;
    short point[2];

    ReadFile(&param, sizeof(param));
    int cnt = 1;

    ReadFile(point, sizeof(point));
    cnt += 2;
    PDF_DEBUG_PRINT(("\tMoveTo(%f, %f)\n", TranslateX(canvas, point[0]),
        TranslateY(canvas, point[1])));
    canvas->MoveTo(TranslateX(canvas, point[0]), TranslateY(canvas, point[1]));
    size -= 3;

    for (unsigned int i = 0; i < size / 2; i++) {
        ReadFile(point, sizeof(point));
        cnt += 2;
        PDF_DEBUG_PRINT(("\tLineTo(%f, %f)\n", TranslateX(canvas, point[0]),
            TranslateY(canvas, point[1])));
        canvas->LineTo(TranslateX(canvas, point[0]), 
            TranslateY(canvas, point[1]));
    }

    if (fPenStyle == WMF_PEN_STYLE_PS_NULL) 
        if (fPolyFillMode == WMF_POLY_FILL_MODE_ALTERNATE)
            canvas->Eofill();
        else
            canvas->Fill();
    else
        if (fPolyFillMode == WMF_POLY_FILL_MODE_ALTERNATE)
            canvas->EofillStroke();
        else
            canvas->FillStroke();

    return cnt;
}

unsigned int 
PdfWMFLoader::PolyLine(PdfContents* canvas, unsigned int size)
{
    unsigned short param;
    short point[2];

    ReadFile(&param, sizeof(param));
    int cnt = 1;

    ReadFile(point, sizeof(point));
    cnt += 2;
    PDF_DEBUG_PRINT(("\tMoveTo(%d, %d)\n", TranslateX(canvas, point[0]),
        TranslateY(canvas, point[1])));
    canvas->MoveTo(TranslateX(canvas, point[0]), TranslateY(canvas, point[1]));
    size -= 3;
    
    for (unsigned int i = 0; i < size / 2; i++) {
        ReadFile(point, sizeof(point));
        cnt += 2;
        PDF_DEBUG_PRINT(("\tLineTo(%d, %d)\n", TranslateX(canvas, point[0]),
            TranslateY(canvas, point[1])));
        canvas->LineTo(TranslateX(canvas, point[0]), 
            TranslateY(canvas, point[1]));
    }

    canvas->Stroke();

    return cnt;
}

unsigned int 
PdfWMFLoader::PolyPolygon(PdfContents* canvas)
{
    unsigned short param;
    short* params;

    ReadFile(&param, sizeof(param));
    int cnt = 1;

    if (param > 32767)
        throw new PdfException(PDF_ERR_INVALID_WMF_FILE, 
                "ERROR: PolyPolygon --invalid num of params");
    
    params = new short[param];
    ReadFile(params, sizeof(short) * param);
    PDF_DEBUG_PRINT(("\tPolyPolygon params-count=%d\n", param));
    cnt += param;

    for (unsigned int i = 0; i < param; i++) {
        for (int j = 0; j < params[i]; j++) {
            short pt[2];

            ReadFile(pt, sizeof(short) * 2);
            PDF_DEBUG_PRINT(("\tPolyPolygon param[%d]=%d,%d\n", pt[0], pt[1]));
            if (j == 0)
               canvas->MoveTo(TranslateX(canvas, pt[0]),
                        TranslateY(canvas, pt[1]));
            else
               canvas->LineTo(TranslateX(canvas, pt[0]),
                        TranslateY(canvas, pt[1]));
            cnt += 2;
        }
    }
    if (fPenStyle == WMF_PEN_STYLE_PS_NULL) 
        if (fPolyFillMode == WMF_POLY_FILL_MODE_ALTERNATE)
            canvas->Eofill();
        else
            canvas->Fill();
    else
        if (fPolyFillMode == WMF_POLY_FILL_MODE_ALTERNATE)
            canvas->EofillStroke();
        else
            canvas->FillStroke();

    return cnt;
}

unsigned int 
PdfWMFLoader::Rectangle(PdfContents* canvas)
{
    short rect[4];

    ReadFile(rect, sizeof(short) * 4);
    PDF_DEBUG_PRINT(("\tRectangle param=[%d,%d,%d,%d]\n",
                rect[0], rect[1], rect[2], rect[3]));

    canvas->Rectangle(rect[0], rect[1], rect[2], rect[3]);
    if (fPolyFillMode == WMF_POLY_FILL_MODE_ALTERNATE)
        canvas->EofillStroke();
    else
        canvas->FillStroke();

    return 4;
}

unsigned int 
PdfWMFLoader::SetPolyFillMode(PdfContents* canvas)
{
    ReadFile(&fPolyFillMode, sizeof(short));
    PDF_DEBUG_PRINT(("\tSetPolyFillMode=%d\n", fPolyFillMode));

    return 1;
}

unsigned int
PdfWMFLoader::TextOut(PdfContents* canvas, unsigned int size)
{
    unsigned short para;
    ReadFile(&para, sizeof(short));

    unsigned int siz = (size - 3) * 2;
    short x;
    short y;
    char* buf = new char[siz + 1];
    ReadFile(buf, siz);
    buf[siz] = 0x00;
    PDF_DEBUG_PRINT(("\tTextOut=%s\n", buf));
    
    ReadFile(&x, sizeof(short));
    ReadFile(&y, sizeof(short));
    PDF_DEBUG_PRINT(("\tx=%f, y=%f\n", TranslateX(canvas, x), 
                TranslateY(canvas, y)));

    canvas->SetRGBFill(0, 0, 0);
    
    canvas->BeginText();
    canvas->SetFontAndSize("Courier", 80 * fYScale);
    canvas->MoveTextPos(TranslateX(canvas, x), TranslateY(canvas, y)
    + 14);
    canvas->ShowText(buf);
    canvas->EndText();

    return size;
}

unsigned int
PdfWMFLoader::SetTextColor(PdfContents* canvas)
{
    unsigned char rgb[4];
    ReadFile(rgb, sizeof(char) * 4);
    fTextColor[0] = rgb[0]; 
    fTextColor[1] = rgb[1]; 
    fTextColor[2] = rgb[2];
    PDF_DEBUG_PRINT(("\tfTextColor=[%d,%d,%d]\n", rgb[1], rgb[2], rgb[3]));
    return 2;
}

unsigned int
PdfWMFLoader::CreateFontIndirect(PdfContents* canvas)
{
    unsigned char siz[2];
    ReadFile(siz, 2);
    PDF_DEBUG_PRINT(("\tCreateFontIndirect size=%d%d\n", siz[0], siz[1]));
    return 1;
}
