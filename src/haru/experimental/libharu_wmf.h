/*
 * << H a r u -- Free PDF Library >> -- libharu_wmf.h
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
 * The opewrators which are supported this program.
 * 1.CreateBrushIndirect
 * 2.CreatePenIndirect
 * 3.SetPolyFillMode
 * 4.Polygon
 * 5.PolyPolygon
 * 6.Polyline
 * 7.SetWindowOrg
 *
 * NOTE:  
 *   Binary raster operation(SetROP2) is not supported. So the WMF file which 
 *   uses SetROP2 operator is not convert collectlly.
 *
 */

#ifndef _LIB_HARU_WMF_H 
#define _LIB_HARU_WMF_H 

#include "libharu.h"

typedef struct _PlaceableMetaHeader
{
    unsigned long Key;
    unsigned short Handle;
    short Left;
    short Top;
    short Right;
    short Bottom;
    unsigned short Inch;
    unsigned long Reserved;
    unsigned short Checksum;
} PLACEABLEMETAHEADER;

typedef struct _WindowsMetaHeader
{
    unsigned short FileType;
    unsigned short HeaderSize;
    unsigned short Version;
    unsigned long FileSize;
    unsigned short NumOfObjects;
    unsigned long MaxRecordSize;
    unsigned short NumOfParams;
} WMFHEAD;

#define WMF_POLY_FILL_MODE_ALTERNATE 1
#define WMF_POLY_FILL_MODE_WINDING   2

#define WMF_BK_MODE_TRANSPARENT 1
#define WMF_BK_MODE_OPAQUE 2

#define WMF_PEN_STYLE_PS_SOLID      0
#define WMF_PEN_STYLE_PS_PS_DASH    1
#define WMF_PEN_STYLE_PS_NULL       5

class PdfWMFLoader : public PdfAutoPtrObject
{
public:
                        PdfWMFLoader();
    PLACEABLEMETAHEADER LoadMetafileHeader(const char* filename);
    void                DrawWMFImage(const char* filename, 
                            PdfContents* canvas, pdf_rect rect);
private:
    FILE*               fFile;
    short               fXOrg;
    short               fYOrg;
    unsigned int        fCnt;
    double              fXScale;
    double              fYScale;
    double              fXOffset;
    double              fYOffset;
    unsigned short      fPenStyle;
    unsigned short      fPolyFillMode;
    unsigned char       fTextColor[3];
    
    PLACEABLEMETAHEADER fHeader;
    WMFHEAD             fWMFHeader;
    void                LoadMetafileHeader();
    bool                DrawMetafileRecord(PdfContents* canvas);
    size_t              ReadFile(void* ptr, size_t siz) {
                            if (fread(ptr, siz, 1, fFile) != 1) {
                                int err = ferror(fFile);
                                throw PdfException(err, "ERROR: cannot "
                                        "read from file(%d).", err);
                            }
                            return siz;
                        };
    double              TranslateY(PdfContents* canvas, short y) {
                            return (fHeader.Bottom - fHeader.Top - 
                                abs(y - fYOrg)) * fYScale + fYOffset;
                        };
    double              TranslateX(PdfContents* canvas, short x) {
                            return (abs(x - fXOrg)) * fXScale + fXOffset;
                        };
    unsigned int        SetWindowOrg(PdfContents* canvas);
    unsigned int        CreateBrushIndirect(PdfContents* canvas);
    unsigned int        CreatePenIndirect(PdfContents* canvas);
    unsigned int        PolyLine(PdfContents* canvas, unsigned int size);
    unsigned int        Polygon(PdfContents* canvas, unsigned int size);
    unsigned int        PolyPolygon(PdfContents* canvas);
    unsigned int        Rectangle(PdfContents* canvas);
    unsigned int        SetPolyFillMode(PdfContents* canvas);
    unsigned int        TextOut(PdfContents* canvas, unsigned int size);
    unsigned int        CreateFontIndirect(PdfContents* canvas);
    unsigned int        SetTextColor(PdfContents* canvas);
};

#endif /* _LIB_HARU_WMF_H */

