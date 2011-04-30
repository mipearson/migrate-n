/*
 * << H a r u --free pdf library >> -- DocMaker.cpp
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
 * 2002.10.21 create.
 *
 *     %T<title> <VERSION> <DATE> 
 *                               --- Add front page with specified parameters.
 *     %S<title> <section-no>    --- Add section named <title> with level n. 
 *     %I<title> <image-file-name> <width> <height>
 *                               --- Insert image file.
 *     %W<title> <WMF-file-name> <width> <height>   
 *                               --- Insert WMF file with specified size.
 *     other                     --- Show text at current position.
 *
 */

#include <assert.h>
#include "libharu.h"
#include "libharu_wmf.h"

#define LEFT_MARGIN   85
#define RIGHT_MARGIN  LEFT_MARGIN
#define TOP_MARGIN    50
#define BOTTOM_MARGIN TOP_MARGIN + 10
#define IMAGE_SPACE   20
#define PAGE_WIDTH    530
#define PAGE_HEIGHT   665

const double BASE_RIGHT = PAGE_WIDTH - RIGHT_MARGIN;
const double BASE_TOP = PAGE_HEIGHT - TOP_MARGIN;

static const int PDF_CHAR_SIZE_TITLE0 = 35;
static const int PDF_CHAR_SIZE_TITLE1 = 18;
static const int PDF_CHAR_SIZE_TITLE2 = 10;
static const int PDF_CHAR_SIZE_TITLE3 = 10;
static const int PDF_CHAR_SIZE_DEFAULT = 10;
static const int PDF_CHAR_SIZE_SMALL = 8;

class DocMaker;

class DocIndex
{
    friend class    DocMaker;
public:
                    DocIndex(int level, const char* title, 
                            PdfDestination* desti, int page, int chapter1,
                            int chapter2);
                    ~DocIndex();
private:
    PdfDestination* fDest;
    char*           fTitle;
    int             fPage;
    int             fLevel;
    int             fChapter[2];
};      

class DocMaker
{
public:
                    DocMaker();
                    ~DocMaker();
    void            WriteBackground(PdfContents* canvas);
    PdfPage*        GetCurrent()        { return fContentsPage; }
    PdfPage*        GetCurrentIdx()     { return fIndexPage; }
    void            AddTitle(const char* buf);
    void            AddSection(const char* title, int level);
    void            ShowText(const char* text);
    void            ShowImage(const char* title, const char* filename,
                            double width = 0, double height = 0);
    void            ShowWMFImage(const char* title, const char* filename,
                            double width = 0, double height = 0);
    void            WriteToFile(const char* filename)
                        { fDoc->WriteToFile(filename); }
    const char*     GetStrParam(const char* str, char* param, int len);
    const char*     GetIntParam(const char* str, int* param);
    const char*     GetFloatParam(const char* str, double* param);
    void            MakeIndex();
    void            MakeOutline();
    void            SetDocName(const char* docname);
    void            SetFont(const char* name);
    void            SetIndent(int level);
    void            ShowLabel(const char* text, int indent);
    int             CurLine()           { return fCurLine; }
    void            IncCurLine()        { fCurLine++; }
    void            IncLine();
private:
    void            DrawPageBack(PdfContents* canvas, double top,
                        bool show_pageno = true);
    void            AddIndexPage();
    void            AddPage();
    PdfDoc*         fDoc;
    PdfContents*    fCanvas;
    PdfPages*       fContentsPages;
    PdfPages*       fFrontPages;
    PdfPage*        fFrontPage;
    PdfPages*       fIndexPages;
    PdfPage*        fContentsPage;
    PdfDestination* fContentsDst;
    PdfPage*        fIndexPage;
    int             fSection[3];
    int             fFigureNo;
    int             fImageNo;
    double           fIndexYPos;
    double           fPageYPos;
    int             fXMargine;
    char            fCurFont[64];
    int             fCurFontSize;
    int             fPage;
    char            fTitle[256];
    char*           fDocName;
    int             fImageIdx;
    double           fTmpXPos;
    double           fTmpYPos;
    double           fLeftMargin;
    bool            fAdjustFlg;
    PdfList*        fIndexes;
    int             fCurLine;
};

DocIndex::DocIndex(int level, const char* title, PdfDestination* dest,
        int page, int chapter_1, int chapter_2)           
{
    int len = strlen(title);
    fTitle = new char[len + 1];
    strncpy(fTitle, title, len);
    fTitle[len] = 0x00;
    fDest = dest;
    fPage = page;
    fLevel = level;
    fChapter[0] = chapter_1;
    fChapter[1] = chapter_2;
}

    
DocIndex::~DocIndex()
{
    if (fTitle != NULL)
        delete[] fTitle;
}

DocMaker::DocMaker()
{
    fDoc = new PdfDoc;
    fDoc->NewDoc();

    PDF_DEBUG_PRINT(("DocMaker::DocMaker()\n"));
    
    /* Set PageMode and PageLayout. */ 
    fDoc->Catalog()->SetPageMode(PDF_USE_OUTLINES);
    fDoc->Catalog()->SetPageLayout(PDF_ONE_COLUMN);
    fDoc->Catalog()->AddPageLabel(0, PDF_PAGE_NUM_LOWER_ROMAN);
    fDoc->RootPages()->SetSize(PAGE_WIDTH, PAGE_HEIGHT);

    PdfInfo* info = fDoc->Info();
    info->SetCreator("DocMaker");

    time_t ct = time(NULL);
    struct tm* lst = localtime(&ct);
    pdf_date date = {lst->tm_year + 1900, lst->tm_mon + 1, lst->tm_mday, 
        lst->tm_hour, lst->tm_min, lst->tm_sec, '-', 0, 0};
    
    info->SetCreationDate(date);

    /* Create Pages objects. */
    PDF_DEBUG_PRINT(("DocMaker::DocMaker() --create pages\n"));
    fFrontPages = fDoc->RootPages()->AddPages();
    fIndexPages = fDoc->RootPages()->AddPages();
    fContentsPages = fDoc->RootPages()->AddPages();

    /* Create Font objects. */
    PDF_DEBUG_PRINT(("DocMaker::DocMaker() --create fonts\n"));
    PdfType1FontDef* fd = new PdfHelveticaFontDef();
    PdfEncodingDef* encoding = new PdfWinAnsiEncoding();
    fDoc->AddType1Font(fd, NULL, encoding);

    fd = new PdfHelveticaBoldFontDef();
    fDoc->AddType1Font(fd, NULL, encoding);

    fd = new PdfHelveticaObliqueFontDef();
    fDoc->AddType1Font(fd, NULL, encoding);

    fd = new PdfTimesRomanFontDef();
    fDoc->AddType1Font(fd, NULL, encoding);

    fd = new PdfTimesItalicFontDef();
    fDoc->AddType1Font(fd, NULL, encoding);

    fd = new PdfTimesBoldFontDef();
    fDoc->AddType1Font(fd, NULL, encoding);
    
    fd = new PdfCourierFontDef();
    fDoc->AddType1Font(fd, NULL, encoding);

    fd = new PdfSymbolFontDef();
    fDoc->AddType1Font(fd);

    strncpy(fCurFont, "Times-Roman", 64);
    fCurFontSize = PDF_CHAR_SIZE_DEFAULT;
    fLeftMargin = LEFT_MARGIN;
    fIndexYPos = 70;
    fSection[0] = 0;
    fSection[1] = 0;
    fSection[2] = 0;
    fPage = 0;

    fImageIdx = 0;
    fDocName = NULL;
    fFrontPage = NULL;
    fIndexes = new PdfList();
    fCurLine = 0;
    PDF_DEBUG_PRINT(("DocMaker::DocMaker() --end\n"));
}

DocMaker::~DocMaker()
{
    for (int i = 0; i < fIndexes->CountItems(); i++)
        delete (DocIndex*)fIndexes->ItemAt(i);

    if (fDocName != NULL)
        delete[] fDocName;

    delete fIndexes;
    delete fDoc;
}

void
DocMaker::MakeIndex()
{
    char buf[255];
    
    PDF_DEBUG_PRINT(("DocMaker::MakeIndex()\n"));
    AddIndexPage();

    fContentsDst = new PdfDestination(fIndexPage);
    fContentsDst->SetFit();
    
    fCanvas = fIndexPage->Canvas();
    fCanvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_TITLE1);
    fCanvas->BeginText();

    double tw = fCanvas->TextWidth("Contents"); 
    double xpos = (fCanvas->Width() - LEFT_MARGIN - RIGHT_MARGIN -
            tw) / 2 + LEFT_MARGIN; 
    fCanvas->MoveTextPos(xpos, BASE_TOP - 60);
    fCanvas->ShowText("Contents");
    fCanvas->EndText();

    fCanvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_TITLE3);    
    fCanvas->BeginText(); 
    fCanvas->MoveTextPos(LEFT_MARGIN + 60, BASE_TOP - 80);
     
    fPageYPos = BASE_TOP - 80;

    for (int i = 0; i < fIndexes->CountItems(); i++) {
        DocIndex* idx = (DocIndex*)fIndexes->ItemAt(i);

        if (idx->fLevel == 0) {
            fCanvas->MoveTextPos(0, -25);
            snprintf(buf, 255, "Chapter %d: %s", idx->fChapter[0],
                    idx->fTitle);
            fCanvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_TITLE3);
            PdfLinkAnnot* annot = fIndexPage->AddLink(PdfRect(
                        fCanvas->TextPos().x,
                        fCanvas->TextPos().y,
                        fCanvas->TextPos().x + fCanvas->TextWidth(buf),
                        fCanvas->TextPos().y + PDF_CHAR_SIZE_TITLE3),
                        idx->fDest);
            annot->SetBorder(0, 0, 0);
            fCanvas->ShowText(buf);
            fCanvas->MoveTextPos(0, -5);

            fPageYPos -= 30;
        } else {
            fCanvas->MoveTextPos(-30, -15);
            snprintf(buf, 255, "%d.%d", idx->fChapter[0], idx->fChapter[1]);
            fCanvas->SetFontAndSize("Times-Roman", PDF_CHAR_SIZE_DEFAULT);
            fCanvas->ShowText(buf);
            fCanvas->MoveTextPos(30, 0);
            PdfLinkAnnot* annot = fIndexPage->AddLink(PdfRect(
                        fCanvas->TextPos().x,
                        fCanvas->TextPos().y,
                        fCanvas->TextPos().x + fCanvas->TextWidth(idx->fTitle),
                        fCanvas->TextPos().y + PDF_CHAR_SIZE_DEFAULT),
                        idx->fDest);
            annot->SetBorder(0, 0, 0);
            fCanvas->ShowText(idx->fTitle);

            fPageYPos -= 15;
        }
        if (fPageYPos < BOTTOM_MARGIN + PDF_CHAR_SIZE_TITLE2 * 2) {
            fCanvas->EndText();
            AddIndexPage();
            fCanvas = fIndexPage->Canvas();
            fCanvas->BeginText();
            fCanvas->MoveTextPos(LEFT_MARGIN + 60, BASE_TOP - 60);
            fPageYPos = BASE_TOP - 60;
        }
    }
    fCanvas->EndText();

    unsigned int idx_count = fIndexPages->GetKidsCount() + 1;
    fDoc->Catalog()->AddPageLabel(idx_count, PDF_PAGE_NUM_DECIMAL);
}


void
DocMaker::MakeOutline()
{
    PDF_DEBUG_PRINT(("DocMaker::MakeOutline()\n"));
    PdfOutlineItem* root = new PdfOutlineItem(fDoc->Outlines());
    if (fDocName != NULL)
        root->SetTitle(fDocName);
    else
        root->SetTitle("no title.");
    root->SetOpened(true);
    if (fFrontPage != NULL) {
        PdfDestination* dst = new PdfDestination(fFrontPage);
        root->SetDestination(dst);
        fDoc->Catalog()->SetOpenAction(dst);
    }

    PdfOutlineItem* item = new PdfOutlineItem(root);
    item->SetTitle("Contents");
    item->SetDestination(fContentsDst);
    
    for (int i = 0; i < fIndexes->CountItems(); i++) {
        DocIndex* idx = (DocIndex*)fIndexes->ItemAt(i);

        if (idx->fLevel == 0) {
            char buf[255];
            item = new PdfOutlineItem(root);
            snprintf(buf, 255, "%d %s", idx->fChapter[0], idx->fTitle); 
            item->SetTitle(buf);
            item->SetDestination(idx->fDest);
        } else {
            char buf[255];
            PdfOutlineItem* tmp_item = new PdfOutlineItem(item);
            snprintf(buf, 255, "%d.%d %s", idx->fChapter[0], 
                    idx->fChapter[1], idx->fTitle);
            tmp_item->SetTitle(buf);    
            tmp_item->SetDestination(idx->fDest);
        }
    }
}

void
DocMaker::DrawPageBack(PdfContents* canvas, double top, bool show_pageno)
{
    PDF_DEBUG_PRINT(("DocMaker::DrawPageBack()\n"));
    try {
        double pw = canvas->Width();

        canvas->SetLineWidth(0.5);
        canvas->MoveTo(LEFT_MARGIN, top);
        canvas->LineTo(BASE_RIGHT, top);
        canvas->MoveTo(LEFT_MARGIN, top + 8);
        canvas->LineTo(LEFT_MARGIN, top - 2);
        canvas->MoveTo(BASE_RIGHT, top + 8);
        canvas->LineTo(BASE_RIGHT, top - 2);
        canvas->MoveTo(pw / 2, top + 4);
        canvas->LineTo(pw / 2, top - 2);
        canvas->Stroke();

        /* print page number */
        if (show_pageno && fPage > 0) {
            char tmp_page[10];

            snprintf(tmp_page, 10, "%d", fPage);

            canvas->BeginText();
            canvas->SetFontAndSize("Helvetica", 8);
            double tw = canvas->TextWidth(tmp_page);
            canvas->MoveTextPos((pw - tw) / 2, top + 6);
            canvas->ShowText(tmp_page);

            canvas->EndText();
        }
    } catch (PdfException& e) {
        fprintf(stderr, "%s\n", e.what());
    }

    return;
}

void
DocMaker::AddIndexPage()
{
    PDF_DEBUG_PRINT(("DocMaker::AddIndexPage()\n"));
    fIndexPage = fIndexPages->AddPage();

    PdfContents* canvas = fIndexPage->Canvas();
    canvas->AddFilter(PDF_FILTER_DEFLATE);

    DrawPageBack(canvas, BASE_TOP, false);
}

void
DocMaker::ShowLabel(const char* text, int indent)
{
    PDF_DEBUG_PRINT(("DocMaker::ShowLabel()\n"));
/*  double sava_pos = fPageYPos; */
/*  double sava_fsize = fCanvas->FontSize(); */
    char save_font[64];
    memset(save_font, 0x00, 64);
    strncpy(save_font, fCanvas->FontName(), 63);

    SetFont("Helvetica-Bold");  
    fCurFontSize -= 1;
    SetIndent(0);
    ShowText(text);
    SetFont("Times-Roman");
    fCurFontSize += 1;
    SetIndent(indent);
/*    fPageYPos = sava_pos - 1; */
    fPageYPos += fCanvas->TextLeading();
}

void
DocMaker::IncLine()
{
    fPageYPos -= PDF_CHAR_SIZE_DEFAULT / 2;
}

void
DocMaker::ShowText(const char* text)
{
    PDF_DEBUG_PRINT(("DocMaker::ShowText()\n"));
    if (strlen(text) == 0) {
        fPageYPos -= fCurFontSize;
        return;
    }

    double w = fCanvas->Width() - LEFT_MARGIN - RIGHT_MARGIN;
    const char* tmp_txt = text;

    fCanvas->BeginText();
    fCanvas->SetFontAndSize(fCurFont, fCurFontSize);
    fCanvas->MoveTextPos(fLeftMargin, fPageYPos + (fCurFontSize - 2));

    unsigned int num_char = 0;
    char buf[1024];

    while (num_char < strlen(text)) {

        unsigned int nc = fCanvas->MeasureText(tmp_txt, w);
        assert(nc < 1024);

        if (fPageYPos < BOTTOM_MARGIN) {
            fCanvas->EndText();
            AddPage();
            fCanvas->BeginText();
            fCanvas->SetFontAndSize(fCurFont, fCurFontSize);
            fCanvas->MoveTextPos(fLeftMargin, fPageYPos +
                    (fCurFontSize - 2));
        }

        memset(buf, 0x00, 1024);
        memcpy(buf, tmp_txt, nc);

        /* adjustting text widths with setting character space. */
        if (num_char + nc < strlen(text)) {
            double tw = fCanvas->TextWidth(buf);
            double cs = (w - tw) / nc;
            fCanvas->SetCharSpace(cs);
        }
        fCanvas->ShowTextNextLine(buf);
        fCanvas->SetCharSpace(0);
        fPageYPos -= fCanvas->TextLeading();

        num_char += nc;
        tmp_txt += nc;
    }
    fCanvas->EndText();
    fTmpXPos = fLeftMargin;
    fTmpYPos = 0;
    fAdjustFlg = false;
}

void
DocMaker::AddSection(const char* title, int level)
{
    char buf[255];

    PDF_DEBUG_PRINT(("DocMaker::AddSection()\n"));
    if (level == 0) {
        double tw;
        double xpos;

        fTitle[255] = 0x00;
        strncpy(fTitle, title, 255);
        fSection[0]++;
        fSection[1] = 0;
        fSection[2] = 0;
        fFigureNo = 0;
        fImageNo = 0;
        AddPage();
        fCanvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_TITLE3);
        fCanvas->BeginText();
        snprintf(buf, 255, "CHAPTER %2d", fSection[0]);

        tw = fCanvas->TextWidth(buf);
        xpos = (fCanvas->Width() - LEFT_MARGIN - RIGHT_MARGIN -
                          tw) / 2 + LEFT_MARGIN;
        fCanvas->MoveTextPos(xpos, BASE_TOP - 60);
        fCanvas->ShowText(buf);
        fCanvas->EndText();

        fCanvas->MoveTo(xpos, BASE_TOP - 64);
        fCanvas->LineTo(xpos + tw, BASE_TOP - 64);
        fCanvas->ClosePathStroke();

        fCanvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_TITLE1);
        fCanvas->BeginText();

        tw = fCanvas->TextWidth(title);
        xpos = (fCanvas->Width() - LEFT_MARGIN - RIGHT_MARGIN -
                        tw) / 2 + LEFT_MARGIN;
        fCanvas->MoveTextPos(xpos, BASE_TOP - 105);
        fCanvas->ShowText(title);
        fCanvas->EndText();

        fPageYPos = BASE_TOP - 140;
        fAdjustFlg = false;
    } else if (level == 1) {
        fSection[1]++;
        fSection[2] = 0;

        if (fPageYPos < (BOTTOM_MARGIN + PDF_CHAR_SIZE_TITLE2 +
                    PDF_CHAR_SIZE_DEFAULT))
            AddPage();

        fCanvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_TITLE2);
        fCanvas->BeginText();

        snprintf(buf, 255, "%d.%d ", fSection[0], fSection[1]);
        if (fAdjustFlg == true)
            fPageYPos -= PDF_CHAR_SIZE_TITLE2;
        else
            fPageYPos -= PDF_CHAR_SIZE_TITLE2 + 10;
        fCanvas->MoveTextPos(LEFT_MARGIN - fCanvas->TextWidth(buf), fPageYPos);
        fCanvas->ShowText(buf);
        fCanvas->ShowText(title);
        fPageYPos -= (PDF_CHAR_SIZE_TITLE2 + 5);

        fCanvas->EndText();
        fAdjustFlg = true; 
    } else {
        fSection[2]++;

        if (fPageYPos < (BOTTOM_MARGIN + PDF_CHAR_SIZE_TITLE3 +
                    PDF_CHAR_SIZE_DEFAULT))
            AddPage();

        fCanvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_TITLE3);
        fCanvas->BeginText();

        snprintf(buf, 255, "%d.%d.%d ", fSection[0], fSection[1], fSection[2]);
        if (fAdjustFlg == true)
            fPageYPos -= PDF_CHAR_SIZE_TITLE3 - 5;
        else
            fPageYPos -= PDF_CHAR_SIZE_TITLE3 + 5;
        fCanvas->MoveTextPos(LEFT_MARGIN - fCanvas->TextWidth(buf), fPageYPos);
        fCanvas->ShowText(buf);
        fCanvas->ShowText(title);
        fPageYPos -= (PDF_CHAR_SIZE_TITLE2 + 5);

        fCanvas->EndText();
        fAdjustFlg = true;
    }

    if (level <= 1) {
        PdfDestination* dst = new PdfDestination(fCanvas->Page());
        if (level != 0)
            dst->SetFitH(fPageYPos + 30);
        DocIndex* di = new DocIndex(level, title, dst, 0, fSection[0], 
                fSection[1]);
        fIndexes->AddItem(di);
    }
    SetIndent(0);
}

void
DocMaker::ShowImage(const char* title, const char* filename,
        double width, double height)
{
    char s[1024];
    double img_width;
    double img_height;

    PDF_DEBUG_PRINT(("DocMaker::ShowImage()\n"));
    PdfPngImage* image = new PdfPngImage(fDoc);
    image->AddFilter(PDF_FILTER_DEFLATE);
    
    try {
        image->LoadFromFile(filename);
    } catch (PdfException& e) {
        fprintf(stderr, "ERROR: failed to load image[%s].\n", filename);
        delete image;
        exit(1);
    }
    try {
        fDoc->AddXObject(image);
    } catch (...) {
        delete image;
        return;
    }

    if (width == 0)
        img_width = image->Width();
    else
        img_width = width;

    if (height == 0)
        img_height = image->Height();
    else
        img_height = height;

    if (fCanvas->Width() - RIGHT_MARGIN < fTmpXPos + img_width ||
            fTmpYPos - img_height - 40 < fPageYPos) {
        if (fPageYPos - img_height < BOTTOM_MARGIN + 10)
            AddPage();
        fTmpYPos = fPageYPos;
        fPageYPos -= (img_height + 40);
        fTmpXPos = LEFT_MARGIN;
        fTmpYPos = 0;
    }

    double base_ypos = fPageYPos + 20 + PDF_CHAR_SIZE_SMALL + 4;

    fCanvas->GSave();
    fCanvas->Concat(img_width, 0, 0, img_height, fTmpXPos, base_ypos);
    fCanvas->ExecuteXObject(image);
    fCanvas->GRestore();

    fImageNo++;
    snprintf(s, 1024, "IMAGE %d.%d  ", fSection[0], fImageNo);

    fCanvas->BeginText();
    fCanvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_SMALL);
    fCanvas->MoveTextPos(fTmpXPos, fPageYPos + 20);
    fCanvas->ShowText(s);
    fCanvas->SetFontAndSize("Times-Roman", PDF_CHAR_SIZE_DEFAULT);
    fCanvas->ShowText(title);
    fCanvas->EndText();
    fCanvas->SetFontAndSize(fCurFont, fCurFontSize);

    fTmpXPos += width + IMAGE_SPACE;
}

void
DocMaker::ShowWMFImage(const char* title, const char* filename,
        double width, double height)
{
    char s[1024];
    double img_width;
    double img_height;

    PDF_DEBUG_PRINT(("DocMaker::ShowWMFImage()\n"));
    img_width = width;
    img_height = height;

    if (fCanvas->Width() - RIGHT_MARGIN < fTmpXPos + img_width ||
            fTmpYPos - img_height - 40 < fPageYPos) {
        if (fPageYPos - img_height < BOTTOM_MARGIN + 10)
            AddPage();
        fTmpYPos = fPageYPos;
        fPageYPos -= (img_height + 40);
        fTmpXPos = LEFT_MARGIN;
        fTmpYPos = 0;
    }

    double base_ypos = fPageYPos + 20 + PDF_CHAR_SIZE_SMALL + 4;

    PdfWMFLoader* loader = new PdfWMFLoader();
    pdf_rect rect = PdfRect(fTmpXPos, base_ypos - img_height, 
            fTmpXPos + img_width, base_ypos);
    loader->DrawWMFImage(filename, fCanvas, rect);

    fCanvas->Rectangle(fTmpXPos, base_ypos, width, height);
    fCanvas->Stroke();
    
    fFigureNo++;
    snprintf(s, 1024, "IMAGE %d.%d  ", fSection[0], fFigureNo);

    
    fCanvas->BeginText();
    fCanvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_SMALL);
    fCanvas->MoveTextPos(fTmpXPos, fPageYPos + 20);
    fCanvas->ShowText(s);
    fCanvas->SetFontAndSize("Times-Roman", PDF_CHAR_SIZE_DEFAULT);
    fCanvas->ShowText(title);
    fCanvas->EndText();
    fCanvas->SetFontAndSize(fCurFont, fCurFontSize);

    fTmpXPos += width + IMAGE_SPACE;
}

void
DocMaker::AddPage()
{
    PDF_DEBUG_PRINT(("DocMaker::AddPage()\n"));
    fContentsPage = fContentsPages->AddPage();
    fCanvas = fContentsPage->Canvas();
    fCanvas->SetFontAndSize(fCurFont, fCurFontSize);
    fCanvas->SetTextLeading(PDF_CHAR_SIZE_DEFAULT + 2);
    fCanvas->AddFilter(PDF_FILTER_DEFLATE);

    fPage++;
    WriteBackground(fCanvas);

    fPageYPos = BASE_TOP - 30;
    fTmpXPos = LEFT_MARGIN;
    fTmpYPos = 0;
}

void
DocMaker::WriteBackground(PdfContents* canvas)
{
    PDF_DEBUG_PRINT(("DocMaker::WriteBackground()\n"));
    DrawPageBack(canvas, BASE_TOP);

    canvas->SetFontAndSize("Helvetica", 8);

    /* print title text */
    int tmp_len;
    if (fTitle != NULL) {
        tmp_len = strlen(fTitle);

        const int FIX_LEN = 20;
        char *tmp_title = new char[tmp_len + FIX_LEN + 1];

        try {
            snprintf(tmp_title, tmp_len + FIX_LEN, "%d.%s ", fSection[0], 
                    fTitle);

            canvas->BeginText();
            canvas->SetFontAndSize("Helvetica-Oblique", 8);
            double tw = canvas->TextWidth(tmp_title);
            canvas->MoveTextPos(BASE_RIGHT - tw - 5, BASE_TOP + 4);
            canvas->ShowText(tmp_title);

            canvas->EndText();
        } catch (PdfException& e) {
            fprintf(stderr, "%s\n", e.what());
        }
        delete[] tmp_title;
    }

    return;
}

const char*
DocMaker::GetStrParam(const char* str, char* param, int len)
{
    if (param == NULL)
        return NULL;
    *param = 0x00;

    if (str == NULL)
        return NULL;

    const char* src = str;
    char* dst = param;
    bool intext = false;

    for (int i = 0; i < len; src++){
        if (*src == 0x00) {
            *dst = *src;
            return NULL;
        } else if (*src == '"')
            intext = !intext;
        else if (!intext && 
                (*src == ' ' || *src == 0x0d || *src == 0x0a || *src == 0x09)) {
            if (dst != param) {
                *dst = 0x00;
                if (*src == 0x0d || *src == 0x0a)
                    return NULL;
                else
                    return src++;
            }
        } else  {
            *dst++ = *src;
        }
    }

    return NULL;
}

const char*
DocMaker::GetIntParam(const char* str, int* param)
{
    char sparam[11];

    str = GetStrParam(str, sparam, 11);
    int i = atoi(sparam);
    *param = i;

    return str;
}

const char*
DocMaker::GetFloatParam(const char* str, double* param)
{
    char sparam[11];

    str = GetStrParam(str, sparam, 11);
    double f = atof(sparam);
    *param = f;

    return str;
}

void
DocMaker::SetDocName(const char* docname)
{
    if (fDocName != NULL)
        delete[] fDocName;

    int len = strlen(docname) + 1;
    fDocName = new char[len];
    strcpy(fDocName, docname);
}

void
DocMaker::AddTitle(const char* buf)
{
    char tmpbuf[255];

    PDF_DEBUG_PRINT(("DocMaker::AddTitle(%s)\n", buf));
    fFrontPage = fFrontPages->AddPage();
    PdfContents* canvas = fFrontPage->Canvas();
    canvas->AddFilter(PDF_FILTER_DEFLATE);

    /* Print the document title. */
    buf = GetStrParam(buf, tmpbuf, 255);
    if (tmpbuf == NULL) 
        return;
    
    SetDocName(tmpbuf);
    canvas->BeginText();
    canvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_TITLE0);
    canvas->MoveTextPos(75, 500);
    int numchars = 0;
    double width = canvas->TextWidth(tmpbuf, &numchars);
    double tmpw = canvas->Width() - 75 * 2 - width;
    canvas->SetCharSpace(tmpw / numchars);
    
    canvas->ShowText(tmpbuf);
    canvas->EndText();
    canvas->SetCharSpace(0);

    /* Print the subtitle. */
    buf = GetStrParam(buf, tmpbuf, 255);
    if (tmpbuf == NULL)
        return;
    
    canvas->BeginText();
    canvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_TITLE1);
    width = canvas->TextWidth(tmpbuf, &numchars);
    double tmpl = canvas->Width() - width - 80;
    canvas->MoveTextPos(tmpl, 470);
    
    canvas->ShowText(tmpbuf);
    canvas->EndText();

    /* Print the version NO. */
    buf = GetStrParam(buf, tmpbuf, 255);
    if (tmpbuf == NULL)
        return;
    
    canvas->BeginText();
    canvas->SetFontAndSize("Helvetica-Bold", PDF_CHAR_SIZE_TITLE2);
    width = canvas->TextWidth(tmpbuf, &numchars);
    
    tmpl = canvas->Width() - width - 80;
    canvas->MoveTextPos(tmpl, 100);
    
    canvas->ShowText(tmpbuf);
    canvas->EndText();

    /* Print the copyright */
    GetStrParam(buf, tmpbuf, 255);
    if (tmpbuf == NULL)
        return;

    canvas->BeginText();
    canvas->SetFontAndSize("Times-Bold", PDF_CHAR_SIZE_TITLE2);
    width = canvas->TextWidth(tmpbuf, &numchars);
    tmpl = canvas->Width() - width - 80;
    canvas->MoveTextPos(tmpl, 80);
    
    canvas->ShowText(tmpbuf);
    canvas->EndText();

}

void
DocMaker::SetFont(const char* name)
{
    double fsize = fCanvas->FontSize();
    const char* fname = fCanvas->FontName();

    if (strcmp(fname, name) != 0) 
        fCanvas->SetFontAndSize(name, fsize);
    memset(fCurFont, 0x00, 64);
    strncpy(fCurFont, name, 64 - 1);
}

void 
DocMaker::SetIndent(int level)
{
    if (level == 0)
        fLeftMargin = LEFT_MARGIN;
    else
        fLeftMargin += level * 10;
}

int main(int argc, char** argv)
{
    FILE* infile;
    const int BUF_LEN = 8192;
    char buf[BUF_LEN];
    char* filename = argv[1];
    char* pdffname = argv[2];

    if (argc < 3) {
        fprintf(stderr, "usage: DocMaker input-file-name output-file-name.\n");
        return -1;
    }
    infile = fopen(filename, "r");
    if (infile == NULL) {
        fprintf(stderr, "error: file open failed.\n");
        return -1;
    }

    DocMaker* doc = NULL;
    try {   
        doc = new DocMaker();
    } catch (PdfException& e) {
        fprintf(stderr, "%s\n", e.what());
        fclose(infile);
        return 1;
    }

    try {
        while (fgets(buf, BUF_LEN, infile) != NULL) {
            unsigned int len = strlen(buf);
            if (buf[len - 1] == 0x0d || buf[len - 1] == 0x0a) {
                buf[len - 1] = 0x00;
                len--;
            }
            if (buf[len - 1] == 0x0d || buf[len - 1] == 0x0a) {
                buf[len - 1] = 0x00;
                len --;
            }
            doc->IncCurLine();

            if (len >= 2) {

                /* insert a brank line */
                if (memcmp(buf, "%;", 2) == 0) 
                    doc->ShowText("");
                else

                /* insert a brank line (harf)*/
                if (memcmp(buf, "%:", 2) == 0) 
                    doc->IncLine();
                else
                    
                /* comment lines */
                if (memcmp(buf, "%#", 2) == 0)
                    ;
                else
                
                /* add front page */
                if (memcmp(buf, "%T", 2) == 0) {
                    char* pbuf = buf + 2;
                    doc->AddTitle(pbuf);
                } else
            
                /* add section. */
                if (memcmp(buf, "%S", 2) == 0) {
                    int level;
                    char title[256];
                    const char* pbuf = buf + 2;
                
                    pbuf = doc->GetStrParam(pbuf, title, 256);
                    doc->GetIntParam(pbuf, &level);
                
                    try {
                        doc->AddSection(title, level);
                    } catch (PdfException& e) {
                        fprintf(stderr, e.what());
                        exit(2);
                    }
                } else

                /* show image. */
                if (memcmp(buf, "%I", 2) == 0) {
                    char title[256];
                    char filename[256];
                    const char* pbuf = buf + 2;

                    pbuf = doc->GetStrParam(pbuf, title, 256);
                    pbuf = doc->GetStrParam(pbuf, filename, 256);

                    int width = 0;
                    int height = 0;
                    if (pbuf != NULL)
                        pbuf = doc->GetIntParam(pbuf, &width);
                    if (pbuf != NULL)
                        doc->GetIntParam(pbuf, &height);
                
                    doc->ShowImage(title, filename, width, height);

                } else
                    
                /*  insert figure */
                if (memcmp(buf, "%W", 2) == 0 && len > 8) {
                    char title[256];
                    char filename[256];
                    double height;
                    double width;
                    const char* pbuf = buf + 2;

                    pbuf = doc->GetStrParam(pbuf, title, 256);
                    pbuf = doc->GetStrParam(pbuf, filename, 256);
                    pbuf = doc->GetFloatParam(pbuf, &width);
                    pbuf = doc->GetFloatParam(pbuf, &height);
                     
                    doc->ShowWMFImage(title, filename, width, height);
                } else 

                /* change default font */

                if (memcmp(buf, "%C", 2) == 0) {
                    char font_name[64];
                    const char* pbuf = buf + 2;
            
                    doc->GetStrParam(pbuf, font_name, 64);
                    doc->SetFont(font_name);
                } else

                /* set indent */
                if (memcmp(buf, "%L", 2) == 0) {
                    int indent;
                    const char* pbuf = buf + 2;
            
                    doc->GetIntParam(pbuf, &indent);
                    doc->SetIndent(indent);
                } else

                /* label */
                if (memcmp(buf, "%A", 2) == 0) {
                    int indent;
                    char label[128];
                    const char* pbuf = buf + 2;

                    pbuf = doc->GetIntParam(pbuf, &indent);
                    doc->GetStrParam(pbuf, label, 128);
                    doc->ShowLabel(label, indent);
                } else
                    
                /* print text */
                    doc->ShowText(buf);
            }
        }
        doc->MakeIndex();
        doc->MakeOutline();
        doc->WriteToFile(pdffname);
    } catch (PdfException& e) {
        fprintf(stderr, "ERROR: %s line:%d.\n", e.what(), doc->CurLine());
        delete doc;
        fclose(infile);
        exit(1);
    }

    fclose(infile);

    delete doc;
}

