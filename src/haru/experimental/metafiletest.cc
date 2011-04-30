#include <stdio.h>
#include "libharu_wmf.h"

int main(int argc, char** argv)
{
    double dheight;
    double dwidth;

        
    PdfDoc* doc = new PdfDoc();
    try {
        doc->NewDoc();
        PdfPage *page = doc->AddPage();
        page->SetSize(500, 500);
        PdfContents *canvas = page->Canvas();
        doc->AddType1Font(new PdfCourierFontDef());
        canvas->SetFontAndSize("Courier", 10);

        PdfWMFLoader* loader = new PdfWMFLoader();
        PLACEABLEMETAHEADER header = loader->LoadMetafileHeader(argv[1]);
        int height = abs(header.Top - header.Bottom);
        int width = abs(header.Right - header.Left);
    
        if (height > width) {
            dheight = 400;
            dwidth = dheight / height * width;
        } else {
            dwidth = 400;
            dheight = dwidth / width * height;
        }
        loader->DrawWMFImage(argv[1], canvas, PdfRect(50, 50, 
                    50 + dwidth, 50 + dheight));

        delete loader;

        doc->WriteToFile("test.pdf");

    } catch (std::exception& e) {
        fprintf(stderr, "%s.\n", e.what());
    }
    delete doc;
}

