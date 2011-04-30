/*
 * << H a r u --free pdf library >> -- PdfCMap.cc
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

#include "libharu.h"

static const unsigned char UNICODE_HEADER[PDF_UNICODE_HEADER_LEN] = {
    0xFE, 0xFF
};

/*------ PdfCMap -------------------------------------------------------------*/

PdfCMap::PdfCMap()
    : PdfAutoPtrObject()
{
    PDF_DEBUG_PRINT(("++ [%x] PdfCMap new[].\n", (int)this));
    fRegistry = NULL;
    fOrdering = NULL;
    fSupplement = 0;

    for (int i = 0; i <= 255; i++) { 
        for (int j = 0; j <= 255; j++) {
            /* undefined charactors are replaced to square */
            fUnicodeArray[i][j] = 0x25A1;
            fCMapArray[i][j] = 0x0000;
        }
    }
}

PdfCMap::~PdfCMap()
{
    PDF_DEBUG_PRINT(("++ [%x] fRegistry delete.\n", (int)fRegistry));
    delete fRegistry;
    
    PDF_DEBUG_PRINT(("++ [%x] fOrdering delete.\n", (int)fOrdering));
    delete fOrdering;

    PDF_DEBUG_PRINT(("++ [%x] PdfCMap delete[].\n", (int)this));
}

unsigned int
PdfCMap::GetUnicode(unsigned int code)
{
    unsigned char l = code;
    unsigned char h = code >> 8;

    return fUnicodeArray[l][h];
}

pdf_cid
PdfCMap::GetCID(unsigned int code)
{
    unsigned char l = code;
    unsigned char h = code >> 8;
          
    return fCMapArray[l][h];
}

void
PdfCMap::AddCMap(const pdf_cid_range* range)
{
    const pdf_cid_range* prange = range;
    
    /* Copy specified pdf_cid_range array to fRangeArray. */
    while (prange->from != 0xffff && prange->to != 0xffff) {
        unsigned short code = prange->from;
        unsigned short cid = prange->cid;

        while (code <= prange->to) {
            unsigned char l = code;
            unsigned char h = code >> 8;

            fCMapArray[l][h] = cid;
            code++;
            cid++;
        }
        prange++;
    }
}

void
PdfCMap::SetCIDSystemInfo(const char* registry,
        const char* ordering, unsigned int supplement)
{
    if (fRegistry != NULL) {
        PDF_DEBUG_PRINT(("++ [%x] fRegistry delete.\n", (int)fRegistry));
        delete fRegistry;
        fRegistry = NULL;
    }
    if (fOrdering != NULL) {
        PDF_DEBUG_PRINT(("++ [%x] fOrdering delete.\n", (int)fOrdering));
        delete fOrdering;
        fOrdering = NULL;
    }
    fSupplement = 0;

    if (registry != NULL) {
        fRegistry = new char[strlen(registry + 1)];
        PDF_DEBUG_PRINT(("++ [%x] fRegistry new.\n", (int)fRegistry));
        strcpy(fRegistry, registry);
    }
    if (ordering != NULL) {
        fOrdering = new char[strlen(ordering + 1)];   
        PDF_DEBUG_PRINT(("++ [%x] fOrdering new.\n", (int)fOrdering));
        strcpy(fOrdering, ordering);
    }
    fSupplement = supplement;
}

void
PdfCMap::AddCIDSystemInfo(PdfCIDFont* font)
{
    if (fRegistry == NULL || fOrdering == NULL) {
        PDF_DEBUG_PRINT(("WARNING PdfCMap::AddCIDSystemInfo "
                    "CIDSystemInfo is not prepared.\n"));
        return;
    }
    
    PdfDictionary* dict = new PdfDictionary(font->GetXref());
    font->AddElement("CIDSystemInfo", dict);
    dict->AddElement("Registry", new PdfText(fRegistry));
    dict->AddElement("Ordering", new PdfText(fOrdering));
    dict->AddElement("Supplement", new PdfNumber(fSupplement));
}

void
PdfCMap::SetUnicodeArray(const pdf_mb_unicode_map1* array1, 
        const pdf_mb_unicode_map2* array2)
{
    if (array1 != NULL)
        while (array1->unicode != 0xffff) {
            unsigned char l = array1->mbchar;
            unsigned char h = array1->mbchar >> 8;
            fUnicodeArray[l][h] = array1->unicode;
            array1++;
        }

    if (array2 != NULL) 
        while (array2->unicode != 0xffff) {
            unsigned short mbchar = array2->from;
            unsigned short unicode = array2->unicode;
            while (mbchar <= array2->to) {
                unsigned char l = mbchar;
                unsigned char h = mbchar >> 8;
                fUnicodeArray[l][h] = unicode;
                mbchar++;
                unicode++;
            }
            array2++;
        }
}

int
PdfCMap::ToUnicode(const char* src, unsigned char* dst, int* len)
{
    int srclen = strlen(src);
    int dstlen;
    int j = 0;

    if (dst != NULL && len != NULL)
        dstlen = *len - 1;
    else
        dstlen = -1;
    
    pdf_byte_type* btype = new pdf_byte_type[srclen];
    try {
        ParseText(src, btype);

        if (dstlen > PDF_UNICODE_HEADER_LEN) {
            memcpy(dst, UNICODE_HEADER, PDF_UNICODE_HEADER_LEN);
            j = PDF_UNICODE_HEADER_LEN;
        } else
            return -1;

        const unsigned char* usrc = (const unsigned char*)src;
        for (int i = 0; i < srclen; i++) {
            unsigned int unicode;
            unsigned char tmp[2];

            if (btype[i] != PDF_BYTE_TYPE_TRIAL) {
                if (btype[i] == PDF_BYTE_TYPE_SINGLE)
                    unicode = GetUnicode(usrc[i]);
                else {
                    unsigned int code = (unsigned int)usrc[i] * 256 + usrc[i+1];
                    unicode = GetUnicode(code);
                }
                tmp[0] = unicode >> 8;
                tmp[1] = unicode;

                if (j < dstlen) {
                    memcpy(dst + j, tmp, 2);
                    j += 2;
                }
            }
        }
    } catch (...) {
        delete btype;
        throw;
    }
    delete btype;
    *len = j;
    return j;
}

/*----------------------------------------------------------------------------*/

