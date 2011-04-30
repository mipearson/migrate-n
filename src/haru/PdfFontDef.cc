/*
 * << H a r u --free pdf library >> -- PdfFontDef.cpp
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
#include <errno.h>
#include "libharu.h"

/*----- PdfType1FontDef class ------------------------------------------------*/

PdfType1FontDef::PdfType1FontDef()
    : PdfAutoPtrObject()
{
    InitAttributes();
    PDF_DEBUG_PRINT(("++ [%x] PdfFontDef new.\n", (int)this));
}

PdfType1FontDef::PdfType1FontDef(const char* basefont)
{
    assert(basefont);

    InitAttributes();
    SetParam(fBaseFont, basefont);
    SetParam(fFontName, basefont);
    PDF_DEBUG_PRINT(("++ [%x] PdfFontDef new.\n", (int)this));
}

PdfType1FontDef::~PdfType1FontDef()
{
    Clear();
    PDF_DEBUG_PRINT(("++ [%x] PdfFontDef delete.\n", (int)this));
}

void
PdfType1FontDef::InitAttributes()
{
    fWidths = NULL;
    fWidthsCount = 0;
    fBaseFont[0] = 0x00;
    fFirstChar = 0;
    fLastChar = 0;
    fAscent = 0;
    fDescent = 0;
	fFlags = PDF_FONT_STD_CHARSET;
    fFontName[0] = 0x00;
    fFontBBox.left = 0;
    fFontBBox.right = 0;
    fFontBBox.top = 0;
    fFontBBox.bottom = 0;
    fItalicAngle = 0;
    fStemV = 0;
    fAvgWidth = 0;
    fLeading = 0;
    fMaxWidth = 0;
    fMissingWidth = 0;
    fStemH = 0;
    fXHeight = 0;
    fCapHeight = 0;
    fCharSet[0] = 0x00;
    fIsBase14Font = false;
    fFontData = NULL;
    fLength1 = 0;
    fLength2 = 0;
    fLength3 = 0;
    fDefaultEncoding = PDF_STANDARD_ENCODING;
    fFontData = NULL;
    fDescriptor = NULL;
}

PdfType1FontDef*
PdfType1FontDef::Duplicate() const
{
	PdfType1FontDef* pNewDef = new PdfType1FontDef();

	try {
		if (fWidths) {
			pNewDef->fWidthsCount = fWidthsCount;
			pNewDef->fWidths = new pdf_char_data[fWidthsCount];
			for (int i = 0; i < fWidthsCount; i++) 
				memset(&pNewDef->fWidths[i], 0x00, sizeof(pdf_char_data));

			for (int j = 0; j < fWidthsCount; j++) {
				pdf_char_data ws = fWidths[j];
				int len = strlen(ws.char_name) + 1;
				char* char_name = new char[len];

				strncpy(char_name, ws.char_name, len);

				pNewDef->fWidths[j].char_cd = ws.char_cd;
				pNewDef->fWidths[j].char_name = char_name;
				pNewDef->fWidths[j].width = ws.width;
			}
		}

		pNewDef->fFirstChar = fFirstChar;
		pNewDef->fLastChar = fLastChar;
		pNewDef->fAscent = fAscent;
		pNewDef->fDescent = fDescent;
		pNewDef->fFlags = fFlags;

		pNewDef->SetParam(pNewDef->fBaseFont, fBaseFont);
		pNewDef->SetParam(pNewDef->fFontName, fFontName);

    	pNewDef->fFontBBox.left = fFontBBox.left;
    	pNewDef->fFontBBox.right = fFontBBox.right;
    	pNewDef->fFontBBox.top = fFontBBox.top;
    	pNewDef->fFontBBox.bottom = fFontBBox.bottom;
    	pNewDef->fItalicAngle = fItalicAngle;
    	pNewDef->fStemV = fStemV;
    	pNewDef->fAvgWidth = fAvgWidth;
    	pNewDef->fLeading = fLeading;
    	pNewDef->fMaxWidth = fMaxWidth;
    	pNewDef->fMissingWidth = fMissingWidth;
    	pNewDef->fStemH = fStemH;
    	pNewDef->fXHeight = fXHeight;
    	pNewDef->fCapHeight = fCapHeight;

		pNewDef->SetParam(pNewDef->fCharSet, fCharSet);

    	pNewDef->fIsBase14Font = fIsBase14Font;

    	pNewDef->fLength1 = fLength1;
    	pNewDef->fLength2 = fLength2;
    	pNewDef->fLength3 = fLength3;
    	pNewDef->fDefaultEncoding = fDefaultEncoding;

    	if (fFontData) {
        	pNewDef->fFontData = fFontData->Duplicate();
    	}

    	pNewDef->fDescriptor = NULL;
	} catch (...) {
		delete pNewDef;
		throw PdfException(errno, "ERROR: cannot duplicate font.");
	}

    return pNewDef;
}

void
PdfType1FontDef::SetWidths(const pdf_char_data_ro array[])
{
    assert(array != NULL);

    FreeWidths();

    int i = 0;
    while (array[i].char_name != NULL)
        i++;
    fWidthsCount = i;

	try {
    	fWidths = new pdf_char_data[fWidthsCount];
		for (i = 0; i < fWidthsCount; i++) 
			memset(&fWidths[i], 0x00, sizeof(pdf_char_data));
		
    	PDF_DEBUG_PRINT(("++ [%x] fWidths new[].\n", (int)fWidths));
    	for (i = 0; i < fWidthsCount; i++) {
        	PDF_DEBUG_PRINT(("PdfType1FontDef::SetWidths [%d]\n", i));
        	pdf_char_data_ro ws = array[i];
        	int len = strlen(ws.char_name) + 1;
        	char* char_name = new char[len];
        	PDF_DEBUG_PRINT(("++ [%x] char_name new[].\n", (int)char_name));

        	strncpy(char_name, ws.char_name, len);

        	fWidths[i].char_cd = ws.char_cd;
        	fWidths[i].char_name = char_name;
        	fWidths[i].width = ws.width;
        	if (strcmp(ws.char_name, "space") == 0) 
            	fMissingWidth = ws.width;
		}
	}	catch (...) {
		PDF_DEBUG_PRINT(("ERROR! PdfType1FontDef::SetWidths [%d]\n", i));
		FreeWidths();
    }
}

void
PdfType1FontDef::FreeWidths()
{
    if (fWidths != NULL) {
        for (int j = 0; j < fWidthsCount; j++)
            if (fWidths[j].char_name != NULL) {
                PDF_DEBUG_PRINT(("++ [%x] char_name delete[].\n",
                            (int)fWidths[j].char_name));
                delete[] fWidths[j].char_name;
            }

        PDF_DEBUG_PRINT(("++ [%x] fWidths delete[].\n", (int)fWidths));
        delete[] fWidths;
        fWidths = NULL;
    }
}

unsigned int
PdfType1FontDef::Widths(const char* char_name)
{
    if (fWidths == NULL)
        return 0;

    if (strcmp(char_name, PDF_CHAR_NOTDEF) == 0)
        return fMissingWidth;

    for (int i = 0; i < fWidthsCount; i++) {
        if (strcmp(fWidths[i].char_name, char_name) == 0) {
            PDF_DEBUG_PRINT(("FontDef::Widths name=%s, width=%d\n",
                char_name, fWidths[i].width));
            return fWidths[i].width;
        }
    }

    PDF_DEBUG_PRINT(("FontDef::Widths --not match [%s]\n", char_name));

    return fMissingWidth;
}

unsigned int
PdfType1FontDef::Widths(unsigned char c)
{
    for (int i = 0; i < fWidthsCount; i++) {

        if (fWidths[i].char_cd == (int)c) {
            PDF_DEBUG_PRINT(("FontDef::Widths char=%d, width=%d\n",
                            (unsigned char)c, fWidths[i].width));
            return fWidths[i].width;
        }
    }

    PDF_DEBUG_PRINT(("FontDef::Widths not match [%d:%d].\n",
                (unsigned char)c, fWidthsCount));
    return fMissingWidth;
}

void
PdfType1FontDef::Clear()
{
    FreeWidths();
    delete fFontData;
    InitAttributes();
}

void
PdfType1FontDef::LoadFontFile(const char* filename)
{
    FILE *f = fopen(filename, "rb");
    if (f == NULL) {
        throw PdfException(errno, "ERROR: file open failed[%s]", filename);
    }
    char buf[PDF_DEF_BUF_SIZE];
    fFontData = new PdfMemStream();
    int len2 = 0;

    char* pbuf = buf + 11;
    fread(buf, 11, 1, f);
    for (;;) {
        len2 = fread(pbuf, 1, PDF_DEF_BUF_SIZE - 11, f);

        if (fLength1 == 0) {
            int ret = FindBuf("eexec", buf, len2);
            if (ret > 0)
                fLength1 = fFontData->GetSize() + ret + 6;
        }
        if (fLength1 > 0 && fLength2 == 0) {
            int ret = FindBuf("cleartomark", buf, len2);
            if (ret > 0)
                fLength2 = fFontData->GetSize() - 520 - fLength1 + ret;
        }
        if (len2 < PDF_DEF_BUF_SIZE - 11) {
            fFontData->Write(buf, len2 + 11);
            break;
        }

        try {
            fFontData->Write(buf, len2);
        } catch (PdfException& e) {
            fclose(f);
            throw PdfException(e.GetCode(), e.what());
        }
        memmove(buf, buf + len2, 11);
        pbuf = buf + 11;
    }
    fLength3 = fFontData->GetSize() - fLength2 - fLength1;

    fclose(f);
}

void
PdfType1FontDef::LoadAfmFile(const char* filename)
{
    fFirstChar = 255;
    fLastChar = 0;
    char buf[PDF_DEF_BUF_SIZE];
    char s[PDF_LIMIT_MAX_NAME];

    FILE* f = fopen(filename, "rb");
    if (f == NULL)
        throw PdfException(errno, "ERROR: file [%s] open failed", filename);

    try {
        if (fgets(buf, PDF_DEF_BUF_SIZE, f) == NULL)
            throw PdfException(errno, "ERROR: file [%s] read failed", filename);

        if (strncmp(buf, "StartFontMetrics", strlen("StartFontMetrics")) != 0)
            throw PdfException(errno, "ERROR: file [%s] is not afm file.",
                filename);

        for (;;) {
            if (fgets(buf, PDF_DEF_BUF_SIZE, f) == NULL) {
                throw PdfException(errno,
                    "ERROR: file [%s] cannot parse AFM file.", filename);
            }
            char* pbuf = GetStrParam(buf, s, PDF_LIMIT_MAX_NAME);
            PDF_DEBUG_PRINT(("DEBUG: afm param=%s\n", s));

            if (strcmp(s, "FontName") == 0) {
                GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
                SetParam(fFontName, s);
                SetParam(fBaseFont, s);
            } else if (strcmp(s, "Weight") == 0) {
                GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
                if (strcmp(s, "Bold") == 0)
                    fFlags |= PDF_FONT_FOURCE_BOLD;
            } else if (strcmp(s, "IsFixedPitch") == 0) {
                GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
                if (strcmp("true", s) == 0)
                    fFlags |= PDF_FONT_FIXED_WIDTH;
            } else if (strcmp(s, "ItalicAngle") == 0) {
                GetIntParam(pbuf, &fItalicAngle);
                if (fItalicAngle != 0)
                    fFlags |= PDF_FONT_ITALIC;
            } else if (strcmp(s, "CharacterSet") == 0) {
                GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
                SetParam(fCharSet, s);
            } else if (strcmp(s, "FontBBox") == 0) {
                pbuf = GetIntParam(pbuf, &fFontBBox.left);
                pbuf = GetIntParam(pbuf, &fFontBBox.bottom);
                pbuf = GetIntParam(pbuf, &fFontBBox.right);
                GetIntParam(pbuf, &fFontBBox.top);
            } else if (strcmp(s, "EncodingScheme") == 0) {
                GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
                fDefaultEncoding = PdfEncodingDef::StringToPdfEncoding(s);
            } else if (strcmp(s, "CapHeight") == 0)
                GetIntParam(pbuf, &fCapHeight);
            else if (strcmp(s, "Ascender") == 0)
                GetIntParam(pbuf, &fAscent);
            else if (strcmp(s, "Descender") == 0)
                GetIntParam(pbuf, &fDescent);
            else if (strcmp(s, "STDHW") == 0)
                GetIntParam(pbuf, &fStemH);
            else if (strcmp(s, "STDHV") == 0)
                GetIntParam(pbuf, &fStemH);
            else if (strcmp(s, "StartCharMetrics") == 0) {
                GetIntParam(pbuf, &fWidthsCount);
                break;
            }
        }

        fWidths = new pdf_char_data[fWidthsCount];
        PDF_DEBUG_PRINT(("++ [%x] fWidths new[].\n", (int)fWidths));
        memset(fWidths, 0x00, sizeof(pdf_char_data) * fWidthsCount);
        pdf_char_data *pcwp = fWidths;

        for (int i = 0; i < fWidthsCount; i++, pcwp++) {
            if (fgets(buf, 1024, f) == NULL || strlen(buf) == 0)
                throw PdfException(errno,
                    "ERROR: reading CharMetrics from AFM file[%s].", filename);
            if (strncmp(buf, "EndCharMetrics", strlen("EndCharMetrics")) == 0)
                break;
            if (buf[0] != 'C')
                throw PdfException(errno,
                    "ERROR: invalid AFM file format[%s].", filename);

            char* pbuf = buf;
            pbuf = GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
            pbuf = GetIntParam(pbuf, &(*pcwp).char_cd);
            if ((*pcwp).char_cd > 0 && (*pcwp).char_cd < fFirstChar)
                fFirstChar = (*pcwp).char_cd;
            if ((*pcwp).char_cd > fLastChar)
                fLastChar = (*pcwp).char_cd;
            pbuf = GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
            pbuf = GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
            pbuf = GetIntParam(pbuf, &(*pcwp).width);
            pbuf = GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
            pbuf = GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
            GetStrParam(pbuf, s, PDF_LIMIT_MAX_NAME);
            int len = strlen(s) + 1;
            char* char_name = new char[len];
            strncpy(char_name, s, len);
            (*pcwp).char_name = char_name;
            if (strcmp(char_name, "space") == 0)
                fMissingWidth = (*pcwp).width;

            PDF_DEBUG_PRINT(("++ [%x] char_name new[].\n", (int)char_name));
            PDF_DEBUG_PRINT(("DEBUG: cd=%d, name=%s, widths=%d.\n",
                        (*pcwp).char_cd, (*pcwp).char_name, (*pcwp).width));
        }
    } catch (PdfException& e) {
        fclose(f);
        throw;
    }

    fclose(f);
}

void
PdfType1FontDef::LoadFromFile(const char* afmfile, const char* fontfile)
{
    Clear();

    PDF_DEBUG_PRINT(("DEBUG: loading file %s.\n", afmfile));
    LoadAfmFile(afmfile);

    if (fontfile != NULL) {
        PDF_DEBUG_PRINT(("DEBUG: loading file %s.\n", fontfile));
        LoadFontFile(fontfile);
    }
}

int
PdfType1FontDef::FindBuf(const char* word, const char* srcbuf, int len)
{
    const char* pos = srcbuf;
    int wlen = strlen(word);

    for (int i = 0; i < len; i++, pos++) {
        if (memcmp(word, pos, wlen) == 0)
            return i;
    }

    return -1;
}

char*
PdfType1FontDef::GetStrParam(char* str, char* param, int len)
{
    if (param == NULL)
        return NULL;
    *param = 0x00;

    if (str == NULL)
        return NULL;

    char* src = str;
    char* dst = param;

    for (int i = 0; i < len; src++){
        if (*src == 0x00)
            return NULL;

        if (*src == ' ' || *src == 0x0d || *src == 0x0a || *src == 0x09) {
            if (dst != param) {
                *dst = 0x00;
                return src++;
            }
        } else
            *dst++ = *src;
    }

    PDF_DEBUG_PRINT(("DEBUG param=%s.\n", param));
    return NULL;
}

char*
PdfType1FontDef::GetIntParam(char* str, int* param)
{
    const int MAX_INT_LEN = 10;
    char sparam[MAX_INT_LEN + 1];

    str = GetStrParam(str, sparam, MAX_INT_LEN + 1);
    int i = atoi(sparam);
    *param = i;

    return str;
}

void
PdfType1FontDef::SetParam(char* param, const char* text)
{
    memset(param, 0x00, PDF_LIMIT_MAX_NAME + 1);
    if (text != NULL) {
        int len = strlen(text);
        if (len <= PDF_LIMIT_MAX_NAME)
            memcpy(param, text, len);
        else
            memcpy(param, text, PDF_LIMIT_MAX_NAME);
    }
}

/*----------------------------------------------------------------------------*/

