/*
 * << H a r u -- Free PDF Library >> -- PdfStreams.cpp
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

#include <errno.h>
#include <locale.h>
#include "libharu.h"

#ifndef NOZLIB
#include <zconf.h>
#include <zlib.h>
#endif /* NOZLIB */

/*----- PdfStreamBase --------------------------------------------------------*/

PdfStreamBase::PdfStreamBase()
{
    fNumBytes = 0;
    fDecimalPoint = 0x00;
    PDF_DEBUG_PRINT(("++ [%x] PdfStreamBase new.\n", (int)this));
}

PdfStreamBase::~PdfStreamBase()
{
    PDF_DEBUG_PRINT(("++ [%x] PdfStreamBase delete.\n", (int)this));
}

PdfStreamBase&
PdfStreamBase::operator<<(const char* s)
{
    int len = strlen(s);

    Write(s, len);
    return *this;
}

PdfStreamBase&
PdfStreamBase::operator<<(unsigned int i)
{
    char Buf[11];
    int len;

#ifdef __WIN32__
    len = _snprintf(Buf, 11, "%u", i);
#else
    len = snprintf(Buf, 11, "%u", i);
#endif
    Write(Buf, len);
    return *this;
}

PdfStreamBase&
PdfStreamBase::operator<<(int i)
{
    char Buf[12];
    int len;

#ifdef __WIN32__
    len = _snprintf(Buf, 12, "%d", i);
#else
    len = snprintf(Buf, 12, "%d", i);
#endif
    Write(Buf, len);
    return *this;
}

PdfStreamBase&
PdfStreamBase::operator<<(double d)
{
    static const int BUF_LEN = 12 + 1;
    char Buf[BUF_LEN];

    if (d < PDF_LIMIT_MIN_REAL || d > PDF_LIMIT_MAX_REAL)
        throw PdfException(PDF_ERR_INVALID_RANGE,
                "error: operator<<(const double d)[%f]\n.", d);
    
    if (fDecimalPoint == 0x00) {
        lconv *lv = localeconv();

        fDecimalPoint = *lv->decimal_point;
    }

#ifdef __WIN32__
    int len = _snprintf(Buf, BUF_LEN, "%.5f", d);
#else
    int len = snprintf(Buf, BUF_LEN, "%.5f", d);
#endif

    if (fDecimalPoint != '.') {
        char* pos = (char*)memchr(Buf, fDecimalPoint, len);
        if (pos != NULL)
            *pos = '.';
    }

    char* pBuf = Buf + len - 1;
    while (pBuf > Buf) {
        if (*pBuf == '0')
            *pBuf = 0x00;
        else {
            if (*pBuf == '.')
                *pBuf = 0x00;
            break;
        }
        pBuf--;
    }

    if (len >= BUF_LEN) 
        throw PdfException(PDF_ERR_INVALID_RANGE, 
            "error: operator<<(const double d)[%d]\n.", len);

    Write(Buf, strlen(Buf));
    return *this;
}

/*----------------------------------------------------------------------------*/
/*----- PdfFileStream class --------------------------------------------------*/

PdfFileStream::PdfFileStream(FILE* file, bool closefile)
        : PdfStreamBase()
{
    fFile = file;
    fCloseFile = closefile;
    fErrorNo = 0;
}

PdfFileStream::PdfFileStream(const char *filename)
        : PdfStreamBase()
{
    fFile = fopen(filename, "wb");
    if (fFile == NULL) {
        fErrorNo = errno;
        fCloseFile = false;
        return;
    }
    fCloseFile = true;
}

PdfFileStream::~PdfFileStream()
{
    if (fFile == NULL)
        return;
    
    fflush(fFile);
    if (fCloseFile)
        fclose(fFile);
}

int
PdfFileStream::Write(const void* ptr, int count)
{
    if (fFile == NULL)
        throw PdfException(PDF_ERR_FILE_OPEN, 
                "ERROR: cannot write to stream(%d).", fErrorNo);

    if (fwrite(ptr, count, 1, fFile) == 1)
        fNumBytes += count;
    else {
        fErrorNo = errno;
        throw PdfException(PDF_ERR_FILE_OPEN, 
                "ERROR: cannot write to stream(%d).", fErrorNo);
    }
    return count;
}

/*----------------------------------------------------------------------------*/
/*----- PdfMemStream class ---------------------------------------------------*/

PdfMemStream::PdfMemStream(const int bufSize)
        : PdfStreamBase()
{
    fBufSize = bufSize;
    fCurrentPos = 0;
    fBuf = NULL;
    fCurrentPtr = NULL;
}

void
PdfMemStream::InitStream()
{
    if (fBuf != NULL)
        return;

    fBuf = new PdfList;
    PDF_DEBUG_PRINT(("++ [%x] fBuf new.\n", (int)fBuf));

    fCurrentPtr = new unsigned char[fBufSize];
    PDF_DEBUG_PRINT(("++ [%x] fCurrentPtr new[].\n", (int)fCurrentPtr));

    fBuf->AddItem(fCurrentPtr);
    fCurrentPos = 0;
}

PdfMemStream::~PdfMemStream()
{
    if (fBuf != NULL) {
        for (int i = 0; i < fBuf->CountItems(); i++) {
            PDF_DEBUG_PRINT(("++ [%x] fBuf->ItemAt(i) delete[].\n",
                        (int)fBuf->ItemAt(i)));
            delete[] (unsigned char*)fBuf->ItemAt(i);
        }
        PDF_DEBUG_PRINT(("++ [%x] fBuf delete.\n", (int)fBuf));
        delete fBuf;
    }
}

PdfMemStream*
PdfMemStream::Duplicate() const
{
    PdfMemStream* pNewStream = new PdfMemStream(fBufSize);

    int numBuf = fBuf->CountItems();

    for (int i=0;i<numBuf-1;i++)
        pNewStream->Write((void*)fBuf->ItemAt(i), fBufSize);

    pNewStream->Write(fBuf->ItemAt(numBuf-1), fCurrentPos);

    return pNewStream;
}

int
PdfMemStream::Write(const void* ptr, int count)
{
    if (fBuf == NULL)
        InitStream();
    
    int writeSize = count;

    while (writeSize > 0)
        InternalWrite((unsigned char**)&ptr, &writeSize);

    fNumBytes += count;
    return count;
}

void
PdfMemStream::InternalWrite(unsigned char** ptr, int *count)
{
    if (*count <= 0) 
        return;
    int rsize = fBufSize - fCurrentPos;

    PDF_DEBUG_PRINT(("PdfMemStream::InternalWrite rsize=%d, count=%d.\n",
                rsize, *count));

    if (rsize >= *count) {
        memcpy(fCurrentPtr, *ptr, *count);
        fCurrentPtr += *count;
        fCurrentPos += *count;
        *count = 0;
    } else {
        if (rsize > 0) {
            memcpy(fCurrentPtr, *ptr, rsize);
            *ptr += rsize;
            *count -= rsize;
            PDF_DEBUG_PRINT(("PdfMemStream::InternalWrite "
                        "rsize=%d count=%d,\n", rsize, *count));
        }
        fCurrentPtr = new unsigned char[fBufSize];
        PDF_DEBUG_PRINT(("++ [%x] fCurrentPtr new[].\n", (int)fCurrentPtr));

        fBuf->AddItem(fCurrentPtr);
        PDF_DEBUG_PRINT(("PdfMemStream::InternalWrite --itemcount=%d.\n",
                    fBuf->CountItems()));
        fCurrentPos = 0;
    }
}

size_t
PdfMemStream::WriteToStream(PdfStreamBase* out, PdfEncryptor* eobj)
{
    size_t numWrote = 0;
    int numBuf = fBuf->CountItems();

    if (eobj == NULL) {
        for (int i = 0; i < numBuf - 1; i++)
            numWrote += out->Write((void*)fBuf->ItemAt(i), fBufSize);
        numWrote += out->Write(fBuf->ItemAt(numBuf-1), fCurrentPos);
    } else {
#ifdef USE_ENCRYPTION
        unsigned char* tmp = new unsigned char[fBufSize];
     
        try {   
            for (int i = 0; i < numBuf - 1; i++) {
                eobj->CryptBuf((unsigned char*)fBuf->ItemAt(i), tmp, fBufSize);
                numWrote += out->Write((void*)tmp, fBufSize);
            }
            eobj->CryptBuf((unsigned char*)fBuf->ItemAt(numBuf - 1), 
                    tmp, fCurrentPos);
            numWrote += out->Write(tmp, fCurrentPos);
        } catch (...) {
            delete[] tmp;
            throw;
        }
        delete[] tmp;
#endif
    }
    return numWrote;
}

size_t
PdfMemStream::WriteToStreamDeflate(PdfStreamBase* out, PdfEncryptor* eobj)
{
#ifndef NOZLIB
    z_stream strm;
    Bytef* otbuf = NULL;
    int ret;
    unsigned int numBytes = 0;
    
    /* check the status of the output stream. */
    if (out == NULL)
        throw PdfException(PDF_INVALID_PARAMETER, 
                "ERROR: write to stream. Output stream is NULL.");

    /* initialize decompression stream. */
    memset(&strm, 0x00, sizeof(z_stream));

    if (eobj == NULL)
        otbuf = new Bytef[fBufSize];
    else
        otbuf = new Bytef[fBufSize * 2];
    PDF_DEBUG_PRINT(("++ [%x] otbuf new.\n", (int)otbuf));
    try {
        strm.next_out = otbuf;
        strm.avail_out = fBufSize;

        if ((ret = deflateInit_(&strm, Z_DEFAULT_COMPRESSION, ZLIB_VERSION,
            sizeof(z_stream))) != Z_OK)
            throw PdfException(PDF_DEFLATEER_ERROR, 
                    "ERROR: PdfMemStream --deflateInit_ returns %d.", ret);

        /* decompress each brocks of the stream. */
        int numBuf = fBuf->CountItems();
        for (int i = 0; i < numBuf - 1; i++) {
            strm.next_in = (Bytef*)fBuf->ItemAt(i);
            strm.avail_in = fBufSize;

            while (strm.avail_in > 0) {
                ret = deflate(&strm, Z_NO_FLUSH);
                if (ret != Z_OK && ret != Z_STREAM_END)
                    throw PdfException(PDF_DEFLATEER_ERROR,
                           "ERROR: PdfMemStream --deflate returns %d.", ret);

                if (strm.avail_out == 0) {
                    if (eobj == NULL)
                        numBytes += out->Write(otbuf, fBufSize);
                    else {
#ifdef USE_ENCRYPTION
                        eobj->CryptBuf(otbuf, otbuf + fBufSize, fBufSize);
                        numBytes += out->Write(otbuf + fBufSize, fBufSize);
#endif
                    }
                    strm.next_out = otbuf;
                    strm.avail_out = fBufSize;
                }
            }
        }

        /* decompress last buffer of the stream. */
        strm.next_in = (Bytef*)fBuf->ItemAt(numBuf - 1);
        strm.avail_in = fCurrentPos;
        while ((ret = deflate(&strm, Z_FINISH)) != Z_STREAM_END) {
            if (ret != Z_OK)
                throw PdfException(PDF_DEFLATEER_ERROR, 
                        "ERROR: PdfMemStream --deflate returns %d.", ret);

            if (strm.avail_out == 0) {
                if (eobj == NULL)
                    numBytes += out->Write(otbuf, fBufSize);
                else {
#ifdef USE_ENCRYPTION
                    eobj->CryptBuf(otbuf, otbuf + fBufSize, fBufSize);
                    numBytes += out->Write(otbuf + fBufSize, fBufSize);
#endif
                }
                strm.next_out = otbuf;
                strm.avail_out = fBufSize;
            }
        }
        if (strm.avail_out < fBufSize)
            if (eobj == NULL)
                numBytes += out->Write(otbuf, fBufSize - strm.avail_out);
            else {
#ifdef USE_ENCRYPTION
                eobj->CryptBuf(otbuf, otbuf + fBufSize, 
                        fBufSize - strm.avail_out);
                numBytes += out->Write(otbuf + fBufSize, 
                        fBufSize - strm.avail_out);
#endif
            }
    } catch (PdfException& e) {
        deflateEnd(&strm);
        delete[] otbuf;
        PDF_DEBUG_PRINT(("++ [%x] otbuf delete.\n", (int)otbuf));
        throw;
    }
    /* finalize */
    deflateEnd(&strm);
    delete[] otbuf;
    PDF_DEBUG_PRINT(("++ [%x] otbuf delete.\n", (int)otbuf));
    return numBytes;
#else 
    throw PdfException(PDF_DEFLATEER_ERROR, 
            "ERROR: flatedecode is not supported.");
#endif /* NOZLIB */
}

void
PdfMemStream::FreeData()
{
    if (fBuf == NULL)
        return;

    for (int i = 0; i < fBuf->CountItems(); i++) {
        PDF_DEBUG_PRINT(("++ [%x] fBuf->ItemAt(i) delete[].\n",
                    (int)fBuf->ItemAt(i)));
        delete[] (unsigned char*)fBuf->ItemAt(i);
    }
    PDF_DEBUG_PRINT(("++ [%x] fBuf delete.\n", (int)fBuf));
    delete fBuf;
    fBuf = NULL;
}

const void*
PdfMemStream::GetBuf(int index, unsigned int* size)
{
    if (fBuf == NULL)
        return NULL;

    int numitem = fBuf->CountItems();

    if (numitem == 0)
        return NULL;
    
    if (index == numitem - 1)
        *size = fCurrentPos;
    else
        *size = fBufSize;
    
    return fBuf->ItemAt(index);
}

/*----------------------------------------------------------------------------*/

