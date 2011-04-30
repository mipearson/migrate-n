/*
 * << H a r u --free pdf library >> -- PdfDoc.cpp
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

/*----------------------------------------------------------------------------*/
/*----- PdfHeader class ------------------------------------------------------*/

void
PdfHeader::WriteToStream(PdfStreamBase* out)
{
    *out << "%PDF-1.2\012%\267\276\255\252\015\012";
}

/*----------------------------------------------------------------------------*/
/*----- PdfTrailer class -----------------------------------------------------*/

PdfTrailer::PdfTrailer(PdfXref *xref)
{
    fXref = xref;
    fXrefAddr = 0;
    fEncryptDict = NULL;
    fAttributes = new PdfDictionary(xref);
}

PdfTrailer::~PdfTrailer()
{
    delete fAttributes;
}

void
PdfTrailer::SetXrefAddr(unsigned int addr)
{
    fXrefAddr = addr;
}

void
PdfTrailer::SetXrefSize(unsigned int size)
{
    fAttributes->AddElement("Size", new PdfNumber(size));
}

void
PdfTrailer::SetRoot(PdfCatalog* root)
{
    fAttributes->AddElement("Root", root);
}

void
PdfTrailer::SetInfo(PdfInfo* info)
{
    fAttributes->AddElement("Info", info);
}

void
PdfTrailer::SetID(const unsigned char* id)
{
    /* the length of ID must be PDF_ID_LEN. */
    memcpy(fID1, id, PDF_ID_LEN);
    memcpy(fID2, id, PDF_ID_LEN);
    PdfArray* a = new PdfArray(fXref);
    fAttributes->AddElement("ID", a);
    PdfBinary* b = new PdfBinary();
    b->SetData(fID1, PDF_ID_LEN);
    a->Add(b);
    b = new PdfBinary();
    b->SetData(fID2, PDF_ID_LEN);
    a->Add(b);
}

void
PdfTrailer::SetEncryptDict(PdfEncryptDict* encrypt)
{
#ifdef USE_ENCRYPTION
    fAttributes->AddElement("Encrypt", encrypt);
    SetID(encrypt->ID());
    fEncryptDict = encrypt;
#endif
}

void
PdfTrailer::WriteToStream(PdfStreamBase* out)
{
    *out << "trailer\015\012";
    fAttributes->WriteToStream(out, NULL);
    *out << "\015\012startxref\015\012"
         << fXrefAddr
         << "\015\012%%EOF\015\012";
}

/*----------------------------------------------------------------------------*/
/*----- PdfDoc class ---------------------------------------------------------*/

PdfDoc::PdfDoc()
{
    Init();
}

void
PdfDoc::Init()
{
    /* initialize PdfDoc object. */
    PDF_DEBUG_PRINT(("PdfDoc::Init(): [%x] \n", (int)this));

    fXref = NULL;
    fHeader = NULL;
    fTrailer = NULL;
    fXObjectMgr = NULL;
    fFontMgr = NULL;
    fAutoPtrMgr = NULL;
    fHasDoc = false;

    FreeDoc();
}

void
PdfDoc::NewDoc()
{
    PDF_DEBUG_PRINT(("PdfDoc::NewDoc() \n"));
    FreeDoc(false);

    PDF_DEBUG_PRINT(("PdfDoc::create xref \n"));
    fXref = new PdfXref(this);

    PDF_DEBUG_PRINT(("PdfDoc::create header \n"));
    fHeader = new PdfHeader();

    PDF_DEBUG_PRINT(("PdfDoc::create trailer \n"));
    fTrailer = new PdfTrailer(fXref);

    PDF_DEBUG_PRINT(("PdfDoc::create xobject-mgr \n"));
    fXObjectMgr = new PdfXObjectMgr(fXref);

    PDF_DEBUG_PRINT(("PdfDoc::create font-mgr \n"));
    fFontMgr = new PdfFontMgr(fXref);

    PDF_DEBUG_PRINT(("PdfDoc::create catalog \n"));
    fCatalog = new PdfCatalog(fXref);
    try {
        fCatalog->Init();
    } catch (...) {
        delete fCatalog;
        throw;
    }

    fRootPages = new PdfPages(fXref, fFontMgr, fXObjectMgr);
    try {
        fRootPages->Init();
    } catch (...) {
        delete fRootPages;
        throw;
    }

    fCatalog->AddElement("Pages", fRootPages);
    fTrailer->SetRoot(fCatalog);
    fCurrentPages = fRootPages;
    fCurrentPage = NULL;

    Info()->SetProducer(LIB_HARU_VERSION_TXT);

    fHasDoc = true;
}

void
PdfDoc::FreeDoc(bool free_all_objects)
{
    if (fHasDoc)
        fHasDoc = false;

    delete fXref;
    delete fHeader;
    delete fTrailer;
    delete fXObjectMgr;
    delete fFontMgr;

    if (free_all_objects) {
        delete fAutoPtrMgr;
        fAutoPtrMgr = NULL;
    }

    fXref = NULL;
    fHeader = NULL;
    fTrailer = NULL;
    fCatalog = NULL;
    fInfo = NULL;
    fFontMgr = NULL;
    fXObjectMgr = NULL;
    fRootPages = NULL;
    fCurrentPages = NULL;
    fCurrentPage = NULL;

    fError = 0;
}

PdfDoc::~PdfDoc()
{
    FreeDoc();
    delete fAutoPtrMgr;
}

void
PdfDoc::SetError(int err)
{
    fError = err;
}

void
PdfDoc::WriteToStream(PdfStreamBase* out)
{
    PdfEncryptor* e = NULL;
    
    if (!fHasDoc)
        throw PdfException(PDF_RUNTIME_ERROR,
            "ERROR: PdfDoc::WriteToStream no document to save.");
#ifdef USE_ENCRYPTION
    if (IsEncrypted()) {
        GetEncryptDict()->EncryptPassword();
        fTrailer->SetEncryptDict(GetEncryptDict());
        e = new PdfEncryptor(GetEncryptDict()->FileKey());
    }

	try {
#endif
    fHeader->WriteToStream(out);
    fXref->WriteToStream(out, e);
    fTrailer->SetXrefAddr(fXref->GetAddr());
    fTrailer->SetXrefSize(fXref->GetCount());
    fTrailer->WriteToStream(out);

#ifdef USE_ENCRYPTION
	} catch (...) {
		delete e;
		throw;
	}
	delete e;
#endif
}

void
PdfDoc::WriteToFile(const char* filename)
{
    PdfFileStream *fs = new PdfFileStream(filename);
    try {
        WriteToStream(fs);
        delete fs;
    } catch (...) {
        delete fs;
        throw;
    }
}

void
PdfDoc::AddType1Font(PdfType1FontDef* fontdef, const char* name,
        PdfEncodingDef* encoding)
{
    if (fontdef == NULL)
        throw PdfException(PDF_RUNTIME_ERROR,
                "PdfDoc::AddType1Font fontdef is null.");

    RegisterObject(fontdef);
    PdfEncodingDef* en_def = encoding;

    if (name == NULL)
        name = fontdef->FontName();

    if (en_def == NULL) {
        pdf_encoding en = fontdef->DefaultEncoding();

        if (en == PDF_ENCODING_EOF)
            en_def = new PdfWinAnsiEncoding();
        else if (en == PDF_STANDARD_ENCODING)
            en_def = new PdfStandardEncoding();
        else if (en == PDF_MAC_ROMAN_ENCODING)
            en_def = new PdfMacRomanEncoding();
        else if (en == PDF_WIN_ANSI_ENCODING)
            en_def = new PdfWinAnsiEncoding();
    }
    if (en_def != NULL) 
        RegisterObject(en_def);

    PdfType1Font* font = new PdfType1Font(fXref);
    if (font->GetObjectType() != PDF_OBJ_TYPE_INDIRECT)
        fXref->AddObject(font);
    font->SetAttributes(name, fontdef, en_def);
    font->AddFilter(PDF_FILTER_DEFLATE);
    FontMgr()->RegisterFont(font);
}

void
PdfDoc::AddType0Font(PdfCIDFontDef* fontdef, const char* name,
        PdfCMap* cmap)
{
    if (fontdef == NULL)
        throw PdfException(PDF_RUNTIME_ERROR,
                "PdfDoc::AddType0Font fontdef is null.");

    if (cmap == NULL)
        throw PdfException(PDF_RUNTIME_ERROR,
                "PdfDoc::AddType0Font cmap is null.");

    RegisterObject(fontdef);
    RegisterObject(cmap);

    PdfType0Font* font = new PdfType0Font(fXref);
    if (font->GetObjectType() != PDF_OBJ_TYPE_INDIRECT)
        fXref->AddObject(font);

    PdfCIDFont* font2 = new PdfCIDType2Font(fXref);
    fXref->AddObject(font2);
    font2->SetAttributes(fontdef);

    font->SetAttributes(name, font2, cmap);
    cmap->AddCIDSystemInfo(font2);

    FontMgr()->RegisterFont(font);
}

void
PdfDoc::AddXObject(PdfXObject* xobject, const char* name)
{
    if (xobject->GetXref() != fXref)
        throw PdfException(PDF_RUNTIME_ERROR,
                "PdfDoc::AddXObject:: cannot add xobject. this object is "
                "owned by other PdfDoc object.");

    if (name == NULL) {
        /* create new name for the xobject. */
        char tmp_name[25];
#ifdef __WIN32__
        _snprintf(tmp_name, 25, "%s%X", "XOBJ", (unsigned int)xobject);
#else
        snprintf(tmp_name, 25, "%s%lX", "XOBJ", (unsigned long)xobject);
#endif
        XObjectMgr()->RegisterXObject(xobject, tmp_name);
    } else
        XObjectMgr()->RegisterXObject(xobject, name);
}

void
PdfDoc::RegisterObject(PdfAutoPtrObject* obj)
{
    if (fAutoPtrMgr == NULL)
        fAutoPtrMgr = new PdfAutoPtrMgr();
    fAutoPtrMgr->RegisterObject(obj);
}

PdfPage*
PdfDoc::AddPage()
{
    return (fCurrentPages == NULL) ? NULL : fCurrentPages->AddPage();
}

PdfInfo*
PdfDoc::Info()
{
    if (fInfo == NULL) {
        fInfo = new PdfInfo(fXref);
        fInfo->Init();
        fTrailer->SetInfo(fInfo);
    }
    return fInfo;
}

void
PdfDoc::SetPassword(const char* owner_passwd, const char* user_passwd)
{
#ifdef USE_ENCRYPTION
    PdfEncryptDict* e = GetEncryptDict();
    if (e == NULL) {
        e = new PdfEncryptDict(fXref, fInfo);
        e->Init();
        fTrailer->SetEncryptDict(e);
    }

    e->SetPassword(owner_passwd, user_passwd);
#else
    throw PdfException(PDF_ERR_NOT_SUPPORTED_FUNCTION, "ERROR: "
            "encryption is not supported.");
#endif
}

void
PdfDoc::SetPermission(int permission)
{
#ifdef USE_ENCRYPTION
    PdfEncryptDict* e = GetEncryptDict();
    if (e == NULL) {
        e = new PdfEncryptDict(fXref, fInfo);
        e->Init();
        fTrailer->SetEncryptDict(e);
    }

    e->SetPermission(permission);
#else
    throw PdfException(PDF_ERR_NOT_SUPPORTED_FUNCTION, "ERROR: "
            "encryption is not supported.");
#endif
}

/*----------------------------------------------------------------------------*/

