/*
 * << H a r u -- Free PDF Library >> -- PdfXref.cpp
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
#include "libharu.h"

#define PDF_FREE_ENTRY          'f'
#define PDF_IN_USE_ENTRY        'n'

#define PDF_MAX_GENERATION_NUM  65535

/*----- PdfXrefEntry class ---------------------------------------------------*/

PdfXrefEntry::~PdfXrefEntry()
{
    if (HasObject())
        delete fObject;
    PDF_DEBUG_PRINT(("++ [%x] PdfXrefEntry delete.\n", (int)this));
}

PdfXrefEntry::PdfXrefEntry(PdfObject *object)
{
    fEntryType = PDF_FREE_ENTRY;
    fByteOffset = 0;
    fObject = NULL;
    SetObject(object);
    PDF_DEBUG_PRINT(("++ [%x] PdfXrefEntry new.\n", (int)this));
}

void
PdfXrefEntry::SetObject(PdfObject *object)
{
    fObject = object;
    if (object != NULL) {
        fGenerationNo = object->GetGenerationNo();
        fEntryType = PDF_IN_USE_ENTRY;
    }
}

/*----------------------------------------------------------------------------*/
/*----- Xref class -----------------------------------------------------------*/

PdfXref::PdfXref(PdfDoc* doc)
{
    fDoc = doc;
    fEntries = NULL;
    PDF_DEBUG_PRINT(("++ [%x] PdfXref new.\n", (int)this));
}

void
PdfXref::Init()
{
    fEntries = new PdfList();
    PdfXrefEntry* newEntry = new PdfXrefEntry((PdfObject*)NULL);
    newEntry->SetEntryType(PDF_FREE_ENTRY);
    newEntry->SetGenerationNo(PDF_MAX_GENERATION_NUM);
    if (fEntries->AddItem(newEntry) == false) {
        delete newEntry;
        throw PdfException(PDF_ERR_MALLOC, "PdfXref::Init()");
    }
}

PdfXref::~PdfXref()
{
    Clear();
    delete fEntries;
    PDF_DEBUG_PRINT(("++ [%x] PdfXref delete.\n", (int)this));
}

void
PdfXref::Clear()
{
    PDF_DEBUG_PRINT(("PdfXref::Clear: item count:[%d]\n", GetCount()));

    if (fEntries == NULL)
        return;
    
    for (int i = 0; i < GetCount(); i++) {
        PdfXrefEntry* entry = GetEntry(i);
        if (entry != NULL) 
            delete entry;
    }
    fEntries->Clear();
    fAddr = 0;
}

PdfObject*
PdfXref::GetObject(PdfOID objectID)
{
    if (fEntries == NULL)
        return NULL;
    
    PdfXrefEntry* entry = GetEntry(objectID);
    if (entry && entry->HasObject())
        return entry->GetObject();
    else
        return NULL;
}

PdfOID
PdfXref::AddObject(PdfObject* object)
{
    PDF_DEBUG_PRINT(("PdfXref::AddObject: [%X]\n", (int)object));

    PdfXrefEntry* newEntry = NULL;
    
    try {
        if (fEntries == NULL)
            Init();
        newEntry = new PdfXrefEntry(object);
    } catch (ALLOC_ERROR& e) {
        delete object;
        throw;
    }
    
    if (fEntries->AddItem((void*)newEntry)) {
        object->SetObjectID(fEntries->CountItems() - 1);
    } else {
        PDF_DEBUG_PRINT(("PdfXref::AddObject: failed to add object at:[%x]\n",
            (int)object));
        delete newEntry;
        throw PdfException(PDF_ERR_MALLOC, "PdfXref::AddObject");   
    }
    return object->GetObjectID();
}

void
PdfXref::WriteToStream(PdfStreamBase* out, PdfEncryptor* e)
{
    char Buf[11];
    PdfXrefEntry* entry;

    if (fEntries == NULL)
        Init();
    
    for (int i = 1; i < GetCount(); i++) {
        entry = GetEntry(i);
        PDF_DEBUG_PRINT(("PdfXref::WriteToStream: objectId = %d\n",
            entry->GetObject()->GetObjectID()));
        entry->SetByteOffset(out->GetPos());
#ifdef USE_ENCRYPTION
        if (e != NULL)
            e->Init(entry->GetObject()->GetObjectID(),
                entry->GetGenerationNo());
#endif
        entry->GetObject()->WriteValueToStream(out, e);
    }

    fAddr = out->GetPos();
    *out << "xref\015\0120 " << (unsigned int)GetCount() << "\015\012";
    for (int j = 0; j < GetCount(); j++) {
        entry = GetEntry(j);
       
#ifdef __WIN32__ 
        _snprintf(Buf, 11, "%010d", entry->GetByteOffset());
#else
        snprintf(Buf, 11, "%010d", entry->GetByteOffset());
#endif
        *out << Buf
             << " ";
#ifdef __WIN32__
        _snprintf(Buf, 6, "%05d", entry->GetGenerationNo());
#else
        snprintf(Buf, 6, "%05d", entry->GetGenerationNo());
#endif
        *out << Buf
             << " "
             << entry->GetEntryType()
             << "\015\012";
    }
}

void
PdfXref::SetError(int err)
{
    fDoc->SetError(err);
}

/*----------------------------------------------------------------------------*/

