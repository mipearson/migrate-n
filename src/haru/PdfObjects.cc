/*
 * << H a r u --free pdf library >> -- PdfObjects.cpp
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

#include <stdarg.h>
#include <ctype.h>
#include <assert.h>
#include <new>

#include "libharu.h"

/*----- PdfObject class ------------------------------------------------------*/

PdfObject::PdfObject()
{
    fObjectID = 0;
    fGenerationNo = 0;
    PDF_DEBUG_PRINT(("++ [%x] PdfObject new.\n", (int)this));
}

PdfObject::PdfObject(PdfOID objectID)
{
    fObjectID = objectID;
    fGenerationNo = 0;
    PDF_DEBUG_PRINT(("++ [%x] PdfObject(%d) new.\n", (int)this, fObjectID));
}

PdfObject::~PdfObject()
{
    PDF_DEBUG_PRINT(("++ [%x] PdfObject(%d) delete.\n", (int)this, fObjectID));
}

void
PdfObject::SetObjectID(PdfOID value)
{
    assert(fObjectID == 0);
    assert(value >= -1);

    PDF_DEBUG_PRINT(("PdfObject[%X] SetObjectID(%d).\n", (int)this, 
                fObjectID));
    fObjectID = value;
}

void
PdfObject::WriteValueToStream(PdfStreamBase* out, PdfEncryptor* e)
{
    // write indirect object to specified stream. this method called by parent
    // object.
    if (GetObjectType() != PDF_OBJ_TYPE_INDIRECT)
        return;

    *out << fObjectID
         << " "
         << fGenerationNo
         << " obj\015\012";
    InternalWriteStream(out, e);
    *out << "\015\012endobj\015\012";
}

#define NEEDS_ESCAPE(c)    (c < 0x21 || \
                            c > 0x7e || \
                            c == '\\' || \
                            c == '%' || \
                            c == '#' || \
                            c == '/' || \
                            c == '(' || \
                            c == ')' || \
                            c == '<' || \
                            c == '>' || \
                            c == '[' || \
                            c == ']' || \
                            c == '{' || \
                            c == '}' )

size_t
PdfObject::WriteEscapeName(PdfStreamBase* stream, const char* text)
{
    char tmpChar[PDF_LIMIT_MAX_NAME * 3 + 2];

    size_t len = strlen(text);
    const unsigned char* pos1 = (unsigned char*)text;
    char* pos2 = tmpChar;

    *pos2++ = '/';
    for (int i = 0; i < (int)len; i++) {
        unsigned char c = *pos1++;
        if (NEEDS_ESCAPE(c)) {
            *pos2 = (char)(*pos1 >> 4);
            if (*pos2 <= 9)
                *pos2 += 0x30;
            else
                *pos2 += 0x41 - 10;
            pos2++;

            *pos2 = (char)(*pos1 & 0x0f);
            if (*pos2 <= 9)
                *pos2 += 0x30;
            else
                *pos2 += 0x41 - 10;
            pos2++;
        } else
            *pos2++ = c;
    }
    *pos2 = 0x00;

    *stream << tmpChar;
    len = strlen(tmpChar);

    return len;
}

size_t
PdfObject::WriteEscapeText(PdfStreamBase* stream, const char* text)
{
    const int WRITE_BUF_SIZ = 256;
    unsigned char buf[WRITE_BUF_SIZ];
    size_t len = strlen(text);
    size_t idx = 0;
    size_t ret = 0;
    const char* pChar = text;

    buf[idx++] = '(';

    for (int i = 0; i < (int)len; i++) {
        unsigned char c = (unsigned char)*pChar++;
        if (NEEDS_ESCAPE(c)) {
            buf[idx++] = '\\';

            buf[idx] = c >> 6;
            buf[idx] += 0x30;
            idx++;
            buf[idx] = (c & 0x38) >> 3;
            buf[idx] += 0x30;
            idx++;
            buf[idx] = (c & 0x07);
            buf[idx] += 0x30;
            idx++;
        }
        else
            buf[idx++] = c;
        
        if (idx > WRITE_BUF_SIZ - 4) {
            stream->Write(buf, idx);
            ret += idx;
            idx = 0;
        }
    }
    buf[idx++] = ')';
    if (idx > 0) {
        stream->Write(buf, idx);
        ret += idx;
    }

    return ret;
}

void
PdfObject::WriteBinary(PdfStreamBase* out, unsigned char* buf, 
        unsigned int len, PdfEncryptor* e)
{
    const unsigned char* pos;
    int buf_len = len * 3 + 1;
    char* tmp = new char[buf_len];
    char* pos2 = tmp;
    unsigned char* tmp2 = NULL;
    PDF_DEBUG_PRINT(("++ [%x] buf new[].\n", (int)tmp));
#ifdef USE_ENCRYPTION    
    if (e == NULL)
#endif
        pos = buf;
#ifdef USE_ENCRYPTION
    else
    try {
        tmp2 = new unsigned char[len];
        e->CryptBuf(buf, tmp2, len);
        pos = tmp2;
    } catch (...) {
        delete[] tmp;
        throw;
    }
#endif

    for (int i = 0; i < (int)len; i++, pos++) {
        *pos2 = (char)(*pos >> 4);
        if (*pos2 <= 9)
            *pos2 += 0x30;
        else
            *pos2 += 0x41 - 10;
        pos2++;

        *pos2 = (char)(*pos & 0x0f);
        if (*pos2 <= 9)
            *pos2 += 0x30;
        else
            *pos2 += 0x41 - 10;
        pos2++;
    }
    *pos2 = 0x00;

    try {
        *out << "<";
        *out << tmp;
        *out << ">";
    } catch (...) {
        PDF_DEBUG_PRINT(("++ [%x] buf delete[].\n", (int)tmp));
        delete[] tmp;
        delete[] tmp2;
        throw;
    }

    PDF_DEBUG_PRINT(("++ [%x] buf delete[].\n", (int)tmp));
    delete[] tmp;
    delete[] tmp2;
}

/*----------------------------------------------------------------------------*/
/*----- PdfName class --------------------------------------------------------*/

PdfName::PdfName(PdfOID objectID, const char* value)
        :PdfObject(objectID)
{
    SetValue(value);
}

PdfName::PdfName(const char* value)
        : PdfObject()
{
    SetValue(value);
}

void
PdfName::SetValue(const char* value)
{
    if (value == NULL) {
        fValue[0] = 0x00;
    } else {
        int len = strlen(value);

        if (len > PDF_LIMIT_MAX_NAME)
            len = PDF_LIMIT_MAX_NAME;

        if (len > 0)
            memcpy(fValue, value, len);
        fValue[len] = 0x00;
    }
}

/*----------------------------------------------------------------------------*/
/*----- PdfText class --------------------------------------------------------*/

PdfText::PdfText(PdfOID objectID, const char* value)
        :PdfObject(objectID)
{
    fValue = NULL;
    SetValue(value);
}

PdfText::PdfText(const char* value)
        : PdfObject()
{
    fValue = NULL;
    SetValue(value);
}

PdfText::~PdfText()
{
    if (fValue!=NULL)
        delete[] fValue;
}

void
PdfText::SetValue(const char* value)
{
    if (fValue != NULL) {
        delete[] fValue;
        fValue = NULL;
    }

    if (value == NULL)
        return;

    int len = strlen(value);
    if (len == 0)
        return;

    if (len > PDF_LIMIT_STRING_MAX)
        len = PDF_LIMIT_STRING_MAX;

    fValue = new char[len + 1];

    memcpy(fValue, value, len);
    fValue[len] = 0x00;
}

void
PdfText::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    if (fValue == NULL)
        *out << "()";
    else if (e == NULL)
        WriteEscapeText(out, fValue);
    else
        WriteBinary(out, (unsigned char*)fValue, strlen(fValue), e);
}

/*----------------------------------------------------------------------------*/
/*----- PdfBinary class ------------------------------------------------------*/

void
PdfBinary::SetData(const void* data, unsigned int length, 
        bool encryptable)
{
    fData = data;
    fLength = length;
    fEncryptable = encryptable;
}

void
PdfBinary::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    if (fEncryptable)
        WriteBinary(out, (unsigned char*)fData, fLength, e);
    else
        WriteBinary(out, (unsigned char*)fData, fLength, NULL);
}

/*----------------------------------------------------------------------------*/
/*----- PdfArray class -------------------------------------------------------*/

PdfArray::PdfArray(PdfXref* xref)
        : PdfObject()
{
    fItems = NULL;
    fXref = xref;
}

PdfArray::PdfArray(PdfOID objectID, PdfXref* xref)
        : PdfObject(objectID)
{
    fItems = NULL;
    fXref = xref;
}

PdfArray::~PdfArray()
{
    PdfObject *obj;

    if (fItems == NULL)
        return;
    for (int i = GetCount()-1; i >= 0; i--) {
        obj = (PdfObject*)fItems->ItemAt(i);

        if (obj != NULL) {
            PDF_DEBUG_PRINT(("obj[%x] class=[%d] locked=[%d] Type=[%d]\n",
                (int)obj, (int)obj->GetClass(),
                (int)obj->IsLockedObject(), (int)obj->GetObjectType()));
            delete obj;
        }
    }
    delete fItems;
}

PdfObject*
PdfArray::GetItem(int index)
{
    CheckList();
    
    if (GetCount() <= index)
        return NULL;

    PdfObject *obj = (PdfObject*)fItems->ItemAt(index);

    if (obj->GetObjectType() == PDF_OBJ_TYPE_DIRECT)
        return obj;
    else if (obj->GetObjectType() == PDF_OBJ_TYPE_VIRTUAL)
        return (fXref ?
            (fXref->GetObject(obj->GetObjectID())) : NULL);
    else
        throw PdfException(PDF_ERR_INVALID_RANGE,
            "ERROR: PdfArray::GetItem[%d] unknown object at[%X]",
                index, (unsigned long)obj);
    return NULL;
}

int
PdfArray::GetAsInteger(int index)
{
    PdfObject *item = GetItem(index);

    if (item == NULL)
        return(0);

    if (item->GetClass() == ocNumber)
        return (((PdfNumber*)item)->GetValue());
    if (item->GetClass() == ocReal)
        return (int)(((PdfNumber*)item)->GetValue());

    PDF_DEBUG_PRINT(("warning: PdfArray::GetAsInteger"
            " --invalid object-type[%d].\n", (int)item->GetClass()));
    return 0;
}

double
PdfArray::GetAsReal(int index)
{
    PdfObject *item = GetItem(index);

    if (item == NULL)
        return(0);

    if (item->GetClass() == ocReal)
        return (((PdfReal*)item)->GetValue());
    if (item->GetClass() == ocNumber)
        return (double)(((PdfNumber*)item)->GetValue());

    PDF_DEBUG_PRINT(("warning: PdfArray::GetAsReal"
            " --invalid object-type[%d].\n", (int)item->GetClass()));
    return 0;
}

void
PdfArray::Add(PdfObject *value)
{
    CheckList();
    
    if (value == NULL)
        return;

    if (fItems->CountItems() >= PDF_LIMIT_MAX_ARRAY) {
        if (value->GetObjectType() == PDF_OBJ_TYPE_DIRECT)
            delete value;
        throw PdfException(PDF_RUNTIME_ERROR, "ERROR: PdfArray::Add "
                "--number of items is over [%d].\n", PDF_LIMIT_MAX_ARRAY);
    }

    if (value->GetObjectType() == PDF_OBJ_TYPE_DIRECT) {
        try {
            fItems->AddItem(value);
            value->SetObjectID(-1);
        } catch (ALLOC_ERROR) {
            delete value;
            throw;
        }
    } else {
        PdfVirtualObject* vo = new PdfVirtualObject(value->GetObjectID());
        try {
            fItems->AddItem(vo);
        } catch (...) {
            delete vo;
            throw;
        }
    }
}

void
PdfArray::Add(const pdf_box value)
{
    Add(new PdfNumber(value.left));
    Add(new PdfNumber(value.bottom));
    Add(new PdfNumber(value.right));
    Add(new PdfNumber(value.top));
}

void
PdfArray::Add(const pdf_rect value)
{
    Add(new PdfReal(value.left));
    Add(new PdfReal(value.bottom));
    Add(new PdfReal(value.right));
    Add(new PdfReal(value.top));
}

void
PdfArray::Insert(PdfObject* value, int index)
{
    CheckList();
    
    int cnt = fItems->CountItems();

    if (index < 0 || cnt < index) {
        if (value->GetObjectType() == PDF_OBJ_TYPE_DIRECT)
            delete value;
        throw PdfException(PDF_RUNTIME_ERROR, "ERROR: PdfArray::Add "
                "--number of items is over [%d].\n", PDF_LIMIT_MAX_ARRAY);
    }
    if (cnt == index)
        Add(value);
    else {
        if (value->GetObjectType() == PDF_OBJ_TYPE_DIRECT) {
            try {
                fItems->AddItem(value, index);
                value->SetObjectID(-1);
            } catch (...) {
                delete value;
                throw;
            }
        } else {
            PdfVirtualObject* vo = new PdfVirtualObject(value->GetObjectID());
            try {
                fItems->AddItem(vo, index);
            } catch (...) {
                delete vo;
                throw;
            }
        }
    }
}

void
PdfArray::Clear()
{
    CheckList();

    for (int i = 0; i < GetCount(); i++) {
        PdfObject* obj = (PdfObject*)fItems->ItemAt(i);
        assert(obj != NULL);
        delete obj;
    }
    fItems->Clear();
}

void
PdfArray::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    CheckList();

    *out << "[";
    for (int i = 0; i < GetCount(); i++) {
        PdfObject* obj = (PdfObject*)fItems->ItemAt(i);
        obj->WriteToStream(out, e);
#ifdef USE_ENCRYPTION
        if (e != NULL)
            e->Reset();
#endif
        *out << " ";
    }
    *out << "]";
}

/*----------------------------------------------------------------------------*/
/*----- PdfDictElement class -------------------------------------------------*/

PdfDictElement::PdfDictElement()
{
    fKey = NULL;
    fValue = NULL;
    PDF_DEBUG_PRINT(("++ [%x] PdfDictElement new.\n", (int)this));
}

PdfDictElement::PdfDictElement(PdfName* key, PdfObject* value)
{
    fKey = key;
    fValue = value;
    PDF_DEBUG_PRINT(("++ [%x] PdfDictElement new.\n", (int)this));
}

PdfDictElement::PdfDictElement(const char* key, PdfObject* value)
{
    fKey = new PdfName(key);
    fValue = value;
    PDF_DEBUG_PRINT(("++ [%x] PdfDictElement new.\n", (int)this));
}

PdfDictElement::~PdfDictElement()
{
    if (fKey!=NULL && fKey->GetObjectType() != PDF_OBJ_TYPE_INDIRECT)
        delete fKey;

    if (fValue!=NULL) {
        PDF_DEBUG_PRINT(("obj[%x] class=[%d] locked=[%d] Type=[%d]\n",
            (int)fValue, (int)fValue->GetClass(),
            (int)fValue->IsLockedObject(), (int)fValue->GetObjectType()));
        delete fValue;
    }
    PDF_DEBUG_PRINT(("++ [%x] PdfDictElement delete.\n", (int)this));
}

void
PdfDictElement::SetValue(PdfObject* value)
{
    if (value!=NULL)
        delete fValue;
    fValue = value;
}

/*----------------------------------------------------------------------------*/
/*----- PdfDictionary class --------------------------------------------------*/

PdfDictionary::PdfDictionary(PdfXref* xref)
        : PdfObject()
{
    fXref = xref;
    fElements = NULL;
}

PdfDictionary::PdfDictionary(PdfOID objectID, PdfXref* xref)
        : PdfObject(objectID)
{
    fXref = xref;
    fElements = NULL;
}

PdfDictionary::~PdfDictionary()
{
    for (int i = GetCount()-1;i>=0;i--) {
        PdfDictElement *element =
            (PdfDictElement*)fElements->ItemAt(i);
        if (element != NULL)
            delete element;
    }
    if (fElements != NULL)
        delete fElements;
}

void
PdfDictionary::SetError(int err)
{
    fXref->SetError(err);
}

int
PdfDictionary::GetCount()
{
    return (fElements == NULL) ? 0 : fElements->CountItems();
}

PdfDictElement*
PdfDictionary::GetElement(const char* key)
{
    if (fElements == NULL) return NULL;
    
    for (int i=0;i<GetCount();i++) {
        PdfDictElement *element = GetElementAt(i);
        if (element->GetKey()->EqualTo(key))
            return element;
    }
    return NULL;
}

void
PdfDictionary::AddElement(const char* key, PdfObject* value)
{
    // adding the specified object to the dictionary. if the object which
    // have same key is exists, replace it to that specified.
    CheckList();
    
    PdfDictElement* element = NULL;
    if (fElements->CountItems() >= PDF_LIMIT_MAX_DICT_ELEMENT) {
        if (value->GetObjectType() == PDF_OBJ_TYPE_DIRECT)
            delete value;
        throw PdfException(PDF_RUNTIME_ERROR,
            "ERROR: PdfDictionary::AddElement --number "
                "of items is over [%d].\n", PDF_LIMIT_MAX_DICT_ELEMENT);
    }

    RemoveElement(key);
    if (value->GetObjectType() == PDF_OBJ_TYPE_DIRECT) {
        try {
            if (value->IsLockedObject()) 
                throw PdfException(PDF_RUNTIME_ERROR, 
                        "this object(id=%d) is owned by another object.",
                        GetObjectID());
            element = new PdfDictElement(key, value);
            value->SetObjectID(-1);
            if (fElements->AddItem(element) == false)
                throw PdfException(PDF_ERR_LIST_ADD, 
                        "ERROR: PdfDictionary::AddElement");
        } catch (...) {
            if (element != NULL) 
                delete element;
            else
                delete value;
            throw;
        }
    } else {
        PdfObject* newobj = NULL;
        try {
            assert(fXref != NULL);
            newobj = new PdfVirtualObject(value->GetObjectID());

            element = new PdfDictElement(key, newobj);
            if (fElements->AddItem(element) == false)
                throw PdfException(PDF_ERR_LIST_ADD,
                        "ERROR: PdfDictionary::AddElement");
        } catch (...) {
            if (element != NULL)
                delete element;
            else
                delete newobj;
            throw;
        }
    }
}

void
PdfDictionary::RemoveElement(const char* key)
{
    PdfDictElement *element = GetElement(key);
    if (element!=NULL) {
        fElements->RemoveItem(element);
        delete element;
    }
}

void
PdfDictionary::RemoveElement(unsigned int index)
{
    CheckList();
    
    PdfDictElement *element =
        (PdfDictElement*)fElements->RemoveItem(index);
    if (element!=NULL)
        delete element;
}

PdfObject*
PdfDictionary::GetValue(const char* key)
{
    PdfDictElement *element = GetElement(key);
    if (element != NULL) {
        PdfObject *obj = element->GetValue();
        if (obj->GetObjectType() != PDF_OBJ_TYPE_VIRTUAL)
            return obj;
        else
            return (fXref ? fXref->GetObject(obj->GetObjectID()) : NULL);
    }
    else
        return NULL;
}

const char*
PdfDictionary::GetKeyValue(unsigned int index)
{
    CheckList();
    PdfDictElement *element = GetElementAt(index);
    if (element != NULL)
        return element->GetKey()->GetValue();
    else
        return NULL;
}

void
PdfDictionary::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    PdfDictElement *element;

    *out << "<<\015\012";
    for (int i = 0; i < GetCount(); i++) {
        element = GetElementAt(i);
        element->GetKey()->WriteToStream(out, e);
        *out << " ";
        element->GetValue()->WriteToStream(out, e);
#ifdef USE_ENCRYPTION
        if (e != NULL)
            e->Reset();
#endif
        *out << "\015\012";
    }
    *out << ">>";
}

bool
PdfDictionary::IsTypeOf(const char* type)
{
    PdfObject* obj;

    obj = GetValue("Type");
    if (obj != NULL && obj->GetClass() == ocName) {
        const char* type_name = ((PdfName*)obj)->GetValue();
        if (strcmp(type, type_name) == 0)
            return true;
    }
    return false;
}

const char*
PdfDictionary::GetTextValue(const char* key)
{
    PdfObject *obj = GetValue(key);

    PDF_DEBUG_PRINT(("PdfDictionary::GetTextValue key=%s\n", key));

    if (obj != NULL && obj->GetClass() == ocText)
        return ((PdfText*)obj)->GetValue();
    else
        return NULL;
}

const char*
PdfDictionary::GetNameValue(const char* key)
{
    PdfObject *obj = GetValue(key);

    PDF_DEBUG_PRINT(("PdfDictionary::GetNameValue key=%s\n", key));

    if (obj != NULL && obj->GetClass() == ocName)
        return ((PdfName*)obj)->GetValue();
    else
        return NULL;
}

/*----------------------------------------------------------------------------*/
/*----- PdfStream class ------------------------------------------------------*/

PdfStream::PdfStream(PdfXref* xref)
        : PdfDictionary(xref)
{
    fStream = NULL;
    fFilter = PDF_FILTER_NONE;
}

PdfStream::PdfStream(PdfOID objectID, PdfXref* xref)
        : PdfDictionary(objectID, xref)
{
    fStream = NULL;
    fFilter = PDF_FILTER_NONE;
}

PdfStream::~PdfStream()
{
    if (fStream != NULL)
        delete fStream;
}

PdfMemStream*
PdfStream::GetStream()
{
    if (fStream == NULL) {
        fStream = new PdfMemStream();
        AddElement("Filter", new PdfArray(GetXref()));
    }
    return fStream;
}

void
PdfStream::AddFilter(pdf_filter filter)
{
    if ((filter & PDF_FILTER_DEFLATE) == PDF_FILTER_DEFLATE) {
        fFilter |= PDF_FILTER_DEFLATE;
    }
    if ((filter & PDF_FILTER_DCT_DECODE) == PDF_FILTER_DCT_DECODE) {
        fFilter |= PDF_FILTER_DCT_DECODE;
    }
}

void
PdfStream::RemoveFilter(pdf_filter filter)
{
    fFilter -= (fFilter & filter);
}

void
PdfStream::InternalWriteStream(PdfStreamBase* out, PdfEncryptor* e)
{
    PdfNumber* length = (PdfNumber*)GetValue("Length");

    PDF_DEBUG_PRINT(("PdfStream::InternalWriteStream..\n"));

    if (fStream == NULL || fStream->GetSize() == 0) {
        /* if the size of the stream is 0, write object as PdfDictionary */
        PDF_DEBUG_PRINT(("PdfStream::InternalWriteStream[%d] size=0.\n",
                    int(this)));

        RemoveElement("Length");
        RemoveElement("Filter");
        PdfDictionary::InternalWriteStream(out, e);
        return;
    } else if (length == NULL) {
        length = new PdfNumber(0);
        GetXref()->AddObject(length);
        AddElement("Length", length);
    } else if (length->GetObjectType() != PDF_OBJ_TYPE_INDIRECT ||
        /* if the lengths element is not indirect object, re-create object to
         * indirect object.
         */
        length->GetObjectID() < GetObjectID()) {
        RemoveElement("Length");

        length = new PdfNumber(0);
        GetXref()->AddObject(length);
        AddElement("Length", length);
    }

#ifdef NOZLIB
	RemoveFilter(PDF_FILTER_DEFLATE);
#endif

	if (fFilter != PDF_FILTER_NONE) {
		PdfArray* filter = new PdfArray(0);
		AddElement("Filter", filter);
		
    	if ((fFilter & PDF_FILTER_DEFLATE) == PDF_FILTER_DEFLATE) 
        	filter->Add(new PdfName("FlateDecode"));

	    if ((fFilter & PDF_FILTER_DCT_DECODE) == PDF_FILTER_DCT_DECODE)
    	    filter->Add(new PdfName("DCTDecode"));

	} else
		RemoveElement("Filter");
		
    PdfDictionary::InternalWriteStream(out, e);

    *out << "\015\012stream\015\012";

    int num_writes;
    if ((fFilter & PDF_FILTER_DEFLATE) == PDF_FILTER_DEFLATE)
        num_writes = fStream->WriteToStreamDeflate(out, e);
    else
        num_writes = fStream->WriteToStream(out, e);

    *out << "\015\012endstream";

    length->SetValue(num_writes);
}

