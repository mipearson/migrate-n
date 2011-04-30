/*
 * << H a r u -- Free PDF Library >> -- PdfAutoPtrObjects.cc
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
/*----- PdfAutoPtrObject class -----------------------------------------------*/

PdfAutoPtrObject::PdfAutoPtrObject()
{
    fMgr = NULL;
    PDF_DEBUG_PRINT(("++ [%x] PdfAutoPtrObject new.\n", (int)this));
}

PdfAutoPtrObject::~PdfAutoPtrObject()
{
    PDF_DEBUG_PRINT(("++ [%x] PdfAutoPtrObject delete.\n", (int)this));
    if (fMgr != NULL)
        fMgr->UnRegisterObject(this);
}

/*----------------------------------------------------------------------------*/
/*----- PdfAutoPtrMgr class --------------------------------------------------*/

PdfAutoPtrMgr::PdfAutoPtrMgr()
{
    fList = NULL;
    PDF_DEBUG_PRINT(("++ [%x] PdfAutoPtrMgr new.\n", (int)this));
}

PdfAutoPtrMgr::~PdfAutoPtrMgr()
{
    if (fList != NULL) {
        for (int i = 0; i < fList->CountItems(); i++) {
            PdfAutoPtrObject* obj = (PdfAutoPtrObject*)fList->ItemAt(i);
            obj->fMgr = NULL;   
            delete obj;
        }
        delete fList;
    }
    PDF_DEBUG_PRINT(("++ [%x] PdfAutoPtrMgr delete.\n", (int)this));
}

void
PdfAutoPtrMgr::RegisterObject(PdfAutoPtrObject* obj)
{
    if (obj->fMgr == NULL) {
        if (fList == NULL) 
            fList = new PdfList();
        if (fList->AddItem(obj)) {
            obj->fMgr = this;
            obj->Init();
        } else {
            delete obj;
            throw new PdfException(PDF_RUNTIME_ERROR, "ERROR: "
                    "PdfAutoPtrMgr::RegisterObject cannot register object.");
        }
    } else if (obj->fMgr != this)
        throw new PdfException(PDF_RUNTIME_ERROR, "ERROR:"
               "PdfAutoPtrMgr::RegisterObject this object owned by other "
               "object.");
}

void
PdfAutoPtrMgr::UnRegisterObject(PdfAutoPtrObject* obj)
{
    if (fList->RemoveItem(obj))
        obj->fMgr = NULL;
}

