/*
 * << H a r u --free pdf library >> -- PdfUtils.cpp
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
#include <new>

/*----- PdfList class ------------------------------------------------------*/

PdfList::PdfList(int itemsPerBlock)
{
    fItemsPerBlock = itemsPerBlock;
    fBlockSize = 0;

    fObjectList = NULL;

    fBlockSize = 0;
    fCount = 0;
    PDF_DEBUG_PRINT(("++ [%x] PdfList new.\n", (int)this));
}

PdfList::~PdfList()
{
    PDF_DEBUG_PRINT(("++ [%x] fObjectList delete[].\n", (int)fObjectList));
    delete[] fObjectList;
    PDF_DEBUG_PRINT(("++ [%x] PdfList delete.\n", (int)this));
}

bool
PdfList::AddItem(void* item)
{
    if (fCount < fBlockSize)
        fObjectList[fCount++] = item;
    else {
        if (Resize(fBlockSize + fItemsPerBlock) == false)
            return false;
        fObjectList[fCount++] = item;
    }
    return true;
}

bool
PdfList::AddItem(void* item, int atIndex)
{
    void** pos;

    if (fCount < atIndex || atIndex < 0)
        return false;
    if (fCount == atIndex)
        AddItem(item);
    else {
        if (fCount == fBlockSize)
            if (Resize(fBlockSize + fItemsPerBlock) == false)
                return false;
        pos = fObjectList;
        pos += atIndex;
        memmove(pos + 1, pos, (fCount - atIndex) * sizeof(void*));
        fObjectList[atIndex] = item;
    }
    return true;
}

void*
PdfList::RemoveItem(int index)
{
    if (fCount <= index || index < 0)
        return NULL;

    void* tmpReturn = fObjectList[index];
    if (index != fCount - 1)
        memmove(fObjectList + index, fObjectList + index + 1,
           (fCount - index - 1) * sizeof(void*));
    fCount--;
    return tmpReturn;
}

int
PdfList::IndexOf(const void* item)
{
    if (fObjectList == NULL) 
        return -1;

    void** fObject = fObjectList;

    for (int i = 0; i < fCount; i++)
        if (*fObject++ == item)
            return i;
    return -1;
}

bool
PdfList::Resize(int count)
{
    // resizing the size of the list to the specified value.
    void** newList;
        
    try {
        newList = (void**)(new void*[count]);
        PDF_DEBUG_PRINT(("++ [%x] fObjectList new[].\n", (int)newList));
    } catch (ALLOC_ERROR& e) {
        return false;
    }
    if (fObjectList != NULL) {
        if (fBlockSize < count)
            memcpy(newList, fObjectList, fBlockSize * sizeof(void*));
        else
            memcpy(newList, fObjectList, count * sizeof(void*));
    }
    fBlockSize = count;
    if (fObjectList != NULL) {
        PDF_DEBUG_PRINT(("++ [%x] fObjectList delete[].\n", (int)fObjectList));
        delete[] fObjectList;
    }
    fObjectList = newList;
    return true;
}

/*----------------------------------------------------------------------------*/
