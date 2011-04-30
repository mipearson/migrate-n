/*
 * << H a r u --free pdf library >> -- PdfLList.cpp
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
#include <stdlib.h>
#include <stdio.h>

/*----- PdfLNode class --------------------------------------------------------*/

PdfLNode::PdfLNode()
{
    fPrev = NULL;
    fNext = NULL;
    fParent = NULL;
}

PdfLNode::~PdfLNode()
{
    if (fParent != NULL)
        fParent->Unlink(this);
}

/*----- PdfLList class --------------------------------------------------------*/

PdfLList::PdfLList()
{
    fFirstNode = NULL;
    fLastNode = NULL;
}

PdfLList::~PdfLList()
{
    while (fFirstNode != NULL)
        delete fFirstNode;
}

bool
PdfLList::Add(PdfLNode* node, PdfLNode* atNode)
{
    if (node == NULL || node->fParent != NULL)
        return false;

    if (atNode != NULL && atNode->fParent != this)
        return false;

    node->fParent = this;

    if (IsEmpty()) {
        fFirstNode = node;
        fLastNode = node;
        return true;
    }

    /* if atNode is NULL, the node add to the last of the list, otherwise the
       node is insert in front of atNode */

    if (atNode != NULL) {
        node->fNext = atNode;
        node->fPrev = atNode->fPrev;
        if (atNode->fPrev != NULL)
            atNode->fPrev->fNext = node;
        atNode->fPrev = node;
        if (atNode == fFirstNode)
            fFirstNode = node;
    } else {
        node->fPrev = fLastNode;
        if (fLastNode != NULL)
            fLastNode->fNext = node;
        fLastNode = node;
    }
    return true;
}

bool
PdfLList::Unlink(PdfLNode* node)
{
    if (node == NULL)
        return false;

    if (node->fPrev != NULL)
        node->fPrev->fNext = node->fNext;
    if (node->fNext != NULL)
        node->fNext->fPrev = node->fPrev;
    if (node == fFirstNode)
        fFirstNode = node->fNext;
    if (node == fLastNode)
        fLastNode = node->fPrev;

    return true;
}

int
PdfLList::CountItems()
{
    int cnt = 0;
    PdfLNode* node = fFirstNode;

    while (node != NULL) {
        node = node->GetNext();
        cnt++;
    }

    return cnt;
}

void
PdfLList::Clear()
{
    while (fFirstNode != NULL)
        delete fFirstNode;
}

/*----------------------------------------------------------------------------*/


