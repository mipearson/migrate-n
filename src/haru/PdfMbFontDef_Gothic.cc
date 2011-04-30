/*
 * << H a r u free pdf library >> -- PdfMbFontDef_Gothic.cc
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

#include "libharu_jpfonts.h"

/*---------------------------------------------------------------------------*/

/* PGothic widths array from CID 1 to 95 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY1[] = {
    305, 219, 500, 500, 500, 500, 594, 203, 305, 305, 
    500, 500, 203, 500, 203, 500, 500, 500, 500, 500, 
    500, 500, 500, 500, 500, 500, 203, 203, 500, 500, 
    500, 453, 668, 633, 637, 664, 648, 566, 551, 680, 
    641, 246, 543, 598, 539, 742, 641, 707, 617, 707, 
    625, 602, 590, 641, 633, 742, 602, 590, 566, 336, 
    504, 336, 414, 305, 414, 477, 496, 500, 496, 500, 
    305, 461, 500, 211, 219, 461, 211, 734, 500, 508, 
    496, 496, 348, 461, 352, 500, 477, 648, 461, 477, 
    457, 234, 234, 234, 414 
};

/* PGothic widths array from CID 326 to 389 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY2[] = {
    500, 441, 441, 441, 441, 441, 547, 523, 445, 480, 
    469, 516, 523, 504, 438, 500, 641, 617, 566, 625, 
    598, 637, 563, 652, 539, 621, 523, 664, 590, 637, 
    645, 555, 527, 602, 602, 602, 461, 645, 598, 578, 
    648, 492, 637, 516, 547, 613, 641, 605, 453, 660, 
    508, 609, 664, 641, 520, 559, 512, 656, 566, 559, 
    590, 563, 250, 230
};

/* PGothic widths array from CID 633 to 695 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY3[] = {
    664, 664, 664, 664, 664, 500, 500, 500, 1000, 1000, 
    500, 500, 500, 500, 500, 500, 1000, 1000, 746, 746, 
    734, 699, 1000, 1000, 1000, 1000, 1000, 961, 1000, 500, 
    1000, 1000, 1000, 1000, 1000, 1000, 1000, 500, 500, 500, 
    500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 
    500, 500, 500, 500, 500, 500, 500, 500, 500, 1000, 
    1000, 1000, 1000
};

/* PGothic widths array from CID 771 to 778 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY4[] = {
    1000, 1000, 1000, 1000, 1000, 500, 500, 500
};

/* PGothic widths array from CID 780 to 789 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY5[] = {
    684, 684, 684, 684, 684, 684, 684, 684, 684, 684
};

/* PGothic widths array from CID 790 to 815 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY6[] = {
    715, 777, 742, 758, 711, 633, 773, 770, 273, 605, 
    754, 629, 934, 770, 805, 711, 805, 758, 742, 617, 
    770, 715, 980, 652, 648, 648
};

/* PGothic widths array from CID 816 to 831 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY7[] = {
    574, 602, 563, 602, 563, 297, 578, 621, 250, 250, 
    594, 250, 938, 621, 605, 605, 602, 379, 570, 336, 
    621, 512, 777, 520, 496, 508
};

/* PGothic widths array from CID 842 to 924 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY8[] = {
    746, 941, 805, 945, 602, 707, 750, 902, 805, 945, 
    1000, 1000, 844, 902, 590, 816, 945, 980, 797, 895, 
    766, 883, 766, 766, 961, 980, 1000, 1000, 922, 961, 
    922, 922, 863, 902, 805, 953, 957, 902, 902, 766, 
    883, 902, 941, 1000, 1000, 1000, 1000, 1000, 1000, 961, 
    961, 961, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 
    1000, 891, 1000, 980, 980, 805, 844, 1000, 844, 980, 
    727, 863, 805, 746, 863, 1000, 844, 863, 1000, 1000, 
    1000, 855, 961
};

/* PGothic widths array from CID 925 to 987 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY9[] = {
    758, 898, 652, 824, 754, 941, 742, 895, 809, 934, 
    824, 922, 961, 965, 805, 941, 930, 961, 797, 891, 
    1000, 1000, 898, 898, 902, 965, 914, 980, 805, 883, 
    766, 922, 910, 961, 734, 863, 922, 887, 961, 648, 
    707, 941, 910, 824, 930, 707, 1000, 1000, 1000, 766, 
    863, 863, 805, 883, 883, 945, 945, 945, 922, 953, 
    953, 902, 668
};

/* PGothic widths array from CID 988 to 1010 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY10[] = {
    977, 719, 898, 805, 980, 813, 961, 629, 727, 809,
    746, 1000, 852, 863, 766, 941, 1000, 1000, 805, 863, 
    961,727,777
};

/*---------------------------------------------------------------------------*/
/*----- Gothic Font ---------------------------------------------------------*/

void
PdfGothicFontDef::Init()
{
    PDF_DEBUG_PRINT(("PdfGothicFontDef %X::Init().\n", (int)this));

    SetBaseFont("Gothic");
    fAscent = 859;
    fDescent = -141;
    fCapHeight = 859;
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH;
    fFontBBox.left = 0;
    fFontBBox.top = 859;
    fFontBBox.right = 1000;
    fFontBBox.bottom = -141;
    fItalicAngle = 0;
    fStemV = 78;
    fMissingWidth = 500;
    AddWidths1(231, 631, 500);
}

/*---------------------------------------------------------------------------*/
/*----- Gothic Font (Bold) --------------------------------------------------*/

void
PdfGothicBoldFontDef::Init()
{
    PdfGothicFontDef::Init();

    SetBaseFont("Gothic,Bold");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_FOURCE_BOLD;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/
/*----- Gothic Font (Italic) ------------------------------------------------*/

void
PdfGothicItalicFontDef::Init()
{
    PdfGothicFontDef::Init();

    SetBaseFont("Gothic,Italic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH;
    fItalicAngle = -11;
}

/*---------------------------------------------------------------------------*/
/*----- Gothic Font (Bold-Italic) -------------------------------------------*/

void
PdfGothicBoldItalicFontDef::Init()
{
    PdfGothicFontDef::Init();

    SetBaseFont("Gothic,BoldItalic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_FOURCE_BOLD;
    fItalicAngle = -11;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/
/*----- Proporsional Gothic Font --------------------------------------------*/

void
PdfPGothicFontDef::Init()
{
    PDF_DEBUG_PRINT(("PdfPGothicFontDef %X::Init().\n", (int)this));

    SetBaseFont("PGothic");
    fAscent = 859;
    fDescent = -141;
    fCapHeight = 859;
    fFlags = PDF_FONT_SYMBOLIC;
    fFontBBox.left = 0;
    fFontBBox.top = 859;
    fFontBBox.right = 1000;
    fFontBBox.bottom = -141;
    fItalicAngle = 0;
    fStemV = 78;
    fMissingWidth = 500;
    fDW = 1000;

    AddWidths2(1, 94, PGOTHIC_FONT_WIDTH_ARRAY1);
    AddWidths2(326, 389, PGOTHIC_FONT_WIDTH_ARRAY2);
    AddWidths2(633, 695, PGOTHIC_FONT_WIDTH_ARRAY3);
    AddWidths2(771, 778, PGOTHIC_FONT_WIDTH_ARRAY4);
    AddWidths2(780, 789, PGOTHIC_FONT_WIDTH_ARRAY5);
    AddWidths2(790, 815, PGOTHIC_FONT_WIDTH_ARRAY6);
    AddWidths2(816, 841, PGOTHIC_FONT_WIDTH_ARRAY7);
    AddWidths2(842, 924, PGOTHIC_FONT_WIDTH_ARRAY8);
    AddWidths2(925, 987, PGOTHIC_FONT_WIDTH_ARRAY9);
    AddWidths2(988, 1010, PGOTHIC_FONT_WIDTH_ARRAY10);
}

/*---------------------------------------------------------------------------*/
/*----- Proporsional Gothic Font (Bold) -------------------------------------*/

void
PdfPGothicBoldFontDef::Init()
{
    PdfPGothicFontDef::Init();

    SetBaseFont("PGothic,Bold");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FOURCE_BOLD;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/
/*----- Proporsional Gothic Font (Italic) -----------------------------------*/

void
PdfPGothicItalicFontDef::Init()
{
    PdfPGothicFontDef::Init();

    SetBaseFont("PGothic,Italic");
    fFlags = PDF_FONT_SYMBOLIC;
    fItalicAngle = -11;
}

/*---------------------------------------------------------------------------*/
/*----- Proporsional Gothic Font (BoldItalic) -------------------------------*/

void
PdfPGothicBoldItalicFontDef::Init()
{
    PdfPGothicFontDef::Init();

    SetBaseFont("PGothic,BoldItalic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FOURCE_BOLD;
    fItalicAngle = -11;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/

