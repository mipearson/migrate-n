/*
 * << H a r u free pdf library >> -- PdfMbFontDef_Mincyo.cc
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

/* PMincyo widths array from CID 1 to 95 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY1[] = {
    305, 305, 461, 500, 500, 500, 613, 305, 305, 305, 
    500, 500, 305, 500, 305, 500, 500, 500, 500, 500, 
    500, 500, 500, 500, 500, 500, 305, 305, 500, 500, 
    500, 500, 727, 664, 621, 699, 691, 598, 598, 711, 
    723, 289, 387, 668, 586, 801, 664, 766, 563, 766, 
    602, 504, 625, 691, 664, 871, 656, 625, 563, 332, 
    500, 332, 305, 305, 305, 453, 500, 465, 500, 473, 
    254, 473, 500, 242, 242, 492, 242, 703, 500, 500, 
    500, 500, 367, 414, 352, 500, 477, 602, 469, 477, 
    453, 242, 219, 242, 500
};

/* PMincyo widths array from CID 326 to 389 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY2[] = {
    500, 441, 441, 441, 441, 441, 547, 523, 445, 480, 
    469, 516, 523, 504, 438, 500, 641, 617, 566, 625, 
    598, 637, 563, 652, 539, 621, 523, 664, 590, 637, 
    645, 555, 527, 602, 602, 602, 461, 645, 598, 578, 
    648, 492, 637, 516, 547, 613, 641, 605, 453, 660, 
    508, 609, 664, 641, 520, 559, 512, 656, 566, 559, 
    590, 563, 250, 230
};

/* PMincyo widths array from CID 633 to 695 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY3[] = {
    664, 664, 664, 664, 664, 500, 500, 500, 1000, 1000, 
    500, 500, 500, 500, 500, 500, 1000, 1000, 648, 801, 
    652, 703, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 500, 
    1000, 1000, 1000, 1000, 1000, 1000, 1000, 500, 500, 500, 
    500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 
    500, 500, 500, 500, 500, 500, 500, 500, 500, 1000, 
    1000, 1000, 1000 
};

/* PMincyo widths array from CID 771 to 778 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY4[] = {
    1000, 1000, 1000, 1000, 1000, 500, 500, 500
};

/* PMincyo widths array from CID 780 to 789 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY5[] = {
    621, 621, 621, 621, 621, 621, 621, 621, 621, 621
};

/* PMincyo widths array from CID 790 to 815 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY6[] = {
    805, 715, 762, 813, 719, 688, 801, 859, 359, 359, 805, 
    676, 1000, 836, 832, 680, 832, 727, 688, 719, 855, 
    770, 977, 730, 777, 656
};

/* PMincyo widths array from CID 816 to 831 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY7[] = {
    531, 574, 531, 574, 543, 387, 559, 613, 293, 293, 
    570, 293, 875, 613, 574, 574, 574, 414, 469, 422, 
    613, 543, 781, 574, 563, 500
};

/* PMincyo widths array from CID 842 to 924 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY8[] = {
    754, 883, 750, 953, 508, 617, 703, 898, 801, 1000, 
    945, 949, 793, 895, 645, 805, 914, 980, 754, 867, 
    754, 883, 777, 777, 1000, 1000, 1000, 1000, 922, 961, 
    906, 949, 902, 902, 855, 1000, 1000, 902, 941, 703, 
    844, 902, 949, 1000, 1000, 949, 969, 1000, 1000, 1000, 
    1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 945, 980, 
    980, 824, 1000, 1000, 953, 758, 875, 1000, 836, 934, 
    688, 785, 766, 641, 793, 984, 863, 801, 953, 945, 
    984, 855, 945
};

/* PMincyo widths array from CID 925 to 987 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY9[] = {
    750, 891, 668, 777, 707, 801, 805, 941, 809, 941, 
    879, 926, 887, 902, 746, 883, 883, 934, 793, 863, 
    953, 961, 902, 902, 820, 902, 930, 949, 754, 855, 
    785, 910, 965, 945, 734, 848, 922, 902, 1000, 590, 
    707, 973, 910, 805, 922, 699, 977, 977, 977, 656, 
    852, 844, 844, 945, 945, 1000, 1000, 1000, 883, 922, 
    922, 930, 609
};

/* PMincyo widths array from CID 988 to 1010 */
static const unsigned int PGOTHIC_FONT_WIDTH_ARRAY10[] = {
    863, 676, 941, 789, 926, 793, 941, 598, 703, 766, 
    609, 980, 832, 785, 699, 805, 965, 961, 785, 863, 
    883, 695, 766
};

/*---------------------------------------------------------------------------*/
/*----- Mincyo Font ---------------------------------------------------------*/

void
PdfMincyoFontDef::Init()
{
    PDF_DEBUG_PRINT(("PdfMincyoFontDef %X::Init().\n", (int)this));

    SetBaseFont("Mincyo");
    fAscent = 859;
    fDescent = -141;
    fCapHeight = 859;
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_SERIF;
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
/*----- Mincyo Font (Bold) --------------------------------------------------*/

void
PdfMincyoBoldFontDef::Init()
{
    PdfMincyoFontDef::Init();

    SetBaseFont("Mincyo,Bold");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + 
            PDF_FONT_SERIF + PDF_FONT_FOURCE_BOLD;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/
/*----- Mincyo Font (Italic) ------------------------------------------------*/

void
PdfMincyoItalicFontDef::Init()
{
    PdfMincyoFontDef::Init();

    SetBaseFont("Mincyo,Italic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + PDF_FONT_SERIF;
    fItalicAngle = -11;
}

/*---------------------------------------------------------------------------*/
/*----- Mincyo Font (Bold-Italic) -------------------------------------------*/

void
PdfMincyoBoldItalicFontDef::Init()
{
    PdfMincyoFontDef::Init();

    SetBaseFont("Mincyo,BoldItalic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_FIXED_WIDTH + 
            PDF_FONT_SERIF + PDF_FONT_FOURCE_BOLD;
    fItalicAngle = -11;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/
/*----- Proporsional Mincyo Font --------------------------------------------*/

void
PdfPMincyoFontDef::Init()
{
    PDF_DEBUG_PRINT(("PdfPMincyoFontDef %X::Init().\n", (int)this));

    SetBaseFont("PMincyo-Medium");
    fAscent = 859;
    fDescent = -141;
    fCapHeight = 859;
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_SERIF;
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
/*----- Proporsional Mincyo Font (Bold) -------------------------------------*/

void
PdfPMincyoBoldFontDef::Init()
{
    PdfPMincyoFontDef::Init();

    SetBaseFont("PMincyo,Bold");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_SERIF + PDF_FONT_FOURCE_BOLD;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/
/*----- Proporsional Mincyo Font (Italic) -----------------------------------*/

void
PdfPMincyoItalicFontDef::Init()
{
    PdfPMincyoFontDef::Init();

    SetBaseFont("PMincyo,Italic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_SERIF;
    fItalicAngle = -11;
}

/*---------------------------------------------------------------------------*/
/*----- Proporsional Mincyo Font (Bold-Italic) ------------------------------*/

void
PdfPMincyoBoldItalicFontDef::Init()
{
    PdfPMincyoFontDef::Init();

    SetBaseFont("PMincyo,BoldItalic");
    fFlags = PDF_FONT_SYMBOLIC + PDF_FONT_SERIF + PDF_FONT_FOURCE_BOLD;
    fItalicAngle = -11;
    fStemV = 156;
}

/*---------------------------------------------------------------------------*/

