/*
 * << H a r u free pdf library >> -- PdfFontDef_Symbol.cc
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

/*----- Symbol Font ---------------------------------------------------------*/

PdfSymbolFontDef::PdfSymbolFontDef()
    : PdfType1FontDef("Symbol")
{
    const pdf_char_data_ro char_data[191] = {
        {32, "space", 250},
        {33, "exclam", 333},
        {34, "universal", 713},
        {35, "numbersign", 500},
        {36, "existential", 549},
        {37, "percent", 833},
        {38, "ampersand", 778},
        {39, "suchthat", 439},
        {40, "parenleft", 333},
        {41, "parenright", 333},
        {42, "asteriskmath", 500},
        {43, "plus", 549},
        {44, "comma", 250},
        {45, "minus", 549},
        {46, "period", 250},
        {47, "slash", 278},
        {48, "zero", 500},
        {49, "one", 500},
        {50, "two", 500},
        {51, "three", 500},
        {52, "four", 500},
        {53, "five", 500},
        {54, "six", 500},
        {55, "seven", 500},
        {56, "eight", 500},
        {57, "nine", 500},
        {58, "colon", 278},
        {59, "semicolon", 278},
        {60, "less", 549},
        {61, "equal", 549},
        {62, "greater", 549},
        {63, "question", 444},
        {64, "congruent", 549},
        {65, "Alpha", 722},
        {66, "Beta", 667},
        {67, "Chi", 722},
        {68, "Delta", 612},
        {69, "Epsilon", 611},
        {70, "Phi", 763},
        {71, "Gamma", 603},
        {72, "Eta", 722},
        {73, "Iota", 333},
        {74, "theta1", 631},
        {75, "Kappa", 722},
        {76, "Lambda", 686},
        {77, "Mu", 889},
        {78, "Nu", 722},
        {79, "Omicron", 722},
        {80, "Pi", 768},
        {81, "Theta", 741},
        {82, "Rho", 556},
        {83, "Sigma", 592},
        {84, "Tau", 611},
        {85, "Upsilon", 690},
        {86, "sigma1", 439},
        {87, "Omega", 768},
        {88, "Xi", 645},
        {89, "Psi", 795},
        {90, "Zeta", 611},
        {91, "bracketleft", 333},
        {92, "therefore", 863},
        {93, "bracketright", 333},
        {94, "perpendicular", 658},
        {95, "underscore", 500},
        {96, "radicalex", 500},
        {97, "alpha", 631},
        {98, "beta", 549},
        {99, "chi", 549},
        {100, "delta", 494},
        {101, "epsilon", 439},
        {102, "phi", 521},
        {103, "gamma", 411},
        {104, "eta", 603},
        {105, "iota", 329},
        {106, "phi1", 603},
        {107, "kappa", 549},
        {108, "lambda", 549},
        {109, "mu", 576},
        {110, "nu", 521},
        {111, "omicron", 549},
        {112, "pi", 549},
        {113, "theta", 521},
        {114, "rho", 549},
        {115, "sigma", 603},
        {116, "tau", 439},
        {117, "upsilon", 576},
        {118, "omega1", 713},
        {119, "omega", 686},
        {120, "xi", 493},
        {121, "psi", 686},
        {122, "zeta", 494},
        {123, "braceleft", 480},
        {124, "bar", 200},
        {125, "braceright", 480},
        {126, "similar", 549},
        {160, "Euro", 750},
        {161, "Upsilon1", 620},
        {162, "minute", 247},
        {163, "lessequal", 549},
        {164, "fraction", 167},
        {165, "infinity", 713},
        {166, "florin", 500},
        {167, "club", 753},
        {168, "diamond", 753},
        {169, "heart", 753},
        {170, "spade", 753},
        {171, "arrowboth", 1042},
        {172, "arrowleft", 987},
        {173, "arrowup", 603},
        {174, "arrowright", 987},
        {175, "arrowdown", 603},
        {176, "degree", 400},
        {177, "plusminus", 549},
        {178, "second", 411},
        {179, "greaterequal", 549},
        {180, "multiply", 549},
        {181, "proportional", 713},
        {182, "partialdiff", 494},
        {183, "bullet", 460},
        {184, "divide", 549},
        {185, "notequal", 549},
        {186, "equivalence", 549},
        {187, "approxequal", 549},
        {188, "ellipsis", 1000},
        {189, "arrowvertex", 603},
        {190, "arrowhorizex", 1000},
        {191, "carriagereturn", 658},
        {192, "aleph", 823},
        {193, "Ifraktur", 686},
        {194, "Rfraktur", 795},
        {195, "weierstrass", 987},
        {196, "circlemultiply", 768},
        {197, "circleplus", 768},
        {198, "emptyset", 823},
        {199, "intersection", 768},
        {200, "union", 768},
        {201, "propersuperset", 713},
        {202, "reflexsuperset", 713},
        {203, "notsubset", 713},
        {204, "propersubset", 713},
        {205, "reflexsubset", 713},
        {206, "element", 713},
        {207, "notelement", 713},
        {208, "angle", 768},
        {209, "gradient", 713},
        {210, "registerserif", 790},
        {211, "copyrightserif", 790},
        {212, "trademarkserif", 890},
        {213, "product", 823},
        {214, "radical", 549},
        {215, "dotmath", 250},
        {216, "logicalnot", 713},
        {217, "logicaland", 603},
        {218, "logicalor", 603},
        {219, "arrowdblboth", 1042},
        {220, "arrowdblleft", 987},
        {221, "arrowdblup", 603},
        {222, "arrowdblright", 987},
        {223, "arrowdbldown", 603},
        {224, "lozenge", 494},
        {225, "angleleft", 329},
        {226, "registersans", 790},
        {227, "copyrightsans", 790},
        {228, "trademarksans", 786},
        {229, "summation", 713},
        {230, "parenlefttp", 384},
        {231, "parenleftex", 384},
        {232, "parenleftbt", 384},
        {233, "bracketlefttp", 384},
        {234, "bracketleftex", 384},
        {235, "bracketleftbt", 384},
        {236, "bracelefttp", 494},
        {237, "braceleftmid", 494},
        {238, "braceleftbt", 494},
        {239, "braceex", 494},
        {241, "angleright", 329},
        {242, "integral", 274},
        {243, "integraltp", 686},
        {244, "integralex", 686},
        {245, "integralbt", 686},
        {246, "parenrighttp", 384},
        {247, "parenrightex", 384},
        {248, "parenrightbt", 384},
        {249, "bracketrighttp", 384},
        {250, "bracketrightex", 384},
        {251, "bracketrightbt", 384},
        {252, "bracerighttp", 494},
        {253, "bracerightmid", 494},
        {254, "bracerightbt", 494},
        {-1, "apple", 790},
        {-1, NULL, 0}
    };
    fFontBBox.left = -180;
    fFontBBox.top = -293;
    fFontBBox.right = 1090;
    fFontBBox.bottom = 1010;
    fAscent = 0;
    fDescent =0;
    fXHeight = 0;
    fIsBase14Font = true;
    fDefaultEncoding = PDF_FONT_SPECIFIC;
    fFirstChar = 32;
    fLastChar = 254;
    
    SetWidths(char_data);
}

/*---------------------------------------------------------------------------*/
/*------ PdfSymbolFontEncoding -----------------------------------------------*/

PdfSymbolFontEncoding::PdfSymbolFontEncoding()
{
    fFirstChar = 32;
    fLastChar = 254;
    fBaseEncoding = PDF_FONT_SPECIFIC;
}

/*----------------------------------------------------------------------------*/

