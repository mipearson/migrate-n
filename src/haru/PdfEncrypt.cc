/*
 * << H a r u --free pdf library >> -- PdfEncryptDict.cpp
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
 *------------------------------------------------------------------------------
 * ARC4 encryption is based on the code of OPEN-BSD. The copyright of an 
 * original code is as follows.
 * 
 * Copyright (c) 2003 Markus Friedl <markus@openbsd.org>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *------------------------------------------------------------------------------
 *
 * The code implements the MD5 message-digest algorithm based on the code 
 * written by Colin Plumb.
 * The copyright of an is as follows.
 *
 * This code implements the MD5 message-digest algorithm.
 * The algorithm is due to Ron Rivest.  This code was
 * written by Colin Plumb in 1993, no copyright is claimed.
 * This code is in the public domain; do with it what you wish.
 *
 * Equivalent code is available from RSA Data Security, Inc.
 * This code has been tested against that, and is equivalent,
 * except that you don't need to include two pages of legalese
 * with every copy.
 *
 * To compute the message digest of a chunk of bytes, declare an
 * MD5Context structure, pass it to MD5Init, call MD5Update as
 * needed on buffers full of bytes, and then call MD5Final, which
 * will fill a supplied 16-byte array with the digest.
 * 
 *------------------------------------------------------------------------------
 */

#ifdef USE_ENCRYPTION

#include "libharu.h"
#include "string.h"

#ifndef pdf_uint32
typedef unsigned long pdf_uint32;
#endif

static const unsigned char PADDING_STRING[] = {
    0x28, 0xBF, 0x4E, 0x5E, 0x4E, 0x75, 0x8A, 0x41,
    0x64, 0x00, 0x4E, 0x56, 0xFF, 0xFA, 0x01, 0x08,
    0x2E, 0x2E, 0x00, 0xB6, 0xD0, 0x68, 0x3E, 0x80,
    0x2F, 0x0C, 0xA9, 0xFE, 0x64, 0x53, 0x69, 0x7A
};

struct MD5Context
{
    pdf_uint32 buf[4];
    pdf_uint32 bits[2];
    unsigned char in[64];
};
typedef struct MD5Context MD5_CTX;

#define RC4STATE 256
#define RC4KEYLEN 16

struct rc4_ctx {
	pdf_uint8 x, y;
	pdf_uint8 state[RC4STATE];
};

#define RC4SWAP(x,y) \
	do { \
		pdf_uint8 t = ctx->state[x];  \
		ctx->state[x] = ctx->state[y]; \
		ctx->state[y] = t; \
	} while(0)


void MD5Init (struct MD5Context *ctx);
void MD5Update (struct MD5Context *ctx, unsigned char const *buf, 
        pdf_uint32 len);
void MD5Final (unsigned char digest[16], struct MD5Context *ctx);
void MD5Transform (pdf_uint32 buf[4], pdf_uint32 const in[16]);
void RC4KeySetup(struct rc4_ctx *, pdf_uint8 *, pdf_uint32);
void RC4Crypt(struct rc4_ctx *, pdf_uint8 *, pdf_uint8 *, pdf_uint32);

void MD5ByteReverse (unsigned char *buf, pdf_uint32 longs);

void MD5ByteReverse (unsigned char *buf, pdf_uint32 longs)
{
    pdf_uint32 t;
    do
    {
        t = (pdf_uint32) ((pdf_uint32) buf[3] << 8 | buf[2]) << 16 |
        ((pdf_uint32) buf[1] << 8 | buf[0]);
        *(pdf_uint32 *) buf = t;
        buf += 4;
    }
    while (--longs);
}

void MD5Init (struct MD5Context *ctx)
{
    ctx->buf[0] = 0x67452301;
    ctx->buf[1] = 0xefcdab89;
    ctx->buf[2] = 0x98badcfe;
    ctx->buf[3] = 0x10325476;

    ctx->bits[0] = 0;
    ctx->bits[1] = 0;
}

void MD5Update (struct MD5Context *ctx, unsigned char const *buf,
     pdf_uint32 len)
{
    pdf_uint32 t;

    /* Update bitcount */

    t = ctx->bits[0];
    if ((ctx->bits[0] = t + ((pdf_uint32) len << 3)) < t)
        ctx->bits[1]++;     /* Carry from low to high */
    ctx->bits[1] += len >> 29;

    t = (t >> 3) & 0x3f; /* Bytes already in shsInfo->data */

    /* Handle any leading odd-sized chunks */

    if (t) {
        unsigned char *p = (unsigned char *) ctx->in + t;

        t = 64 - t;
        if (len < t)
        {
            memcpy (p, buf, len);
            return;
        }
        memcpy (p, buf, t);
        MD5ByteReverse (ctx->in, 16);
        MD5Transform (ctx->buf, (pdf_uint32 *) ctx->in);
        buf += t;
        len -= t;
    }
    /* Process data in 64-byte chunks */

    while (len >= 64) {
        memcpy (ctx->in, buf, 64);
        MD5ByteReverse (ctx->in, 16);
        MD5Transform (ctx->buf, (pdf_uint32 *) ctx->in);
        buf += 64;
        len -= 64;
    }

    /* Handle any remaining bytes of data. */

    memcpy (ctx->in, buf, len);
}

/*
 * Final wrapup - pad to 64-byte boundary with the bit pattern 
 * 1 0* (64-bit count of bits processed, MSB-first)
 */
void MD5Final (unsigned char digest[16], struct MD5Context *ctx)
{
    pdf_uint32 count;
    unsigned char *p;

    /* Compute number of bytes mod 64 */
    count = (ctx->bits[0] >> 3) & 0x3F;

    /* Set the first char of padding to 0x80.  This is safe since there is
       always at least one byte free */
    p = ctx->in + count;
    *p++ = 0x80;

    /* Bytes of padding needed to make 64 bytes */
    count = 64 - 1 - count;

    /* Pad out to 56 mod 64 */
    if (count < 8) {
        /* Two lots of padding:  Pad the first block to 64 bytes */
        memset (p, 0, count);
        MD5ByteReverse (ctx->in, 16);
        MD5Transform (ctx->buf, (pdf_uint32 *) ctx->in);

        /* Now fill the next block with 56 bytes */
        memset (ctx->in, 0, 56);
    } else {
        /* Pad block to 56 bytes */
        memset (p, 0, count - 8);
    }
    MD5ByteReverse (ctx->in, 14);

    /* Append length in bits and transform */
    ((pdf_uint32 *) ctx->in)[14] = ctx->bits[0];
    ((pdf_uint32 *) ctx->in)[15] = ctx->bits[1];

    MD5Transform (ctx->buf, (pdf_uint32 *) ctx->in);
    MD5ByteReverse ((unsigned char *) ctx->buf, 4);
    memcpy (digest, ctx->buf, 16);
    memset (ctx, 0, sizeof (ctx));   /* In case it's sensitive */
}

/* The four core functions - F1 is optimized somewhat */

/* #define F1(x, y, z) (x & y | ~x & z) */
#define F1(x, y, z) (z ^ (x & (y ^ z)))
#define F2(x, y, z) F1(z, x, y)
#define F3(x, y, z) (x ^ y ^ z)
#define F4(x, y, z) (y ^ (x | ~z))

/* This is the central step in the MD5 algorithm. */
#define MD5STEP(f, w, x, y, z, data, s) \
 ( w += f(x, y, z) + data,  w = w<<s | w>>(32-s),  w += x )

/*
 * The core of the MD5 algorithm, this alters an existing MD5 hash to
 * reflect the addition of 16 longwords of new data.  MD5Update blocks
 * the data and converts bytes into longwords for this routine.
 */
void MD5Transform (pdf_uint32 buf[4], pdf_uint32 const in[16])
{
    register pdf_uint32 a, b, c, d;

    a = buf[0];
    b = buf[1];
    c = buf[2];
    d = buf[3];

    MD5STEP (F1, a, b, c, d, in[0] + 0xd76aa478, 7);
    MD5STEP (F1, d, a, b, c, in[1] + 0xe8c7b756, 12);
    MD5STEP (F1, c, d, a, b, in[2] + 0x242070db, 17);
    MD5STEP (F1, b, c, d, a, in[3] + 0xc1bdceee, 22);
    MD5STEP (F1, a, b, c, d, in[4] + 0xf57c0faf, 7);
    MD5STEP (F1, d, a, b, c, in[5] + 0x4787c62a, 12);
    MD5STEP (F1, c, d, a, b, in[6] + 0xa8304613, 17);
    MD5STEP (F1, b, c, d, a, in[7] + 0xfd469501, 22);
    MD5STEP (F1, a, b, c, d, in[8] + 0x698098d8, 7);
    MD5STEP (F1, d, a, b, c, in[9] + 0x8b44f7af, 12);
    MD5STEP (F1, c, d, a, b, in[10] + 0xffff5bb1, 17);
    MD5STEP (F1, b, c, d, a, in[11] + 0x895cd7be, 22);
    MD5STEP (F1, a, b, c, d, in[12] + 0x6b901122, 7);
    MD5STEP (F1, d, a, b, c, in[13] + 0xfd987193, 12);
    MD5STEP (F1, c, d, a, b, in[14] + 0xa679438e, 17);
    MD5STEP (F1, b, c, d, a, in[15] + 0x49b40821, 22);

    MD5STEP (F2, a, b, c, d, in[1] + 0xf61e2562, 5);
    MD5STEP (F2, d, a, b, c, in[6] + 0xc040b340, 9);
    MD5STEP (F2, c, d, a, b, in[11] + 0x265e5a51, 14);
    MD5STEP (F2, b, c, d, a, in[0] + 0xe9b6c7aa, 20);
    MD5STEP (F2, a, b, c, d, in[5] + 0xd62f105d, 5);
    MD5STEP (F2, d, a, b, c, in[10] + 0x02441453, 9);
    MD5STEP (F2, c, d, a, b, in[15] + 0xd8a1e681, 14);
    MD5STEP (F2, b, c, d, a, in[4] + 0xe7d3fbc8, 20);
    MD5STEP (F2, a, b, c, d, in[9] + 0x21e1cde6, 5);
    MD5STEP (F2, d, a, b, c, in[14] + 0xc33707d6, 9);
    MD5STEP (F2, c, d, a, b, in[3] + 0xf4d50d87, 14);
    MD5STEP (F2, b, c, d, a, in[8] + 0x455a14ed, 20);
    MD5STEP (F2, a, b, c, d, in[13] + 0xa9e3e905, 5);
    MD5STEP (F2, d, a, b, c, in[2] + 0xfcefa3f8, 9);
    MD5STEP (F2, c, d, a, b, in[7] + 0x676f02d9, 14);
    MD5STEP (F2, b, c, d, a, in[12] + 0x8d2a4c8a, 20);

    MD5STEP (F3, a, b, c, d, in[5] + 0xfffa3942, 4);
    MD5STEP (F3, d, a, b, c, in[8] + 0x8771f681, 11);
    MD5STEP (F3, c, d, a, b, in[11] + 0x6d9d6122, 16);
    MD5STEP (F3, b, c, d, a, in[14] + 0xfde5380c, 23);
    MD5STEP (F3, a, b, c, d, in[1] + 0xa4beea44, 4);
    MD5STEP (F3, d, a, b, c, in[4] + 0x4bdecfa9, 11);
    MD5STEP (F3, c, d, a, b, in[7] + 0xf6bb4b60, 16);
    MD5STEP (F3, b, c, d, a, in[10] + 0xbebfbc70, 23);
    MD5STEP (F3, a, b, c, d, in[13] + 0x289b7ec6, 4);
    MD5STEP (F3, d, a, b, c, in[0] + 0xeaa127fa, 11);
    MD5STEP (F3, c, d, a, b, in[3] + 0xd4ef3085, 16);
    MD5STEP (F3, b, c, d, a, in[6] + 0x04881d05, 23);
    MD5STEP (F3, a, b, c, d, in[9] + 0xd9d4d039, 4);
    MD5STEP (F3, d, a, b, c, in[12] + 0xe6db99e5, 11);
    MD5STEP (F3, c, d, a, b, in[15] + 0x1fa27cf8, 16);
    MD5STEP (F3, b, c, d, a, in[2] + 0xc4ac5665, 23);

    MD5STEP (F4, a, b, c, d, in[0] + 0xf4292244, 6);
    MD5STEP (F4, d, a, b, c, in[7] + 0x432aff97, 10);
    MD5STEP (F4, c, d, a, b, in[14] + 0xab9423a7, 15);
    MD5STEP (F4, b, c, d, a, in[5] + 0xfc93a039, 21);
    MD5STEP (F4, a, b, c, d, in[12] + 0x655b59c3, 6);
    MD5STEP (F4, d, a, b, c, in[3] + 0x8f0ccc92, 10);
    MD5STEP (F4, c, d, a, b, in[10] + 0xffeff47d, 15);
    MD5STEP (F4, b, c, d, a, in[1] + 0x85845dd1, 21);
    MD5STEP (F4, a, b, c, d, in[8] + 0x6fa87e4f, 6);
    MD5STEP (F4, d, a, b, c, in[15] + 0xfe2ce6e0, 10);
    MD5STEP (F4, c, d, a, b, in[6] + 0xa3014314, 15);
    MD5STEP (F4, b, c, d, a, in[13] + 0x4e0811a1, 21);
    MD5STEP (F4, a, b, c, d, in[4] + 0xf7537e82, 6);
    MD5STEP (F4, d, a, b, c, in[11] + 0xbd3af235, 10);
    MD5STEP (F4, c, d, a, b, in[2] + 0x2ad7d2bb, 15);
    MD5STEP (F4, b, c, d, a, in[9] + 0xeb86d391, 21);

    buf[0] += a;
    buf[1] += b;
    buf[2] += c;
    buf[3] += d;
}

void
RC4KeySetup(struct rc4_ctx *ctx, pdf_uint8 *key, pdf_uint32 klen)
{
	pdf_uint8 x, y;
	pdf_uint32 i;

	x = y = 0;
	for (i = 0; i < RC4STATE; i++)
		ctx->state[i] = i;
	for (i = 0; i < RC4STATE; i++) {
		y = (key[x] + ctx->state[i] + y) % RC4STATE;
		RC4SWAP(i, y);
		x = (x + 1) % klen;
	}
	ctx->x = ctx->y = 0;
}

void
RC4Crypt(struct rc4_ctx *ctx, pdf_uint8 *src, pdf_uint8 *dst,
    pdf_uint32 len)
{
	pdf_uint32 i;

	for (i = 0; i < len; i++) {
		ctx->x = (ctx->x + 1) % RC4STATE;
		ctx->y = (ctx->state[ctx->x] + ctx->y) % RC4STATE;
		RC4SWAP(ctx->x, ctx->y);
		dst[i] = src[i] ^ ctx->state[
		   (ctx->state[ctx->x] + ctx->state[ctx->y]) % RC4STATE];
	}
}

/*------ PdfEncryptDict ------------------------------------------------------*/

PdfEncryptDict::PdfEncryptDict(PdfXref* xref, PdfInfo* info)
        : PdfDictionary(xref)
{
    memcpy(fOwnerPasswd, PADDING_STRING, PDF_PASSWD_LEN);
    memcpy(fUserPasswd, PADDING_STRING, PDF_PASSWD_LEN);
    fPermission = PDF_ENABLE_PRINT | PDF_ENABLE_EDIT_ALL |
        PDF_ENABLE_COPY | PDF_ENABLE_EDIT | PDF_PERMISSION_PAD;
    fInfo = info;
    memset(fFileKey, 0x00, PDF_ENCRYPT_KEY_LEN);
    memset(fOwnerPasswdValue, 0x00, PDF_PASSWD_LEN);
    memset(fUserPasswdValue, 0x00, PDF_PASSWD_LEN);
    memset(fID, 0x00, PDF_ID_LEN);
}

void
PdfEncryptDict::SetPassword(const char* owner_passwd, const char* user_passwd)
{
    if (owner_passwd == NULL)
        throw PdfException(PDF_ERR_INVALID_OPERATION, "ERROR: "
            "PdfEncrypt::SetPassword -- owner_passwd cannot be set to NULL");

    if (user_passwd == NULL)
        throw PdfException(PDF_ERR_INVALID_OPERATION, "ERROR: "
            "PdfEncrypt::SetPassword -- user_passwd cannot be set to NULL");

    if (strcmp(owner_passwd, user_passwd) == 0)
        throw PdfException(PDF_ERR_INVALID_OPERATION, 
            "ERROR: PdfEncrypt::SetPassword -- user and owner password "
            "must be diiferent.");
    
    PadOrTruncatePasswd(owner_passwd, fOwnerPasswd);
    PadOrTruncatePasswd(user_passwd, fUserPasswd);
}

void
PdfEncryptDict::CreateID()
{
    MD5Context ctx;

    MD5Init(&ctx);
    
    /* create File Identifier from various elements. */
    time_t t = time(NULL);
    MD5Update(&ctx, (unsigned char*)&t, sizeof(t));

    if (fInfo != NULL) {
        int len;
        const char* author = fInfo->Author();
        if (author != NULL && (len = strlen(author)) > 0)
            MD5Update(&ctx, (const unsigned char*)author, len);
        
        const char* creator = fInfo->Creator();
        if (creator != NULL && (len = strlen(creator)) > 0)
            MD5Update(&ctx, (const unsigned char*)creator, len);

        const char* subject = fInfo->Subject();
        if (subject != NULL && (len = strlen(subject)) > 0)
            MD5Update(&ctx, (const unsigned char*)subject, len);

        const char* keywords = fInfo->Keywords();
        if (keywords != NULL && (len = strlen(keywords)) > 0)
            MD5Update(&ctx, (const unsigned char*)keywords, len);

        const char* producer = fInfo->Producer();
        if (producer != NULL && (len = strlen(producer)) > 0)
            MD5Update(&ctx, (const unsigned char*)producer, len);
    }

    int cnt = GetXref()->GetCount();
    MD5Update(&ctx, (const unsigned char*)&cnt, sizeof(cnt));

    MD5Final(fID, &ctx);
}

void
PdfEncryptDict::CreateUserKey()
{
    rc4_ctx key;

    /* create U(User passwd) key */
    RC4KeySetup(&key, fFileKey, PDF_ENCRYPT_KEY_LEN);
    RC4Crypt(&key, (pdf_uint8*)PADDING_STRING, fUserPasswdValue, PDF_PASSWD_LEN);
    
    PDF_PRINT_BINARY(fUserPasswd, PDF_PASSWD_LEN, "UserPasswd(Padded)");

    PdfBinary* b = new PdfBinary();
    AddElement("U", b);
    b->SetData(fUserPasswdValue, PDF_PASSWD_LEN, false);
}

void
PdfEncryptDict::CreateOwnerKey()
{
    unsigned char digest[PDF_MD5_KEY_LEN];
    rc4_ctx key;
    MD5_CTX c;
   
    MD5Init(&c);
    MD5Update(&c, fOwnerPasswd, PDF_PASSWD_LEN);
    PDF_PRINT_BINARY(fOwnerPasswd, PDF_PASSWD_LEN, "Owner password(Padded)");
    MD5Final(digest, &c);
    PDF_PRINT_BINARY(digest, PDF_ENCRYPT_KEY_LEN, "Owner password digest");

    RC4KeySetup(&key, digest, PDF_ENCRYPT_KEY_LEN);
    RC4Crypt(&key, fUserPasswd, fOwnerPasswdValue, PDF_PASSWD_LEN);
    PDF_PRINT_BINARY(fOwnerPasswdValue, PDF_PASSWD_LEN, "OwnerPasswdValue");

    PdfBinary* b = new PdfBinary();
    AddElement("O", b);
    b->SetData(fOwnerPasswdValue, PDF_PASSWD_LEN, false);
}

void
PdfEncryptDict::CreateFileKey()
{
    MD5_CTX c;

    /* Algorithm3.2 step2 */
    MD5Init(&c);
    MD5Update(&c, fUserPasswd, PDF_PASSWD_LEN);

    /* Algorithm3.2 step3 */
    MD5Update(&c, fOwnerPasswdValue, PDF_PASSWD_LEN);

    /* Algorithm3.2 step4 */
    unsigned char tmp_flg[4];
   
    tmp_flg[0] = fPermission;
    tmp_flg[1] = (fPermission >> 8);
    tmp_flg[2] = (fPermission >> 16);
    tmp_flg[3] = (fPermission >> 24);
    
    MD5Update(&c, tmp_flg, 4);

    /* Algorithm3.2 step5 */
    MD5Update(&c, fID, PDF_ID_LEN);
    PDF_PRINT_BINARY(fID, PDF_ID_LEN, "step1 id");

    MD5Final(fFileKey, &c);
    PDF_PRINT_BINARY(fFileKey, PDF_ENCRYPT_KEY_LEN, "User password digest");
}

void
PdfEncryptDict::Init()
{
    GetXref()->AddObject(this);
    AddElement("Filter", new PdfName("Standard"));
    AddElement("V", new PdfNumber(1));
    AddElement("R", new PdfNumber(2));
}

void
PdfEncryptDict::PadOrTruncatePasswd(const char* passwd, unsigned char* out)
{
    int len = strlen(passwd);

    memset(out, 0x00, PDF_PASSWD_LEN);
    if (len >= PDF_PASSWD_LEN)
        memcpy(out, passwd, PDF_PASSWD_LEN);
    else {
        memcpy(out, passwd, len);
        memcpy(out + len, PADDING_STRING, PDF_PASSWD_LEN - len);
    }
}

void
PdfEncryptDict::EncryptPassword()
{
    CreateID();
    CreateOwnerKey();
    CreateFileKey();
    CreateUserKey();
    AddElement("P", new PdfNumber(fPermission));
}

/*------ PdfEncryptor --------------------------------------------------------*/

PdfEncryptor::PdfEncryptor(const unsigned char* key)
{
    memcpy(fKey, key, PDF_ENCRYPT_KEY_LEN);
    memset(fKey + PDF_ENCRYPT_KEY_LEN, 0x00, 5);
    memset(fDigest, 0x00, PDF_MD5_KEY_LEN);
    fARC4Key = (unsigned char*)new rc4_ctx;

    PDF_PRINT_BINARY(fKey, PDF_ENCRYPT_KEY_LEN + 5, "PdfEncryptor Key: ");
}

PdfEncryptor::~PdfEncryptor()
{
    delete fARC4Key;
}

void
PdfEncryptor::Init(PdfOID id, unsigned int gen_no)
{
    MD5Context ctx;
   
    fKey[5] = id;
    fKey[6] = (id >> 8);
    fKey[7] = (id >> 16);
    fKey[8] = gen_no;
    fKey[9] = (gen_no >> 8);

    MD5Init(&ctx);
    MD5Update(&ctx, fKey, PDF_ENCRYPT_KEY_LEN + 5);
    MD5Final(fDigest, &ctx);
    RC4KeySetup((rc4_ctx*)fARC4Key, fDigest, PDF_ENCRYPT_KEY_LEN + 5);

    PDF_PRINT_BINARY(fKey, PDF_ENCRYPT_KEY_LEN + 5, "PdfEncryptor Init: ");
}

void
PdfEncryptor::Reset()
{
    RC4KeySetup((rc4_ctx*)fARC4Key, fDigest, PDF_ENCRYPT_KEY_LEN + 5);
}

void
PdfEncryptor::CryptBuf(const unsigned char* src,
        unsigned char* dst, unsigned int len)
{
    RC4Crypt((rc4_ctx*)fARC4Key, (pdf_uint8*)src, (pdf_uint8*)dst, len);
    PDF_DEBUG_PRINT(("PdfEncryptor CryptBuf %d bytes\n", len));
}

/*----------------------------------------------------------------------------*/

#endif /* USE_ENCRYPTION */

