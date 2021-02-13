/*
  Copyright (c) 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014 Jouni Sirén
  Copyright (c) 2016 Genome Research Ltd.

  Author: Jouni Sirén <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#ifndef _RLCSA_DEFINITIONS_H
#define _RLCSA_DEFINITIONS_H

#include <algorithm>
#include <climits>

namespace CSA
{

//--------------------------------------------------------------------------

typedef uint64_t usint;
typedef int64_t  sint;

inline usint popcount(usint field)
{
  return __builtin_popcountl(field);
}

#ifndef uchar
typedef unsigned char uchar;
#endif

#ifndef uint
typedef uint32_t uint;
#endif

typedef std::pair<usint, usint> pair_type;

//--------------------------------------------------------------------------

inline usint length(usint n)
{
  usint b = 0;
  while(n > 0) { b++; n >>= 1; }
  return b;
}

inline bool isEmpty(const pair_type data)
{
  return (data.first > data.second);
}

inline usint length(const pair_type data)
{
  return data.second + 1 - data.first;
}

inline usint bound(usint value, usint low, usint high)
{
  return std::max(std::min(value, high), low);
}

inline usint nextMultipleOf(usint multiplier, usint value)
{
  return multiplier * ((value / multiplier) + 1);
}

const pair_type EMPTY_PAIR = pair_type(1, 0);

//--------------------------------------------------------------------------

const usint MILLION          = 1000000;
const usint BILLION          = 1000 * MILLION;

const usint KILOBYTE         = 1024;
const usint MEGABYTE         = KILOBYTE * KILOBYTE;
const usint GIGABYTE         = KILOBYTE * MEGABYTE;

const double MILLION_DOUBLE  = 1000000.0;
const double BILLION_DOUBLE  = 1000.0 * MILLION_DOUBLE;

const double KILOBYTE_DOUBLE = 1024.0;
const double MEGABYTE_DOUBLE = KILOBYTE_DOUBLE * KILOBYTE_DOUBLE;
const double GIGABYTE_DOUBLE = KILOBYTE_DOUBLE * MEGABYTE_DOUBLE;

inline double
inMegabytes(usint bytes)
{
  return bytes / MEGABYTE_DOUBLE;
}

inline double
inBPC(usint bytes, usint size)
{
  return (CHAR_BIT * bytes) / (double)size;
}

inline double
inMicroseconds(double seconds)
{
  return seconds * MILLION_DOUBLE;
}

inline double
inNanoseconds(double seconds)
{
  return seconds * BILLION_DOUBLE;
}

//--------------------------------------------------------------------------

const usint CHARS = ((usint)1 << CHAR_BIT);
const usint WORD_BITS = CHAR_BIT * sizeof(usint);
const usint WORD_MAX = ~((usint)0);

// Previous GET was broken when BITS == WORD_BITS
// Current version works for usints and less
//#define GET(FIELD, BITS) ((FIELD) & ((1 << (BITS)) - 1))
#define GET(FIELD, BITS) ((FIELD) & (WORD_MAX >> (WORD_BITS - (BITS))))
#define LOWER(FIELD, N)  ((FIELD) >> (N))
#define HIGHER(FIELD, N) ((FIELD) << (N))

#define BITS_TO_BYTES(BITS) (((BITS) + CHAR_BIT - 1) / CHAR_BIT)
#define BYTES_TO_WORDS(BYTES) (((BYTES) + sizeof(usint) - 1) / sizeof(usint))
#define BITS_TO_WORDS(BITS) (((BITS) + WORD_BITS - 1) / WORD_BITS)

//--------------------------------------------------------------------------

} // namespace CSA


#endif // _RLCSA_DEFINITIONS_H
