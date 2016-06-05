#ifndef _sarray_h_
#define _sarray_h_

typedef unsigned char uchar;

int lcpa(const int *a, const int *s, int *b, int n, int *inv);

/** @returns the character in @param seq before the position @param pos */
int pred_char(int const* seq, const int pos);



void suffixArray(int* s, int* SA, int n, int K) ;
void suffixsort(int *x, int *p, int n, int k, int l);

#endif // _sarray_h_
