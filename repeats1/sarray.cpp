/*
Hybrid suffix-array builder, written by Sean Quinlan and Sean Doward,
distributed under the Plan 9 license, which reads in part

3.3 With respect to Your distribution of Licensed Software (or any
portion thereof), You must include the following information in a
conspicuous location governing such distribution (e.g., a separate
file) and on all copies of any Source Code version of Licensed
Software You distribute:

    "The contents herein includes software initially developed by
    Lucent Technologies Inc. and others, and is subject to the terms
    of the Lucent Technologies Inc. Plan 9 Open Source License
    Agreement.  A copy of the Plan 9 Open Source License Agreement is
    available at: http://plan9.bell-labs.com/plan9dist/download.html
    or by contacting Lucent Technologies at http: //www.lucent.com.
    All software distributed under such Agreement is distributed on,
    obligations and limitations under such Agreement.  Portions of
    the software developed by Lucent Technologies Inc. and others are
    Copyright (c) 2002.  All rights reserved.
    Contributor(s):___________________________"
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sarray.h"


int pred_char(int const* seq, const int position){
	if(position>0) return seq[position-1];
	return -1;
}

/*
   int *lcp(const int *a, const char *s, int n)
   Precondition: a is suffix array for string s of
      length n *including* the terminating '\0'.
   Return value: longest-common-prefix array; 0 on error.
   Reference: T. Kasai, G. Lee, H. Arimura, S.Arikawa
   and K. Park, "Linear-time longest-common-prefix
   computation in suffix arrays and its applications",
   Proc 12th Annual Conference on Combinatorial Pattern
   Matching, Springer, LNCS 2089 (2001) 181-192.

   lcp[x] is the length of the longest common prefix of
   suffixes s[a[x-1]..] and s[a[x]..].

   The algorithm determines the elements of lcp in the
   order that the suffixes occur in s.  It uses this fact:
   If the lcp for suffix s[i..] has length h, where h>0,
   then the lcp for suffix s[i+1..] is at least h-1.

   Proof.
   Let the immediate lexicographic predecessor of suffix
      s[i..] be s[j..], i.e. lex[i]=lex[j]+1.
   If s[i..] and s[j..] have a common prefix of length h,
      where h>0, then s[i+1..] and s[j+1..] have a common
      prefix of length h-1.
   Since s[i+1..] and s[j+1..] differ from s[i..] and s[j..]
      respectively only by the deletion of a common first
      letter, the two pairs must be similarly ordered.
      Hence s[j+1..] lexicographically precedes s[i+1..].
   Since s[i+1..] shares a common prefix of length h-1 with
      some lexicographic predecessor, namely s[j+1..], it
      must share a common prefix of length at least h-1 with
      its immediate predecessor.  Otherwise the suffix array
      would be out of order.

   Running time is O(n).

   Proof.
   h is bounded by n; and h is decreased by 1 at most
   n times. Hence h is increased at most 2n times.
   This bounds the number of executions of the inner loop.
*/

/*
   inv is the inverse of a: if inv[i]=x then a[x]=i.
   In other words, inv[i] is the index x of the
   pointer (in array a) to suffix s[i..].
*/

/* int* */
/* lcp(const int *a, const char *s, int n)  */
/* { */
/* 	int *lcp = (int*)malloc(n*sizeof(int)); */

/* 	if(lcp == 0) */
/* 		return 0; */
/* 	if(lcpa(a, s, lcp, n) == 0) { */
/* 		free(lcp); */
/* 		return 0; */
/* 	} */
/* 	return lcp; */
/* } */

/* lcpa is used by the java native method */

// mgalle modification: changed int -> char in declaration of s0. The same?

int
lcpa(const int *a, const int *s0, int *lcp, int n, int *inv)
{
	int i, h;
	int *s = (int*) s0;

	h = 0;			/* visit in string order */
	for(i=0; i<n-1; i++) {  /* omit last, least suff */
	  //	  if(i%1000==0) printf(".");
		int x = inv[i]; /* i,j,x,h as in intro */
		int j = a[x-1];

		int *p1 = s + i + h;
		int *p0 = s + j + h;
		while(*p1++ == *p0++)
			h++;


		lcp[x] = h;
		if(h > 0)
		  h--;
	}

	lcp[0] = 0;	/* least suffix has no predecessor */
	//free(inv);
	return 1;
}




////////////////////////////////////
/* AUTHORS:
 * Juha Karkaainen (University of Helsinki, Helsinki, Finland)
 * Peter Sanders (University of Karlsruhe, Karlsruhe, Germany)
 * Stefan Burkhardt (Google, Inc., Zurich, Switzerland)

 downloaded from http://www.mpi-inf.mpg.de/~sanders/programs/suffix/

*/

inline bool leq(int a1, int a2,   int b1, int b2) { // lexic. order for pairs
  return(a1 < b1 || a1 == b1 && a2 <= b2);
}                                                   // and triples
inline bool leq(int a1, int a2, int a3,   int b1, int b2, int b3) {
  return(a1 < b1 || a1 == b1 && leq(a2,a3, b2,b3));
}
// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
static void radixPass(int* a, int* b, int* r, int n, int K)
{ // count occurrences
  int* c = new int[K + 1];                          // counter array
  for (int i = 0;  i <= K;  i++) c[i] = 0;         // reset counters
  for (int i = 0;  i < n;  i++) c[r[a[i]]]++;    // count occurences
  for (int i = 0, sum = 0;  i <= K;  i++) { // exclusive prefix sums
     int t = c[i];  c[i] = sum;  sum += t;
  }
  for (int i = 0;  i < n;  i++) b[c[r[a[i]]]++] = a[i];      // sort
  delete [] c;
}

// find the suffix array SA of s[0..n-1] in {1..K}^n
// require s[n]=s[n+1]=s[n+2]=0, n>=2
void suffixArray(int* s, int* SA, int n, int K) {
  int n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2;
  int* s12  = new int[n02 + 3];  s12[n02]= s12[n02+1]= s12[n02+2]=0;
  int* SA12 = new int[n02 + 3]; SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;
  int* s0   = new int[n0];
  int* SA0  = new int[n0];

  // generate positions of mod 1 and mod  2 suffixes
  // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
  for (int i=0, j=0;  i < n+(n0-n1);  i++) if (i%3 != 0) s12[j++] = i;

  // lsb radix sort the mod 1 and mod 2 triples
  radixPass(s12 , SA12, s+2, n02, K);
  radixPass(SA12, s12 , s+1, n02, K);
  radixPass(s12 , SA12, s  , n02, K);

  // find lexicographic names of triples
  int name = 0, c0 = -1, c1 = -1, c2 = -1;
  for (int i = 0;  i < n02;  i++) {
    if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) {
      name++;  c0 = s[SA12[i]];  c1 = s[SA12[i]+1];  c2 = s[SA12[i]+2];
    }
    if (SA12[i] % 3 == 1) { s12[SA12[i]/3]      = name; } // left half
    else                  { s12[SA12[i]/3 + n0] = name; } // right half
  }

  // recurse if names are not yet unique
  if (name < n02) {
    suffixArray(s12, SA12, n02, name);
    // store unique names in s12 using the suffix array
    for (int i = 0;  i < n02;  i++) s12[SA12[i]] = i + 1;
  } else // generate the suffix array of s12 directly
    for (int i = 0;  i < n02;  i++) SA12[s12[i] - 1] = i;

  // stably sort the mod 0 suffixes from SA12 by their first character
  for (int i=0, j=0;  i < n02;  i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
  radixPass(s0, SA0, s, n0, K);

  // merge sorted SA0 suffixes and sorted SA12 suffixes
  for (int p=0,  t=n0-n1,  k=0;  k < n;  k++) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
    int i = GetI(); // pos of current offset 12 suffix
    int j = SA0[p]; // pos of current offset 0  suffix
    if (SA12[t] < n0 ?
        leq(s[i],       s12[SA12[t] + n0], s[j],       s12[j/3]) :
        leq(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]))
    { // suffix from SA12 is smaller
      SA[k] = i;  t++;
      if (t == n02) { // done --- only SA0 suffixes left
        for (k++;  p < n0;  p++, k++) SA[k] = SA0[p];
      }
    } else {
      SA[k] = j;  p++;
      if (p == n0)  { // done --- only SA12 suffixes left
        for (k++;  t < n02;  t++, k++) SA[k] = GetI();
      }
    }
  }
  delete [] s12; delete [] SA12; delete [] SA0; delete [] s0;
}











/* qsufsort.c
   Copyright 1999, N. Jesper Larsson, all rights reserved.

   This file contains an implementation of the algorithm presented in "Faster
   Suffix Sorting" by N. Jesper Larsson (jesper@cs.lth.se) and Kunihiko
   Sadakane (sada@is.s.u-tokyo.ac.jp).

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.*/

#include <limits.h>

static int *I,                  /* group array, ultimately suffix array.*/
   *V,                          /* inverse array, ultimately inverse of I.*/
   r,                           /* number of symbols aggregated by transform.*/
   h;                           /* length of already-sorted prefixes.*/

#define KEY(p)          (V[*(p)+(h)])
#define SWAP(p, q)      (tmp=*(p), *(p)=*(q), *(q)=tmp)
#define MED3(a, b, c)   (KEY(a)<KEY(b) ?                        \
        (KEY(b)<KEY(c) ? (b) : KEY(a)<KEY(c) ? (c) : (a))       \
        : (KEY(b)>KEY(c) ? (b) : KEY(a)>KEY(c) ? (c) : (a)))

/* Subroutine for select_sort_split and sort_split. Sets group numbers for a
   group whose lowest position in I is pl and highest position is pm.*/

static void update_group(int *pl, int *pm)
{
   int g;

   g=pm-I;                      /* group number.*/
   V[*pl]=g;                    /* update group number of first position.*/
   if (pl==pm)
      *pl=-1;                   /* one element, sorted group.*/
   else
      do                        /* more than one element, unsorted group.*/
         V[*++pl]=g;            /* update group numbers.*/
      while (pl<pm);
}

/* Quadratic sorting method to use for small subarrays. To be able to update
   group numbers consistently, a variant of selection sorting is used.*/

static void select_sort_split(int *p, int n) {
   int *pa, *pb, *pi, *pn;
   int f, v, tmp;

   pa=p;                        /* pa is start of group being picked out.*/
   pn=p+n-1;                    /* pn is last position of subarray.*/
   while (pa<pn) {
      for (pi=pb=pa+1, f=KEY(pa); pi<=pn; ++pi)
         if ((v=KEY(pi))<f) {
            f=v;                /* f is smallest key found.*/
            SWAP(pi, pa);       /* place smallest element at beginning.*/
            pb=pa+1;            /* pb is position for elements equal to f.*/
         } else if (v==f) {     /* if equal to smallest key.*/
            SWAP(pi, pb);       /* place next to other smallest elements.*/
            ++pb;
         }
      update_group(pa, pb-1);   /* update group values for new group.*/
      pa=pb;                    /* continue sorting rest of the subarray.*/
   }
   if (pa==pn) {                /* check if last part is single element.*/
      V[*pa]=pa-I;
      *pa=-1;                   /* sorted group.*/
   }
}

/* Subroutine for sort_split, algorithm by Bentley & McIlroy.*/

static int choose_pivot(int *p, int n) {
   int *pl, *pm, *pn;
   int s;

   pm=p+(n>>1);                 /* small arrays, middle element.*/
   if (n>7) {
      pl=p;
      pn=p+n-1;
      if (n>40) {               /* big arrays, pseudomedian of 9.*/
         s=n>>3;
         pl=MED3(pl, pl+s, pl+s+s);
         pm=MED3(pm-s, pm, pm+s);
         pn=MED3(pn-s-s, pn-s, pn);
      }
      pm=MED3(pl, pm, pn);      /* midsize arrays, median of 3.*/
   }
   return KEY(pm);
}

/* Sorting routine called for each unsorted group. Sorts the array of integers
   (suffix numbers) of length n starting at p. The algorithm is a ternary-split
   quicksort taken from Bentley & McIlroy, "Engineering a Sort Function",
   Software -- Practice and Experience 23(11), 1249-1265 (November 1993). This
   function is based on Program 7.*/

static void sort_split(int *p, int n)
{
   int *pa, *pb, *pc, *pd, *pl, *pm, *pn;
   int f, v, s, t, tmp;

   if (n<7) {                   /* multi-selection sort smallest arrays.*/
      select_sort_split(p, n);
      return;
   }

   v=choose_pivot(p, n);
   pa=pb=p;
   pc=pd=p+n-1;
   while (1) {                  /* split-end partition.*/
      while (pb<=pc && (f=KEY(pb))<=v) {
         if (f==v) {
            SWAP(pa, pb);
            ++pa;
         }
         ++pb;
      }
      while (pc>=pb && (f=KEY(pc))>=v) {
         if (f==v) {
            SWAP(pc, pd);
            --pd;
         }
         --pc;
      }
      if (pb>pc)
         break;
      SWAP(pb, pc);
      ++pb;
      --pc;
   }
   pn=p+n;
   if ((s=pa-p)>(t=pb-pa))
      s=t;
   for (pl=p, pm=pb-s; s; --s, ++pl, ++pm)
      SWAP(pl, pm);
   if ((s=pd-pc)>(t=pn-pd-1))
      s=t;
   for (pl=pb, pm=pn-s; s; --s, ++pl, ++pm)
      SWAP(pl, pm);

   s=pb-pa;
   t=pd-pc;
   if (s>0)
      sort_split(p, s);
   update_group(p+s, p+n-t-1);
   if (t>0)
      sort_split(p+n-t, t);
}

/* Bucketsort for first iteration.

   Input: x[0...n-1] holds integers in the range 1...k-1, all of which appear
   at least once. x[n] is 0. (This is the corresponding output of transform.) k
   must be at most n+1. p is array of size n+1 whose contents are disregarded.

   Output: x is V and p is I after the initial sorting stage of the refined
   suffix sorting algorithm.*/

static void bucketsort(int *x, int *p, int n, int k)
{
   int *pi, i, c, d, g;

   for (pi=p; pi<p+k; ++pi)
      *pi=-1;                   /* mark linked lists empty.*/
   for (i=0; i<=n; ++i) {
      x[i]=p[c=x[i]];           /* insert in linked list.*/
      p[c]=i;
   }
   for (pi=p+k-1, i=n; pi>=p; --pi) {
      d=x[c=*pi];               /* c is position, d is next in list.*/
      x[c]=g=i;                 /* last position equals group number.*/
      if (d>=0) {               /* if more than one element in group.*/
         p[i--]=c;              /* p is permutation for the sorted x.*/
         do {
            d=x[c=d];           /* next in linked list.*/
            x[c]=g;             /* group number in x.*/
            p[i--]=c;           /* permutation in p.*/
         } while (d>=0);
      } else
         p[i--]=-1;             /* one element, sorted group.*/
   }
}

/* Transforms the alphabet of x by attempting to aggregate several symbols into
   one, while preserving the suffix order of x. The alphabet may also be
   compacted, so that x on output comprises all integers of the new alphabet
   with no skipped numbers.

   Input: x is an array of size n+1 whose first n elements are positive
   integers in the range l...k-1. p is array of size n+1, used for temporary
   storage. q controls aggregation and compaction by defining the maximum value
   for any symbol during transformation: q must be at least k-l; if q<=n,
   compaction is guaranteed; if k-l>n, compaction is never done; if q is
   INT_MAX, the maximum number of symbols are aggregated into one.

   Output: Returns an integer j in the range 1...q representing the size of the
   new alphabet. If j<=n+1, the alphabet is compacted. The global variable r is
   set to the number of old symbols grouped into one. Only x[n] is 0.*/

static int transform(int *x, int *p, int n, int k, int l, int q)
{
   int b, c, d, e, i, j, m, s;
   int *pi, *pj;

   for (s=0, i=k-l; i; i>>=1)
      ++s;                      /* s is number of bits in old symbol.*/
   e=INT_MAX>>s;                /* e is for overflow checking.*/
   for (b=d=r=0; r<n && d<=e && (c=d<<s|(k-l))<=q; ++r) {
      b=b<<s|(x[r]-l+1);        /* b is start of x in chunk alphabet.*/
      d=c;                      /* d is max symbol in chunk alphabet.*/
   }
   m=(1<<(r-1)*s)-1;            /* m masks off top old symbol from chunk.*/
   x[n]=l-1;                    /* emulate zero terminator.*/
   if (d<=n) {                  /* if bucketing possible, compact alphabet.*/
      for (pi=p; pi<=p+d; ++pi)
         *pi=0;                 /* zero transformation table.*/
      for (pi=x+r, c=b; pi<=x+n; ++pi) {
         p[c]=1;                /* mark used chunk symbol.*/
         c=(c&m)<<s|(*pi-l+1);  /* shift in next old symbol in chunk.*/
      }
      for (i=1; i<r; ++i) {     /* handle last r-1 positions.*/
         p[c]=1;                /* mark used chunk symbol.*/
         c=(c&m)<<s;            /* shift in next old symbol in chunk.*/
      }
      for (pi=p, j=1; pi<=p+d; ++pi)
         if (*pi)
            *pi=j++;            /* j is new alphabet size.*/
      for (pi=x, pj=x+r, c=b; pj<=x+n; ++pi, ++pj) {
         *pi=p[c];              /* transform to new alphabet.*/
         c=(c&m)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
      }
      while (pi<x+n) {          /* handle last r-1 positions.*/
         *pi++=p[c];            /* transform to new alphabet.*/
         c=(c&m)<<s;            /* shift right-end zero in chunk.*/
      }
   } else {                     /* bucketing not possible, don't compact.*/
      for (pi=x, pj=x+r, c=b; pj<=x+n; ++pi, ++pj) {
         *pi=c;                 /* transform to new alphabet.*/
         c=(c&m)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
      }
      while (pi<x+n) {          /* handle last r-1 positions.*/
         *pi++=c;               /* transform to new alphabet.*/
         c=(c&m)<<s;            /* shift right-end zero in chunk.*/
      }
      j=d+1;                    /* new alphabet size.*/
   }
   x[n]=0;                      /* end-of-string symbol is zero.*/
   return j;                    /* return new alphabet size.*/
}

/* Makes suffix array p of x. x becomes inverse of p. p and x are both of size
   n+1. Contents of x[0...n-1] are integers in the range l...k-1. Original
   contents of x[n] is disregarded, the n-th symbol being regarded as
   end-of-string smaller than all other symbols.*/

void suffixsort(int *x, int *p, int n, int k, int l)
{
   int *pi, *pk;
   int i, j, s, sl;

   V=x;                         /* set global values.*/
   I=p;

   if (n>=k-l) {                /* if bucketing possible,*/
      j=transform(V, I, n, k, l, n);
      bucketsort(V, I, n, j);   /* bucketsort on first r positions.*/
   } else {
      transform(V, I, n, k, l, INT_MAX);
      for (i=0; i<=n; ++i)
         I[i]=i;                /* initialize I with suffix numbers.*/
      h=0;
      sort_split(I, n+1);       /* quicksort on first r positions.*/
   }
   h=r;                         /* number of symbols aggregated by transform.*/

   while (*I>=-n) {
      pi=I;                     /* pi is first position of group.*/
      sl=0;                     /* sl is negated length of sorted groups.*/
      do {
         if ((s=*pi)<0) {
            pi-=s;              /* skip over sorted group.*/
            sl+=s;              /* add negated length to sl.*/
         } else {
            if (sl) {
               *(pi+sl)=sl;     /* combine sorted groups before pi.*/
               sl=0;
            }
            pk=I+V[s]+1;        /* pk-1 is last position of unsorted group.*/
            sort_split(pi, pk-pi);
            pi=pk;              /* next group.*/
         }
      } while (pi<=I+n);
      if (sl)                   /* if the array ends with a sorted group.*/
         *(pi+sl)=sl;           /* combine sorted groups at end of I.*/
      h=2*h;                    /* double sorted-depth.*/
   }

   for (i=0; i<=n; ++i)         /* reconstruct suffix array from inverse.*/
      I[V[i]]=i;
}
