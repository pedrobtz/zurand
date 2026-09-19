/* SIMD across chunks, not within one: 4 consecutive chunks, one per AVX2
 * lane, each running the SHIPPED sequential recurrence. Output must equal
 * seq() word for word -- if it does, SIMD needs no stream change at all.
 * Four steps are buffered and transposed so each chunk gets contiguous
 * vector stores instead of a 4 KiB-strided scatter. */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include "Random123/philox.h"
#define WORDS 512
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);
  return t.tv_sec+t.tv_nsec*1e-9;}
static philox4x64_key_t PK={{0x243f6a8885a308d3ULL,0x13198a2e03707344ULL}};
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline double u01(uint64_t b){uint64_t s=(b>>12)|UINT64_C(0x3ff0000000000000);
  double d; memcpy(&d,&s,sizeof d); return d-(1.0-0x1.0p-53);}
static inline void seed(uint64_t c,uint64_t*s){
  philox4x64_ctr_t k={{c,0,1,0}}; philox4x64_ctr_t r=philox4x64(k,PK);
  s[0]=r.v[0];s[1]=r.v[1];s[2]=r.v[2];s[3]=r.v[3]; if(!(s[0]|s[1]|s[2]|s[3]))s[0]=1;}
/* the shipped engine: one chunk, sequential */
static void seq(uint64_t*buf,uint64_t c){ uint64_t s[4]; seed(c,s);
  for(int j=0;j<WORDS;j++){ uint64_t r=rotl(s[0]+s[3],23)+s[0]; uint64_t t=s[1]<<17;
    s[2]^=s[0];s[3]^=s[1];s[1]^=s[2];s[0]^=s[3];s[2]^=t;s[3]=rotl(s[3],45); buf[j]=r; } }
static inline __m256i vrotl(__m256i x,int k){return _mm256_or_si256(_mm256_slli_epi64(x,k),_mm256_srli_epi64(x,64-k));}
/* 4 chunks c..c+3 -> buf[0..2047], chunk l at buf+l*WORDS. Same stream as seq. */
static void seq4c(uint64_t*buf,uint64_t c){ uint64_t st[4][4];
  for(int l=0;l<4;l++) seed(c+l,st[l]);
  __m256i s0=_mm256_set_epi64x(st[3][0],st[2][0],st[1][0],st[0][0]);
  __m256i s1=_mm256_set_epi64x(st[3][1],st[2][1],st[1][1],st[0][1]);
  __m256i s2=_mm256_set_epi64x(st[3][2],st[2][2],st[1][2],st[0][2]);
  __m256i s3=_mm256_set_epi64x(st[3][3],st[2][3],st[1][3],st[0][3]);
  for(int j=0;j<WORDS;j+=4){ __m256i v[4];
    for(int k=0;k<4;k++){ v[k]=_mm256_add_epi64(vrotl(_mm256_add_epi64(s0,s3),23),s0);
      __m256i t=_mm256_slli_epi64(s1,17);
      s2=_mm256_xor_si256(s2,s0); s3=_mm256_xor_si256(s3,s1);
      s1=_mm256_xor_si256(s1,s2); s0=_mm256_xor_si256(s0,s3);
      s2=_mm256_xor_si256(s2,t);  s3=vrotl(s3,45); }
    /* 4x4 transpose of 64-bit lanes: rows = steps, cols = chunks */
    __m256i t0=_mm256_unpacklo_epi64(v[0],v[1]), t1=_mm256_unpackhi_epi64(v[0],v[1]);
    __m256i t2=_mm256_unpacklo_epi64(v[2],v[3]), t3=_mm256_unpackhi_epi64(v[2],v[3]);
    _mm256_storeu_si256((__m256i*)(buf+0*WORDS+j),_mm256_permute2x128_si256(t0,t2,0x20));
    _mm256_storeu_si256((__m256i*)(buf+1*WORDS+j),_mm256_permute2x128_si256(t1,t3,0x20));
    _mm256_storeu_si256((__m256i*)(buf+2*WORDS+j),_mm256_permute2x128_si256(t0,t2,0x31));
    _mm256_storeu_si256((__m256i*)(buf+3*WORDS+j),_mm256_permute2x128_si256(t1,t3,0x31)); } }
static void fill_seq(double*out,long n){ for(long c=0;c<n/WORDS;c++){uint64_t b[WORDS];seq(b,c);
  double*o=out+c*WORDS; for(int j=0;j<WORDS;j++)o[j]=u01(b[j]);}}
static void fill_seq4c(double*out,long n){ long nc=n/WORDS; long c=0;
  for(;c+4<=nc;c+=4){uint64_t b[4*WORDS];seq4c(b,c);
    double*o=out+c*WORDS; for(int j=0;j<4*WORDS;j++)o[j]=u01(b[j]);}
  for(;c<nc;c++){uint64_t b[WORDS];seq(b,c);double*o=out+c*WORDS;for(int j=0;j<WORDS;j++)o[j]=u01(b[j]);}}
int main(int argc,char**argv){ long n=argc>1?atol(argv[1]):10000000L; n-=n%WORDS; int reps=20;
  { uint64_t a[WORDS],b[4*WORDS]; int ok=1; seq4c(b,100);
    for(int l=0;l<4;l++){ seq(a,100+l); if(memcmp(a,b+l*WORDS,sizeof a)) ok=0; }
    printf("seq4c == seq for chunks 100..103, word for word: %s\n", ok?"yes":"NO"); }
  double*out=aligned_alloc(64,(size_t)n*8); double acc=0,e1,e2;
  fill_seq(out,n); double t=now(); for(int r=0;r<reps;r++)fill_seq(out,n); e1=(now()-t)/reps; acc+=out[n-1];
  fill_seq4c(out,n); t=now(); for(int r=0;r<reps;r++)fill_seq4c(out,n); e2=(now()-t)/reps; acc+=out[n-1];
  printf("  seq   (shipped, scalar)      %6.1f M/s\n", n/e1/1e6);
  printf("  seq4c (4 chunks per AVX2 reg) %6.1f M/s   %.2fx  -- SAME stream\n", n/e2/1e6, e1/e2);
  free(out); return acc<-1e9; }
