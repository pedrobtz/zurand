/* Lane layout for xoshiro256pp: what does interleaving cost in scalar, and
 * what does SIMD actually pay? zurand's two-pass uniform fill, reused buffer.
 *   seq   : 1 lane, sequential (the engine as shipped)
 *   il4   : 4 lanes interleaved, scalar   -- the proposed stream definition
 *   il8   : 8 lanes interleaved, scalar   -- if AVX-512 were the target
 *   avx4  : 4 lanes in one AVX2 register  -- identical output to il4          */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "Random123/philox.h"
#ifdef __AVX2__
#include <immintrin.h>
#endif
#define WORDS 512
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);
  return t.tv_sec+t.tv_nsec*1e-9;}
static philox4x64_key_t PK={{0x243f6a8885a308d3ULL,0x13198a2e03707344ULL}};
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline double u01(uint64_t b){uint64_t s=(b>>12)|UINT64_C(0x3ff0000000000000);
  double d; memcpy(&d,&s,sizeof d); return d-(1.0-0x1.0p-53);}
static inline void seed(uint64_t c,uint64_t lane,uint64_t*s){
  philox4x64_ctr_t k={{c,lane,1,0}}; philox4x64_ctr_t r=philox4x64(k,PK);
  s[0]=r.v[0];s[1]=r.v[1];s[2]=r.v[2];s[3]=r.v[3]; if(!(s[0]|s[1]|s[2]|s[3]))s[0]=1;}
#define STEP(s0,s1,s2,s3,out) do{ uint64_t r_=rotl(s0+s3,23)+s0; uint64_t t_=s1<<17; \
  s2^=s0;s3^=s1;s1^=s2;s0^=s3;s2^=t_;s3=rotl(s3,45); out=r_; }while(0)

static void seq(uint64_t*buf,uint64_t c){ uint64_t s[4]; seed(c,0,s);
  for(int j=0;j<WORDS;j++) STEP(s[0],s[1],s[2],s[3],buf[j]); }
#define ILN(NAME,L) \
static void NAME(uint64_t*buf,uint64_t c){ uint64_t s0[L],s1[L],s2[L],s3[L]; \
  for(int l=0;l<L;l++){uint64_t s[4];seed(c,l,s);s0[l]=s[0];s1[l]=s[1];s2[l]=s[2];s3[l]=s[3];} \
  for(int j=0;j<WORDS/L;j++) for(int l=0;l<L;l++) STEP(s0[l],s1[l],s2[l],s3[l],buf[j*L+l]); }
ILN(il2,2) ILN(il4,4) ILN(il8,8)
static void seq8(uint64_t*buf,uint64_t c){ uint64_t s[4],d[4]; seed(c,0,s);
  for(int l=1;l<8;l++){seed(c,l,d); s[0]^=d[0]&0;}   /* 7 more seeds, discarded */
  for(int j=0;j<WORDS;j++) STEP(s[0],s[1],s[2],s[3],buf[j]); }
#ifdef __AVX2__
static inline __m256i vrotl(__m256i x,int k){return _mm256_or_si256(_mm256_slli_epi64(x,k),_mm256_srli_epi64(x,64-k));}
static void avx4(uint64_t*buf,uint64_t c){ uint64_t st[4][4];
  for(int l=0;l<4;l++) seed(c,l,st[l]);
  __m256i s0=_mm256_set_epi64x(st[3][0],st[2][0],st[1][0],st[0][0]);
  __m256i s1=_mm256_set_epi64x(st[3][1],st[2][1],st[1][1],st[0][1]);
  __m256i s2=_mm256_set_epi64x(st[3][2],st[2][2],st[1][2],st[0][2]);
  __m256i s3=_mm256_set_epi64x(st[3][3],st[2][3],st[1][3],st[0][3]);
  for(int j=0;j<WORDS/4;j++){
    __m256i r=_mm256_add_epi64(vrotl(_mm256_add_epi64(s0,s3),23),s0);
    __m256i t=_mm256_slli_epi64(s1,17);
    s2=_mm256_xor_si256(s2,s0); s3=_mm256_xor_si256(s3,s1);
    s1=_mm256_xor_si256(s1,s2); s0=_mm256_xor_si256(s0,s3);
    s2=_mm256_xor_si256(s2,t);  s3=vrotl(s3,45);
    _mm256_storeu_si256((__m256i*)(buf+4*j),r); } }
#endif
static void fill(double*out,long n,void(*src)(uint64_t*,uint64_t)){
  for(long c=0;c<n/WORDS;c++){uint64_t buf[WORDS];src(buf,c);
    double*o=out+c*WORDS; for(int j=0;j<WORDS;j++)o[j]=u01(buf[j]);}}
int main(int argc,char**argv){ long n=argc>1?atol(argv[1]):10000000L; n-=n%WORDS;
  int reps=n>1000000?20:400; double*out=aligned_alloc(64,(size_t)n*8); double acc=0;
  struct{const char*nm;void(*f)(uint64_t*,uint64_t);}t[]={{"seq (shipped)",seq},{"seq8 (8 seeds)",seq8},{"il2 scalar   ",il2},{"il4 scalar   ",il4},{"il8 scalar   ",il8}
#ifdef __AVX2__
  ,{"avx4 (4 lanes)",avx4}
#endif
  }; int nt=sizeof t/sizeof*t; double base=0;
  /* output identity: il4 and avx4 must agree word for word */
  { uint64_t a[WORDS],b[WORDS]; il4(a,7);
#ifdef __AVX2__
    avx4(b,7); printf("il4 == avx4 word-for-word: %s\n", memcmp(a,b,sizeof a)?"NO":"yes");
#endif
  }
  for(int k=0;k<nt;k++){ fill(out,n,t[k].f);
    double b=now(); for(int r=0;r<reps;r++) fill(out,n,t[k].f); double e=(now()-b)/reps;
    if(!base)base=e; acc+=out[n-1];
    printf("  %s %6.1f M/s  %5.2fx vs seq  %.2f ns/value\n",t[k].nm,n/e/1e6,base/e,e/n*1e9);}
  free(out); return acc<-1e9; }
