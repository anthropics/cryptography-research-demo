// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0

/* ===================================================================
 * Key-recovery attack on 13-round LEA-128 (chosen-plaintext encryption oracle)
 *
 * Technique: conditional differential-linear cryptanalysis.
 *  - A 1-round "conditional differential" (6 key-dependent bit conditions on
 *    the chosen plaintexts) forces Delta x^(1) = (0x80000000,)*4 with prob 1.
 *  - That MSB difference then propagates with probability 1 for two more
 *    rounds:  -> (0,0,0,2^31) -> (0,0,2^28,0)  (MSB differences are free
 *    through modular addition).
 *  - The resulting 13-round differential-linear distinguisher gives huge
 *    biases on Delta x1^(11) (up to 4e-2).
 *  - bit j of Delta x1^(11) depends only on bits 0..j-1 of four 32-bit
 *    unknowns derived from the last round keys, enabling a bit-by-bit
 *    beam search key recovery.
 *
 * The beam search scores a candidate by recomputing, per ciphertext pair,
 * the word chain
 *      a1 = ror9(c0) - (c3 ^ A0)          (Delta x1^(12))
 *      a2 = rol5(c1) - (a1 ^ T12)         (Delta x2^(12))
 *      a3 = rol3(c2) - (a2 ^ T13)         (Delta x0^(11))
 *      b1 = ror9(c3) - (a3 ^ U)           (Delta x1^(11))
 * and correlating bit b+1 of the four pair-deltas (with candidate bits
 * 0..b assimilated, bit b+1 of each delta is exact: unknown higher key
 * bits cancel in the XOR across the pair).  Scoring is stateless - no
 * per-candidate data survives between bit levels - so memory is
 * independent of the beam width:  32 bytes per pair, period.  Wide
 * rescue searches that previously needed >100 GiB now fit in a few.
 *
 * Because some keys have runs of bit positions where the distinguisher
 * signal vanishes (key-intrinsic "dead zones"), the beam can adapt its
 * width (experimental, ADAPTIVE=1; OFF by default - every number in the
 * README was measured with it off): a level with no statistical
 * separation keeps ALL children (capped) instead of cutting on noise,
 * and a weakly-separated level keeps 4x the normal width.  See README.
 * =================================================================== */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include "lea.h"
#include "oracle.h"

typedef struct { uint64_t s[4]; } rng_t;
static inline uint64_t rotl64(const uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nextr(rng_t*r){ const uint64_t res=rotl64(r->s[1]*5,7)*9; const uint64_t t=r->s[1]<<17;
 r->s[2]^=r->s[0];r->s[3]^=r->s[1];r->s[1]^=r->s[2];r->s[0]^=r->s[3];r->s[2]^=t;r->s[3]=rotl64(r->s[3],45); return res;}
static void rseed(rng_t*r,uint64_t x){for(int i=0;i<4;i++){x+=0x9E3779B97f4A7C15ULL;uint64_t z=x;z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;r->s[i]=z^(z>>31);} }

const uint32_t DP[4] = {0x80000000u, 0x80400000u, 0x80400010u, 0x80400014u};

static void build_pt(uint32_t p[4], uint32_t r0,uint32_t r1,uint32_t r2,uint32_t r3, int combo){
    uint32_t p0b22=(combo>>0)&1, p1b4=(combo>>1)&1, p2b2=(combo>>2)&1;
    uint32_t X=(combo>>3)&1, Y=(combo>>4)&1, Z=(combo>>5)&1;
    p[0] = (p0b22<<22) | (r0 & 0xFF800000u);
    uint32_t p1b22 = (r1>>22)&1;
    p[1] = (p1b4<<4) | (p1b22<<22) | (r1 & 0xFF800000u);
    uint32_t p2b4 = (r2>>4)&1;
    uint32_t p2b22 = p1b22 ^ X;
    p[2] = (p2b2<<2) | (p2b4<<4) | (r2 & 0x003FFFE0u) | (p2b22<<22) | (r2 & 0xFF800000u);
    uint32_t p3b4 = p2b4 ^ Z;
    uint32_t p3b22 = p2b22 ^ Y;
    p[3] = (r3 & 0x0000000Cu) | (p3b4<<4) | (r3 & 0x003FFFE0u) | (p3b22<<22) | (r3 & 0xFF800000u);
}

/* ---------- overflow-checked allocation: every buffer here is sized from
   runtime parameters, so every size computation is checked ---------- */
static void* xalloc(size_t nmemb, size_t size){
    if(nmemb && size && nmemb > (size_t)-1/size){ fprintf(stderr,"memory allocation failed\n"); exit(3); }
    void *p = malloc(nmemb*size);
    if(!p){ fprintf(stderr,"memory allocation failed\n"); exit(3); }
    return p;
}

/* ---------- ciphertext pair data: 8 structure-of-arrays columns ----------
   Per pair member we store the 4 rotated ciphertext words the scoring chain
   consumes: ror9(c0), rol5(c1), rol3(c2), c3.  32 bytes per pair total -
   the same bytes the attack has always needed, just word-addressed. */
static uint32_t *COL[8];     /* [m*4+w], member m=0,1 */
static long long N;
static int NTHREADS=8;

static void alloc_cols(long long n){
    /* N multiple of 512 keeps every aligned_alloc size a multiple of 64 */
    N=(n/512)*512;
    for(int i=0;i<8;i++){
        COL[i]=aligned_alloc(64,(size_t)N*4);
        if(!COL[i]){fprintf(stderr,"memory allocation failed\n");exit(3);}
    }
}
static void free_cols(void){ for(int i=0;i<8;i++){ free(COL[i]); COL[i]=NULL; } }
static void insert_pair(long long i,const uint32_t c0[4],const uint32_t c1[4]){
    COL[0][i]=ror32(c0[0],9); COL[1][i]=rol32(c0[1],5); COL[2][i]=rol32(c0[2],3); COL[3][i]=c0[3];
    COL[4][i]=ror32(c1[0],9); COL[5][i]=rol32(c1[1],5); COL[6][i]=rol32(c1[2],3); COL[7][i]=c1[3];
}

/* ---------- oracle batch: generate N pairs for a combo, fetch ciphertexts ---------- */
static void gather(long long npairs, int combo, rng_t *rng){
    static uint8_t *ptbuf=NULL,*ctbuf=NULL;
    if(!ptbuf){ ptbuf=xalloc(16,MAX_BLKS); ctbuf=xalloc(16,MAX_BLKS); }
    long long done=0;
    int PAIRS_PER_REQ = MAX_BLKS/2;     /* 2 blocks per pair */
    while(done<npairs){
        long long chunk = npairs-done; if(chunk>PAIRS_PER_REQ) chunk=PAIRS_PER_REQ;
        for(long long i=0;i<chunk;i++){
            uint32_t p[4],q[4];
            build_pt(p,(uint32_t)nextr(rng),(uint32_t)nextr(rng),(uint32_t)nextr(rng),(uint32_t)nextr(rng),combo);
            for(int k=0;k<4;k++)q[k]=p[k]^DP[k];
            for(int k=0;k<4;k++){
                uint32_t v=p[k]; uint8_t*d=ptbuf+ (2*i)*16 + 4*k;
                d[0]=v; d[1]=v>>8; d[2]=v>>16; d[3]=v>>24;
                v=q[k]; d=ptbuf+(2*i+1)*16 + 4*k;
                d[0]=v; d[1]=v>>8; d[2]=v>>16; d[3]=v>>24;
            }
        }
        ora_encrypt(ptbuf,ctbuf,(int)chunk*2);
        for(long long i=0;i<chunk;i++){
            uint32_t c0[4],c1[4];
            for(int k=0;k<4;k++){
                const uint8_t*s=ctbuf+(2*i)*16+4*k;
                c0[k]=(uint32_t)s[0]|((uint32_t)s[1]<<8)|((uint32_t)s[2]<<16)|((uint32_t)s[3]<<24);
                s=ctbuf+(2*i+1)*16+4*k;
                c1[k]=(uint32_t)s[0]|((uint32_t)s[1]<<8)|((uint32_t)s[2]<<16)|((uint32_t)s[3]<<24);
            }
            insert_pair(done+i,c0,c1);
        }
        done+=chunk;
    }
}

/* ---------- stateless tiled scoring ----------
   For each candidate, recompute the chain over a tile of pairs held in L2
   and count set delta bits at the measured position.  One read of the pair
   data serves every candidate in the caller's block; counts are integers,
   so any split over pairs or candidates gives identical totals. */
typedef struct { uint32_t A0,T12,T13,U; double score,lsc; } cand_t;

#define SCORE_TILE 2048   /* pairs per tile: 8 cols x 8 KB = 64 KB (L2) */

static void score_tiles(const cand_t*cands,int nc,int bit,long long lo,long long hi,long long (*cnt)[4]){
    for(long long s=lo;s<hi;s+=SCORE_TILE){
        long long e=s+SCORE_TILE; if(e>hi)e=hi; int len=(int)(e-s);
        const uint32_t *p0=COL[0]+s,*p1=COL[1]+s,*p2=COL[2]+s,*p3=COL[3]+s,
                       *p4=COL[4]+s,*p5=COL[5]+s,*p6=COL[6]+s,*p7=COL[7]+s;
        for(int ci=0;ci<nc;ci++){
            const uint32_t A0=cands[ci].A0,T12=cands[ci].T12,T13=cands[ci].T13,U=cands[ci].U;
            uint32_t t0=0,t1=0,t2=0,t3=0;
            for(int j=0;j<len;j++){
                uint32_t a1 =p0[j]-(p3[j]^A0),      a1b=p4[j]-(p7[j]^A0);
                uint32_t a2 =p1[j]-(a1^T12),        a2b=p5[j]-(a1b^T12);
                uint32_t a3 =p2[j]-(a2^T13),        a3b=p6[j]-(a2b^T13);
                uint32_t b1 =ror32(p3[j],9)-(a3^U), b1b=ror32(p7[j],9)-(a3b^U);
                t0+=((a1^a1b)>>bit)&1; t1+=((a2^a2b)>>bit)&1;
                t2+=((a3^a3b)>>bit)&1; t3+=((b1^b1b)>>bit)&1;
            }
            cnt[ci][0]+=t0; cnt[ci][1]+=t1; cnt[ci][2]+=t2; cnt[ci][3]+=t3;
        }
    }
}

typedef struct { int tid; const cand_t*cands; int nc,bit; long long lo,hi; long long (*cnt)[4]; } sc_job_t;
static void* score_worker_cands(void*arg){   /* split candidates across threads */
    sc_job_t*J=(sc_job_t*)arg;
    int c0=(int)((long long)J->nc*J->tid/NTHREADS), c1=(int)((long long)J->nc*(J->tid+1)/NTHREADS);
    score_tiles(J->cands+c0,c1-c0,J->bit,J->lo,J->hi,J->cnt+c0);
    return NULL;
}
static void* score_worker_pairs(void*arg){   /* split pairs across threads */
    sc_job_t*J=(sc_job_t*)arg;
    /* chunk boundaries are tile-aligned RELATIVE to the range start, so the
       chunks seam exactly at tile edges and cover [lo,hi) with no overlap */
    long long d=J->hi-J->lo;
    long long lo=J->lo + (d*J->tid/NTHREADS/SCORE_TILE)*SCORE_TILE;
    long long hi=(J->tid==NTHREADS-1)? J->hi : J->lo + (d*(J->tid+1)/NTHREADS/SCORE_TILE)*SCORE_TILE;
    score_tiles(J->cands,J->nc,J->bit,lo,hi,J->cnt);
    return NULL;
}

/* Score the level contribution of cands[0..nc) at guess-bit level m (delta
   bits measured at m+1) over pairs [lo,hi).  Writes cand.lsc. */
static void score_cands(cand_t*cands,int nc,int m,long long lo,long long hi){
    int bit=m+1;
    pthread_t *th=xalloc(NTHREADS,sizeof *th); sc_job_t *J=xalloc(NTHREADS,sizeof *J);
    long long (*cnt)[4];
    if(nc>=NTHREADS*4){           /* each thread owns a candidate slice */
        cnt=xalloc(nc,sizeof *cnt); memset(cnt,0,(size_t)nc*sizeof *cnt);
        for(int i=0;i<NTHREADS;i++){ J[i]=(sc_job_t){i,cands,nc,bit,lo,hi,cnt};
            if(pthread_create(&th[i],NULL,score_worker_cands,&J[i])){ fprintf(stderr,"thread creation failed\n"); exit(3); } }
        for(int i=0;i<NTHREADS;i++) pthread_join(th[i],NULL);
    }else{                        /* each thread owns a pair range; merge integer counts */
        cnt=xalloc((size_t)nc*NTHREADS,sizeof *cnt); memset(cnt,0,(size_t)nc*NTHREADS*sizeof *cnt);
        for(int i=0;i<NTHREADS;i++){ J[i]=(sc_job_t){i,cands,nc,bit,lo,hi,cnt+(size_t)nc*i};
            if(pthread_create(&th[i],NULL,score_worker_pairs,&J[i])){ fprintf(stderr,"thread creation failed\n"); exit(3); } }
        for(int i=0;i<NTHREADS;i++) pthread_join(th[i],NULL);
        for(int t=1;t<NTHREADS;t++)for(int c=0;c<nc;c++)for(int k=0;k<4;k++) cnt[c][k]+=cnt[(size_t)nc*t+c][k];
    }
    long long np=hi-lo;
    for(int ci=0;ci<nc;ci++){ double s=0;
        for(int k=0;k<4;k++){ double corr=2.0*((double)cnt[ci][k]/np-0.5); s += corr*corr*(double)np; }
        cands[ci].lsc=s;
    }
    free(cnt); free(th); free(J);
}

static int cmp_cum(const void*a,const void*b){   /* by cumulative score desc;
    ties broken on the candidate words so the order is total and the search
    is reproducible (dead levels make exact score ties systematic) */
    const cand_t*p=a,*q=b;
    double x=p->score+p->lsc, y=q->score+q->lsc;
    if(x<y)return 1;
    if(x>y)return -1;
    if(p->A0 !=q->A0 )return p->A0 <q->A0 ?-1:1;
    if(p->T12!=q->T12)return p->T12<q->T12?-1:1;
    if(p->T13!=q->T13)return p->T13<q->T13?-1:1;
    if(p->U  !=q->U  )return p->U  <q->U  ?-1:1;
    return 0;
}
static int cmp_dbl(const void*a,const void*b){ double x=*(const double*)a,y=*(const double*)b; return x<y?-1:(x>y?1:0); }

/* ---------- adaptive-width knobs (env-overridable; see README) ----------
   Thresholds are in ABSOLUTE score units: under the null the level score is
   chi^2(4)-like independent of N (corr^2*N ~ chi^2(1)), so thresholds do not
   scale with the data size. */
static double env_f(const char*name,double dflt){
    const char*e=getenv(name); if(!e) return dflt;
    char*end; double v=strtod(e,&end);
    if(*end || !(v>=0.0) || v>1e12){ fprintf(stderr,"%s out of range [0..1e12]\n",name); exit(2); }
    return v;
}
static long long env_i(const char*name,long long dflt,long long lo,long long hi){
    const char*e=getenv(name); long long v=dflt;
    if(e){ char*end; v=strtoll(e,&end,0);
        if(*end){ fprintf(stderr,"%s is not a number\n",name); exit(2); } }
    if(v<lo||v>hi){ fprintf(stderr,"%s out of range [%lld..%lld]\n",name,lo,hi); exit(2); }
    return v;
}

/* Diagnostic (TRUTH=32hexkey in the environment): track the rank of the true
   candidate prefix through the beam.  Prints only; never affects the search.
   Useful for calibrating the adaptive thresholds and for understanding any
   failed key: it shows exactly where and how the true prefix was lost. */
static int TRUTH_SET=0; static uint32_t TRUTHW[4];
static void truth_init(void){
    const char*e=getenv("TRUTH");
    if(!e||strlen(e)!=32||strspn(e,"0123456789abcdefABCDEF")!=32) return;
    uint32_t k[4];
    for(int w=0;w<4;w++){uint32_t v=0;for(int b=0;b<4;b++){unsigned x;sscanf(e+2*(4*w+b),"%2x",&x);v|=(uint32_t)x<<(8*b);}k[w]=v;}
    lea_ctx c; lea_setkey(&c,k,13);
    TRUTHW[0]=c.rk[12][0]; TRUTHW[1]=c.rk[12][1]^c.rk[12][2];
    TRUTHW[2]=c.rk[12][1]^c.rk[12][3]; TRUTHW[3]=c.rk[12][1]^c.rk[11][0];
    TRUTH_SET=1;
    printf("TRUTH diagnostic on: A0=%08x T12=%08x T13=%08x U=%08x\n",TRUTHW[0],TRUTHW[1],TRUTHW[2],TRUTHW[3]);
}
static void truth_rank(const cand_t*c,long long n,int m,long long kept,const char*tag){
    if(!TRUTH_SET) return;
    uint32_t mask=(m>=31)?0xffffffffu:((2u<<m)-1);
    long long r=-1;
    for(long long i=0;i<n;i++)
        if(((c[i].A0^TRUTHW[0])&mask)==0&&((c[i].T12^TRUTHW[1])&mask)==0&&
           ((c[i].T13^TRUTHW[2])&mask)==0&&((c[i].U^TRUTHW[3])&mask)==0){r=i;break;}
    if(r<0) printf("    TRUTH %s: LOST at bit %d\n",tag,m);
    else    printf("    TRUTH %s: rank %lld of %lld%s\n",tag,r,n,r<kept?" (kept)":" (DROPPED)");
    fflush(stdout);
}

/* ---------- the beam search ----------
   Scores are CUMULATIVE across levels (child = parent + level contribution):
   a level with no signal adds ~equal noise to every candidate, so earlier
   separation persists through it.  The adaptive tests below consume the
   PER-LEVEL contribution (cand.lsc), never the cumulative score - on
   cumulative scores a dead level would be masked by previously accumulated
   separation and the detector would never fire.

   Adaptive width on top of the cumulative beam:
    - DEAD level (level contributions show no separation; a subsampled verdict
      is re-tested at full N before it counts): keep ALL children, stratified-
      capped at DEADCAP; the level contributes no score.  Never taken at the
      final level - an uncut pool must not reach phase 3, whose cost is
      ~2^(4*(32-MAXBIT)) trial encryptions PER candidate.
    - FLAT level (separation present but below threshold): keep FLATX x BEAM.
   Wide fronts are ranked in two stages (subsample rank, full rescore of the
   top 8 x BEAM; interior rescores capped at 2^26 pairs, the final level uses
   everything).  The stage-A subsample shrinks when the front is huge so one
   level costs a bounded number of chain evaluations (SAMPBUDGET), with a
   2^20-pair floor to keep the ranking statistically meaningful.
   These rules are an experimental alternative to the fixed-width beam and
   are OFF by default (ADAPTIVE=1 enables them). The level score is
   chi^2(4)-like under the null INDEPENDENT of N, so the threshold defaults
   are constants, not N-scaled. */
#define DEADTHR_DFLT 48.0    /* full-N dead gap threshold (experimental path) */
#define FLATTHR_DFLT 60.0    /* flat spread threshold (experimental path) */

static cand_t* beam_search(int BEAM,int maxbit,int verbose,int adaptive,int *out_n){
    long long DEADCAP = env_i("DEADCAP", 1LL<<20, 1LL<<10, 1LL<<22);
    double    GAPTHR  = env_f("GAPTHR",  0);   /* 0 = built-in default */
    double    FLATTHR = env_f("FLATTHR", 0);
    int       FLATX   = (int)env_i("FLATX", 4, 1, 16);
    long long SAMPBUDGET = env_i("SAMPBUDGET", 1LL<<41, 1LL<<34, 1LL<<48);
    /* widthcap bounds every kept-front size: plain BEAM, the DEAD keep-all
       cap, AND the FLAT keep of FLATX*BEAM -- all three, so no documented
       parameter combination can write past cur[] or children[] */
    long long widthcap = BEAM;
    if(adaptive){
        if(DEADCAP>widthcap) widthcap=DEADCAP;
        if((long long)FLATX*BEAM>widthcap) widthcap=(long long)FLATX*BEAM;
    }
    long long chcap = widthcap*16;
    cand_t *cur=xalloc(widthcap,sizeof(cand_t));
    cand_t *children=xalloc(chcap,sizeof(cand_t));
    double *tmpd=xalloc(chcap>8192?8192:chcap,sizeof(double));
    long long nb=1; cur[0]=(cand_t){0,0,0,0,0.0,0.0};
    for(int m=0;m<maxbit;m++){
        long long nch=nb*16, nn=0;
        for(long long b=0;b<nb;b++) for(int g=0;g<16;g++){
            children[nn]=(cand_t){
                cur[b].A0 |((uint32_t)((g>>0)&1)<<m),
                cur[b].T12|((uint32_t)((g>>1)&1)<<m),
                cur[b].T13|((uint32_t)((g>>2)&1)<<m),
                cur[b].U  |((uint32_t)((g>>3)&1)<<m),
                cur[b].score, 0.0 };
            nn++;
        }
        int dead=0; double gap=0;
        double thrF = GAPTHR>0 ? GAPTHR : DEADTHR_DFLT;
        long long nres;
        int final_level = (m==maxbit-1);
        if(nch<=4096 || !adaptive){
            /* narrow front: score every child on all pairs (this is the
               entire fast path; identical decisions to the classic engine).
               The dead test runs directly on the full-N level contributions,
               so no extra scoring pass and no subsample escalation needed. */
            score_cands(children,(int)nch,m,0,N);
            if(adaptive && !final_level){
                double mx=-1e300; for(long long i=0;i<nch;i++) if(children[i].lsc>mx)mx=children[i].lsc;
                for(long long i=0;i<nch;i++) tmpd[i]=children[i].lsc;
                qsort(tmpd,nch,sizeof(double),cmp_dbl);
                gap = mx - tmpd[nch/2];
                dead = gap < thrF;
            }
            nres=nch;
        }else{
            /* ---- sampled dead test: head children + strided sample, gap =
               max - median(strided) of the level contributions.  A below-
               threshold gap at the subsample is re-tested at full N before it
               counts (subsample noise must not trigger keep-alls).  Skipped
               at the final level, which always cuts. ---- */
            if(!final_level){
                long long ndt = N<(1LL<<24)?N:(1LL<<24);
                long long nhead = 1024, nstr = 4096, nsamp = nhead+nstr;
                cand_t *samp=xalloc(nsamp,sizeof(cand_t));
                for(long long i=0;i<nsamp;i++) samp[i]=children[i<nhead? i : (i-nhead)*nch/nstr];
                score_cands(samp,(int)nsamp,m,0,ndt);
                double thr0 = thrF*1.2;
                double mx=-1e300; for(long long i=0;i<nsamp;i++) if(samp[i].lsc>mx)mx=samp[i].lsc;
                for(long long i=0;i<nstr;i++) tmpd[i]=samp[nhead+i].lsc;
                qsort(tmpd,nstr,sizeof(double),cmp_dbl);
                gap = mx - tmpd[nstr/2];
                dead = gap < thr0;
                if(dead && N>ndt){
                    for(long long i=0;i<nsamp;i++) samp[i]=children[i<nhead? i : (i-nhead)*nch/nstr];
                    score_cands(samp,(int)nsamp,m,0,N);
                    mx=-1e300; for(long long i=0;i<nsamp;i++) if(samp[i].lsc>mx)mx=samp[i].lsc;
                    for(long long i=0;i<nstr;i++) tmpd[i]=samp[nhead+i].lsc;
                    qsort(tmpd,nstr,sizeof(double),cmp_dbl);
                    gap = mx - tmpd[nstr/2];
                    dead = gap < thrF;
                    if(!dead && verbose){ printf("  bit %2d: revived at full N (gap=%.1f)\n",m,gap); fflush(stdout); }
                }
                free(samp);
            }
            nres=0;
            if(!dead){
                /* stage A: rank everything on a deterministic, budgeted subsample */
                long long naa = nch>1000000 ? (1LL<<23) : (1LL<<24);
                if(naa*nch>SAMPBUDGET){ naa=SAMPBUDGET/nch; if(naa<(1LL<<20))naa=1LL<<20; }  /* the 2^20 floor
                   (statistical power) can override the budget: the true per-level bound is
                   max(SAMPBUDGET, 16*DEADCAP*2^20) chain evaluations */
                if(naa>N)naa=N;
                long long off = N>naa ? (long long)((((uint64_t)m*2654435761ULL)^0x9E3779B9ULL)%(uint64_t)(N-naa)) : 0;
                score_cands(children,(int)nch,m,off,off+naa);
                qsort(children,nch,sizeof(cand_t),cmp_cum);
                /* stage B: rescore the shortlist; interior levels capped at
                   2^26 pairs for cost, the final level uses everything - its
                   order decides what phase 3 gets to try */
                long long nfull = final_level ? N : (N<(1LL<<26)?N:(1LL<<26));
                nres = nch<8LL*BEAM?nch:8LL*BEAM;
                score_cands(children,(int)nres,m,0,nfull);
            }
        }
        if(dead){
            /* keep all children (stratified cap); no score contribution */
            long long keep;
            if(nch>DEADCAP){ long long step=nch/DEADCAP;
                /* Children are laid out 16-per-parent, so a bare stride that
                   shares a factor with 16 aliases onto fixed child slots
                   (step=16 would keep ONLY the g=0 child of each parent -
                   all four new bits forced to 0, truth lost w.p. 15/16).
                   The i-dependent offset cycles through the child slots.
                   For DEADCAP<nch<2*DEADCAP (step=1, offset 0) this is a
                   head-keep of top-ranked parents' children; acceptable at a
                   dead level, where order carries no signal. */
                for(long long i=0;i<DEADCAP;i++){
                    long long off = step>=16 ? (i&15) : (step>1 ? i%step : 0);
                    cur[i]=children[i*step+off];
                }
                keep=DEADCAP;
            }else{ memcpy(cur,children,sizeof(cand_t)*nch); keep=nch; }
            nb=keep;
            if(verbose){ printf("  bit %2d: DEAD (gap=%.1f < %.1f) keeping %lld of %lld\n",m,gap,thrF,keep,nch); fflush(stdout); }
            if(verbose) truth_rank(cur,nb,m,nb,"deadkeep");
            continue;
        }
        qsort(children,nres,sizeof(cand_t),cmp_cum);
        long long keep = nres<BEAM?nres:BEAM;
        if(adaptive){
            /* ---- flat test on the rescored shortlist's level contributions:
               separation present (not dead) but max - median below threshold
               means this level's order is mostly noise, and cutting to BEAM
               would prune on noise; keep FLATX x BEAM instead ---- */
            long long nf = nres<8192?nres:8192;
            for(long long i=0;i<nf;i++) tmpd[i]=children[i].lsc;
            qsort(tmpd,nf,sizeof(double),cmp_dbl);
            double spread = tmpd[nf-1]-tmpd[nf/2];
            double fthr = FLATTHR>0?FLATTHR:FLATTHR_DFLT;
            if(spread<fthr){
                long long kw=(long long)FLATX*BEAM; if(kw>nres)kw=nres;
                if(kw>widthcap)kw=widthcap;   /* cannot bind (widthcap >= FLATX*BEAM); kept as an explicit invariant */
                if(kw>keep){ keep=kw;
                    if(verbose){ printf("  bit %2d: FLAT (spread=%.1f < %.1f) keeping %lld\n",m,spread,fthr,keep); fflush(stdout); } }
            }
        }
        if(verbose) truth_rank(children,nres,m,keep,"level");
        for(long long i=0;i<keep;i++){ cur[i]=children[i]; cur[i].score+=cur[i].lsc; cur[i].lsc=0; }
        nb=keep;
        if(verbose){
            if(adaptive) printf("  bit %2d: top=(A0=%08x T12=%08x T13=%08x U=%08x) score=%.1f gap=%.1f (kept %lld)\n",
                m,cur[0].A0,cur[0].T12,cur[0].T13,cur[0].U,cur[0].score,gap,nb);
            else printf("  bit %2d: top=(A0=%08x T12=%08x T13=%08x U=%08x) score=%.1f\n",
                m,cur[0].A0,cur[0].T12,cur[0].T13,cur[0].U,cur[0].score);
            fflush(stdout); }
    }
    free(children); free(tmpd);
    *out_n=(int)nb;
    return cur;
}

static const int SH[4]={1,3,6,11};
/* B0 = rk[11][0] from A0 = rk[12][0]:  A0 = rol(B0 + rol(DELTA[0],12),1) */
static uint32_t B0_from_A0(uint32_t A0){ return ror32(A0,1) - rol32(LEA_DELTA[12&3],12); }

/* ---------- final brute force over the top bits ----------
   Loop nesting hoists the per-chain key-schedule work:
     hA0 -> A0 -> k0 -> t0[r]      (outermost)
     hU  -> A1 = U^B0(A0) -> k1 -> t1[r]
     hT12-> A2 = A1^T12 -> k2 -> t2[r]
     hT13-> A3 = A1^T13 -> k3 -> t3[r]  then encrypt
*/
typedef struct { int tid; uint32_t kA0,kT12,kT13,kU; int knownbits;
                 const uint8_t*kp; const uint8_t*kc; int found; uint32_t key[4]; } bf_job_t;
/* invert a single chain: given t_j at round index 12, recover k_j */
static uint32_t invert_chain(uint32_t v,int j,int rounds){
    for(int i=rounds-1;i>=0;i--) v = ror32(v,SH[j]) - rol32(LEA_DELTA[i&3], i+j);
    return v;
}
static void expand_chain(uint32_t k,int j,int rounds,uint32_t*out){
    uint32_t v=k;
    for(int i=0;i<rounds;i++){ v = rol32(v + rol32(LEA_DELTA[i&3], i+j), SH[j]); out[i]=v; }
}
static void* bfw(void*arg){
    bf_job_t*J=(bf_job_t*)arg;
    int kb=J->knownbits, hb=32-kb;
    uint32_t HM=(hb>=32)?0xffffffffu:((1u<<hb)-1);
    uint32_t lowmask=(kb==32)?0xffffffffu:((1u<<kb)-1);
    uint32_t pt[4],ct[4];
    for(int k=0;k<4;k++){
        pt[k]=(uint32_t)J->kp[4*k]|((uint32_t)J->kp[4*k+1]<<8)|((uint32_t)J->kp[4*k+2]<<16)|((uint32_t)J->kp[4*k+3]<<24);
        ct[k]=(uint32_t)J->kc[4*k]|((uint32_t)J->kc[4*k+1]<<8)|((uint32_t)J->kc[4*k+2]<<16)|((uint32_t)J->kc[4*k+3]<<24);
    }
    uint32_t T0[13],T1[13],T2[13],T3[13];
    uint32_t nspan = HM+1;
    /* split hA0 across threads */
    for(uint32_t hA0=J->tid; hA0<nspan; hA0+=NTHREADS){
        uint32_t A0=(J->kA0&lowmask)|(hA0<<kb);
        uint32_t k0=invert_chain(A0,0,13); expand_chain(k0,0,13,T0);
        uint32_t B0=B0_from_A0(A0);
        for(uint32_t hU=0;hU<nspan;hU++){
            uint32_t U=(J->kU&lowmask)|(hU<<kb);
            uint32_t A1=U^B0;
            uint32_t k1=invert_chain(A1,1,13); expand_chain(k1,1,13,T1);
            for(uint32_t h12=0;h12<nspan;h12++){
                uint32_t T12=(J->kT12&lowmask)|(h12<<kb);
                uint32_t A2=A1^T12;
                uint32_t k2=invert_chain(A2,2,13); expand_chain(k2,2,13,T2);
                for(uint32_t h13=0;h13<nspan;h13++){
                    uint32_t T13=(J->kT13&lowmask)|(h13<<kb);
                    uint32_t A3=A1^T13;
                    uint32_t k3=invert_chain(A3,3,13); expand_chain(k3,3,13,T3);
                    uint32_t x0=pt[0],x1=pt[1],x2=pt[2],x3=pt[3];
                    for(int r=0;r<13;r++){
                        uint32_t y0=rol32((x0^T0[r])+(x1^T1[r]),9);
                        uint32_t y1=ror32((x1^T2[r])+(x2^T1[r]),5);
                        uint32_t y2=ror32((x2^T3[r])+(x3^T1[r]),3);
                        x3=x0; x0=y0; x1=y1; x2=y2;
                    }
                    if(x0==ct[0]&&x1==ct[1]&&x2==ct[2]&&x3==ct[3]){
                        J->found=1; J->key[0]=k0;J->key[1]=k1;J->key[2]=k2;J->key[3]=k3;
                        return NULL;
                    }
                }
            }
        }
    }
    return NULL;
}

static void usage(const char*prog){
    fprintf(stderr,
        "usage: %s [N1 [N2 [BEAM [THREADS [MAXBIT [SEED [SIMKEY]]]]]]]\n"
        "  N1      pairs per combo in phase 1          [2^10 .. 2^36]  (default 2^20)\n"
        "  N2      pairs in phase 2                    [2^10 .. 2^36]  (default 2^27)\n"
        "  BEAM    beam width in phase 2               [1 .. 1048576]  (default 64)\n"
        "  THREADS worker threads                      [1 .. 1024]     (default 8)\n"
        "  MAXBIT  bits/word recovered by the beam     [8 .. 31]       (default 24;\n"
        "          phase 3 then brute-forces 2^(4*(32-MAXBIT)) per candidate)\n"
        "  SEED    PRNG seed for plaintext generation  (number; 0x.. accepted)\n"
        "  SIMKEY  32 hex chars: run against a built-in local oracle with this key\n"
        "          instead of the HTTP oracle at ORACLE_HOST:ORACLE_PORT\n"
        "environment: ADAPTIVE=1 enables the adaptive beam (default 0 = classic\n"
        "  search semantics at every width); with it: DEADCAP [2^10..2^22],\n"
        "  GAPTHR, FLATTHR (absolute chi^2-scale units; the null is N-independent),\n"
        "  FLATX [1..16], SAMPBUDGET [2^34..2^48].\n"
        "  PH3MAX [1..2^20] (default 1024) caps phase-3 candidate walks\n", prog);
    exit(2);
}

int main(int argc,char**argv){
    long long N1 = argc>1? atoll(argv[1]) : (1LL<<20);
    long long N2 = argc>2? atoll(argv[2]) : (1LL<<27);
    long long BEAML = argc>3? atoll(argv[3]) : 64;
    NTHREADS     = argc>4? atoi(argv[4]) : 8;
    int MAXBIT   = argc>5? atoi(argv[5]) : 24;
    uint64_t SEED= argc>6? strtoull(argv[6],NULL,0) : 0xC0FFEE;
    if(N1<(1LL<<10)||N1>(1LL<<36)){ fprintf(stderr,"N1 out of range\n"); usage(argv[0]); }
    if(N2<(1LL<<10)||N2>(1LL<<36)){ fprintf(stderr,"N2 out of range\n"); usage(argv[0]); }
    if(BEAML<1||BEAML>(1LL<<20)){ fprintf(stderr,"BEAM out of range\n"); usage(argv[0]); }
    int BEAM=(int)BEAML;
    if(NTHREADS<1||NTHREADS>1024){ fprintf(stderr,"THREADS out of range\n"); usage(argv[0]); }
    if(MAXBIT<8||MAXBIT>31){ fprintf(stderr,"MAXBIT out of range\n"); usage(argv[0]); }
    long long PH3MAX = env_i("PH3MAX", 1024, 1, 1LL<<20);
    /* argv[7]: if present, run in LOCAL SIMULATION mode with this hex key */
    if(argc>7){
        const char*ks=argv[7];
        if(strlen(ks)!=32 || strspn(ks,"0123456789abcdefABCDEF")!=32){
            fprintf(stderr,"SIMKEY must be exactly 32 hex characters\n"); usage(argv[0]); }
        ORA_SIM=1;
        for(int k=0;k<4;k++){
            uint32_t v=0;
            for(int b=0;b<4;b++){ unsigned x; sscanf(ks+2*(4*k+b),"%2x",&x); v |= (uint32_t)x<<(8*b); }
            ORA_SIM_KEY[k]=v;
        }
        printf("### LOCAL SIMULATION MODE (self-test): oracle key = %s\n", ks);
    }

    truth_init();
    struct timespec T0,T1; clock_gettime(CLOCK_MONOTONIC,&T0);
    rng_t rng; rseed(&rng,SEED);

    /* grab two known plaintext/ciphertext blocks: block 0 drives the phase-3
       search, block 1 is only used to re-verify the recovered key at the end */
    uint8_t kp[32],kc[32];
    for(int i=0;i<32;i++) kp[i]=(uint8_t)(i*37+11);
    ora_encrypt(kp,kc,2);
    printf("known pt = "); for(int i=0;i<32;i++) printf("%02x",kp[i]); printf("\n");
    printf("known ct = "); for(int i=0;i<32;i++) printf("%02x",kc[i]); printf("\n\n");

    /* ============ PHASE 1: find the 6-bit combo ============
       For each of the 64 candidate plaintext structures, run a short beam search
       over key bit positions 0..6 and record the best cumulative score.  The
       correct combo reaches the huge bit-6 bias (~4e-2) and scores ~1000x higher. */
    printf("=== PHASE 1: conditional-differential combo search (N1=%lld pairs/combo) ===\n",N1);
    fflush(stdout);
    alloc_cols(N1);
    double allsc[64];
    int bestcombo=-1; double bestsc=-1;
    for(int combo=0;combo<64;combo++){
        gather(N,combo,&rng);
        int tn; cand_t *tiny=beam_search(8,7,0,0,&tn);   /* bits 0..6 -> scores Delta b1[1..7] */
        allsc[combo]=tiny[0].score; free(tiny);
        if(allsc[combo]>bestsc){bestsc=allsc[combo];bestcombo=combo;}
        printf("  combo %2d: score=%10.1f\n",combo,allsc[combo]); fflush(stdout);
    }
    free_cols();
    double second=0; for(int i=0;i<64;i++) if(i!=bestcombo && allsc[i]>second) second=allsc[i];
    double ratio=bestsc/(second+1e-9);
    printf("\n>>> best combo=%d  score=%.1f   second best=%.1f   (ratio %.1f)\n",
        bestcombo,bestsc,second,ratio);
    if(ratio<10.0)
        printf("WARNING: weak combo separation (expected >50x at N1=2^20); "
               "phase 2 may be built on the wrong plaintext structure\n");
    printf("\n"); fflush(stdout);

    /* ============ PHASE 2: bit-by-bit beam search ============ */
    int ADAPTIVE0 = (int)env_i("ADAPTIVE", 0, 0, 1);
    long long DEADCAP = ADAPTIVE0 ? env_i("DEADCAP", 1LL<<20, 1LL<<10, 1LL<<22) : 0;
    long long FLATXW  = ADAPTIVE0 ? env_i("FLATX", 4, 1, 16)*(long long)BEAM : 0;
    long long widthcap = BEAM;
    if(DEADCAP>widthcap) widthcap=DEADCAP;
    if(FLATXW>widthcap) widthcap=FLATXW;
    printf("=== PHASE 2: bit-by-bit key recovery (N2=%lld pairs, beam=%d%s) ===\n",N2,BEAM,ADAPTIVE0?", adaptive":"");
    /* transient scoring counters (32 B/candidate scored) can add up to
       ~0.5 GiB on top of this in adaptive mode after a dead keep-all */
    printf("estimated memory: %.1f GiB (pair data %.1f + candidate pools %.1f)\n",
        ((double)N2*32.0 + (double)widthcap*17.0*sizeof(cand_t))/(1024.0*1024.0*1024.0),
        (double)N2*32.0/(1024.0*1024.0*1024.0),
        (double)widthcap*17.0*sizeof(cand_t)/(1024.0*1024.0*1024.0));
    fflush(stdout);
    alloc_cols(N2);
    gather(N,bestcombo,&rng);
    printf("data gathered: N=%lld pairs, noise(bias)=%.3e\n",N,0.5/sqrt((double)N)); fflush(stdout);

    int nb;
    cand_t *beam=beam_search(BEAM,MAXBIT,1,ADAPTIVE0,&nb);
    clock_gettime(CLOCK_MONOTONIC,&T1);
    printf("\nPhase 2 done (%.0f s). Recovered low %d bits.\n",
        (T1.tv_sec-T0.tv_sec)+1e-9*(T1.tv_nsec-T0.tv_nsec), MAXBIT);

    /* ============ PHASE 3: brute force the top bits, try the top beam candidates ============ */
    long long NTRY = nb;
    if(NTRY>PH3MAX){
        printf("\n(final beam holds %d candidates; walking the top %lld by score - raise PH3MAX to walk more)\n",nb,PH3MAX);
        NTRY=PH3MAX;
    }
    printf("\n=== PHASE 3: brute force top %d bits x4 = 2^%d (over %lld of %d beam candidates) ===\n",
        32-MAXBIT,4*(32-MAXBIT),NTRY,nb);
    fflush(stdout);
    for(long long t=0;t<NTRY;t++){
        printf("  trying beam candidate %lld: A0=%08x T12=%08x T13=%08x U=%08x\n",
            t,beam[t].A0,beam[t].T12,beam[t].T13,beam[t].U); fflush(stdout);
        pthread_t *th=xalloc(NTHREADS,sizeof *th); bf_job_t *jobs=xalloc(NTHREADS,sizeof *jobs);
        for(int i=0;i<NTHREADS;i++){
            jobs[i]=(bf_job_t){i,beam[t].A0,beam[t].T12,beam[t].T13,beam[t].U,MAXBIT,kp,kc,0,{0,0,0,0}};
            if(pthread_create(&th[i],NULL,bfw,&jobs[i])){ fprintf(stderr,"thread creation failed\n"); exit(3); }
        }
        int found=0; uint32_t key[4];
        for(int i=0;i<NTHREADS;i++){ pthread_join(th[i],NULL); if(jobs[i].found){found=1;memcpy(key,jobs[i].key,16);} }
        free(th); free(jobs);
        if(found){
            printf("\n*** KEY FOUND: %08x %08x %08x %08x ***\n",key[0],key[1],key[2],key[3]);
            printf("KEY = ");
            for(int k=0;k<4;k++) for(int b=0;b<4;b++) printf("%02x",(key[k]>>(8*b))&0xff);
            printf("\n");
            /* independent check: re-encrypt the second known block with lea.h */
            lea_ctx c; lea_setkey(&c,key,13);
            uint32_t x[4];
            for(int k=0;k<4;k++){ const uint8_t*s=kp+16+4*k;
                x[k]=(uint32_t)s[0]|((uint32_t)s[1]<<8)|((uint32_t)s[2]<<16)|((uint32_t)s[3]<<24); }
            lea_enc_nr(&c,x,13);
            int ok=1;
            for(int k=0;k<4;k++){ const uint8_t*s=kc+16+4*k;
                uint32_t w=(uint32_t)s[0]|((uint32_t)s[1]<<8)|((uint32_t)s[2]<<16)|((uint32_t)s[3]<<24);
                if(x[k]!=w) ok=0; }
            printf("verification on a second known block: %s\n", ok?"PASS":"FAIL");
            printf("oracle usage: %lld queries, %lld blocks (2^%.2f)\n",
                ORA_QUERIES,ORA_BLOCKS,log2((double)ORA_BLOCKS));
            if(ORA_SIM){
                int simok=(memcmp(key,ORA_SIM_KEY,16)==0);
                printf("simulation self-check: recovered key %s the oracle key\n",
                    simok?"matches":"DOES NOT MATCH");
                ok&=simok;
            }
            clock_gettime(CLOCK_MONOTONIC,&T1);
            printf("\ntotal time: %.0f s\n",(T1.tv_sec-T0.tv_sec)+1e-9*(T1.tv_nsec-T0.tv_nsec));
            free(beam); free_cols();
            return ok?0:1;
        }
    }
    clock_gettime(CLOCK_MONOTONIC,&T1);
    printf("oracle usage: %lld queries, %lld blocks (2^%.2f)\n",
        ORA_QUERIES,ORA_BLOCKS,log2((double)ORA_BLOCKS));
    printf("\nATTACK FAILED: key not found in the top %lld beam candidates "
           "(N1=%lld N2=%lld BEAM=%d MAXBIT=%d SEED=%llu).\n"
           "Retry with a wider beam and/or more data (run.sh escalates automatically).\n",
           NTRY,N1,N2,BEAM,MAXBIT,(unsigned long long)SEED);
    printf("total time: %.0f s\n",(T1.tv_sec-T0.tv_sec)+1e-9*(T1.tv_nsec-T0.tv_nsec));
    free(beam); free_cols();
    return 1;
}
