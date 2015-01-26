/*
 * bwape.h
 *
 *  Created on: 2014Äê7ÔÂ17ÈÕ
 *      Author: Administrator
 */

#ifndef BWAPE_H_
#define BWAPE_H_
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include "bwtaln.h"
#include "kvec.h"
#include "bntseq.h"
#include "utils.h"
#include "stdaln.h"

 inline uint64_t hash_64(uint64_t key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}
// here v>=u. When ii is set, we check insert size with ii; otherwise with opt->max_isize
#define __pairing_aux(u,v) do {		\
	bwtint_t l = ((v) >> 32) + p[(v)&1] ->len -((u)>>32);	\
	if((u) != (uint64_t) -1 &&(v) >>32 > (u)>>32 && l >= max_len \
			&& ((ii->high && l <=ii -> high_bayesian) || (ii->high == 0 && l <= opt->max_isize))) \
			{	\
				uint64_t s = d->aln[(v)&1].a [(uint32_t)(v)>>1].score + d->aln[(u)&1].a[(uint32_t)(u)>>1].score; \
				s *=10;\
				if(ii->high) s+=(int) (-4.343*log(0.5 * erfc(M_SQRT1_2*fabs(l - ii->avg)/ii ->std)) + 0.499); \
				s = s<<32 | (uint32_t)hash_64((u)>>32<<32 | (v)>>32);		\
				if (s>>32 == o_score>>32) ++o_n;							\
				else if ( s>>32 < o_score<<32) {subo_n += o_n; o_n=1;}			\
				else ++subo_n;				\
				if(s<o_score) subo_score = o_score, o_score = s, o_pos[(u)&1] = (u), o_pos[(v)&1] =(v); 		\
				else if(s < subo_score) subo_score =s;			\
			}				\
}while(0)			\


#define __pairing_aux2(q, w) do {										\
		const bwt_aln1_t *r = d->aln[(w)&1].a + ((uint32_t)(w)>>1);		\
		(q)->extra_flag |= SAM_FPP;										\
		if ((q)->pos != (w)>>32 || (q)->strand != r->a) {				\
			(q)->n_mm = r->n_mm; (q)->n_gapo = r->n_gapo; (q)->n_gape = r->n_gape; (q)->strand = r->a; \
			(q)->score = r->score;										\
			(q)->pos = (w)>>32;											\
			if ((q)->mapQ > 0) ++cnt_chg;								\
		}																\
	} while (0)

#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
	int n;
	bwtint_t *a;
} poslist_t;

typedef struct {
	double avg, std, ap_prior;
	bwtint_t low, high, high_bayesian;
} isize_info_t;



typedef struct {
	kvec_t(uint64_t) arr;
	kvec_t(uint64_t) pos[2];
	kvec_t(bwt_aln1_t) aln[2];
} pe_data_t;

#define MIN_HASH_WIDTH 1000

extern int g_log_n[256]; // in bwase.c


void bwase_initialize();
pe_opt_t *bwa_init_pe_opt();
void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);
void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);
ubyte_t * bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq, bntseq_t *ntbns);
int bwa_approx_mapQ(const bwa_seq_t *p, int mm);
void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2);
bntseq_t *bwa_open_nt(const char *prefix);
ubyte_t *bwa_paired_sw(const bntseq_t *bns, const ubyte_t *_pacseq, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii,int mode);
void bwa_print_sam_SQ(const bntseq_t *bns);
void bwa_print_sam_PG();
 int pairing(bwa_seq_t *p[2], pe_data_t *d, const pe_opt_t *opt, int s_mm, const isize_info_t *ii);
 int infer_isize(int n_seqs, bwa_seq_t *seqs[2], isize_info_t *ii, double ap_prior, int64_t L);

#ifdef __cplusplus
}
#endif


#endif /* BWAPE_H_ */
