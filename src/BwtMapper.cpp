/*
 * BwtMapper.cpp
 *
 *  Created on: 2014Äê7ÔÂ9ÈÕ
 *      Author: Administrator
 */

#include "BwtMapper.h"
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <sstream>
#include "../libbwa/bwtgap.h"
#include "../libbwa/bwase.h"
#include "../libbwa/bwape.h"
#include "../libbwa/khash.h"
#include "../libmpu/Error.h"
#include <algorithm>
#include <gperftools/profiler.h>


using namespace std;
//extern string Prefix;
extern void notice(const char*, ...);
extern void warning(const char*, ...);
extern void error(const char*, ...);
KHASH_MAP_INIT_INT64(64, poslist_t)
kh_64_t *g_hash;
//#define DEBUG 0
#ifdef DEBUG
#define DBG(CODE) CODE
#else
#define DBG(CODE)
#endif

#define N_OCC 3

#ifdef HAVE_PTHREAD
#define THREAD_BLOCK_SIZE 1024
#include <pthread.h>
static pthread_mutex_t g_seq_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

extern ubyte_t * bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs,
	ubyte_t *_pacseq, bntseq_t *ntbns);
extern void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p,
	const bwa_seq_t *mate, int mode, int max_top2);
extern void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s,
	int set_main, int n_multi);
extern int bwa_approx_mapQ(const bwa_seq_t *p, int mm);
extern int bwt_cal_width(const bwt_t *rbwt, int len, const ubyte_t *str,
	bwt_width_t *width);
extern bwa_seq_t *bwa_read_bam(bwa_seqio_t *bs, int n_needed, int *n, int is_comp, int trim_qual);
extern int bwa_trim_read(int trim_qual, bwa_seq_t *p);
static void bwa_cal_sa_reg_gap(int tid, bwt_t * const bwt[2], int n_seqs,
	bwa_seq_t *seqs, const gap_opt_t *opt, const BwtIndexer * Indexer)
{
	int i, max_l = 0, max_len;
	gap_stack_t *stack;
	bwt_width_t *w[2], *seed_w[2];
	const ubyte_t *seq[2];
	gap_opt_t local_opt = *opt;

	// initiate priority stack
	for (i = max_len = 0; i != n_seqs; ++i)
		if (seqs[i].len > max_len)
			max_len = seqs[i].len;
	if (opt->fnr > 0.0)
		local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo)
		local_opt.max_gapo = local_opt.max_diff;
	stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo,
		local_opt.max_gape, &local_opt);

	seed_w[0] = (bwt_width_t*)calloc(opt->seed_len + 1, sizeof(bwt_width_t));
	seed_w[1] = (bwt_width_t*)calloc(opt->seed_len + 1, sizeof(bwt_width_t));
	w[0] = w[1] = 0;
	uint32_t unmapped_num = 0;
	for (i = 0; i != n_seqs; ++i)
	{
		bwa_seq_t *p = seqs + i;
		//notice("Read %s",p->name);
		/*decoupling multithread*/
		/*
  #ifdef HAVE_PTHREAD
  if (opt->n_threads > 1)
  {
  pthread_mutex_lock(&g_seq_lock);
  if (p->tid < 0)
  { // unassigned
  int j;
  for (j = i; j < n_seqs && j < i + THREAD_BLOCK_SIZE; ++j)
  seqs[j].tid = tid;
  }
  else if (p->tid != tid)
  {
  pthread_mutex_unlock(&g_seq_lock);
  continue;
  }
  pthread_mutex_unlock(&g_seq_lock);
  }
  #endif*/
		//if (strcmp(p->name, "ERR018525.148353") == 0)
		//{
		// fprintf(stderr, "ERR018525.148353	is coming\n");
		//}
		p->sa = 0;
		p->type = BWA_TYPE_NO_MATCH;
		p->c1 = p->c2 = 0;
		p->n_aln = 0;
		p->aln = 0;
		if (/* drand48() >opt->frac ||*/p->filtered)//|| Indexer->IsReadFiltered(p->seq, p->qual, p->len))//here I use BWA_TYPE_UNIQUE as a flag of being downsampled
		{
			//if(strcmp(p->name,"ERR018525.148353") ==0)
			//{
			//fprintf(stderr,"ERR018525.148353	is filtered\n" );
			//exit(1);
			//}

			unmapped_num++;
			continue;
		}

		seq[0] = p->seq;
		seq[1] = p->rseq;
		if (max_l < p->len)
		{
			max_l = p->len;
			w[0] = (bwt_width_t*)realloc(w[0],
				(max_l + 1) * sizeof(bwt_width_t));
			w[1] = (bwt_width_t*)realloc(w[1],
				(max_l + 1) * sizeof(bwt_width_t));
			memset(w[0], 0, (max_l + 1) * sizeof(bwt_width_t));
			memset(w[1], 0, (max_l + 1) * sizeof(bwt_width_t));
		}
		bwt_cal_width(bwt[0], p->len, seq[0], w[0]);
		bwt_cal_width(bwt[1], p->len, seq[1], w[1]);
		if (opt->fnr > 0.0)
			local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		local_opt.seed_len =
			opt->seed_len < p->len ? opt->seed_len : 0x7fffffff;
		if (p->len > opt->seed_len)
		{
			bwt_cal_width(bwt[0], opt->seed_len,
				seq[0] + (p->len - opt->seed_len), seed_w[0]);
			bwt_cal_width(bwt[1], opt->seed_len,
				seq[1] + (p->len - opt->seed_len), seed_w[1]);
		}
		// core function

		p->aln = bwt_match_gap(bwt, p->len, seq, w,
			p->len <= opt->seed_len ? 0 : seed_w, &local_opt, &p->n_aln,
			stack);
		/*if(strcmp(p->name,"ERR013170.1716") ==0)
		 {
		 fwrite(w[0], sizeof(bwt_width_t), max_l, stdout);
		 fwrite(w[1], sizeof(bwt_width_t), max_l, stdout);
		 fwrite(seed_w[0], sizeof(bwt_width_t),opt->seed_len, stdout);
		 fwrite(seed_w[1], sizeof(bwt_width_t),opt->seed_len, stdout);
		 //fwrite(&local_opt, sizeof(gap_opt_t), 1, stdout);
		 //fwrite(stack, sizeof(gap_stack_t), 1, stdout);
		 fprintf(stderr,"n_stacks:%d\tbest:%d\tn_entries:%d\tstacks_n_entries:%d\tstacks_last_diff_pos:%d\n",stack->n_stacks,stack->best,stack->n_entries,stack->stacks->n_entries,stack->stacks->stack->last_diff_pos);
		 for(int iter=0;iter<p->n_aln;++iter)
		 fprintf(stderr,"l:%d\tk:%d\t strand:%d\t score:%d\n",(p->aln+iter)->l,(p->aln+iter)->k,(p->aln+iter)->a, (p->aln+iter)->score);
		 }*/
		// store the alignment
		//free(p->name); free(p->seq); free(p->rseq); free(p->qual);
		//p->name = 0; p->seq = p->rseq = p->qual = 0;
	}
	notice("RollingHash filtered %d reads...", unmapped_num);
	free(seed_w[0]);
	free(seed_w[1]);
	free(w[0]);
	free(w[1]);
	gap_destroy_stack(stack);
}
BwtMapper::BwtMapper()
{}
int BwtMapper::bwa_cal_pac_pos(BwtIndexer& BwtIndex, int n_seqs,bwa_seq_t *seqs, int max_mm, float fnr)
{
	int i, j;
	//char str[1024];
	bwt_t *bwt;
	// load forward SA
	//strcpy(str, prefix); strcat(str, ".bwt");
	bwt = BwtIndex.bwt_d; //bwt_restore_bwt(str);
	//strcpy(str, prefix); strcat(str, ".sa");
	//bwt_restore_sa(str, bwt);
	for (i = 0; i != n_seqs; ++i)
	{
		if (seqs[i].strand)
			bwa_cal_pac_pos_core(bwt, 0, &seqs[i], max_mm, fnr);
		for (j = 0; j < seqs[i].n_multi; ++j)
		{
			bwt_multi1_t *p = seqs[i].multi + j;
			if (p->strand)
				p->pos = bwt_sa(bwt, p->pos); //transform pos into actual position
		}
	}
	//bwt_destroy(bwt);
	// load reverse BWT and SA
	//strcpy(str, prefix); strcat(str, ".rbwt");
	bwt = BwtIndex.rbwt_d; //bwt_restore_bwt(str);
	//strcpy(str, prefix); strcat(str, ".rsa");
	//bwt_restore_sa(str, bwt);
	for (i = 0; i != n_seqs; ++i)
	{
		if (!seqs[i].strand)
			bwa_cal_pac_pos_core(0, bwt, &seqs[i], max_mm, fnr);
		for (j = 0; j < seqs[i].n_multi; ++j)
		{
			bwt_multi1_t *p = seqs[i].multi + j;
			if (!p->strand)
				p->pos = bwt->seq_len - (bwt_sa(bwt, p->pos) + seqs[i].len);
		}
	}
	//bwt_destroy(bwt);
	return 1;
}
//bool BwtMapper::Skip(BwtIndexer& BwtIndex, int mode, const char* seq1, const char* qual1, const char* seq2 = 0, const char* qual2 = 0, int len = 0, double frac = 1.)
//{
//	if (drand48() >frac)
//	{
//		//PassFlag[index] = 1;
//		return true;
//	}
//	else
//	{
//		ubyte_t* Seq = new ubyte_t[len];
//		ubyte_t* Qual = new ubyte_t[len];
//		ubyte_t* Seq2 = 0;
//		ubyte_t* Qual2 = 0;
//		for (int i = 0; i != len; ++i)
//			Seq[i] = nst_nt4_table[(int)seq1[i]];
//		if (qual1) { // copy quality
//			if (mode | BWA_MODE_IL13)
//				for (int i = 0; i < len; ++i) Qual[i] = qual1[i] - 31;
//			else
//				for (int i = 0; i < len; ++i) Qual[i] = qual1[i];
//		}
//		if (seq2 != 0 && qual2 != 0)
//		{
//			Seq2 = new ubyte_t[len];
//			Qual2 = new ubyte_t[len];
//			for (int i = 0; i != len; ++i)
//				Seq2[i] = nst_nt4_table[(int)seq2[i]];
//			if (qual2) { // copy quality
//				if (mode | BWA_MODE_IL13)
//					for (int i = 0; i < len; ++i) Qual2[i] = qual2[i] - 31;
//				else
//					for (int i = 0; i < len; ++i) Qual2[i] = qual2[i];
//			}
//		}
//		if (BwtIndex.IsReadFiltered(Seq, Qual, len) && BwtIndex.IsReadFiltered(Seq2, Qual2, len))
//		{
//			delete[] Seq, Qual;
//			if (Seq2 != 0)
//			{
//				delete[] Seq2, Qual2;
//			}
//			return true;
//		}
//		else
//		{
//			delete[] Seq, Qual;
//			if (Seq2 != 0)
//			{
//				delete[] Seq2, Qual2;
//			}
//			return false;
//		}
//	}
//	return false;
//}
typedef struct
{
	kvec_t(bwt_aln1_t)
		aln;
}
aln_buf_t;
//specific wrapper for bwa_read_seq[begin]
#include "../libbwa/kseq.h"
#include "../libbwa/bamlite.h"
#include "zlib.h"
KSEQ_INIT_FPC(gzFile, gzread)
struct __bwa_seqio_t {
	// for BAM input
	int is_bam, which; // 1st bit: read1, 2nd bit: read2, 3rd: SE
	bamFile fp;
	// for fastq input
	kseq_t *ks;
};
#define BARCODE_LOW_QUAL 13

static bwa_seq_t* bwa_read_seq_with_hash(BwtIndexer* BwtIndex, bwa_seqio_t *bs, int n_needed, int *n, int mode, int trim_qual, double frac, uint32_t seed)
{
	
	struct drand48_data randBuffer;
	srand48_r(seed,&randBuffer);
	bwa_seq_t *seqs, *p;
	kseq_t *seq = bs->ks;
	int n_seqs, l, i, is_comp = mode&BWA_MODE_COMPREAD, is_64 = mode&BWA_MODE_IL13, l_bc = mode >> 24;
	long n_trimmed = 0, n_tot = 0;

	if (l_bc > 15) {
		fprintf(stderr, "[%s] the maximum barcode length is 15.\n", __func__);
		return 0;
	}
	if (bs->is_bam) return bwa_read_bam(bs, n_needed, n, is_comp, trim_qual); // l_bc has no effect for BAM input
	n_seqs = 0;
	//seqs = (bwa_seq_t*)calloc(n_needed, sizeof(bwa_seq_t));
	seqs = new bwa_seq_t[n_needed];
	//ProfilerStart("FastPopCon.prof");
	/*while ((l = kseq_read(seq)) >= 0) {*/

	/*if (Skip(BwtIndex, mode, seq->seq.s, seq->qual.s, 0, 0, seq->seq.l, frac)) continue;*/
	//while (1)
	//{
	//	double rand_num = 0;
	//	drand48_r(&randBuffer,&rand_num);
	//	if ( rand_num> frac)
	//	{
	//		//notice("before length:%d", seq->seq.l);
	//		if ((l = kseq_read4_fpc(seq)) < 0) break;
	//		//notice("after length:%d", seq->seq.l);
	//		continue;
	//	}
	//	if ((l = kseq_read3_fpc(seq)) < 0) break;

	//	if (seq->seq.l <= l_bc) continue; // sequence length equals or smaller than the barcode length
	//	p = &seqs[n_seqs++];

	//	if (l_bc) { // then trim barcode
	//		for (i = 0; i < l_bc; ++i)
	//			p->bc[i] = (seq->qual.l && seq->qual.s[i] - 33 < BARCODE_LOW_QUAL) ? tolower(seq->seq.s[i]) : toupper(seq->seq.s[i]);
	//		p->bc[i] = 0;
	//		for (; i < seq->seq.l; ++i)
	//			seq->seq.s[i - l_bc] = seq->seq.s[i];
	//		seq->seq.l -= l_bc; seq->seq.s[seq->seq.l] = 0;
	//		if (seq->qual.l) {
	//			for (i = l_bc; i < seq->qual.l; ++i)
	//				seq->qual.s[i - l_bc] = seq->qual.s[i];
	//			seq->qual.l -= l_bc; seq->qual.s[seq->qual.l] = 0;
	//		}
	//		l = seq->seq.l;
	//	}
	//	else p->bc[0] = 0;
	//	p->tid = -1; // no assigned to a thread
	//	p->qual = 0;
	//	//p->count=0;
	//	p->full_len = p->clip_len = p->len = l;
	//	n_tot += p->full_len;
	//	p->seq = (ubyte_t*)calloc(p->len, 1);
	//	for (i = 0; i != p->full_len; ++i)
	//		p->seq[i] = nst_nt4_table[(int)seq->seq.s[i]];
	//	if (seq->qual.l) { // copy quality
	//		if (is_64)
	//			for (i = 0; i < seq->qual.l; ++i) seq->qual.s[i] -= 31;
	//		p->qual = (ubyte_t*)strdup((char*)seq->qual.s);
	//		if (trim_qual >= 1) n_trimmed += bwa_trim_read(trim_qual, p);
	//	}

	//	if (BwtIndex->IsReadFiltered(p->seq, p->qual, p->len))
	//	{
	//		p->filtered |= 1;
	//		if (n_seqs == n_needed) break;
	//		continue;
	//	}

	//	p->rseq = (ubyte_t*)calloc(p->full_len, 1);
	//	memcpy(p->rseq, p->seq, p->len);
	//	//fprintf(stderr, "I have been here: %d times!\n",i);
	//	seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()//reversing here might affect hash filtering result comparing to old version that put hash after this 
	//	seq_reverse(p->len, p->rseq, is_comp);

	//	p->name = strdup((const char*)seq->name.s);
	//	{ // trim /[12]$
	//		int t = strlen(p->name);
	//		if (t > 2 && p->name[t - 2] == '/' && (p->name[t - 1] == '1' || p->name[t - 1] == '2')) p->name[t - 2] = '\0';
	//	}

	//	if (n_seqs == n_needed) break;
	//}
	while (1)
		{
			double rand_num = 0;
			drand48_r(&randBuffer,&rand_num);
			if ( rand_num> frac)
			{
				//notice("before length:%d", seq->seq.l);
				if ((l = kseq_read4_fpc(seq)) < 0) break;
				//notice("after length:%d", seq->seq.l);
				continue;
			}
			if ((l = kseq_read3_fpc(seq)) < 0) break;

			if (seq->seq.l <= l_bc) continue; // sequence length equals or smaller than the barcode length
			p = &seqs[n_seqs++];

			if (l_bc) { // then trim barcode
				for (i = 0; i < l_bc; ++i)
					p->bc[i] = (seq->qual.l && seq->qual.s[i] - 33 < BARCODE_LOW_QUAL) ? tolower(seq->seq.s[i]) : toupper(seq->seq.s[i]);
				p->bc[i] = 0;
				for (; i < seq->seq.l; ++i)
					seq->seq.s[i - l_bc] = seq->seq.s[i];
				seq->seq.l -= l_bc; seq->seq.s[seq->seq.l] = 0;
				if (seq->qual.l) {
					for (i = l_bc; i < seq->qual.l; ++i)
						seq->qual.s[i - l_bc] = seq->qual.s[i];
					seq->qual.l -= l_bc; seq->qual.s[seq->qual.l] = 0;
				}
				l = seq->seq.l;
			}
			else p->bc[0] = 0;
			p->full_len = p->clip_len = p->len = l;
			n_tot += p->full_len;
			//if (p->len>SEQ_INIT_Len)
			{
				delete[] p->seq;
				delete[] p->qual;
				p->seq = new ubyte_t[p->len];
				p->qual = new ubyte_t[p->len];
			}
			for (i = 0; i != p->full_len; ++i)
				p->seq[i] = nst_nt4_table[(int)seq->seq.s[i]];
			if (seq->qual.l) { // copy quality
				if (is_64)
					for (i = 0; i < seq->qual.l; ++i)
					{
						seq->qual.s[i] -= 31;
						p->qual[i] = seq->qual.s[i];
					}
				else
				{
					for (i = 0; i < seq->qual.l; ++i)
					{
						p->qual[i] = seq->qual.s[i];
					}
				}
				if (trim_qual >= 1) n_trimmed += bwa_trim_read(trim_qual, p);
			}

			if (BwtIndex->IsReadFiltered(p->seq, p->qual, p->len))
			{
				p->filtered |= 1;
				if (n_seqs == n_needed) break;
				continue;
			}

			p->rseq = (ubyte_t*)calloc(p->full_len, 1);
			memcpy(p->rseq, p->seq, p->len);
			//fprintf(stderr, "I have been here: %d times!\n",i);
			seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()//reversing here might affect hash filtering result comparing to old version that put hash after this 
			seq_reverse(p->len, p->rseq, is_comp);

			//p->name = strdup((const char*)seq->name.s);
			strncpy(p->name,seq->name.s,seq->name.l);
			{ // trim /[12]$
				int t = strlen(p->name);
				if (t > 2 && p->name[t - 2] == '/' && (p->name[t - 1] == '1' || p->name[t - 1] == '2')) p->name[t - 2] = '\0';
			}

			if (n_seqs == n_needed) break;
		}
	*n = n_seqs;
	if (n_seqs && trim_qual >= 1)
		fprintf(stderr, "[bwa_read_seq] %.1f%% bases are trimmed.\n", 100.0f * n_trimmed / n_tot);
	if (n_seqs == 0) {
		//free(seqs);
		delete[] seqs;
		return 0;
	}
	//ProfilerStop();
	return seqs;
}
static int bwa_read_seq_with_hash_dev(BwtIndexer* BwtIndex, bwa_seqio_t *bs, int n_needed, int *n, int mode, int trim_qual, double frac, uint32_t seed,bwa_seq_t* seqs, int read_len)
{
	
	struct drand48_data randBuffer;
	srand48_r(seed, &randBuffer);
	bwa_seq_t /**seqs,*/ *p;
	kseq_t *seq = bs->ks;
	int n_seqs, l, i, is_comp = mode&BWA_MODE_COMPREAD, is_64 = mode&BWA_MODE_IL13, l_bc = mode >> 24;
	long n_trimmed = 0, n_tot = 0;

	if (l_bc > 15) {
		fprintf(stderr, "[%s] the maximum barcode length is 15.\n", __func__);
		return 0;
	}
	if (bs->is_bam) return bwa_read_bam(bs, n_needed, n, is_comp, trim_qual); // l_bc has no effect for BAM input
	n_seqs = 0;
	//ProfilerStart("FastPopCon.prof");
	while (1)
	{
		double rand_num = 0;
		drand48_r(&randBuffer, &rand_num);
		if (rand_num> frac)
		{
			if ((l = kseq_read4_fpc(seq)) < 0) {
			break;
			}
			continue;
		}
		if ((l = kseq_read3_fpc(seq)) < 0) {
			break;
		}
		if (seq->seq.l <= l_bc) continue; // sequence length equals or smaller than the barcode length
		p = &seqs[n_seqs++];
		if (l_bc) { // then trim barcode
			for (i = 0; i < l_bc; ++i)
				p->bc[i] = (seq->qual.l && seq->qual.s[i] - 33 < BARCODE_LOW_QUAL) ? tolower(seq->seq.s[i]) : toupper(seq->seq.s[i]);
			p->bc[i] = 0;
			for (; i < seq->seq.l; ++i)
				seq->seq.s[i - l_bc] = seq->seq.s[i];
			seq->seq.l -= l_bc; seq->seq.s[seq->seq.l] = 0;
			if (seq->qual.l) {
				for (i = l_bc; i < seq->qual.l; ++i)
					seq->qual.s[i - l_bc] = seq->qual.s[i];
				seq->qual.l -= l_bc; seq->qual.s[seq->qual.l] = 0;
			}
			l = seq->seq.l;
		}
		else p->bc[0] = 0;
		p->full_len = p->clip_len = p->len = l;
		n_tot += p->full_len;
		p->filtered = 0;
		if (p->len>read_len)
		{
			//fprintf(stderr, "the length is weird:%d\np->seq:%x\n%s\n",p->len,p->seq,(char*)p->seq);
			free(p->seq);
			free(p->qual);
			p->seq = (ubyte_t*)calloc(p->len, 1);
			p->qual = (ubyte_t*)calloc(p->len, 1);
			p->rseq = (ubyte_t*)calloc(p->len, 1);
		}
		for (i = 0; i != p->full_len; ++i)
			p->seq[i] = nst_nt4_table[(int)seq->seq.s[i]];
		if (seq->qual.l) { // copy quality
			if (is_64)
				for (i = 0; i < seq->qual.l; ++i)
				{
					seq->qual.s[i] -= 31;
					p->qual[i] = seq->qual.s[i];
				}
			else
			{
				for (i = 0; i < seq->qual.l; ++i)
				{
					p->qual[i] = seq->qual.s[i];
				}
			}
			if (trim_qual >= 1) n_trimmed += bwa_trim_read(trim_qual, p);
		}
		if (BwtIndex->IsReadFiltered(p->seq, p->qual, p->len))
		{
			p->filtered |= 1;
			if (n_seqs == n_needed) break;
			continue;
		}
		//p->rseq = (ubyte_t*)calloc(p->full_len, 1);
		memcpy(p->rseq, p->seq, p->len);
		seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()//reversing here might affect hash filtering result comparing to old version that put hash after this 
		seq_reverse(p->len, p->rseq, is_comp);
		//p->name = strdup((const char*)seq->name.s);
		strncpy(p->name, seq->name.s, seq->name.l);
		{ // trim /[12]$
			int t = strlen(p->name);
			if (t > 2 && p->name[t - 2] == '/' && (p->name[t - 1] == '1' || p->name[t - 1] == '2')) p->name[t - 2] = '\0';
		}

		if (n_seqs == n_needed) break;
	}
	*n = n_seqs;
	if (n_seqs && trim_qual >= 1)
		fprintf(stderr, "[bwa_read_seq] %.1f%% bases are trimmed.\n", 100.0f * n_trimmed / n_tot);
	if (n_seqs == 0) {
		//free(seqs);
		//delete[] seqs;
		return 0;
	}
	//ProfilerStop();
	//return seqs;
	//ProfilerStop();
	return n_seqs;//return success
}
//int BwtMapper::bwa_read_seq_pair_with_hash(BwtIndexer& BwtIndex, bwa_seqio_t *bs, bwa_seqio_t *bs2, int n_needed, int *n, int mode, int trim_qual, double frac, bwa_seq_t * seqs, bwa_seq_t * seqs2)
//{
//	bwa_seq_t /**seqs,*seqs2,*/ *p, *p2;
//	kseq_t *seq = bs->ks;
//	kseq_t *seq2 = bs2->ks;
//	int n_seqs, l, i, is_comp = mode&BWA_MODE_COMPREAD, is_64 = mode&BWA_MODE_IL13, l_bc = mode >> 24;
//	long n_trimmed = 0, n_tot = 0;
//
//	if (l_bc > 15) {
//		fprintf(stderr, "[%s] the maximum barcode length is 15.\n", __func__);
//		return 0;
//	}
//	//if (bs->is_bam) return bwa_read_bam(bs, n_needed, n, is_comp, trim_qual); // l_bc has no effect for BAM input
//	n_seqs = 0;
//	seqs = (bwa_seq_t*)calloc(n_needed, sizeof(bwa_seq_t));
//	seqs2 = (bwa_seq_t*)calloc(n_needed, sizeof(bwa_seq_t));
//	while ((l = kseq_read2(seq, seq2)) >= 0) {
//
//		if (Skip(BwtIndex, mode, seq->seq.s, seq->qual.s, seq2->seq.s, seq2->qual.s, seq2->seq.l, frac)) continue;
//
//		if (is_64 && seq->qual.l)
//		{
//			for (i = 0; i < seq->qual.l; ++i) seq->qual.s[i] -= 31;
//			for (i = 0; i < seq2->qual.l; ++i) seq2->qual.s[i] -= 31;
//		}
//		if (seq->seq.l <= l_bc) continue; // sequence length equals or smaller than the barcode length
//		p2 = &seqs2[n_seqs];
//		p = &seqs[n_seqs++];
//
//		if (l_bc) { // then trim barcode
//			for (i = 0; i < l_bc; ++i)
//			{
//				p->bc[i] = (seq->qual.l && seq->qual.s[i] - 33 < BARCODE_LOW_QUAL) ? tolower(seq->seq.s[i]) : toupper(seq->seq.s[i]);
//				p2->bc[i] = (seq2->qual.l && seq2->qual.s[i] - 33 < BARCODE_LOW_QUAL) ? tolower(seq2->seq.s[i]) : toupper(seq2->seq.s[i]);
//			}
//			p->bc[i] = 0;
//			p2->bc[i] = 0;
//			for (; i < seq->seq.l; ++i)
//			{
//				seq->seq.s[i - l_bc] = seq->seq.s[i];
//				seq2->seq.s[i - l_bc] = seq2->seq.s[i];
//			}
//			seq->seq.l -= l_bc; seq->seq.s[seq->seq.l] = 0;
//			seq2->seq.l -= l_bc; seq2->seq.s[seq->seq.l] = 0;
//			if (seq->qual.l) {
//				for (i = l_bc; i < seq->qual.l; ++i)
//					seq->qual.s[i - l_bc] = seq->qual.s[i];
//				seq->qual.l -= l_bc; seq->qual.s[seq->qual.l] = 0;
//			}
//			if (seq2->qual.l) {
//				for (i = l_bc; i < seq2->qual.l; ++i)
//					seq2->qual.s[i - l_bc] = seq2->qual.s[i];
//				seq2->qual.l -= l_bc; seq2->qual.s[seq2->qual.l] = 0;
//			}
//			l = seq->seq.l;
//		}
//		else
//		{
//			p->bc[0] = 0;
//			p2->bc[0] = 0;
//		}
//		p->tid = -1; // no assigned to a thread
//		p2->tid = -1; // no assigned to a thread
//		p->qual = 0;
//		p2->qual = 0;
//		//p->count=0;
//		p->full_len = p->clip_len = p->len = l;
//		p2->full_len = p2->clip_len = p2->len = l;
//		n_tot += p->full_len;
//		n_tot += p2->full_len;
//
//		p->seq = (ubyte_t*)calloc(p->len, 1);
//		p2->seq = (ubyte_t*)calloc(p2->len, 1);
//		/*
//		p->original_seq =(char*)calloc(p->len, 1);
//		strcpy(p->original_seq,seq->seq.s);
//		*/
//		for (i = 0; i != p->full_len; ++i)
//			p->seq[i] = nst_nt4_table[(int)seq->seq.s[i]];
//		for (i = 0; i != p2->full_len; ++i)
//			p2->seq[i] = nst_nt4_table[(int)seq2->seq.s[i]];
//		if (seq->qual.l) { // copy quality
//			p->qual = (ubyte_t*)strdup((char*)seq->qual.s);
//			if (trim_qual >= 1) n_trimmed += bwa_trim_read(trim_qual, p);
//		}
//		if (seq2->qual.l) { // copy quality
//			p2->qual = (ubyte_t*)strdup((char*)seq2->qual.s);
//			if (trim_qual >= 1) n_trimmed += bwa_trim_read(trim_qual, p);
//		}
//		p->rseq = (ubyte_t*)calloc(p->full_len, 1);
//		memcpy(p->rseq, p->seq, p->len);
//		p2->rseq = (ubyte_t*)calloc(p2->full_len, 1);
//		memcpy(p2->rseq, p2->seq, p2->len);
//		//fprintf(stderr, "I have been here: %d times!\n",i);
//		seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
//		seq_reverse(p->len, p->rseq, is_comp);
//		seq_reverse(p2->len, p2->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
//		seq_reverse(p2->len, p2->rseq, is_comp);
//		p->name = strdup((const char*)seq->name.s);
//		{ // trim /[12]$
//			int t = strlen(p->name);
//			if (t > 2 && p->name[t - 2] == '/' && (p->name[t - 1] == '1' || p->name[t - 1] == '2')) p->name[t - 2] = '\0';
//		}
//		p2->name = strdup((const char*)seq2->name.s);
//		{ // trim /[12]$
//			int t = strlen(p2->name);
//			if (t > 2 && p2->name[t - 2] == '/' && (p2->name[t - 1] == '1' || p2->name[t - 1] == '2')) p2->name[t - 2] = '\0';
//		}
//		if (n_seqs == n_needed) break;
//	}
//	*n = n_seqs;//now it's pairs
//	if (n_seqs && trim_qual >= 1)
//		fprintf(stderr, "[bwa_read_seq] %.1f%% bases are trimmed.\n", 100.0f * n_trimmed / n_tot);
//	if (n_seqs == 0) {
//		free(seqs);
//		free(seqs2)
//			return 0;
//	}
//	return 1;
//}
//specific wrapper for bwa_read_seq[end]
#ifdef HAVE_PTHREAD
static void bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt[2], int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, const BwtIndexer * Indexer);
typedef struct
{
	int tid;
	bwt_t *bwt[2];
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
	const BwtIndexer* Indexer_Ptr;
}
thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	bwa_cal_sa_reg_gap(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt, d->Indexer_Ptr);
	return 0;
}

typedef struct// if((seqs[0] = bwa_read_seq(ks[0], 10000000, &n_seqs, opt->mode, opt->trim_qual))==0) ReadIsGood=1;
{
	BwtIndexer * BwtIndex;
	bwa_seq_t * seqAddress;
	bwa_seqio_t * ksAddress;
	int* n_seqs;
	int mode;
	int trim_qual;
	double frac;
	uint32_t round;
	int* ret;
	int read_len;
}
thread_IO_t;
static void *IOworker(void *data)
{
	thread_IO_t *d = (thread_IO_t*)data;
	//d->seqAddress = bwa_read_seq_with_hash(d->BwtIndex, d->ksAddress, READ_BUFFER_SIZE, d->n_seqs, d->mode, d->trim_qual, d->frac, d->round);
	*(d->ret) = bwa_read_seq_with_hash_dev(d->BwtIndex, d->ksAddress, READ_BUFFER_SIZE, d->n_seqs, d->mode, d->trim_qual, d->frac, d->round, d->seqAddress,d->read_len);
	return 0;
}
#endif
int BwtMapper::bwa_cal_pac_pos_pe(bwt_t * const _bwt[2], const int n_seqs,
	bwa_seq_t *seqs[2], isize_info_t *ii, const pe_opt_t *opt,
	const gap_opt_t *gopt, const isize_info_t *last_ii)
{
	int i, j, cnt_chg = 0;
	//char str[1024];
	bwt_t *bwt[2];
	pe_data_t *d;
	aln_buf_t *buf[2];

	d = (pe_data_t*)calloc(1, sizeof(pe_data_t));
	buf[0] = (aln_buf_t*)calloc(n_seqs, sizeof(aln_buf_t));
	buf[1] = (aln_buf_t*)calloc(n_seqs, sizeof(aln_buf_t));
	/*
	 if (_bwt[0] == 0) { // load forward SA
	 strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
	 strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt[0]);
	 strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
	 strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt[1]);
	 } else
	 */
	bwt[0] = _bwt[0], bwt[1] = _bwt[1];
	DBG(cerr << "Arrived check point 1....\n";)
		// SE
		for (i = 0; i != n_seqs; ++i)
		{ //for each seqs
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j)
		{ // for each end
			uint32_t n_aln;
			p[j] = seqs[j] + i;
			p[j]->n_multi = 0;
			p[j]->extra_flag |= SAM_FPD | (j == 0 ? SAM_FR1 : SAM_FR2);
			if (p[j]->filtered) continue;
			//read(&n_aln, 4, 1, fp_sa[j]);// read in total number of aln
			n_aln = p[j]->n_aln;
			if (n_aln > kv_max(d->aln[j]))
				kv_resize(bwt_aln1_t, d->aln[j], n_aln);
			d->aln[j].n = n_aln; // update total number
			//fread(d->aln[j].a, sizeof(bwt_aln1_t), n_aln, fp_sa[j]);// read in aln of one end
			//d->aln[j].a=p[j]->aln;
			memcpy(d->aln[j].a, p[j]->aln, n_aln * sizeof(bwt_aln1_t));
			kv_copy(bwt_aln1_t, buf[j][i].aln, d->aln[j]); // backup d->aln[j]
			// generate SE alignment and mapping quality
			bwa_aln2seq(n_aln, d->aln[j].a, p[j]);
			if (p[j]->type == BWA_TYPE_UNIQUE || p[j]->type == BWA_TYPE_REPEAT)
			{
				int max_diff =
					gopt->fnr > 0.0 ?
					bwa_cal_maxdiff(p[j]->len, BWA_AVG_ERR,
					gopt->fnr) :
					gopt->max_diff;
				p[j]->pos =
					p[j]->strand ?
					bwt_sa(bwt[0], p[j]->sa) :
					bwt[1]->seq_len
					- (bwt_sa(bwt[1], p[j]->sa) + p[j]->len);
				p[j]->seQ = p[j]->mapQ = bwa_approx_mapQ(p[j], max_diff);
			}
		}
		}
	DBG(cerr << "Arrived check point 2....\n";)
		// infer isize
		infer_isize(n_seqs, seqs, ii, opt->ap_prior, bwt[0]->seq_len);
	if (ii->avg < 0.0 && last_ii->avg > 0.0) *ii = *last_ii;
	if (opt->force_isize) {
		notice("[%s] discard insert size estimate as user's request.\n", __func__);
		ii->low = ii->high = 0;
		ii->avg = ii->std = -1.0;
	}

	// PE
	for (i = 0; i != n_seqs; ++i)
	{
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j)
		{
			p[j] = seqs[j] + i;
			if (p[j]->filtered) continue;
			kv_copy(bwt_aln1_t, d->aln[j], buf[j][i].aln); //copy aln back to d
		}
		if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
			&& (p[1]->type == BWA_TYPE_UNIQUE
			|| p[1]->type == BWA_TYPE_REPEAT))
		{ // only when both ends mapped
			uint64_t x;
			uint32_t j, k, n_occ[2];
			for (j = 0; j < 2; ++j)
			{
				n_occ[j] = 0;
				for (k = 0; k < d->aln[j].n; ++k) // for each aln
					n_occ[j] += d->aln[j].a[k].l - d->aln[j].a[k].k + 1;
			}
			if (n_occ[0] > opt->max_occ || n_occ[1] > opt->max_occ)
				continue; //if any end of the pair exceeded max occ  then process next pair of sequence
			d->arr.n = 0;
			for (j = 0; j < 2; ++j)
			{
				for (k = 0; k < d->aln[j].n; ++k)
				{ // for each alignment
					bwt_aln1_t *r = d->aln[j].a + k;
					bwtint_t l;
					if (r->l - r->k + 1 >= MIN_HASH_WIDTH)
					{ // then check hash table
						uint64_t key = (uint64_t)r->k << 32 | r->l; // key is formed by lower and upper bound
						int ret;
						khint_t iter = kh_put(64, g_hash, key, &ret);
						if (ret)
						{ // if this key is not in the hash table; ret must equal 1 as we never remove elements
							poslist_t *z = &kh_val(g_hash, iter); // return the bwtint_ts pointed by this iter
							z->n = r->l - r->k + 1;
							z->a = (bwtint_t*)malloc(sizeof(bwtint_t) * z->n);
							for (l = r->k; l <= r->l; ++l)
								z->a[l - r->k] =
								r->a ? bwt_sa(bwt[0], l) : bwt[1]->seq_len
								- (bwt_sa(bwt[1], l)
								+ p[j]->len); //call forward / reverse bwt respectively
						}
						for (l = 0; l < kh_val(g_hash, iter).n; ++l)
						{ // ret will surelly show this key in hash, just get its value
							x = kh_val(g_hash, iter).a[l];
							x = x << 32 | k << 1 | j; //packed by  bwtint, lower bound k and pair end j
							kv_push(uint64_t, d->arr, x);
						}
					}
					else
					{ // then calculate on the fly
						for (l = r->k; l <= r->l; ++l)
						{
							x = r->a ?
								bwt_sa(bwt[0], l) :
								bwt[1]->seq_len
								- (bwt_sa(bwt[1], l) + p[j]->len);
							x = x << 32 | k << 1 | j;
							kv_push(uint64_t, d->arr, x);
						}
					}
				}
			}
			cnt_chg += pairing(p, d, opt, gopt->s_mm, ii);
		}
		DBG(cerr << "Arrived check point 3....\n";)
			if (opt->N_multi || opt->n_multi)
			{
			DBG(cerr << "Arrived check point 3-a....\n";)
				for (j = 0; j < 2; ++j)
				{
				DBG(cerr << "Arrived check point 3-0....\n";)
					if (p[j]->type != BWA_TYPE_NO_MATCH)
					{
					int k;
					if (!(p[j]->extra_flag & SAM_FPP)
						&& p[1 - j]->type != BWA_TYPE_NO_MATCH)
					{
						DBG(cerr << "Arrived check point 3-1....\n";)
							bwa_aln2seq_core(d->aln[j].n, d->aln[j].a, p[j], 0,
							p[j]->c1 + p[j]->c2 - 1 > opt->N_multi ?
							opt->n_multi : opt->N_multi);
						DBG(cerr << "Arrived check point 3-2....\n";)
					}
					else
					{
						DBG(cerr << "Arrived check point 3-3....\n";)
							bwa_aln2seq_core(d->aln[j].n, d->aln[j].a, p[j], 0,
							opt->n_multi);
						DBG(cerr << "Arrived check point 3-4....\n";)
					}

					for (k = 0; k < p[j]->n_multi; ++k)
					{
						bwt_multi1_t *q = p[j]->multi + k;
						q->pos =
							q->strand ?
							bwt_sa(bwt[0], q->pos) :
							bwt[1]->seq_len
							- (bwt_sa(bwt[1], q->pos) + p[j]->len);
					}
					}
				}
			}
	}
	DBG(cerr << "Arrived check point 4....\n";)
		// free

		for (i = 0; i < n_seqs; ++i)
		{
		kv_destroy(buf[0][i].aln);
		kv_destroy(buf[1][i].aln);
		}
	free(buf[0]);
	free(buf[1]);
	/*
	if (_bwt[0] == 0) {
	bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	}*/
	kv_destroy(d->arr);
	kv_destroy(d->pos[0]);
	kv_destroy(d->pos[1]);
	kv_destroy(d->aln[0]);
	kv_destroy(d->aln[1]);
	free(d);
	return cnt_chg;
}

static int64_t pos_5(const bwa_seq_t *p)
{
	if (p->type != BWA_TYPE_NO_MATCH)
		return p->strand ? pos_end(p) : p->pos;
	return -1;
}
extern char *bwa_escape(char *s);

int BwtMapper::bwa_set_rg(const char *s)
{
	char *p, *q, *r;
	if (strstr(s, "@RG") != s)
		return -1;
	if (bwa_rg_line)
		free(bwa_rg_line);
	if (bwa_rg_id)
		free(bwa_rg_id);
	bwa_rg_line = strdup(s);
	bwa_rg_id = 0;
	bwa_escape(bwa_rg_line);
	p = strstr(bwa_rg_line, "\tID:");
	if (p == 0)
		return -1;
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q)
		;
	bwa_rg_id = static_cast<char*> (calloc(q - p + 1, 1));
	for (q = p, r = bwa_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	return 0;
}

extern int64_t pos_end_multi(const bwt_multi1_t *p, int len); // analogy to pos_end()
bool BwtMapper::SetSamFileHeader(SamFileHeader& SFH, const bntseq_t * bns)
{

	/*PACKAGE_VERSION*/
	if (!SFH.setPGTag("VN", "1.0.0", "FastqA"))
		std::cerr << "WARNING:SetPGTag failed" << endl;


	if (bwa_rg_line&&strstr(bwa_rg_line, "@RG") == bwa_rg_line)
	{
		stringstream RGline(bwa_rg_line);
		string token, value;
		char tag[3];
		while (RGline >> token)
		{
			tag[0] = token[0];
			tag[1] = token[1];
			tag[2] = 0;
			value = token.substr(3, token.size() - 3);
			if (!SFH.setRGTag(tag, value.c_str(), bwa_rg_id))
				std::cerr << "WARNING:SetRGTag failed" << endl;
		}
	}

	for (int i = 0; i < bns->n_seqs; ++i)
	{

		if (!SFH.setSQTag("LN", std::to_string(bns->anns[i].len).c_str(), bns->anns[i].name))
			std::cerr << "WARNING:SetSQTag failed" << endl;
		//printf("@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
	}

	return 0;
}
bool BwtMapper::SetSamRecord(const bntseq_t *bns, bwa_seq_t *p,
	const bwa_seq_t *mate, int mode, int max_top2, SamFileHeader& SFH,
	SamRecord& SR)
{
	int j;
	int is_64 = mode&BWA_MODE_IL13;
	if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH))
	{
		int seqid, nn, am = 0, flag = p->extra_flag;
		char XT[2];

		if (p->type == BWA_TYPE_NO_MATCH)
		{
			p->pos = mate->pos;
			p->strand = mate->strand;
			flag |= SAM_FSU;
			j = 1;
		}
		else
			j = pos_end(p) - p->pos; // j is the length of the reference in the alignment

		// get seqid
		nn = bns_coor_pac2real(bns, p->pos, j, &seqid);
		if (p->type != BWA_TYPE_NO_MATCH
			&& p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len)
			flag |= SAM_FSU; // flag UNMAP as this alignment bridges two adjacent reference sequences

		// update flag and print it
		if (p->strand)
			flag |= SAM_FSR;
		if (mate)
		{
			if (mate->type != BWA_TYPE_NO_MATCH)
			{
				if (mate->strand)
					flag |= SAM_FMR;
			}
			else
				flag |= SAM_FMU;
		}
		//printf("%s\t%d\t%s\t", p->name, flag, bns->anns[seqid].name);
		//ss<<p->name<<"\t"<<flag<<"\t"<<bns->anns[seqid].name<<"\t";
		SR.setReadName(p->name);
		SR.setFlag(flag);
		SR.setReferenceName(SFH, bns->anns[seqid].name);
		//printf("%d\t%d\t", (int) (p->pos - bns->anns[seqid].offset + 1), p->mapQ);
		//ss<<(int) (p->pos - bns->anns[seqid].offset + 1)<<"\t"<<p->mapQ;
		SR.set1BasedPosition((int)(p->pos - bns->anns[seqid].offset + 1));
		//SR.set0BasedPosition((int) (p->pos - bns->anns[seqid].offset ));
		SR.setMapQuality(p->mapQ);
		// print CIGAR
		ostringstream ss;
		if (p->cigar)
		{
			for (j = 0; j != p->n_cigar; ++j)
				//printf("%d%c", __cigar_len(p->cigar[j]),"MIDS"[__cigar_op(p->cigar[j])]);
				ss << __cigar_len(p->cigar[j]) << "MIDS"[__cigar_op(p->cigar[j])];
		}
		else if (p->type == BWA_TYPE_NO_MATCH)
			//printf("*");
			ss << "*";
		else
			//printf("%dM", p->len);
			ss << p->len << "M";
		SR.setCigar(ss.str().c_str());
		// ss.clear();
		//ss.str("");
		// print mate coordinate
		if (mate && mate->type != BWA_TYPE_NO_MATCH)
		{
			int m_seqid, m_is_N;
			long long isize;
			am = mate->seQ < p->seQ ? mate->seQ : p->seQ; // smaller single-end mapping quality
			// redundant calculation here, but should not matter too much
			m_is_N = bns_coor_pac2real(bns, mate->pos, mate->len, &m_seqid);
			//printf("\t%s\t", (seqid == m_seqid) ? "=" : bns->anns[m_seqid].name);
			// ss << (seqid == m_seqid) ? "=" : bns->anns[m_seqid].name;
			(seqid == m_seqid) ? SR.setMateReferenceName(SFH, "=") : SR.setMateReferenceName(SFH, bns->anns[m_seqid].name);
			isize = (seqid == m_seqid) ? pos_5(mate) - pos_5(p) : 0;
			if (p->type == BWA_TYPE_NO_MATCH)
				isize = 0;
			//printf("%d\t%lld\t", (int) (mate->pos - bns->anns[m_seqid].offset + 1), isize);
			//ss<<(int) (mate->pos - bns->anns[m_seqid].offset + 1)<<"\t"<<isize;
			SR.set1BasedMatePosition((int)(mate->pos - bns->anns[m_seqid].offset + 1));
			SR.setInsertSize(isize);
		}
		else if (mate)
			//printf("\t=\t%d\t0\t", (int) (p->pos - bns->anns[seqid].offset + 1));
			//ss<<"\t=\t"<<(int) (p->pos - bns->anns[seqid].offset + 1)<<"\t0\t";
		{
			SR.setMateReferenceName(SFH, "=");
			SR.set1BasedMatePosition((int)(p->pos - bns->anns[seqid].offset + 1));
			SR.setInsertSize(0);
		}
		else
			//printf("\t*\t0\t0\t");
			//ss<<"\t*0\t0\t";
		{
			SR.setMateReferenceName(SFH, "*");
			SR.set1BasedMatePosition(0);
			SR.setInsertSize(0);
		}
		ss.clear();
		ss.str("");
		// print sequence and quality
		if (p->strand == 0)
		{
			for (j = 0; j != p->full_len; ++j)
				//putchar("ACGTN"[(int) p->seq[j]]);
				ss << "ACGTN"[(int)p->seq[j]];
			SR.setSequence(ss.str().c_str());
		}
		else
		{
			for (j = 0; j != p->full_len; ++j)
				//putchar("TGCAN"[p->seq[p->full_len - 1 - j]]);
				ss << "TGCAN"[p->seq[p->full_len - 1 - j]];
			SR.setSequence(ss.str().c_str());
		}
		//putchar('\t');
		//ss<<"\t";
		ss.clear();
		ss.str("");
		if (p->qual)
		{
			if (is_64)
			{
				for (int i = 0; i < p->len; ++i) p->qual[i] += 31;
			}
			if (p->strand)
				seq_reverse(p->len, p->qual, 0); // reverse quality
			//printf("%s", p->qual);
			//ss << p->qual;
			SR.setQuality((char*)p->qual);
		}
		else
			//printf("*");
			//ss<<"*";
			SR.setQuality("*");

		if (bwa_rg_id)
			//printf("\tRG:Z:%s", bwa_rg_id);
			//ss<<"\tRG:Z:"<<bwa_rg_id;
			SR.addTag("RG", 'Z', bwa_rg_id);
		if (p->bc[0])
			//printf("\tBC:Z:%s", p->bc);
			//ss<<"\tBC:Z:"<<p->bc;
			SR.addTag("BC", 'Z', p->bc);
		if (p->clip_len < p->full_len)
			//printf("\tXC:i:%d", p->clip_len);
			//ss<<"\tXC:i:"<<p->clip_len;
			SR.addTag("XC", 'i', std::to_string(p->clip_len).c_str());
		if (p->type != BWA_TYPE_NO_MATCH)
		{
			int i;
			// calculate XT tag
			XT[0] = "NURM"[p->type];
			if (nn > 10)
				XT[0] = 'N';
			// print tags
			//printf("\tXT:A:%c\t%s:i:%d", XT, (mode & BWA_MODE_COMPREAD) ? "NM" : "CM",p->nm);
			//ss<<"\tXT:A:"<<XT<<"\t";
			XT[1] = '\0';
			SR.addTag("XT", 'A', XT);
			ss.clear();
			ss.str("");
			if (mode & BWA_MODE_COMPREAD)
			{
				ss << "NM";
			}
			else
			{
				ss << "CM";
			}
			//ss<<":i:"<<p->nm;
			SR.addTag(ss.str().c_str(), 'i', to_string(p->nm).c_str());
			if (nn)
				//printf("\tXN:i:%d", nn);
				//ss<<"\tXN:i:"<<nn;
				SR.addTag("XN", 'i', to_string(nn).c_str());
			if (mate)
				//printf("\tSM:i:%d\tAM:i:%d", p->seQ, am);
				//ss<<"\tSM:i:"<<p->seQ<<"\tAM:i:"<<am;
			{
				SR.addTag("SM", 'i', to_string(p->seQ).c_str());
				SR.addTag("AM", 'i', to_string(am).c_str());
			}
			if (p->type != BWA_TYPE_MATESW)
			{ // X0 and X1 are not available for this type of alignment
				//printf("\tX0:i:%d", p->c1);
				//ss<<"\tX0:i:"<<p->c1;
				SR.addTag("X0", 'i', to_string(p->c1).c_str());
				if (p->c1 <= max_top2)
					//printf("\tX1:i:%d", p->c2);
					//ss<<"\tX1:i:"<<p->c2;
					SR.addTag("X1", 'i', to_string(p->c2).c_str());
			}
			//printf("\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo,p->n_gapo + p->n_gape);
			//ss<<"\tXM:i:"<<p->n_mm<<"\tXO:i:"<<p->n_gapo<<"\tXG:i:"<<p->n_gapo + p->n_gape;
			SR.addTag("XM", 'i', to_string(p->n_mm).c_str());
			SR.addTag("XO", 'i', to_string(p->n_gapo).c_str());
			SR.addTag("XG", 'i', to_string(p->n_gapo + p->n_gape).c_str());
			if (p->md)
				//printf("\tMD:Z:%s", p->md);
				//ss<<"\tMD:Z:"<<p->md;
				SR.addTag("MD", 'Z', p->md);
			// print multiple hits
			ss.clear();
			ss.str("");
			if (p->n_multi)
			{
				//printf("\tXA:Z:");
				//ss<<"\tXA:Z:";

				for (i = 0; i < p->n_multi; ++i)
				{
					bwt_multi1_t *q = p->multi + i;
					int k;
					j = pos_end_multi(q, p->len) - q->pos;
					nn = bns_coor_pac2real(bns, q->pos, j, &seqid);
					//printf("%s,%c%d,", bns->anns[seqid].name, q->strand ? '-' : '+',(int) (q->pos - bns->anns[seqid].offset + 1));
					ss << bns->anns[seqid].name << ",";
					if (q->strand)
						ss << '-';
					else
						ss << '+';
					ss << (int)(q->pos - bns->anns[seqid].offset + 1) << ",";
					if (q->cigar)
					{
						for (k = 0; k < q->n_cigar; ++k)
							//printf("%d%c", __cigar_len(q->cigar[k]),"MIDS"[__cigar_op(q->cigar[k])]);
							ss << __cigar_len(q->cigar[k])
							<< "MIDS"[__cigar_op(q->cigar[k])];
					}
					else
						//printf("%dM", p->len);
						ss << p->len << "M";
					//printf(",%d;", q->gap + q->mm);
					ss << "," << q->gap + q->mm << ";";
				}
				SR.addTag("XA", 'Z', ss.str().c_str());
			}
		}
		//putchar('\n');
		//ss<<endl;
	}
	else
	{ // this read has no match
		ostringstream ss;
		ubyte_t *s = p->strand ? p->rseq : p->seq;
		int flag = p->extra_flag | SAM_FSU;
		if (mate && mate->type == BWA_TYPE_NO_MATCH)
			flag |= SAM_FMU;
		//printf("%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", p->name, flag);
		//ss << p->name << "\t" << flag << "\t*\t0\t0\t*\t*\t0\t0\t";
		SR.setReadName(p->name);
		SR.setFlag(flag);
		SR.setReferenceName(SFH, "*");
		SR.set1BasedPosition(0);
		SR.setMapQuality(0);
		SR.setCigar("*");
		SR.setMateReferenceName(SFH, "*");
		SR.set1BasedMatePosition(0);
		SR.setInsertSize(0);
		ss.clear();
		ss.str("");
		for (j = 0; j != p->len; ++j)
			//putchar("ACGTN"[(int) s[j]]);
			ss << "ACGTN"[(int)s[j]];
		//putchar('\t');
		//ss << "\t";
		SR.setSequence(ss.str().c_str());
		ss.clear();
		ss.str("");
		if (p->qual)
		{
			if (p->strand)
				seq_reverse(p->len, p->qual, 0); // reverse quality
			//printf("%s", p->qual);
			//ss << p->qual;
			SR.setQuality((char*)p->qual);
		}
		else
			//printf("*");
			//ss < "*";
			SR.setQuality("*");
		//SR.setQuality(ss.str().c_str());
		if (bwa_rg_id)
			//printf("\tRG:Z:%s", bwa_rg_id);
			//ss << "\tRG:Z:" << bwa_rg_id;
			SR.addTag("RG", 'Z', bwa_rg_id);
		if (p->bc[0])
			//printf("\tBC:Z:%s", p->bc);
			//ss << "\tBC:Z:", p->bc;
			SR.addTag("BC", 'Z', p->bc);
		if (p->clip_len < p->full_len)
			//printf("\tXC:i:%d", p->clip_len);
			//ss << "\tXC:i:" << p->clip_len;
			SR.addTag("XC", 'i', to_string(p->clip_len).c_str());
		//putchar('\n');
		//ss << "\n";
	}
	return 0;
}

/*********bam format*********************************************************/
//ToDo:add in hash filtering in single versoin
bool BwtMapper::SingleEndMapper(BwtIndexer& BwtIndex, const char *fn_fa,
	const gap_opt_t * opt, SamFileHeader& SFH, BamInterface & BamIO, IFILE BamFile, StatGenStatus& StatusTracker, std::ofstream& fout, int &total_add)
{
	int i, n_seqs, tot_seqs = 0; //,m_aln;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	bwt_t *bwt[2];
	bntseq_t *ntbns = 0;
	// initialization
	bwase_initialize();
	for (i = 1; i != 256; ++i)
		g_log_n[i] = (int)(4.343 * log(i) + 0.5);
	srand48(BwtIndex.bns->seed);
	ks = bwa_seq_open(fn_fa);

	bwt[0] = BwtIndex.bwt_d;
	bwt[1] = BwtIndex.rbwt_d;

	ubyte_t *pacseq = 0;
	t = clock();
	FileStatCollector FSC(fn_fa);
	while ((seqs = bwa_read_seq_with_hash(&BwtIndex, ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual, opt->frac,t))
		!= 0)
	{
		tot_seqs += n_seqs;
		FSC.NumRead += n_seqs;
		fprintf(stderr, "NOTICE - Reading in %d sequences into buffer...", n_seqs);
		//      fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
		//      t = clock();
		//t = clock();

#ifdef HAVE_PTHREAD

		if (opt->n_threads <= 1)
		{ // no multi-threading at all
			bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt, &BwtIndex);
		}
		else
		{
			pthread_t *tid;
			pthread_attr_t attr;
			thread_aux_t *data;
			int j;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
			tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
			for (j = 0; j < opt->n_threads; ++j)
			{
				data[j].tid = j;
				data[j].bwt[0] = bwt[0];
				data[j].bwt[1] = bwt[1];
				data[j].n_seqs = n_seqs;
				data[j].seqs = seqs;
				data[j].opt = opt;
				data[j].Indexer_Ptr=&BwtIndex;
				pthread_create(&tid[j], &attr, worker, data + j);
			}
			for (j = 0; j < opt->n_threads; ++j)
				pthread_join(tid[j], 0);
			free(data);
			free(tid);
	}
#else
		//DBG(fprintf(stderr,"Before come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
		bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt, &BwtIndex);
		//DBG(fprintf(stderr,"After come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
#endif

		fprintf(stderr, "NOTICE - Calculate SA coordinate and ");
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();

		for (i = 0; i < n_seqs; ++i)
		{
			bwa_seq_t *p = seqs + i;
			FSC.NumBase += p->full_len;
			bwa_aln2seq_core(p->n_aln, p->aln, p, 1, N_OCC);

		}

		fprintf(stderr, "NOTICE - Convert to sequence coordinate... ");
		bwa_cal_pac_pos(BwtIndex, n_seqs, seqs, opt->max_diff, opt->fnr); // forward bwt will be destroyed here
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();

		fprintf(stderr, "NOTICE - Refine gapped alignments... ");

		//DBG(fprintf(stderr,"Before come into refined gap...%s\n%s\n",seqs->name,seqs->seq);)
		if (pacseq == 0)
			pacseq = bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs, BwtIndex.pac_buf, ntbns);
		else
			pacseq = bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs, pacseq, ntbns);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();

		fprintf(stderr, "NOTICE - Print alignments... ");

		if (!opt->out_bam)
		{
			for (i = 0; i < n_seqs; ++i)
			{
				bwa_seq_t* p = seqs + i;
				if (p->filtered) continue;
				if ((p == 0 || p->type == BWA_TYPE_NO_MATCH)) continue;
				//collector.addAlignment(string(BwtIndex.bns->anns[seqid].name),(seqs+i)->seq,(seqs+i)->qual,(seqs+i)->n_cigar,(seqs+i)->cigar,(seqs+i)->md,(int)((seqs+i)->pos - BwtIndex.bns->anns[seqid].offset + 1),opt);
				//if (
				collector.addAlignment(BwtIndex.bns, seqs + i, 0, opt, fout, total_add);
				//==0)
				//continue; //failed
				bwa_print_sam1(BwtIndex.bns, seqs + i, 0, opt->mode, opt->max_top2);
				//exit(1);
			}
		}
		else
		{
			for (i = 0; i < n_seqs; ++i)
			{
				bwa_seq_t* p = seqs + i;
				if (p->filtered) continue;
				if ((p == 0 || p->type == BWA_TYPE_NO_MATCH)) continue;
				//if (
				collector.addAlignment(BwtIndex.bns, seqs + i, 0, opt, fout, total_add);
				//  ==0)
				//continue; //failed
				SamRecord SR;
				SetSamRecord(BwtIndex.bns, seqs + i, 0, opt->mode, opt->max_top2, SFH, SR);
				//std::cerr<<"\nPassed Read:"<<(seqs+i)->name<<"\t"<<SR.getCigar()<<"\t"<<SR.getSequence()<<"\t"<<SR.getQuality()<<endl;//<<std::for_each((seqs+i),(seqs+i)+(seqs+i)->len-1, [](bwa_seq_t* s, int j ){return  "ACGTN"[(int) *s+j];})<<endl;
				BamIO.writeRecord(BamFile, SFH, SR, SamRecord::SequenceTranslation::NONE);
			}
		}

		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();
		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "NOTICE - %d sequences have been processed.\n", tot_seqs);
		//t = clock();
} //end while
	collector.addFSC(FSC);
	//bam_destroy1(b);
	if (pacseq)
		free(pacseq);
	// destroy
	//bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	bwa_seq_close(ks);
	return 0;
}
bool BwtMapper::PairEndMapper(BwtIndexer& BwtIndex, const char *fn_fa1,
	const char * fn_fa2, const pe_opt_t *popt,
	const gap_opt_t* opt, SamFileHeader& SFH, BamInterface & BamIO, IFILE BamFile, StatGenStatus& StatusTracker, std::ofstream& fout, int &total_add)
{

	int i, j, n_seqs, n_seqs_buff, tot_seqs = 0; //,m_aln;
	bwa_seq_t *seqs[2];
	bwa_seq_t *seqs_buff[2];
	bwa_seqio_t *ks[2];
	clock_t t;
	bwt_t *bwt[2];
	bntseq_t *ntbns = 0;
	khint_t iter;
	isize_info_t last_ii; // this is for the last batch of reads
	// initialization
	bwase_initialize(); // initialize g_log_n[] in bwase.c
	for (i = 1; i != 256; ++i)
		g_log_n[i] = (int)(4.343 * log(i) + 0.5);
	srand48(BwtIndex.bns->seed);
	g_hash = kh_init(64);
	last_ii.avg = -1.0;
	ks[0] = bwa_seq_open(fn_fa1);
	ks[1] = bwa_seq_open(fn_fa2);

	bwt[0] = BwtIndex.bwt_d;
	bwt[1] = BwtIndex.rbwt_d;

	ubyte_t *pacseq = 0;
	FileStatCollector FSC(fn_fa1, fn_fa2);
	int ReadIsGood = 2;
	/*while ((seqs[0] = bwa_read_seq(ks[0], 10000000, &n_seqs, opt->mode,
									opt->trim_qual)) != 0)*/
	uint32_t round = 0;
	//notice("Read 1 length: %d\n", kseq_get_seq_len_fpc(ks[0]->ks));
	//notice("Read 2 length: %d\n", kseq_get_seq_len_fpc(ks[1]->ks));
	while (ReadIsGood)
	{ // opt should be different for two fa files theoretically
		if (ReadIsGood == 2)
		{
			if ((seqs[0] = bwa_read_seq_with_hash(&BwtIndex, ks[0], READ_BUFFER_SIZE, &n_seqs, opt->mode, opt->trim_qual, opt->frac,round)) != 0) ReadIsGood = 1;
			else ReadIsGood = 0;
			FSC.NumRead += n_seqs;
		}
		//FSC.NumRead+=n_seqs;
		int cnt_chg;
		isize_info_t ii;

		t = clock();
		// seqs[1] = bwa_read_seq(ks[1], READ_BUFFER_SIZE, &n_seqs, opt->mode, opt->trim_qual);

		tot_seqs += n_seqs;

		//fprintf(stderr, "NOTICE - Reading in %d sequences into buffer...%fsecs\n", n_seqs, (float)(clock() - t) / CLOCKS_PER_SEC);

		//#pragma omp parallel for 
#ifdef HAVE_PTHREAD

		if (opt->n_threads <= 1)
		{ // no multi-threading at all
			seqs[1] = bwa_read_seq_with_hash(&BwtIndex, ks[1], READ_BUFFER_SIZE, &n_seqs, opt->mode, opt->trim_qual, opt->frac,round);
			FSC.NumRead += n_seqs;
			for (int pair_idx = 0; pair_idx < 2; ++pair_idx)
			{
				bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs[pair_idx], opt, &BwtIndex);
			}
			round++;
			if ((seqs_buff[0] = bwa_read_seq_with_hash(&BwtIndex,ks[0], READ_BUFFER_SIZE, &n_seqs_buff, opt->mode, opt->trim_qual, opt->frac,round)) != 0) 
			{
				ReadIsGood = 1;
				FSC.NumRead+=n_seqs;
			}
			else ReadIsGood = 0;
			FSC.NumRead += n_seqs_buff;
		}
		else
		{
			for (int pair_idx = 0; pair_idx < 2; ++pair_idx)
			{
				pthread_t *tid;
				pthread_attr_t attr;
				thread_aux_t *data;
				int j;
				/*added for IO begin*/
				int n_align_thread=opt->n_threads-1;
				thread_IO_t* IO_param=(thread_IO_t*)calloc(1,sizeof(thread_IO_t));
				/*added for IO end*/
				pthread_attr_init(&attr);
				pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
				data = (thread_aux_t*)calloc(n_align_thread, sizeof(thread_aux_t));
				tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
				/*decopling of the multithread*/
				size_t grain_size=n_seqs/n_align_thread;
				for (j = 0; j < n_align_thread; ++j)
				{
					data[j].tid = j;//+pair_idx*opt->n_threads;
					data[j].bwt[0] = bwt[0];
					data[j].bwt[1] = bwt[1];
					//data[j].n_seqs = n_seqs;
					//data[j].seqs = seqs[pair_idx];
					if(j==n_align_thread-1) data[j].n_seqs = n_seqs-grain_size*(n_align_thread-1);
					else data[j].n_seqs = grain_size;
					data[j].seqs = seqs[pair_idx]+j*grain_size;
					data[j].opt = opt;
					data[j].Indexer_Ptr=&BwtIndex;
					pthread_create(&tid[j], &attr, worker, data + j);
				}
				  {
					  IO_param->BwtIndex = &BwtIndex;
					  IO_param->ksAddress = ks[1-pair_idx];
					  IO_param->n_seqs = &n_seqs_buff;
					  IO_param->mode = opt->mode;
					  IO_param->trim_qual = opt->trim_qual;
					  IO_param->frac=opt->frac;
					  IO_param->round = round;
					  pthread_create(&tid[opt->n_threads-1],&attr,IOworker,IO_param);
				  }

				for (j = 0; j < opt->n_threads; ++j)
					pthread_join(tid[j], 0);

				/*IO thread*/
				if(pair_idx==0)
				{
					seqs[1]=IO_param->seqAddress;
					round++;
				}
				else
				{  
					seqs_buff[0]=IO_param->seqAddress;
				}
				FSC.NumRead += n_seqs;
				if(IO_param->seqAddress==0) ReadIsGood=0;
				free(IO_param);
				free(data);
				free(tid);
			}
	}
#else


		seqs[1] = bwa_read_seq_with_hash(&BwtIndex, ks[1], READ_BUFFER_SIZE, &n_seqs, opt->mode, opt->trim_qual, opt->frac,round);
		FSC.NumRead += n_seqs;
		for (int pair_idx = 0; pair_idx < 2; ++pair_idx)
		{
			//DBG(fprintf(stderr,"Before come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
			bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs[pair_idx], opt, &BwtIndex);
			//DBG(fprintf(stderr,"After come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
		}
		round++;
		if ((seqs[0] = bwa_read_seq_with_hash(&BwtIndex, ks[0], READ_BUFFER_SIZE, &n_seqs, opt->mode, opt->trim_qual, opt->frac,round)) != 0) ReadIsGood = 1;
		else ReadIsGood = 0;
		FSC.NumRead += n_seqs;
#endif

		//}
		fprintf(stderr, "NOTICE - Calculate SA coordinate... ");
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();

		//t = clock();
		//fprintf(stderr, "[bwa_aln_core] write to the disk... ");
		/*for(int iter=0;iter!=2;++iter)
		 {
		 for (i = 0; i < n_seqs; ++i) {
		 bwa_seq_t *p = seqs[iter] + i;
		 //fwrite(&p->n_aln, 4, 1, stdout);
		 //if (p->n_aln) fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, stdout);
		 }
		 }*/

		fprintf(stderr, "NOTICE - convert to sequence coordinate... \n");
		cnt_chg = bwa_cal_pac_pos_pe(bwt, n_seqs, seqs, &ii, popt, opt, &last_ii);
		fprintf(stderr, "NOTICE - time elapses: %.2f sec\n",
			(float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();
		fprintf(stderr, "NOTICE - changing coordinates of %d alignments.\n", cnt_chg);

		fprintf(stderr, "NOTICE - align unmapped mate...\n");
		if (pacseq == 0) //indexing path
			/*pacseq = */
			pacseq = bwa_paired_sw(BwtIndex.bns, BwtIndex.pac_buf, n_seqs, seqs, popt, &ii, opt->mode);
		else
			/*pacseq = */
			pacseq = bwa_paired_sw(BwtIndex.bns, pacseq, n_seqs, seqs, popt, &ii, opt->mode);
		fprintf(stderr, "NOTICE - time elapses: %.2f sec\n",
			(float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();

		fprintf(stderr, "NOTICE - refine gapped alignments... ");
		for (j = 0; j < 2; ++j)
			/*  if (BwtIndex.bns->fp_pac == 0) //indexing path
				bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs[j], BwtIndex.pac_buf,
				ntbns);
				else*/
				bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs[j], pacseq, ntbns);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		t = clock();
		//if (pacseq!= 0) free(pacseq);

		fprintf(stderr, "NOTICE - print alignments... ");
		if (!opt->out_bam)
		{
			for (i = 0; i < n_seqs; ++i)
			{
				bwa_seq_t *p[2];
				p[0] = seqs[0] + i;
				p[1] = seqs[1] + i;

				FSC.NumBase += p[0]->full_len;
				FSC.NumBase += p[1]->full_len;
				//if (p[0]->filtered || p[1]->filtered) continue;
				if ((p[0] == 0 || p[0]->type == BWA_TYPE_NO_MATCH) && (p[1] == 0 || p[1]->type == BWA_TYPE_NO_MATCH)) continue;
				if (p[0]->bc[0] || p[1]->bc[0])
				{
					strcat(p[0]->bc, p[1]->bc);
					strcpy(p[1]->bc, p[0]->bc);
				}

				// if (
				collector.addAlignment(BwtIndex.bns, seqs[0] + i, seqs[1] + i, opt,
					fout, total_add);
				//							  ==0)
				// continue;
				bwa_print_sam1(BwtIndex.bns, p[0], p[1], opt->mode, opt->max_top2);
				bwa_print_sam1(BwtIndex.bns, p[1], p[0], opt->mode, opt->max_top2);
			}
		}
		else
		{
			for (i = 0; i < n_seqs; ++i)
			{
				bwa_seq_t *p[2];
				p[0] = seqs[0] + i;
				p[1] = seqs[1] + i;

				FSC.NumBase += p[0]->full_len;
				FSC.NumBase += p[1]->full_len;
				//if (p[0]->filtered || p[1]->filtered) continue;
				if ((p[0] == 0 || p[0]->type == BWA_TYPE_NO_MATCH) && (p[1] == 0 || p[1]->type == BWA_TYPE_NO_MATCH)) continue;
				if (p[0]->bc[0] || p[1]->bc[0])
				{
					strcat(p[0]->bc, p[1]->bc);
					strcpy(p[1]->bc, p[0]->bc);
				}

				// if (
				collector.addAlignment(BwtIndex.bns, seqs[0] + i, seqs[1] + i, opt,
					fout, total_add);
				// == 0)// means no successfully added reads
				//   if((seqs[0]+i)->type== BWA_TYPE_NO_MATCH)
				//continue;
				SamRecord SR[2];
				SetSamRecord(BwtIndex.bns, seqs[0] + i, seqs[1] + i, opt->mode, opt->max_top2, SFH, SR[0]);
				BamIO.writeRecord(BamFile, SFH, SR[0], SamRecord::SequenceTranslation::NONE);
				SetSamRecord(BwtIndex.bns, seqs[1] + i, seqs[0] + i, opt->mode, opt->max_top2, SFH, SR[1]);
				//std::cerr<<"\nPassed Read:"<<(seqs+i)->name<<"\t"<<SR.getCigar()<<"\t"<<SR.getSequence()<<"\t"<<SR.getQuality()<<endl;//<<std::for_each((seqs+i),(seqs+i)+(seqs+i)->len-1, [](bwa_seq_t* s, int j ){return  "ACGTN"[(int) *s+j];})<<endl;
				BamIO.writeRecord(BamFile, SFH, SR[1], SamRecord::SequenceTranslation::NONE);
			}
		}

		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		//cerr<<"In total "<<total_add<<" reads were calculated!"<<endl;
		t = clock();

		for (j = 0; j < 2; ++j)
		{
			bwa_free_read_seq(n_seqs, seqs[j]);
			delete[] seqs[j];
			n_seqs = n_seqs_buff;
			seqs[j] = seqs_buff[j];
		}
		fprintf(stderr, "NOTICE - %d sequences have been processed.\n", tot_seqs);

		last_ii = ii;
} //end while

	collector.addFSC(FSC);
	if (pacseq)
		free(pacseq);
	//bns_destroy(bns);
	if (ntbns)
		bns_destroy(ntbns);

	for (i = 0; i < 2; ++i)
	{
		bwa_seq_close(ks[i]);
	}
	for (iter = kh_begin(g_hash); iter != kh_end(g_hash); ++iter)
		if (kh_exist(g_hash, iter))
			free(kh_val(g_hash, iter).a);
	kh_destroy(64, g_hash);
	return 0;
	}

	bool BwtMapper::PairEndMapper_dev(BwtIndexer& BwtIndex, const char *fn_fa1,
		const char * fn_fa2, const pe_opt_t *popt,
		const gap_opt_t* opt, SamFileHeader& SFH, BamInterface & BamIO, IFILE BamFile, StatGenStatus& StatusTracker, std::ofstream& fout, int &total_add)
	{

		int i, j, n_seqs, n_seqs_buff, tot_seqs = 0; //,m_aln;
		bwa_seq_t *seqs[2];
		bwa_seq_t *seqs_buff[2];
		bwa_seqio_t *ks[2];
		clock_t t;
		bwt_t *bwt[2];
		bntseq_t *ntbns = 0;
		khint_t iter;
		isize_info_t last_ii; // this is for the last batch of reads
		// initialization
		bwase_initialize(); // initialize g_log_n[] in bwase.c
		for (i = 1; i != 256; ++i)
			g_log_n[i] = (int)(4.343 * log(i) + 0.5);
		srand48(BwtIndex.bns->seed);
		g_hash = kh_init(64);
		last_ii.avg = -1.0;
		ks[0] = bwa_seq_open(fn_fa1);
		ks[1] = bwa_seq_open(fn_fa2);

		bwt[0] = BwtIndex.bwt_d;
		bwt[1] = BwtIndex.rbwt_d;

		ubyte_t *pacseq = 0;
		FileStatCollector FSC(fn_fa1, fn_fa2);
		//try{
			seqs[0] = (bwa_seq_t*) calloc(READ_BUFFER_SIZE,sizeof(bwa_seq_t));
			seqs[1] = (bwa_seq_t*)calloc(READ_BUFFER_SIZE, sizeof(bwa_seq_t));
			seqs_buff[0] = (bwa_seq_t*)calloc(READ_BUFFER_SIZE, sizeof(bwa_seq_t));
			seqs_buff[1] = (bwa_seq_t*)calloc(READ_BUFFER_SIZE, sizeof(bwa_seq_t));
			bwa_init_read_seq(READ_BUFFER_SIZE, seqs[0],opt);
			bwa_init_read_seq(READ_BUFFER_SIZE, seqs[1],opt);
			bwa_init_read_seq(READ_BUFFER_SIZE, seqs_buff[0],opt);
			bwa_init_read_seq(READ_BUFFER_SIZE, seqs_buff[1],opt);
		//}
		//catch (std::bad_alloc& exc)
		//{
		//	warning("Allocating memory failed!\n");
		//	exit(EXIT_FAILURE);
		//}
		int ReadIsGood = 2;
		uint32_t round = 0;
		int ret(-1);
		while (ReadIsGood)
		{ // opt should be different for two fa files theoretically
			if (ReadIsGood == 2)
			{
				if ((ret = bwa_read_seq_with_hash_dev(&BwtIndex, ks[0], READ_BUFFER_SIZE, &n_seqs, opt->mode, opt->trim_qual, opt->frac, round, seqs[0],opt->read_len)) != 0) ReadIsGood = 1;
				else ReadIsGood = 0;
				FSC.NumRead += n_seqs;
			}
			//FSC.NumRead+=n_seqs;
			int cnt_chg;
			isize_info_t ii;

			t = clock();
			// seqs[1] = bwa_read_seq(ks[1], READ_BUFFER_SIZE, &n_seqs, opt->mode, opt->trim_qual);

			tot_seqs += n_seqs;

			//fprintf(stderr, "NOTICE - Reading in %d sequences into buffer...%fsecs\n", n_seqs, (float)(clock() - t) / CLOCKS_PER_SEC);

			//#pragma omp parallel for 
#ifdef HAVE_PTHREAD

			if (opt->n_threads <= 1)
			{ // no multi-threading at all
				ret = bwa_read_seq_with_hash_dev(&BwtIndex, ks[1], READ_BUFFER_SIZE, &n_seqs, opt->mode, opt->trim_qual, opt->frac, round, seqs[1],opt->read_len);
				FSC.NumRead += n_seqs;
				for (int pair_idx = 0; pair_idx < 2; ++pair_idx)
				{
					bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs[pair_idx], opt, &BwtIndex);
				}
				round++;
				if ((ret = bwa_read_seq_with_hash_dev(&BwtIndex, ks[0], READ_BUFFER_SIZE, &n_seqs_buff, opt->mode, opt->trim_qual, opt->frac, round, seqs_buff[0],opt->read_len)) != 0)
				{
					ReadIsGood = 1;
					FSC.NumRead += n_seqs;
				}
				else ReadIsGood = 0;
				FSC.NumRead += n_seqs_buff;
			}
			else
			{
				for (int pair_idx = 0; pair_idx < 2; ++pair_idx)
				{
					pthread_t *tid;
					pthread_attr_t attr;
					thread_aux_t *data;
					int j;
					/*added for IO begin*/
					int n_align_thread = opt->n_threads - 1;
					thread_IO_t* IO_param = (thread_IO_t*)calloc(1, sizeof(thread_IO_t));
					/*added for IO end*/
					pthread_attr_init(&attr);
					pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
					data = (thread_aux_t*)calloc(n_align_thread, sizeof(thread_aux_t));
					tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
					/*decopling of the multithread*/
					size_t grain_size = n_seqs / n_align_thread;
					for (j = 0; j < n_align_thread; ++j)
					{
						data[j].tid = j;//+pair_idx*opt->n_threads;
						data[j].bwt[0] = bwt[0];
						data[j].bwt[1] = bwt[1];
						//data[j].n_seqs = n_seqs;
						//data[j].seqs = seqs[pair_idx];
						if (j == n_align_thread - 1) data[j].n_seqs = n_seqs - grain_size*(n_align_thread - 1);
						else data[j].n_seqs = grain_size;
						data[j].seqs = seqs[pair_idx] + j*grain_size;
						data[j].opt = opt;
						data[j].Indexer_Ptr = &BwtIndex;
						pthread_create(&tid[j], &attr, worker, data + j);
					}
				  {
					  IO_param->BwtIndex = &BwtIndex;
					  IO_param->ksAddress = ks[1 - pair_idx];
					  IO_param->n_seqs = &n_seqs_buff;
					  IO_param->mode = opt->mode;
					  IO_param->trim_qual = opt->trim_qual;
					  IO_param->frac = opt->frac;
					  IO_param->round = round;
					  IO_param->seqAddress = seqs_buff[1 - pair_idx];
					  IO_param->ret = &ret;
					  pthread_create(&tid[opt->n_threads - 1], &attr, IOworker, IO_param);
				  }

				  for (j = 0; j < opt->n_threads; ++j)
					  pthread_join(tid[j], 0);

				  /*IO thread*/
				  if (pair_idx == 0)
				  {
					  bwa_seq_t * tmp = seqs[1];
					  seqs[1] = seqs_buff[1 - pair_idx];
					  seqs_buff[1 - pair_idx] = tmp;		 
				  }
				  else
				  {
					  //seqs_buff[0] = IO_param->seqAddress;
					  round++;
				  }
				  FSC.NumRead += n_seqs;
				  //if (IO_param->seqAddress == 0) ReadIsGood = 0;
				  if (ret == 0) ReadIsGood = 0;
				  free(IO_param);
				  free(data);
				  free(tid);
				}
			}
#else


			ret = bwa_read_seq_with_hash_dev(&BwtIndex, ks[1], READ_BUFFER_SIZE, &n_seqs, opt->mode, opt->trim_qual, opt->frac, round,seqs[1],opt->read_len);
			FSC.NumRead += n_seqs;
			for (int pair_idx = 0; pair_idx < 2; ++pair_idx)
			{
				//DBG(fprintf(stderr,"Before come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
				bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs[pair_idx], opt, &BwtIndex);
				//DBG(fprintf(stderr,"After come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
			}
			round++;
			if ((ret = bwa_read_seq_with_hash_dev(&BwtIndex, ks[0], READ_BUFFER_SIZE, &n_seqs, opt->mode, opt->trim_qual, opt->frac, round,seqs[0],opt->read_len)) != 0) ReadIsGood = 1;
			else ReadIsGood = 0;
			FSC.NumRead += n_seqs;
#endif

			//}
			fprintf(stderr, "NOTICE - Calculate SA coordinate... ");
			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
			t = clock();

			//t = clock();
			//fprintf(stderr, "[bwa_aln_core] write to the disk... ");
			/*for(int iter=0;iter!=2;++iter)
			{
			for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs[iter] + i;
			//fwrite(&p->n_aln, 4, 1, stdout);
			//if (p->n_aln) fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, stdout);
			}
			}*/

			fprintf(stderr, "NOTICE - convert to sequence coordinate... \n");
			cnt_chg = bwa_cal_pac_pos_pe(bwt, n_seqs, seqs, &ii, popt, opt, &last_ii);
			fprintf(stderr, "NOTICE - time elapses: %.2f sec\n",
				(float)(clock() - t) / CLOCKS_PER_SEC);
			t = clock();
			fprintf(stderr, "NOTICE - changing coordinates of %d alignments.\n", cnt_chg);

			fprintf(stderr, "NOTICE - align unmapped mate...\n");
			if (pacseq == 0) //indexing path
				/*pacseq = */
				pacseq = bwa_paired_sw(BwtIndex.bns, BwtIndex.pac_buf, n_seqs, seqs, popt, &ii, opt->mode);
			else
				/*pacseq = */
				pacseq = bwa_paired_sw(BwtIndex.bns, pacseq, n_seqs, seqs, popt, &ii, opt->mode);
			fprintf(stderr, "NOTICE - time elapses: %.2f sec\n",
				(float)(clock() - t) / CLOCKS_PER_SEC);
			t = clock();

			fprintf(stderr, "NOTICE - refine gapped alignments... ");
			for (j = 0; j < 2; ++j)
				/*  if (BwtIndex.bns->fp_pac == 0) //indexing path
				bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs[j], BwtIndex.pac_buf,
				ntbns);
				else*/
				bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs[j], pacseq, ntbns);
			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
			t = clock();
			//if (pacseq!= 0) free(pacseq);

			fprintf(stderr, "NOTICE - print alignments... ");
			if (!opt->out_bam)
			{
				for (i = 0; i < n_seqs; ++i)
				{
					bwa_seq_t *p[2];
					p[0] = seqs[0] + i;
					p[1] = seqs[1] + i;

					FSC.NumBase += p[0]->full_len;
					FSC.NumBase += p[1]->full_len;
					//if (p[0]->filtered || p[1]->filtered) continue;
					if ((p[0] == 0 || p[0]->type == BWA_TYPE_NO_MATCH) && (p[1] == 0 || p[1]->type == BWA_TYPE_NO_MATCH)) continue;
					if (p[0]->bc[0] || p[1]->bc[0])
					{
						strcat(p[0]->bc, p[1]->bc);
						strcpy(p[1]->bc, p[0]->bc);
					}

					// if (
					collector.addAlignment(BwtIndex.bns, seqs[0] + i, seqs[1] + i, opt,
						fout, total_add);
					//							  ==0)
					// continue;
					bwa_print_sam1(BwtIndex.bns, p[0], p[1], opt->mode, opt->max_top2);
					bwa_print_sam1(BwtIndex.bns, p[1], p[0], opt->mode, opt->max_top2);
				}
			}
			else
			{
				for (i = 0; i < n_seqs; ++i)
				{
					bwa_seq_t *p[2];
					p[0] = seqs[0] + i;
					p[1] = seqs[1] + i;

					FSC.NumBase += p[0]->full_len;
					FSC.NumBase += p[1]->full_len;
					//if (p[0]->filtered || p[1]->filtered) continue;
					if ((p[0] == 0 || p[0]->type == BWA_TYPE_NO_MATCH) && (p[1] == 0 || p[1]->type == BWA_TYPE_NO_MATCH)) continue;
					if (p[0]->bc[0] || p[1]->bc[0])
					{
						strcat(p[0]->bc, p[1]->bc);
						strcpy(p[1]->bc, p[0]->bc);
					}

					// if (
					collector.addAlignment(BwtIndex.bns, seqs[0] + i, seqs[1] + i, opt,
						fout, total_add);
					// == 0)// means no successfully added reads
					//   if((seqs[0]+i)->type== BWA_TYPE_NO_MATCH)
					//continue;
					SamRecord SR[2];
					SetSamRecord(BwtIndex.bns, seqs[0] + i, seqs[1] + i, opt->mode, opt->max_top2, SFH, SR[0]);
					BamIO.writeRecord(BamFile, SFH, SR[0], SamRecord::SequenceTranslation::NONE);
					SetSamRecord(BwtIndex.bns, seqs[1] + i, seqs[0] + i, opt->mode, opt->max_top2, SFH, SR[1]);
					//std::cerr<<"\nPassed Read:"<<(seqs+i)->name<<"\t"<<SR.getCigar()<<"\t"<<SR.getSequence()<<"\t"<<SR.getQuality()<<endl;//<<std::for_each((seqs+i),(seqs+i)+(seqs+i)->len-1, [](bwa_seq_t* s, int j ){return  "ACGTN"[(int) *s+j];})<<endl;
					BamIO.writeRecord(BamFile, SFH, SR[1], SamRecord::SequenceTranslation::NONE);
				}
			}

			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
			//cerr<<"In total "<<total_add<<" reads were calculated!"<<endl;
			t = clock();

			for (j = 0; j < 2; ++j)
			{
				//bwa_free_read_seq(n_seqs, seqs[j]);
				//delete[] seqs[j];
				n_seqs = n_seqs_buff;
				bwa_seq_t * tmp = seqs[j];
				seqs[j] = seqs_buff[j];
				seqs_buff[j] = tmp;
			}
			fprintf(stderr, "NOTICE - %d sequences have been processed.\n", tot_seqs);

			last_ii = ii;
		} //end while

		for (j = 0; j < 2; ++j)
		{
			bwa_free_read_seq(READ_BUFFER_SIZE, seqs[j]);
			bwa_free_read_seq(READ_BUFFER_SIZE, seqs_buff[j]);
		}
		free(seqs[0]);
		free(seqs[1]);
		free(seqs_buff[0]);
		free(seqs_buff[1]);
		collector.addFSC(FSC);
		if (pacseq)
			free(pacseq);
		//bns_destroy(bns);
		if (ntbns)
			bns_destroy(ntbns);

		for (i = 0; i < 2; ++i)
		{
			bwa_seq_close(ks[i]);
		}
		for (iter = kh_begin(g_hash); iter != kh_end(g_hash); ++iter)
			if (kh_exist(g_hash, iter))
				free(kh_val(g_hash, iter).a);
		kh_destroy(64, g_hash);
		return 0;
	}

BwtMapper::BwtMapper(BwtIndexer& BwtIndex, const string & Fastq_1,
	const string & Fastq_2, const string & Prefix, const string & RefPath, const pe_opt_t* popt,
	const gap_opt_t * opt)
{
	//std::cerr<<"Open Fastq  ... "<<endl;
	bwa_set_rg(opt->RG);
	if (opt->in_bam != 0)
	{
		notice(" Input alignments from Bam file...\n");
		/*
		SamFileHeader SFH;
		SamFile SFIO;*/
		notice("Restore Variant Site Info...\n");
		collector.restoreVcfSites(RefPath, opt); collector.getGenomeSize(BwtIndex.RefPath);
		ofstream fout(Prefix + ".InsertSizeTable");
		int total_add = 0;
		collector.ReadAlignmentFromBam(opt, /*SFH, SFIO,*/ opt->in_bam, fout, total_add);
		notice("%d reads were calculated...\n", total_add);
		fout.close();
		// BamFile->ifclose();
		// destroy
		notice("Calculate distributions...\n ");
		collector.processCore(Prefix, opt);
	}
	else if (Fastq_2 != "Empty")
	{
		notice("Using Pair End mapping...\n");
		SamFileHeader SFH;
		BamInterface BamIO;
		IFILE BamFile = new InputFile((Prefix + ".bam").c_str(), "w", InputFile::ifileCompression::BGZF);
		StatGenStatus StatusTracker;
		StatusTracker.setStatus(StatGenStatus::Status::SUCCESS, "Initialization when start.\n");
		/*********bam header init end***********/

		if (!opt->out_bam)
		{
			bwa_print_sam_SQ(BwtIndex.bns);
			bwa_print_sam_PG();
		}
		else
		{
			if (!BamFile->isOpen())
			{
				warning("Open Bam file for writing failed, abort!\n");
				exit(1);
			}
			SetSamFileHeader(SFH, BwtIndex.bns);
			BamIO.writeHeader(BamFile, SFH, StatusTracker);
		}
		notice("Restore Variant Site Info...\n");
		collector.restoreVcfSites(RefPath, opt);
		ofstream fout(Prefix + ".InsertSizeTable");
		int total_add = 0;
		PairEndMapper(BwtIndex, Fastq_1.c_str(), Fastq_2.c_str(), popt, opt, SFH, BamIO, BamFile, StatusTracker, fout, total_add);
		notice(" %d reads were calculated...\n", total_add);
		fout.close();
		BamFile->ifclose();
		// destroy
		notice("Calculate distributions... \n");
		collector.processCore(Prefix, opt);

	}
	else
	{
		notice("Using Single End mapping...\n");
		SamFileHeader SFH;
		BamInterface BamIO;
		IFILE BamFile = new InputFile((Prefix + ".bam").c_str(), "w", InputFile::ifileCompression::BGZF);
		StatGenStatus StatusTracker;
		StatusTracker.setStatus(StatGenStatus::Status::SUCCESS, "Initialization when start.\n");
		/*********bam header init end***********/

		if (!opt->out_bam)
		{
			bwa_print_sam_SQ(BwtIndex.bns);
			bwa_print_sam_PG();
		}
		else
		{
			if (!BamFile->isOpen())
			{
				warning("Open Bam file for writing failed, abort!\n");
				exit(1);
			}
			SetSamFileHeader(SFH, BwtIndex.bns);
			BamIO.writeHeader(BamFile, SFH, StatusTracker);
		}
		notice("Restore Variant Site Info...\n");
		collector.restoreVcfSites(RefPath, opt); collector.getGenomeSize(BwtIndex.RefPath);
		ofstream fout(Prefix + ".InsertSizeTable");
		int total_add = 0;
		SingleEndMapper(BwtIndex, Fastq_1.c_str(), opt, SFH, BamIO, BamFile, StatusTracker, fout, total_add);
		fout.close();
		BamFile->ifclose();
		delete BamFile;
		// destroy
		notice("Calculate distributions...\n ");
		collector.processCore(Prefix, opt);
	}
}
BwtMapper::BwtMapper(BwtIndexer& BwtIndex, const string & FaList,
	const string & Prefix, const string & RefPath, const pe_opt_t* popt,
	const gap_opt_t * opt)
{
	bwa_set_rg(opt->RG);
	if (opt->in_bam != 0)
	{
		notice("Input alignments from Bam file...\n");
		/*
		SamFileHeader SFH;
		SamFile SFIO;*/
		notice("Restore Variant Site Info...\n");
		collector.restoreVcfSites(RefPath, opt); collector.getGenomeSize(BwtIndex.RefPath);
		ofstream fout(Prefix + ".InsertSizeTable");
		int total_add = 0;
		collector.ReadAlignmentFromBam(opt, /*SFH, SFIO,*/ opt->in_bam, fout, total_add);
		notice("%d reads were calculated!\n", total_add);
		fout.close();
		// BamFile->ifclose();
		notice("Calculate distributions...\n");
		collector.processCore(Prefix, opt);
	}
	else
	{
		notice("Open Fastq List ...\n");
		SamFileHeader SFH;
		BamInterface BamIO;
		IFILE BamFile = new InputFile((Prefix + ".bam").c_str(), "w", InputFile::ifileCompression::BGZF);
		StatGenStatus StatusTracker;
		StatusTracker.setStatus(StatGenStatus::Status::SUCCESS, "Initialization when start.\n");
		if (!opt->out_bam)
		{
			bwa_print_sam_SQ(BwtIndex.bns);
			bwa_print_sam_PG();
		}
		else
		{
			if (!BamFile->isOpen())
			{
				warning("Open Bam file for writing failed, abort!\n");
				exit(EXIT_FAILURE);
			}
			SetSamFileHeader(SFH, BwtIndex.bns);
			BamIO.writeHeader(BamFile, SFH, StatusTracker);
		}
		double t_tmp = realtime();

		collector.restoreVcfSites(RefPath, opt); collector.getGenomeSize(BwtIndex.RefPath);
		notice("Restore Variant Site Info...%f sec\n", realtime() - t_tmp);
		ofstream fout(Prefix + ".InsertSizeTable");
		int total_add = 0;
		ifstream fin(FaList);
		if (!fin.is_open()) error("Open file %s failed", FaList.c_str());
		std::string line;
		int i(0);
		while (getline(fin, line))
		{
			++i;
			std::string Fastq_1(""), Fastq_2("");
			stringstream ss(line);
			ss >> Fastq_1 >> Fastq_2;
			if (Fastq_2 != "")
			{
				notice("Processing Pair End mapping\t%d\n%s\n%s\n", i, Fastq_1.c_str(), Fastq_2.c_str());
				t_tmp = realtime();
				/*PairEndMapper*/PairEndMapper_dev(BwtIndex, Fastq_1.c_str(), Fastq_2.c_str(), popt, opt, SFH, BamIO, BamFile, StatusTracker, fout, total_add);
				notice("Processed Pair End mapping in %f sec\n", realtime() - t_tmp);
			}
			else
			{
				notice("Processing Single End mapping\t%d\n%s\n", i, Fastq_1.c_str());
				SingleEndMapper(BwtIndex, Fastq_1.c_str(), opt, SFH, BamIO, BamFile, StatusTracker, fout, total_add);
			}
		}
		notice("%d reads were calculated!", total_add);
		fout.close();
		BamFile->ifclose();
		delete BamFile;
		// destroy
		t_tmp = realtime();
		collector.processCore(Prefix, opt);
		notice("Calculate distributions... %f sec\n",realtime()-t_tmp);
	}
}

BwtMapper::~BwtMapper()
{
	// TODO Auto-generated destructor stub
}

