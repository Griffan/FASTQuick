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
#include "./libbwa/bwtgap.h"
#include "./libbwa/bwase.h"
#include "./libbwa/bwape.h"
#include "./libbwa/khash.h"
KHASH_MAP_INIT_INT64(64, poslist_t)

using namespace std;
extern BwtIndexer Indexer;


kh_64_t *g_hash;
#define DEBUG 0
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

#ifdef HAVE_PTHREAD
static void bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt[2], int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);
typedef struct {
	int tid;
	bwt_t *bwt[2];
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
} thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	bwa_cal_sa_reg_gap(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt);
	return 0;
}
#endif


extern void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq, bntseq_t *ntbns);
extern void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2);
extern void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);
extern int bwa_approx_mapQ(const bwa_seq_t *p, int mm);
extern int bwt_cal_width(const bwt_t *rbwt, int len, const ubyte_t *str, bwt_width_t *width);

static void bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt[2], int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt)
{
	int i, max_l = 0, max_len;
	gap_stack_t *stack;
	bwt_width_t *w[2], *seed_w[2];
	const ubyte_t *seq[2];
	gap_opt_t local_opt = *opt;

	// initiate priority stack
	for (i = max_len = 0; i != n_seqs; ++i)
		if (seqs[i].len > max_len) max_len = seqs[i].len;
	if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;
	stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);

	seed_w[0] = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	seed_w[1] = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	w[0] = w[1] = 0;
	uint32_t unmapped_num=0;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;

		if(Indexer.IsReadFiltered(p->seq, p->qual, p->len))
		{
			/*if(strcmp(p->name,"ERR013170.1716") ==0)
				{
				fprintf(stdout,"ERR013170.1716	is filtered\n" );
				exit(1);
				}*/
			p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = (bwt_aln1_t*)calloc(4, sizeof(bwt_aln1_t));
			unmapped_num++;
			continue;
		}


#ifdef HAVE_PTHREAD
		if (opt->n_threads > 1) {
			pthread_mutex_lock(&g_seq_lock);
			if (p->tid < 0) { // unassigned
				int j;
				for (j = i; j < n_seqs && j < i + THREAD_BLOCK_SIZE; ++j)
					seqs[j].tid = tid;
			} else if (p->tid != tid) {
				pthread_mutex_unlock(&g_seq_lock);
				continue;
			}
			pthread_mutex_unlock(&g_seq_lock);
		}
#endif
		p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
		seq[0] = p->seq; seq[1] = p->rseq;
		if (max_l < p->len) {
			max_l = p->len;
			w[0] = (bwt_width_t*)realloc(w[0], (max_l + 1) * sizeof(bwt_width_t));
			w[1] = (bwt_width_t*)realloc(w[1], (max_l + 1) * sizeof(bwt_width_t));
			memset(w[0], 0, (max_l + 1) * sizeof(bwt_width_t));
			memset(w[1], 0, (max_l + 1) * sizeof(bwt_width_t));
		}
		bwt_cal_width(bwt[0], p->len, seq[0], w[0]);
		bwt_cal_width(bwt[1], p->len, seq[1], w[1]);
		if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		local_opt.seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;
		if (p->len > opt->seed_len) {
			bwt_cal_width(bwt[0], opt->seed_len, seq[0] + (p->len - opt->seed_len), seed_w[0]);
			bwt_cal_width(bwt[1], opt->seed_len, seq[1] + (p->len - opt->seed_len), seed_w[1]);
		}
		// core function

		p->aln = bwt_match_gap(bwt, p->len, seq, w, p->len <= opt->seed_len? 0 : seed_w, &local_opt, &p->n_aln, stack);
		/*if(strcmp(p->name,"ERR013170.1716") ==0)
		{
			/*fwrite(w[0], sizeof(bwt_width_t), max_l, stdout);
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
	fprintf(stderr,"RollingHash filtered %d reads...",unmapped_num);
	free(seed_w[0]); free(seed_w[1]);
	free(w[0]); free(w[1]);
	gap_destroy_stack(stack);
}

int BwtMapper::bwa_cal_pac_pos(BwtIndexer& BwtIndex, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr)
{
	int i,j;
		//char str[1024];
		bwt_t *bwt;
		// load forward SA
		//strcpy(str, prefix); strcat(str, ".bwt");
		bwt =BwtIndex.bwt_d; //bwt_restore_bwt(str);
		//strcpy(str, prefix); strcat(str, ".sa");
		//bwt_restore_sa(str, bwt);
		for (i = 0; i != n_seqs; ++i) {
				if (seqs[i].strand) bwa_cal_pac_pos_core(bwt, 0, &seqs[i], max_mm, fnr);
				for (j = 0; j < seqs[i].n_multi; ++j) {
					bwt_multi1_t *p = seqs[i].multi + j;
					if (p->strand) p->pos = bwt_sa(bwt, p->pos);//transform pos into actual position
				}
			}
		//bwt_destroy(bwt);
		// load reverse BWT and SA
		//strcpy(str, prefix); strcat(str, ".rbwt");
		bwt = BwtIndex.rbwt_d; //bwt_restore_bwt(str);
		//strcpy(str, prefix); strcat(str, ".rsa");
		//bwt_restore_sa(str, bwt);
		for (i = 0; i != n_seqs; ++i) {
				if (!seqs[i].strand) bwa_cal_pac_pos_core(0, bwt, &seqs[i], max_mm, fnr);
				for (j = 0; j < seqs[i].n_multi; ++j) {
					bwt_multi1_t *p = seqs[i].multi + j;
					if (!p->strand) p->pos = bwt->seq_len - (bwt_sa(bwt, p->pos) + seqs[i].len);
				}
			}
		//bwt_destroy(bwt);
		return 1;
}

typedef struct {
	kvec_t(bwt_aln1_t) aln;
} aln_buf_t;
int BwtMapper::bwa_cal_pac_pos_pe(bwt_t *const _bwt[2], int n_seqs, bwa_seq_t *seqs[2],  isize_info_t *ii,
					   const pe_opt_t *opt, const gap_opt_t *gopt, const isize_info_t *last_ii)
{
	int i, j, cnt_chg = 0;
	char str[1024];
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
		DBG(cerr<<"Arrived check point 1....\n";)
	// SE
	for (i = 0; i != n_seqs; ++i) {//for each seqs
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {// for each end
			int n_aln;
			p[j] = seqs[j] + i;
			p[j]->n_multi = 0;
			p[j]->extra_flag |= SAM_FPD | (j == 0? SAM_FR1 : SAM_FR2);
			//read(&n_aln, 4, 1, fp_sa[j]);// read in total number of aln
			n_aln=p[j]->n_aln;
			if (n_aln > kv_max(d->aln[j]))
				kv_resize(bwt_aln1_t, d->aln[j], n_aln);
			d->aln[j].n = n_aln;// update total number
			//fread(d->aln[j].a, sizeof(bwt_aln1_t), n_aln, fp_sa[j]);// read in aln of one end
			//d->aln[j].a=p[j]->aln;
			memcpy(d->aln[j].a, p[j]->aln, n_aln*sizeof(bwt_aln1_t));
			kv_copy(bwt_aln1_t, buf[j][i].aln, d->aln[j]); // backup d->aln[j]
			// generate SE alignment and mapping quality
			bwa_aln2seq(n_aln, d->aln[j].a, p[j]);
			if (p[j]->type == BWA_TYPE_UNIQUE || p[j]->type == BWA_TYPE_REPEAT) {
				int max_diff = gopt->fnr > 0.0? bwa_cal_maxdiff(p[j]->len, BWA_AVG_ERR, gopt->fnr) : gopt->max_diff;
				p[j]->pos = p[j]->strand? bwt_sa(bwt[0], p[j]->sa)
					: bwt[1]->seq_len - (bwt_sa(bwt[1], p[j]->sa) + p[j]->len);
				p[j]->seQ = p[j]->mapQ = bwa_approx_mapQ(p[j], max_diff);
			}
		}
	}
		DBG(cerr<<"Arrived check point 2....\n";)
	// infer isize
	//infer_isize(n_seqs, seqs, ii, opt->ap_prior, bwt[0]->seq_len);
	//if (ii->avg < 0.0 && last_ii->avg > 0.0) *ii = *last_ii;
	//if (opt->force_isize) {
		//fprintf(stderr, "[%s] discard insert size estimate as user's request.\n", __func__);
		ii->low = ii->high = 0; ii->avg = ii->std = -1.0;
	//}

	// PE
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {
			p[j] = seqs[j] + i;
			kv_copy(bwt_aln1_t, d->aln[j], buf[j][i].aln);//copy aln back to d
		}
		if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
			&& (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT))
		{ // only when both ends mapped
			uint64_t x;
			int j, k, n_occ[2];
			for (j = 0; j < 2; ++j) {
				n_occ[j] = 0;
				for (k = 0; k < d->aln[j].n; ++k)// for each aln
					n_occ[j] += d->aln[j].a[k].l - d->aln[j].a[k].k + 1;
			}
			if (n_occ[0] > opt->max_occ || n_occ[1] > opt->max_occ) continue;//if any end of the pair exceeded max occ  then process next pair of sequence
			d->arr.n = 0;
			for (j = 0; j < 2; ++j) {
				for (k = 0; k < d->aln[j].n; ++k) {// for each alignment
					bwt_aln1_t *r = d->aln[j].a + k;
					bwtint_t l;
					if (r->l - r->k + 1 >= MIN_HASH_WIDTH) { // then check hash table
						uint64_t key = (uint64_t)r->k<<32 | r->l;// key is formed by lower and upper bound
						int ret;
						khint_t iter = kh_put(64, g_hash, key, &ret);
						if (ret) { // if this key is not in the hash table; ret must equal 1 as we never remove elements
							poslist_t *z = &kh_val(g_hash, iter);// return the bwtint_ts pointed by this iter
							z->n = r->l - r->k + 1;
							z->a = (bwtint_t*)malloc(sizeof(bwtint_t) * z->n);
							for (l = r->k; l <= r->l; ++l)
								z->a[l - r->k] = r->a? bwt_sa(bwt[0], l) : bwt[1]->seq_len - (bwt_sa(bwt[1], l) + p[j]->len);//call forward / reverse bwt respectively
						}
						for (l = 0; l < kh_val(g_hash, iter).n; ++l) {// ret will surelly show this key in hash, just get its value
							x = kh_val(g_hash, iter).a[l];
							x = x<<32 | k<<1 | j;//packed by  bwtint, lower bound k and pair end j
							kv_push(uint64_t, d->arr, x);
						}
					} else { // then calculate on the fly
						for (l = r->k; l <= r->l; ++l) {
							x = r->a? bwt_sa(bwt[0], l) : bwt[1]->seq_len - (bwt_sa(bwt[1], l) + p[j]->len);
							x = x<<32 | k<<1 | j;
							kv_push(uint64_t, d->arr, x);
						}
					}
				}
			}
			cnt_chg += pairing(p, d, opt, gopt->s_mm, ii);
		}
		DBG(cerr<<"Arrived check point 3....\n";)
		if (opt->N_multi || opt->n_multi) {
			DBG(cerr<<"Arrived check point 3-a....\n";)
			for (j = 0; j < 2; ++j) {
				DBG(cerr<<"Arrived check point 3-0....\n";)
				if (p[j]->type != BWA_TYPE_NO_MATCH) {
					int k;
					if (!(p[j]->extra_flag&SAM_FPP) && p[1-j]->type != BWA_TYPE_NO_MATCH) {
						DBG(cerr<<"Arrived check point 3-1....\n";)
						bwa_aln2seq_core(d->aln[j].n, d->aln[j].a, p[j], 0, p[j]->c1+p[j]->c2-1 > opt->N_multi? opt->n_multi : opt->N_multi);
						DBG(cerr<<"Arrived check point 3-2....\n";)
					} else
						{
						DBG(cerr<<"Arrived check point 3-3....\n";)
						bwa_aln2seq_core(d->aln[j].n, d->aln[j].a, p[j], 0, opt->n_multi);
						DBG(cerr<<"Arrived check point 3-4....\n";)
						}

					for (k = 0; k < p[j]->n_multi; ++k) {
						bwt_multi1_t *q = p[j]->multi + k;
						q->pos = q->strand? bwt_sa(bwt[0], q->pos) : bwt[1]->seq_len - (bwt_sa(bwt[1], q->pos) + p[j]->len);
					}
				}
			}
		}
	}
	DBG(cerr<<"Arrived check point 4....\n";)
	// free
	/*
	for (i = 0; i < n_seqs; ++i) {
		kv_destroy(buf[0][i].aln);
		kv_destroy(buf[1][i].aln);
	}
	free(buf[0]); free(buf[1]);
	if (_bwt[0] == 0) {
		bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	}*/
	kv_destroy(d->arr);
	kv_destroy(d->pos[0]); kv_destroy(d->pos[1]);
	kv_destroy(d->aln[0]); kv_destroy(d->aln[1]);
	free(d);
	return cnt_chg;
}

BwtMapper::BwtMapper()
{

}



bool BwtMapper::SingleEndMapper(BwtIndexer& BwtIndex,const char *fn_fa, const gap_opt_t * opt)
{
	//void bwa_aln_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt)
//	{
		int i, n_seqs, tot_seqs = 0;//,m_aln;
		bwa_seq_t *seqs;
		bwa_seqio_t *ks;
		clock_t t;
		bwt_t *bwt[2];
		bntseq_t  *ntbns = 0;
		// initialization
		bwase_initialize();
		for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);
		srand48(BwtIndex.bns->seed);
		ks = bwa_seq_open(fn_fa);

		/*{ // load BWT
			char *str = (char*)calloc(strlen(prefix) + 10, 1);
			strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
			strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
			free(str);
		}*/
		bwt[0]=BwtIndex.bwt_d;
		bwt[1]=BwtIndex.rbwt_d;

		//bwt_dump_bwt("DEBUG.bwt", bwt[0]);

		//cerr<<"The Mapper hs rbwt_d->bwt_size is:"<<rbwt_d->bwt_size<<endl;
		//bwt_dump_bwt("DEBUG.rbwt", bwt[1]);
		//DBG(fprintf(stderr,"checkpoint 1....bwt->seq_len:%d\n",bwt[0]->seq_len);)
		// core loop
		//m_aln=0;
		//fwrite(opt, sizeof(gap_opt_t), 1, stdout);
		//revealopt(opt);
		//if (!(opt.mode & BWA_MODE_COMPREAD)) // in color space; initialize ntpac
		//bntseq_t* ntbns = bwa_open_nt(prefix)
		bwa_print_sam_SQ(BwtIndex.bns);
		bwa_print_sam_PG();
		//exit(1);
		t = clock();
		while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode , opt->trim_qual)) != 0) {
			tot_seqs += n_seqs;

			fprintf(stderr,"Reading in %d sequences into buffer...",n_seqs);
			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
			t = clock();
			fprintf(stderr, "Calculate SA coordinate... ");



	#ifdef HAVE_PTHREAD
			if (opt->n_threads <= 1) { // no multi-threading at all
				bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
			} else {
				pthread_t *tid;
				pthread_attr_t attr;
				thread_aux_t *data;
				int j;
				pthread_attr_init(&attr);
				pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
				data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
				tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
				for (j = 0; j < opt->n_threads; ++j) {
					data[j].tid = j; data[j].bwt[0] = bwt[0]; data[j].bwt[1] = bwt[1];
					data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
					pthread_create(&tid[j], &attr, worker, data + j);
				}
				for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
				free(data); free(tid);
			}
	#else
			//DBG(fprintf(stderr,"Before come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
			bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
			//DBG(fprintf(stderr,"After come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
	#endif

			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

			t = clock();
			//fprintf(stderr, "[bwa_aln_core] write to the disk... ");
			for (i = 0; i < n_seqs; ++i) {
				bwa_seq_t *p = seqs + i;
				//fwrite(&p->n_aln, 4, 1, stdout);
				//if (p->n_aln) fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, stdout);
				//bwa_aln2seq(p->n_aln, p->aln, p);
				bwa_aln2seq_core(p->n_aln, p->aln, p, 1, N_OCC);

			}

		/*	for (i = 0; i < n_seqs; ++i) {
					bwa_seq_t *p = seqs + i;
					fwrite(&p->n_multi, 4, 1, stdout);
					if (p->n_multi) fwrite(p->multi, sizeof(bwt_multi1_t), p->n_multi, stdout);
					//fwrite(&p->sa, sizeof(bwtint_t), 1, stdout);
					//int c1=p->c1;
					//int c2=p->c2;
					//fwrite(&c1, sizeof(int), 1, stdout);
					//fwrite(&c2, sizeof(int), 1, stdout);
					/*fprintf(stdout,"%s\tn_aln:%d\t",p->name,p->n_aln);
					if (p->n_aln)
						{
						int iter=0;
						for(;iter<p->n_aln;++iter)
						fprintf(stdout,"aln.k:%d-aln.l:%d\t",p->aln[iter].k, p->aln[iter].l);
						fprintf(stdout,"\n");
						}
					else
						fprintf(stdout,"\n");





		}
		exit(1);
		{*/
			fprintf(stderr, "Convert to sequence coordinate... ");
			bwa_cal_pac_pos(BwtIndex, n_seqs, seqs, opt->max_diff, opt->fnr); // forward bwt will be destroyed here
			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

/*
			for (i = 0; i < n_seqs; ++i) {
					bwa_seq_t *p = seqs + i;
					fwrite(&p->n_multi, 4, 1, stdout);
					if (p->n_multi) fwrite(p->multi, sizeof(bwt_multi1_t), p->n_multi, stdout);
			}
*/

			fprintf(stderr, "Refine gapped alignments... ");

			//DBG(fprintf(stderr,"Before come into refined gap...%s\n%s\n",seqs->name,seqs->seq);)
			if(BwtIndex.bns->fp_pac==0)
			bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs, BwtIndex.pac_buf, ntbns);
			else
				bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs, 0, ntbns);
			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();



			fprintf(stderr, "Print alignments... ");
			for (i = 0; i < n_seqs; ++i)
				bwa_print_sam1(BwtIndex.bns, seqs + i, 0, opt->mode, opt->max_top2);

			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
			bwa_free_read_seq(n_seqs, seqs);
			fprintf(stderr, " %d sequences have been processed.\n", tot_seqs);
			t = clock();
		}//end while

		// destroy
		//bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
		bwa_seq_close(ks);
}
bool BwtMapper::PairEndMapper(BwtIndexer& BwtIndex,const char *fn_fa1,const char * fn_fa2,const  pe_opt_t *popt,const gap_opt_t* opt)
{

		int i,j, n_seqs, tot_seqs = 0;//,m_aln;
		bwa_seq_t *seqs[2];
		bwa_seqio_t *ks[2];
		clock_t t;
		bwt_t *bwt[2];
		bntseq_t  *ntbns = 0;
		khint_t iter;
		isize_info_t last_ii; // this is for the last batch of reads
		// initialization
		bwase_initialize(); // initialize g_log_n[] in bwase.c
		for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);
		srand48(BwtIndex.bns->seed);
		g_hash = kh_init(64);
		last_ii.avg = -1.0;
		ks[0] = bwa_seq_open( fn_fa1);
		ks[1] = bwa_seq_open( fn_fa2);

		bwt[0]=BwtIndex.bwt_d;
		bwt[1]=BwtIndex.rbwt_d;

		bwa_print_sam_SQ(BwtIndex.bns);
		bwa_print_sam_PG();

		while ((seqs[0] = bwa_read_seq(ks[0], 0x40000, &n_seqs, opt->mode , opt->trim_qual)) != 0) {// opt should be different for two fa files theoretically
			int cnt_chg;
			isize_info_t ii;
			ubyte_t *pacseq=0;

			seqs[1] = bwa_read_seq(ks[1], 0x40000, &n_seqs, opt->mode, opt->trim_qual);
			tot_seqs += n_seqs;
			t = clock();
			fprintf(stderr,"Reading in %d sequences into buffer...\n",n_seqs);
			fprintf(stderr, "Calculate SA coordinate... ");


for(int iter=0;iter!=2;++iter)
{
	#ifdef HAVE_PTHREAD
			if (opt->n_threads <= 1) { // no multi-threading at all
				bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs[iter], opt);
			} else {
				pthread_t *tid;
				pthread_attr_t attr;
				thread_aux_t *data;
				int j;
				pthread_attr_init(&attr);
				pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
				data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
				tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
				for (j = 0; j < opt->n_threads; ++j) {
					data[j].tid = j; data[j].bwt[0] = bwt[0]; data[j].bwt[1] = bwt[1];
					data[j].n_seqs = n_seqs; data[j].seqs = seqs[iter]; data[j].opt = opt;
					pthread_create(&tid[j], &attr, worker, data + j);
				}
				for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
				free(data); free(tid);
			}
	#else
			//DBG(fprintf(stderr,"Before come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
			bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs[iter], opt);
			//DBG(fprintf(stderr,"After come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
	#endif
}
			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

			t = clock();
			//fprintf(stderr, "[bwa_aln_core] write to the disk... ");
			/*for(int iter=0;iter!=2;++iter)
			{
				for (i = 0; i < n_seqs; ++i) {
					bwa_seq_t *p = seqs[iter] + i;
				//fwrite(&p->n_aln, 4, 1, stdout);
				//if (p->n_aln) fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, stdout);
				}
			}*/

			fprintf(stderr, "[bwa_sai2sam_pe_core] convert to sequence coordinate... \n");
			cnt_chg = bwa_cal_pac_pos_pe( bwt, n_seqs, seqs,  &ii, popt, opt, &last_ii);
			fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
			fprintf(stderr, "[bwa_sai2sam_pe_core] changing coordinates of %d alignments.\n", cnt_chg);

			fprintf(stderr, "[bwa_sai2sam_pe_core] align unmapped mate...\n");
			pacseq = bwa_paired_sw(BwtIndex.bns,pacseq, n_seqs, seqs, popt, &ii);
			fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

			fprintf(stderr, "[bwa_sai2sam_pe_core] refine gapped alignments... ");
			for (j = 0; j < 2; ++j)
				bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs[j], pacseq, ntbns);
			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
			if (pacseq!= 0) free(pacseq);

			fprintf(stderr, "[bwa_sai2sam_pe_core] print alignments... ");
			for (i = 0; i < n_seqs; ++i) {
				bwa_seq_t *p[2];
				p[0] = seqs[0] + i; p[1] = seqs[1] + i;
				if (p[0]->bc[0] || p[1]->bc[0]) {
					strcat(p[0]->bc, p[1]->bc);
					strcpy(p[1]->bc, p[0]->bc);
				}
				bwa_print_sam1(BwtIndex.bns, p[0], p[1], opt->mode, opt->max_top2);
				bwa_print_sam1(BwtIndex.bns, p[1], p[0], opt->mode, opt->max_top2);
			}
			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

			for (j = 0; j < 2; ++j)
				bwa_free_read_seq(n_seqs, seqs[j]);
			fprintf(stderr, "[bwa_sai2sam_pe_core] %d sequences have been processed.\n", tot_seqs);
			last_ii = ii;
		}

		// destroy
		//bns_destroy(bns);
		//if (ntbns) bns_destroy(ntbns);
		for (i = 0; i < 2; ++i) {
			bwa_seq_close(ks[i]);
		}
		for (iter = kh_begin(g_hash); iter != kh_end(g_hash); ++iter)
			if (kh_exist(g_hash, iter)) free(kh_val(g_hash, iter).a);
		kh_destroy(64, g_hash);
		return 1;
}

BwtMapper::BwtMapper(BwtIndexer& BwtIndex,const string & Fastq_1, const string & Fastq_2, const pe_opt_t* popt, const gap_opt_t * opt )
{
	if(Fastq_2 != "Empty")
	{
		cerr<<"Now processing Pair End mapping..."<<endl;
		PairEndMapper(BwtIndex,Fastq_1.c_str(), Fastq_2.c_str(),  popt, opt);

	}
	else
	{
		cerr<<"Now processing Single End mapping..."<<endl;
		SingleEndMapper(BwtIndex,Fastq_1.c_str(),  opt);
	}
}



BwtMapper::~BwtMapper()
{
	// TODO Auto-generated destructor stub
}

