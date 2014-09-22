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
#include "./libbwa/bwtgap.h"
#include "./libbwa/bwase.h"
#include "./libbwa/bwape.h"
#include "./libbwa/khash.h"
#include "BamInterface.h"
#include <algorithm>
KHASH_MAP_INIT_INT64(64, poslist_t)
using namespace std;

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
#endif

extern void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs,
                                ubyte_t *_pacseq, bntseq_t *ntbns);
extern void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p,
                             const bwa_seq_t *mate, int mode, int max_top2);
extern void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s,
                               int set_main, int n_multi);
extern int bwa_approx_mapQ(const bwa_seq_t *p, int mm);
extern int bwt_cal_width(const bwt_t *rbwt, int len, const ubyte_t *str,
                           bwt_width_t *width);

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

  seed_w[0] = (bwt_width_t*) calloc(opt->seed_len + 1, sizeof(bwt_width_t));
  seed_w[1] = (bwt_width_t*) calloc(opt->seed_len + 1, sizeof(bwt_width_t));
  w[0] = w[1] = 0;
  uint32_t unmapped_num = 0;
  for (i = 0; i != n_seqs; ++i)
    {
      bwa_seq_t *p = seqs + i;

      if (Indexer->IsReadFiltered(p->seq, p->qual, p->len))
        {
          /*if(strcmp(p->name,"ERR013170.1716") ==0)
           {
           fprintf(stdout,"ERR013170.1716	is filtered\n" );
           exit(1);
           }*/
          p->sa = 0;
          p->type = BWA_TYPE_NO_MATCH;
          p->c1 = p->c2 = 0;
          p->n_aln = 0;
          p->aln = (bwt_aln1_t*) calloc(4, sizeof(bwt_aln1_t));
          unmapped_num++;
          continue;
        }

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
#endif
      p->sa = 0;
      p->type = BWA_TYPE_NO_MATCH;
      p->c1 = p->c2 = 0;
      p->n_aln = 0;
      p->aln = 0;
      seq[0] = p->seq;
      seq[1] = p->rseq;
      if (max_l < p->len)
        {
          max_l = p->len;
          w[0] = (bwt_width_t*) realloc(w[0],
                                        (max_l + 1) * sizeof(bwt_width_t));
          w[1] = (bwt_width_t*) realloc(w[1],
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
  fprintf(stderr, "RollingHash filtered %d reads...\n", unmapped_num);
  free(seed_w[0]);
  free(seed_w[1]);
  free(w[0]);
  free(w[1]);
  gap_destroy_stack(stack);
}
BwtMapper::BwtMapper()
{}
int BwtMapper::bwa_cal_pac_pos(BwtIndexer& BwtIndex, int n_seqs,
                               bwa_seq_t *seqs, int max_mm, float fnr)
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

typedef struct
  {
    kvec_t(bwt_aln1_t)
    aln;
  }
aln_buf_t;
int BwtMapper::bwa_cal_pac_pos_pe(bwt_t * const _bwt[2], int n_seqs,
                                  bwa_seq_t *seqs[2], isize_info_t *ii, const pe_opt_t *opt,
                                  const gap_opt_t *gopt, const isize_info_t *last_ii)
{
  int i, j, cnt_chg = 0;
  //char str[1024];
  bwt_t *bwt[2];
  pe_data_t *d;
  aln_buf_t *buf[2];

  d = (pe_data_t*) calloc(1, sizeof(pe_data_t));
  buf[0] = (aln_buf_t*) calloc(n_seqs, sizeof(aln_buf_t));
  buf[1] = (aln_buf_t*) calloc(n_seqs, sizeof(aln_buf_t));
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
  for (i = 0; i != n_seqs; ++i)
    { //for each seqs
      bwa_seq_t *p[2];
      for (j = 0; j < 2; ++j)
        { // for each end
          int n_aln;
          p[j] = seqs[j] + i;
          p[j]->n_multi = 0;
          p[j]->extra_flag |= SAM_FPD | (j == 0 ? SAM_FR1 : SAM_FR2);
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
  DBG(cerr<<"Arrived check point 2....\n";)
  // infer isize
  //infer_isize(n_seqs, seqs, ii, opt->ap_prior, bwt[0]->seq_len);
  //if (ii->avg < 0.0 && last_ii->avg > 0.0) *ii = *last_ii;
  //if (opt->force_isize) {
  //fprintf(stderr, "[%s] discard insert size estimate as user's request.\n", __func__);
  ii->low = ii->high = 0;
  ii->avg = ii->std = -1.0;
  //}

  // PE
  for (i = 0; i != n_seqs; ++i)
    {
      bwa_seq_t *p[2];
      for (j = 0; j < 2; ++j)
        {
          p[j] = seqs[j] + i;
          kv_copy(bwt_aln1_t, d->aln[j], buf[j][i].aln); //copy aln back to d
        }
      if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
          && (p[1]->type == BWA_TYPE_UNIQUE
              || p[1]->type == BWA_TYPE_REPEAT))
        { // only when both ends mapped
          uint64_t x;
          int j, k, n_occ[2];
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
                      uint64_t key = (uint64_t) r->k << 32 | r->l; // key is formed by lower and upper bound
                      int ret;
                      khint_t iter = kh_put(64, g_hash, key, &ret);
                      if (ret)
                        { // if this key is not in the hash table; ret must equal 1 as we never remove elements
                          poslist_t *z = &kh_val(g_hash, iter); // return the bwtint_ts pointed by this iter
                          z->n = r->l - r->k + 1;
                          z->a = (bwtint_t*) malloc(sizeof(bwtint_t) * z->n);
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
      DBG(cerr<<"Arrived check point 3....\n";)
      if (opt->N_multi || opt->n_multi)
        {
          DBG(cerr<<"Arrived check point 3-a....\n";)
          for (j = 0; j < 2; ++j)
            {
              DBG(cerr<<"Arrived check point 3-0....\n";)
              if (p[j]->type != BWA_TYPE_NO_MATCH)
                {
                  int k;
                  if (!(p[j]->extra_flag & SAM_FPP)
                      && p[1 - j]->type != BWA_TYPE_NO_MATCH)
                    {
                      DBG(cerr<<"Arrived check point 3-1....\n";)
                      bwa_aln2seq_core(d->aln[j].n, d->aln[j].a, p[j], 0,
                                       p[j]->c1 + p[j]->c2 - 1 > opt->N_multi ?
                                       opt->n_multi : opt->N_multi);
                      DBG(cerr<<"Arrived check point 3-2....\n";)
                    }
                  else
                    {
                      DBG(cerr<<"Arrived check point 3-3....\n";)
                      bwa_aln2seq_core(d->aln[j].n, d->aln[j].a, p[j], 0,
                                       opt->n_multi);
                      DBG(cerr<<"Arrived check point 3-4....\n";)
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
  bwa_rg_id = calloc(q - p + 1, 1);
  for (q = p, r = bwa_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
    *r++ = *q;
  return 0;
}

extern int64_t pos_end_multi(const bwt_multi1_t *p, int len); // analogy to pos_end()
bool BwtMapper::SetSamFileHeader(SamFileHeader& SFH, const bntseq_t * bns)
{
  //cerr<<"Setting  Bam file..."<<endl;
  //ostringstream handler;
  //handler << "@PG\tID:FastqA\tPN:FastqA\tVN:" << PACKAGE_VERSION << endl;
  //string PG(handler.str());

  if(!SFH.setPGTag("Version", PACKAGE_VERSION, "FastqA"))
    std::cerr<<"SetPGTag failed"<<endl;
  if (bwa_rg_line&&strstr(bwa_rg_line, "@RG") == bwa_rg_line)
    if(!SFH.setRGTag("RG", bwa_rg_line, bwa_rg_id))
      std::cerr<<"SetRGTag failed"<<endl;

  for (int i = 0; i < bns->n_seqs; ++i)
    {

      if(!SFH.setSQTag("LN", std::to_string(bns->anns[i].len).c_str(), bns->anns[i].name))
        std::cerr<<"SetSQTag failed"<<endl;
      //printf("@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
    }
  /*
   h->l_text = PG.length();
   h->text = (char*) malloc(h->l_text + 1);
   h->text[h->l_text] = 0;
   strcpy(h->text, PG.c_str());
   h->n_targets = BwtIndex.bns->n_seqs;
   h->target_name = (char**) calloc(h->n_targets, sizeof(char*));
   h->target_len = (uint32_t*) calloc(h->n_targets, 4);
   for (i = 0; i != h->n_targets; ++i)
   {
   int name_len = strlen(BwtIndex.bns->anns[i].name);
   //if (fp->is_be) ed_swap_4p(&name_len);
   h->target_name[i] = (char*) calloc(name_len, 1);
   strcpy(h->target_name[i], BwtIndex.bns->anns[i].name);
   h->target_len[i] = BwtIndex.bns->anns[i].len;
   //if (fp->is_be) ed_swap_4p(&h->target_len[i]);
   }*/
  return 0;
}
bool BwtMapper::SetSamRecord(const bntseq_t *bns, bwa_seq_t *p,
                             const bwa_seq_t *mate, int mode, int max_top2, SamFileHeader& SFH,
                             SamRecord& SR)
{
  int j;
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
      SR.setReferenceName(SFH,bns->anns[seqid].name);
      //printf("%d\t%d\t", (int) (p->pos - bns->anns[seqid].offset + 1), p->mapQ);
      //ss<<(int) (p->pos - bns->anns[seqid].offset + 1)<<"\t"<<p->mapQ;
      SR.set1BasedPosition((int) (p->pos - bns->anns[seqid].offset + 1));
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
          (seqid == m_seqid) ? SR.setMateReferenceName(SFH, "="):SR.setMateReferenceName(SFH, bns->anns[m_seqid].name);
          isize = (seqid == m_seqid) ? pos_5(mate) - pos_5(p) : 0;
          if (p->type == BWA_TYPE_NO_MATCH)
            isize = 0;
          //printf("%d\t%lld\t", (int) (mate->pos - bns->anns[m_seqid].offset + 1), isize);
          //ss<<(int) (mate->pos - bns->anns[m_seqid].offset + 1)<<"\t"<<isize;
          SR.set1BasedMatePosition((int) (mate->pos - bns->anns[m_seqid].offset + 1));
          SR.setInsertSize(isize);
        }
      else if (mate)
        //printf("\t=\t%d\t0\t", (int) (p->pos - bns->anns[seqid].offset + 1));
        //ss<<"\t=\t"<<(int) (p->pos - bns->anns[seqid].offset + 1)<<"\t0\t";
        {
          SR.setMateReferenceName(SFH, "=");
          SR.set1BasedMatePosition((int) (p->pos - bns->anns[seqid].offset + 1));
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
            ss << "ACGTN"[(int) p->seq[j]];
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
                  ss << (int) (q->pos - bns->anns[seqid].offset + 1)<<",";
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
      SR.setReferenceName(SFH,"*");
      SR.set1BasedPosition(0);
      SR.setMapQuality(0);
      SR.setCigar("*");
      SR.setMateReferenceName(SFH,"*");
      SR.set1BasedMatePosition(0);
      SR.setInsertSize(0);
      ss.clear();
      ss.str("");
      for (j = 0; j != p->len; ++j)
        //putchar("ACGTN"[(int) s[j]]);
        ss << "ACGTN"[(int) s[j]];
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
        SR.addTag("RG",'Z',bwa_rg_id);
      if (p->bc[0])
        //printf("\tBC:Z:%s", p->bc);
        //ss << "\tBC:Z:", p->bc;
        SR.addTag("BC",'Z',p->bc);
      if (p->clip_len < p->full_len)
        //printf("\tXC:i:%d", p->clip_len);
        //ss << "\tXC:i:" << p->clip_len;
        SR.addTag("XC",'i',to_string(p->clip_len).c_str());
      //putchar('\n');
      //ss << "\n";
    }
  return 0;
}

/*********bam format*********************************************************/
/*
 static inline int hts_reg2bin(int64_t beg, int64_t end, int min_shift, int n_lvls)
 {
 int l, s = min_shift, t = ((1<<((n_lvls<<1) + n_lvls)) - 1) / 7;
 for (--end, l = n_lvls; l > 0; --l, s += 3, t -= 1<<((l<<1)+l))
 if (beg>>s == end>>s) return t + (beg>>s);
 return 0;
 }
 static inline int bam_reg2bin(uint32_t beg, uint32_t end)
 {
 return hts_reg2bin(beg, end, 14, 5);
 }
 #define BAM_CIGAR_STR   "MIDNSHP=XB"
 #define BAM_CIGAR_SHIFT 4
 #define BAM_CIGAR_MASK  0xf
 #define BAM_CIGAR_TYPE  0x3C1A7
 #define bam_get_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
 #define bam1_cigar(b) (bam_get_cigar((b)))
 #define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference
 #define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
 #define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
 static int bam_cigar2rlen(int n_cigar, const uint32_t *cigar)
 {
 int k, l;
 for (k = l = 0; k < n_cigar; ++k)
 if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
 l += bam_cigar_oplen(cigar[k]);
 return l;
 }
 static inline uint32_t bam_calend(const bam1_core_t *c, const uint32_t *cigar) { return c->pos + (c->n_cigar? bam_cigar2rlen(c->n_cigar, cigar) : 1); }

 */

/*********bam format*********************************************************/
bool BwtMapper::SingleEndMapper(BwtIndexer& BwtIndex, const char *fn_fa,
                                const string & VcfPath, const gap_opt_t * opt)
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
    g_log_n[i] = (int) (4.343 * log(i) + 0.5);
  srand48(BwtIndex.bns->seed);
  ks = bwa_seq_open(fn_fa);

  bwt[0] = BwtIndex.bwt_d;
  bwt[1] = BwtIndex.rbwt_d;

  /*********bam header init***********/
  /*ostringstream handler;
   BGZF *fp;
   fp =
   strcmp(opt->bam_name, "-") ?
   bgzf_open(opt->bam_name, "r") : bgzf_dopen(fileno(stdin), "r");
   bam1_t *b;
   b = bam_init1();
   bam_hdr_t *h;
   h = bam_hdr_init();*/
  SamFileHeader SFH;
  BamInterface BamIO;
  IFILE BamFile=new InputFile(opt->bam_name,"w", InputFile::ifileCompression::BGZF);
  StatGenStatus StatusTracker;
  /*********bam header init end***********/

  if (!opt->bam_name)
    {
      bwa_print_sam_SQ(BwtIndex.bns);
      bwa_print_sam_PG();
    }
  else
    {
      if(!BamFile->isOpen())
        {
          cerr<<"Open Bam file for writing failed, abort!"<<endl;
          exit(1);
        }
      SetSamFileHeader(SFH,BwtIndex.bns);
      BamIO.writeHeader(BamFile,SFH,StatusTracker);
    }

  t = clock();
  fprintf(stderr, "Restore Variant Site Info...\n ");
  collector.restoreVcfSites(VcfPath, opt);
  while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual))
         != 0)
    {
      tot_seqs += n_seqs;

      fprintf(stderr, "Reading in %d sequences into buffer...", n_seqs);
      fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
      t = clock();
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

      fprintf(stderr, "Calculate SA coordinate and ");
      fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
      t = clock();

      for (i = 0; i < n_seqs; ++i)
        {
          bwa_seq_t *p = seqs + i;
          bwa_aln2seq_core(p->n_aln, p->aln, p, 1, N_OCC);

        }

      fprintf(stderr, "Convert to sequence coordinate... ");
      bwa_cal_pac_pos(BwtIndex, n_seqs, seqs, opt->max_diff, opt->fnr); // forward bwt will be destroyed here
      fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
      t = clock();

      fprintf(stderr, "Refine gapped alignments... ");

      //DBG(fprintf(stderr,"Before come into refined gap...%s\n%s\n",seqs->name,seqs->seq);)
      if (BwtIndex.bns->fp_pac == 0)
        bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs, BwtIndex.pac_buf, ntbns);
      else
        bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs, 0, ntbns);
      fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
      t = clock();

      fprintf(stderr, "Print alignments... ");

      if (!opt->bam_name)
        {
          for (i = 0; i < n_seqs; ++i)
            {
              //collector.addAlignment(string(BwtIndex.bns->anns[seqid].name),(seqs+i)->seq,(seqs+i)->qual,(seqs+i)->n_cigar,(seqs+i)->cigar,(seqs+i)->md,(int)((seqs+i)->pos - BwtIndex.bns->anns[seqid].offset + 1),opt);
              if (!collector.addAlignment(BwtIndex.bns, seqs + i, opt))
                continue; //failed
              bwa_print_sam1(BwtIndex.bns, seqs + i, 0, opt->mode, opt->max_top2);
              //exit(1);
            }
        }
      else
        {
          for (i = 0; i < n_seqs; ++i)
            {
              if (!collector.addAlignment(BwtIndex.bns, seqs + i, opt))
                continue; //failed
              SamRecord SR;
              SetSamRecord(BwtIndex.bns, seqs + i, 0, opt->mode, opt->max_top2, SFH,SR);
              //std::cerr<<"\nPassed Read:"<<(seqs+i)->name<<"\t"<<SR.getCigar()<<"\t"<<SR.getSequence()<<"\t"<<SR.getQuality()<<endl;//<<std::for_each((seqs+i),(seqs+i)+(seqs+i)->len-1, [](bwa_seq_t* s, int j ){return  "ACGTN"[(int) *s+j];})<<endl;
              BamIO.writeRecord(BamFile,SFH,SR,SamRecord::SequenceTranslation::NONE);
            }
        }

      fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
      t = clock();
      bwa_free_read_seq(n_seqs, seqs);
      fprintf(stderr, " %d sequences have been processed.\n", tot_seqs);
      //t = clock();
    } //end while
  //bam_destroy1(b);
  if(BamFile->ifclose())
    std::cerr<<"Bam file close failed!"<<endl;
  delete BamFile;
  fprintf(stderr, "Calculate distributions... ");
  collector.processCore(VcfPath);
  // destroy
  //bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
  bwa_seq_close(ks);
}
bool BwtMapper::PairEndMapper(BwtIndexer& BwtIndex, const char *fn_fa1,
                              const char * fn_fa2, const string & VcfPath, const pe_opt_t *popt,
                              const gap_opt_t* opt)
{

  int i, j, n_seqs, tot_seqs = 0; //,m_aln;
  bwa_seq_t *seqs[2];
  bwa_seqio_t *ks[2];
  clock_t t;
  bwt_t *bwt[2];
  bntseq_t *ntbns = 0;
  khint_t iter;
  isize_info_t last_ii; // this is for the last batch of reads
  // initialization
  bwase_initialize(); // initialize g_log_n[] in bwase.c
  for (i = 1; i != 256; ++i)
    g_log_n[i] = (int) (4.343 * log(i) + 0.5);
  srand48(BwtIndex.bns->seed);
  g_hash = kh_init(64);
  last_ii.avg = -1.0;
  ks[0] = bwa_seq_open(fn_fa1);
  ks[1] = bwa_seq_open(fn_fa2);

  bwt[0] = BwtIndex.bwt_d;
  bwt[1] = BwtIndex.rbwt_d;

  SamFileHeader SFH;

  BamInterface BamIO;
  IFILE BamFile=new InputFile(opt->bam_name,"w", InputFile::ifileCompression::BGZF);
  StatGenStatus StatusTracker;
  /*********bam header init end***********/

  if (!opt->bam_name)
    {
      bwa_print_sam_SQ(BwtIndex.bns);
      bwa_print_sam_PG();
    }
  else
    {
      if(!BamFile->isOpen())
        {
          cerr<<"Open Bam file for writing failed, abort!"<<endl;
          exit(1);
        }
      SetSamFileHeader(SFH,BwtIndex.bns);
      BamIO.writeHeader(BamFile,SFH,StatusTracker);
    }
  fprintf(stderr, "Restore Variant Site Info...\n");
  collector.restoreVcfSites(VcfPath, opt);
  ofstream fout(VcfPath + ".InsertSizeTable");
  int total_add = 0;
  while ((seqs[0] = bwa_read_seq(ks[0], 0x40000, &n_seqs, opt->mode,
                                 opt->trim_qual)) != 0)
    { // opt should be different for two fa files theoretically
      int cnt_chg;
      isize_info_t ii;
      ubyte_t *pacseq = 0;

      seqs[1] = bwa_read_seq(ks[1], 0x40000, &n_seqs, opt->mode, opt->trim_qual);
      tot_seqs += n_seqs;
      t = clock();
      fprintf(stderr, "Reading in %d sequences into buffer...\n", n_seqs);

      for (int iter = 0; iter != 2; ++iter)
        {
#ifdef HAVE_PTHREAD
          if (opt->n_threads <= 1)
            { // no multi-threading at all
              bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs[iter], opt, &BwtIndex);
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
                  data[j].seqs = seqs[iter];
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
          bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs[iter], opt, &BwtIndex);
          //DBG(fprintf(stderr,"After come into cal sa reg gap...%s\n%d\nlength:%d\n",seqs->name,seqs->seq[0],seqs->len);)
#endif

        }
      fprintf(stderr, "Calculate SA coordinate... ");
      fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
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

      fprintf(stderr, "convert to sequence coordinate... \n");
      cnt_chg = bwa_cal_pac_pos_pe(bwt, n_seqs, seqs, &ii, popt, opt, &last_ii);
      fprintf(stderr, "time elapses: %.2f sec\n",
              (float) (clock() - t) / CLOCKS_PER_SEC);
      t = clock();
      fprintf(stderr, "changing coordinates of %d alignments.\n", cnt_chg);

      fprintf(stderr, "align unmapped mate...\n");
      if (BwtIndex.bns->fp_pac == 0) //indexing path
        /*pacseq = */
        bwa_paired_sw(BwtIndex.bns, BwtIndex.pac_buf, n_seqs, seqs,
                      popt, &ii);
      else
        /*pacseq = */
        bwa_paired_sw(BwtIndex.bns, pacseq, n_seqs, seqs, popt, &ii);
      fprintf(stderr, "time elapses: %.2f sec\n",
              (float) (clock() - t) / CLOCKS_PER_SEC);
      t = clock();

      fprintf(stderr, "refine gapped alignments... ");
      for (j = 0; j < 2; ++j)
        if (BwtIndex.bns->fp_pac == 0) //indexing path
          bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs[j], BwtIndex.pac_buf,
                            ntbns);
        else
          bwa_refine_gapped(BwtIndex.bns, n_seqs, seqs[j], 0, ntbns);
      fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
      t = clock();
      //if (pacseq!= 0) free(pacseq);

      fprintf(stderr, "print alignments... ");
      if (!opt->bam_name)
        {
          for (i = 0; i < n_seqs; ++i)
            {
              bwa_seq_t *p[2];
              p[0] = seqs[0] + i;
              p[1] = seqs[1] + i;
              if (p[0]->bc[0] || p[1]->bc[0])
                {
                  strcat(p[0]->bc, p[1]->bc);
                  strcpy(p[1]->bc, p[0]->bc);
                }

              if (!collector.addPairAlignment(BwtIndex.bns, seqs[0] + i, seqs[1] + i, opt,
                                              fout, total_add))
                continue;
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
              if (p[0]->bc[0] || p[1]->bc[0])
                {
                  strcat(p[0]->bc, p[1]->bc);
                  strcpy(p[1]->bc, p[0]->bc);
                }

              if (!collector.addPairAlignment(BwtIndex.bns, seqs[0] + i, seqs[1] + i, opt,
                                              fout, total_add))
                continue;
              SamRecord SR[2];
              SetSamRecord(BwtIndex.bns, seqs[0] + i, seqs[1]+i, opt->mode, opt->max_top2, SFH,SR[0]);
              BamIO.writeRecord(BamFile,SFH,SR[0],SamRecord::SequenceTranslation::NONE);
              SetSamRecord(BwtIndex.bns, seqs[1] + i, seqs[0]+i, opt->mode, opt->max_top2, SFH,SR[1]);
              //std::cerr<<"\nPassed Read:"<<(seqs+i)->name<<"\t"<<SR.getCigar()<<"\t"<<SR.getSequence()<<"\t"<<SR.getQuality()<<endl;//<<std::for_each((seqs+i),(seqs+i)+(seqs+i)->len-1, [](bwa_seq_t* s, int j ){return  "ACGTN"[(int) *s+j];})<<endl;
              BamIO.writeRecord(BamFile,SFH,SR[1],SamRecord::SequenceTranslation::NONE);
            }
        }

      fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
      //cerr<<"In total "<<total_add<<" reads were calculated!"<<endl;
      t = clock();

      for (j = 0; j < 2; ++j)
        bwa_free_read_seq(n_seqs, seqs[j]);
      fprintf(stderr, " %d sequences have been processed.\n", tot_seqs);

      last_ii = ii;
    } //end while
  cerr << "In total " << total_add << " reads were calculated!" << endl;
  fout.close();
  BamFile->ifclose();
  // destroy
  fprintf(stderr, "Calculate distributions... ");
  collector.processCore(VcfPath);
  //bns_destroy(bns);
  //if (ntbns) bns_destroy(ntbns);
  for (i = 0; i < 2; ++i)
    {
      bwa_seq_close(ks[i]);
    }
  for (iter = kh_begin(g_hash); iter != kh_end(g_hash); ++iter)
    if (kh_exist(g_hash, iter))
      free(kh_val(g_hash, iter).a);
  kh_destroy(64, g_hash);
  return 1;
}

BwtMapper::BwtMapper(BwtIndexer& BwtIndex, const string & Fastq_1,
                     const string & Fastq_2, const string & VcfPath, const pe_opt_t* popt,
                     const gap_opt_t * opt)
{
	//std::cerr<<"Open Fastq  ... "<<endl;
  bwa_set_rg(opt->RG);
  if (Fastq_2 != "Empty")
    {
      cerr << "Now processing Pair End mapping..." << endl;

      PairEndMapper(BwtIndex, Fastq_1.c_str(), Fastq_2.c_str(), VcfPath, popt, opt);

    }
  else
    {
      cerr << "Now processing Single End mapping..." << endl;

      SingleEndMapper(BwtIndex, Fastq_1.c_str(), VcfPath, opt);
    }
}
BwtMapper::BwtMapper(BwtIndexer& BwtIndex, const string & FaList,
                     const string & VcfPath, const pe_opt_t* popt,
                     const gap_opt_t * opt)
{
	std::cerr<<"Open Fastq List ... "<<endl;
  bwa_set_rg(opt->RG);
  ifstream fin(FaList);
  std::string line,Fastq_1,Fastq_2;
  int i(0);
  while(getline(fin,line))
    {
	  ++i;
	  stringstream ss(line);
      ss>>Fastq_1>>Fastq_2;
      if (Fastq_2 != "")
        {
          cerr << "Now processing Pair End mapping\t"<<i<<"\t..." << endl;
          cerr<<Fastq_1<<endl;
          cerr<<Fastq_2<<endl;
          PairEndMapper(BwtIndex, Fastq_1.c_str(), Fastq_2.c_str(), VcfPath, popt, opt);

        }
      else
        {
          cerr << "Now processing Single End mapping\t"<<i<<"\t..." << endl;
          cerr<<Fastq_1<<endl;
          SingleEndMapper(BwtIndex, Fastq_1.c_str(), VcfPath, opt);
        }

    }
}

BwtMapper::~BwtMapper()
{
  // TODO Auto-generated destructor stub
}

