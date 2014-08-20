/*
 * BwtMapper.h
 *
 *  Created on: 2014Äê7ÔÂ9ÈÕ
 *      Author: Administrator
 */

#ifndef BWTMAPPER_H_
#define BWTMAPPER_H_
#include "BwtIndexer.h"
#include "./libbwa/bwtaln.h"
#include <iostream>
#ifndef BWAPE_H_
#include "./libbwa/bwape.h"
#endif
class BwtMapper
{
public:

	BwtMapper();
	BwtMapper(BwtIndexer& BwtIndex,const string & Fastq_1, const string & Fastq_2, const pe_opt_t* popt, const gap_opt_t * opt );
	int bwa_cal_pac_pos(BwtIndexer& BwtIndex, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr);
	int bwa_cal_pac_pos_pe(bwt_t *const _bwt[2], int n_seqs, bwa_seq_t *seqs[2],  isize_info_t *ii, const pe_opt_t *opt, const gap_opt_t *gopt, const isize_info_t *last_ii);
	//static void bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt[2], int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);
	bool SingleEndMapper(BwtIndexer& BwtIndex,const char *fn_fa, const gap_opt_t * opt);
	bool PairEndMapper(BwtIndexer& BwtIndex,const char *fn_fa1,const char * fn_fa2,const  pe_opt_t *popt,const gap_opt_t* opt);
	virtual ~BwtMapper();

	inline void  revealopt(const gap_opt_t* opt)
	{
		cerr<<opt->s_mm<<endl;
		cerr<<opt->s_gapo<<endl;
		cerr<<opt->s_gape<<endl;
		cerr<<opt->mode<<endl
				<<opt->indel_end_skip<<endl
				<<opt->max_del_occ<<endl
				<<opt->max_entries<<endl
				<<opt->fnr<<endl
				<<opt->max_diff<<endl
				<<opt->max_gapo<<endl
				<<opt->max_gape<<endl
				<<opt->max_seed_diff<<endl
				<<opt->seed_len<<endl
				<<opt->n_threads<<endl
				<<opt->max_top2<<endl
				<<opt->trim_qual<<endl
				<<"opt reveal end"<<endl;

	}
};

#endif /* BWTMAPPER_H_ */
