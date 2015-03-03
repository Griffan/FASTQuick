/*
 * BwtIndexer.h
 *
 *  Created on: 2014Äê7ÔÂ9ÈÕ
 *      Author: Administrator
 */

#ifndef BWTINDEXER_H_
#define BWTINDEXER_H_
#include "../libbwa/bntseq.h"
#include "../libbwa/bwt.h"
#include <cstdio>
#include <cstdlib>
#include "RefBuilder.h"
#include <sys/stat.h>
#include <assert.h>
#include <cmath>
#include <time.h>
#include <iostream>
//Bloom filter parameters
#define VCF_SITES 10000
#define WINDOW_SIZE 1000
#define KMER_SIZE  32
//#define BLOOM_FPP 0.1
#define PROCESS_RATIO 0.5
#define READ_STEP_SIZE 1
static bool debug_flag=false;
struct _RollParam
{
	int kmer_size;
	int read_step_size;
	int thresh;
};
//extern struct _RollParam RollParam;


class BwtIndexer
{
public:


	std::string RefPath;
	ubyte_t * pac_buf;
	ubyte_t * rpac_buf;
	unsigned char* bwt_buf;
	unsigned char* rbwt_buf;
	bntseq_t *bns;
	int l_buf;
	int32_t m_pac_buf, m_rpac_buf, m_bwt_buf, m_rbwt_buf;

	bwt_t *bwt_d;
	bwt_t *rbwt_d;

	std::unordered_map<std::string, bool> longRefTable;

	BwtIndexer();

	BwtIndexer(std::string & NewRef);

	BwtIndexer(RefBuilder & ArtiRef, std::string & NewRef);

	bool LoadIndex(std::string & NewRef);

	bool BuildIndex(RefBuilder & ArtiRef, string & OldRef, std::string & NewRef,
			const gap_opt_t * opt);

	bool Fa2Pac(RefBuilder & ArtiRef, const char *prefix, const gap_opt_t* op);

	bool Fa2RevPac(const char * prefix);

	bwt_t* Pac2Bwt(unsigned char * pac);

	void bwt_bwtupdate_core(bwt_t *bwt);

	void bwt_cal_sa(bwt_t *bwt, int intv);

	bool IsReadHighQ(const ubyte_t *Q, int len)const;

	bool IsKmerInHash(uint64_t kmer)const;

	bool IsReadInHash(const ubyte_t * S, int len)const;

	//typedef uint32_t v4si __attribute__ ((vector_size (16)));
	bool IsReadFiltered(const ubyte_t * S, const ubyte_t * Q, int len)const;

#ifdef BLOOM_FPP
#include "bloom_filter.hpp"
	bloom_parameters parameters;
	bloom_filter filter;
	ubyte_t BloomReadBuff[KMER_SIZE+1];
	inline void AddSeq2BloomFliter(const std::string & Seq)
	{
		//std::string tmp;
		for(int i=0;i!=Seq.length()-KMER_SIZE/*last one considered*/;++i)
		{
			//tmp=Seq.substr(i,KMER_SIZE);
			//for(int j=0;j!=KMER_SIZE;++j)
			//tmp[i]=nst_nt4_table[(int)tmp[i]];
			//printf("/*******now insert: %s**************/\n",Seq.substr(i,KMER_SIZE).c_str());
			filter.insert(Seq.substr(i,KMER_SIZE));

		}
	}

	inline bool IsReadInHash( const char * S, int len)
	{
		//std::string Seq(S);
		int processed(0);
		//int total= len-KMER_SIZE+1;

		for(int i=0;i<len-KMER_SIZE/*last one considered*/;i+=READ_STEP_SIZE,++processed)
		{
			// if(filter.contains(Seq.substr(i,KMER_SIZE)))
			//assert(BloomReadBuff!=0);
			//printf("/*********************/\n");
			//for(int j=0;j!=KMER_SIZE;++j)
			//printf("number %d : %d\n",i,S[i]);
			//printf("/*********************/\n");
			//memcpy(BloomReadBuff,&S[i],KMER_SIZE*sizeof(char));
			//BloomReadBuff[KMER_SIZE]='\0';
			//printf("/*******now test: %s**************/\n",std::string(BloomReadBuff).c_str());
			if(filter.contains(std::string(S).substr(i,KMER_SIZE)))
			{
				return true;
			}
			//else if (total*READ_STEP_SIZE>PROCESS_RATIO*total)
			//return false;
		}
//		if(double(hit)/total > HIT_RATIO)// is in hash
//		{

//			return true;
//		}
//		else
		{
			//printf("the total is: %d\n the hit is : %d\n the ratio is :%f\n",total,hit, double(hit)/total);
			return false;
		}
	}

#else
	_RollParam RollParam;
	//Fast Rolling Hash parameters
#define OVERFLOWED_KMER_SIZE 32
#define LSIZE(k)     (1L << (2*(k)))    /* language size, 4096, 16384          */
#define LOMEGA(k)    (LSIZE(k)-1)      /* last (all-one) word, 4095, 16383,   */
#define MERMASK(k)   ((1 << (k))-1)    /* 6->63 and 7->127 (K-bits-all-one)   */

#define min_only(X,Y) ((X) < (Y) ? (X) : (Y))

#define KMER_64_BIT 1
#ifdef KMER_8_BIT
	unsigned int mask[6]=
	{	240, //11110000
		15,//00001111
		195,//11000011
		60,//00111100
		204,//11001100
		51}; //00110011
#endif
#ifdef KMER_16_BIT
	unsigned int mask[6] =
	{	65280, //1111111100000000
		255,//0000000011111111
		61455,//1111000000001111
		4080,//00001111111110000
		12528,//1111000011110000
		3855}; //0000111100001111
#endif
#ifdef KMER_32_BIT
	unsigned int mask[6] =
	{ 0xffff0000, 0xffff, 0xff0000ff, 0xffff00, 0xff00ff00, 0xff00ff };
#endif
#ifdef KMER_64_BIT
	uint64_t mask[6] =
	{ 0xffffffff00000000, 0xffffffff, 0xffff00000000ffff, 0xffffffff0000, 0xffff0000ffff0000, 0xffff0000ffff };
#endif

	unsigned char * roll_hash_table[6]; //11110000

	long long int hash_table_size;

	void InitializeRollHashTable();

	void ReadRollHashTable(const std::string& prefix);

	void DumpRollHashTable(const std::string& prefix);

	void DestroyRollHashTable();

	std::string ReverseComplement(const std::string & a);

	uint32_t KmerShrinkage(uint64_t a, unsigned int mask)const;

	void AddSeq2HashCore(const std::string & Seq, int iter);

	void AddSeq2Hash(const std::string & Seq, const gap_opt_t * opt);

#endif
	virtual ~BwtIndexer();
};
static double QualQuickTable[46] =
{ 1, 0.7943, 0.6310, 0.5012, 0.3981, 0.3162, 0.2512, 0.1995, 0.1585, 0.1259,
		0.1000, 0.0794, 0.0631, 0.0501, 0.0398, 0.0316, 0.0251, 0.0200, 0.0158,
		0.0126, 0.0100, 0.0079, 0.0063, 0.0050, 0.0040, 0.0032, 0.0025, 0.0020,
		0.0016, 0.0013, 0.0010, 0.0008, 0.0006, 0.0005, 0.0004, 0.0003, 0.0003,
		0.0002, 0.0002, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001 };
static std::unordered_map<char, char> match_table =
	{
	{ 'a', 'T' },
	{ 'A', 'T' },
	{ 'c', 'G' },
	{ 'C', 'G' },
	{ 'G', 'C' },
	{ 'g', 'C' },
	{ 't', 'A' },
	{ 'T', 'A' } };

inline bool BwtIndexer::IsReadHighQ(const ubyte_t *Q, int len)const//check if first 25 bases are of high quality
{
	double tmp(1);
	for (int i = 0; i != 25; ++i)
	{
		//std::cerr << int(Q[i]) << "\t";
		tmp *= (1 - QualQuickTable[Q[i]]);
	}
	//std::cerr << std::endl;
	if (1 - tmp < 0.1)
		return true;
	else
		return false;
}

inline uint32_t BwtIndexer::KmerShrinkage(uint64_t kmer, unsigned int iter)const
{
	//uint32_t tmp(0), i(1);
	////printf("the incoming kmer:%016llx\t",a);
	//while (mask) //32 circle
	//{
	//	if (mask & 1)
	//	{
	//		tmp |= ((a & 3) << (2 * i - 2));
	//		++i;
	//	}

	//	mask >>= 1;
	//	a >>= 2;
	//}
	////printf("the resulted i is :%d\n",i);
	//return tmp;

	uint32_t shrinked(0);
	switch (iter)
	{
	case 0:
	/*0xffffffff00000000*/
	shrinked |= (kmer & 0xffffffff00000000) >> 32;
	return shrinked;
	case 1:
	/*0xffffffff*/
	shrinked |= kmer & 0xffffffff;
	return shrinked;
	case 2:
	/*0xffff00000000ffff*/
	shrinked |= (kmer & 0xffff000000000000) >> 32;
	shrinked |= (kmer & 0xffff);
	return shrinked;
	case 3:
	/*0xffffffff0000*/
	shrinked |= (kmer & 0xffffffff0000) >> 16;
	return shrinked;
	case 4:
	/*0xffff0000ffff0000*/
	shrinked |= (kmer & 0xffff000000000000) >> 32;
	shrinked |= (kmer & 0xffff0000) >> 16;
	return shrinked;
	case 5:
	/*0xffff0000ffff*/
	shrinked |= (kmer & 0xffff00000000) >> 16;
	shrinked |= (kmer & 0xffff);
	return shrinked;
	default:
		fprintf(stderr,"Fatal error: unknown mask.\n");
		exit(EXIT_FAILURE);
	}
	
}

inline void BwtIndexer::AddSeq2Hash(const std::string & Seq,
		const gap_opt_t * opt)
{
	for (int i = 0; i != 6; ++i)
	{
		AddSeq2HashCore(Seq, i);
	}
}


#endif /* BWTINDEXER_H_ */
