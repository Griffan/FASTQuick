/*
 * BwtIndexer.cpp
 *
 *  Created on: 2014Äê7ÔÂ9ÈÕ
 *      Author: Administrator
 */

#include "BwtIndexer.h"
#include "RefBuilder.h"
#include "../libbwa/bwt.h"
#include "stdint.h"
#include "../libbwa/utils.h"
#include "stdio.h"
#include "../libbwa/bntseq.h"
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "../libbwa/kseq.h"
#include <iostream>
#include "../misc/Error.h"
using namespace std;
extern void notice(const char*,...);
extern void warning(const char*,...);
extern void error(const char*,...);
typedef uint32_t u_int32_t;

#ifdef DEBUG
#define DBG(CODE) CODE
#else
#define DBG(CODE)
#endif
struct stat sb;

extern int is_bwt(ubyte_t *T, int n);
//extern int64_t bwa_seq_len(const char *fn_pac);
extern bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is);
extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);
BwtIndexer::BwtIndexer() :RefPath(),
		pac_buf(0), rpac_buf(0), bwt_buf(0), rbwt_buf(0), bns(0), l_buf(0), m_pac_buf(
				0), m_rpac_buf(0), m_bwt_buf(0), m_rbwt_buf(0), bwt_d(0), rbwt_d(
				0)
{
	//***bloom filter initialization
#ifdef BLOOM_FPP
	// How many elements roughly do we expect to insert?
	parameters.projected_element_count = VCF_SITES *( WINDOW_SIZE-KMER_SIZE+1+KMER_SIZE) * 2;

	// Maximum tolerable false positive probability? (0,1)
	parameters.false_positive_probability = BLOOM_FPP;// 1 in 10000

	// Simple randomizer (optional)
	parameters.random_seed = 0xA5A5A5A5;

	if (!parameters)
	{
		std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
		exit(1);
	}
	parameters.compute_optimal_parameters();

	filter.set_param(parameters);

	//***bloom filter initialization end
#endif

	InitializeRollHashTable();

}
BwtIndexer::BwtIndexer(RefBuilder & ArtiRef, string & prefix) :RefPath(prefix),
		pac_buf(0), rpac_buf(0), bwt_buf(0), rbwt_buf(0), bns(0), l_buf(0), m_pac_buf(
				0), m_rpac_buf(0), m_bwt_buf(0), m_rbwt_buf(0), bwt_d(0), rbwt_d(
				0)
{
#ifdef BLOOM_FPP
	//***bloom filter initialization

	// How many elements roughly do we expect to insert?
	parameters.projected_element_count = VCF_SITES *( WINDOW_SIZE-KMER_SIZE+1+KMER_SIZE) * 2;

	// Maximum tolerable false positive probability? (0,1)
	parameters.false_positive_probability = BLOOM_FPP;// 1 in 10000

	// Simple randomizer (optional)
	parameters.random_seed = 0xA5A5A5A5;

	if (!parameters)
	{
		std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
		exit(1);
	}
	parameters.compute_optimal_parameters();

	filter.set_param(parameters);

	//***bloom filter initialization end
#endif
	InitializeRollHashTable();
}

BwtIndexer::BwtIndexer(string & prefix) :RefPath(prefix),
		pac_buf(0), rpac_buf(0), bwt_buf(0), rbwt_buf(0), bns(0), l_buf(0), m_pac_buf(
				0), m_rpac_buf(0), m_bwt_buf(0), m_rbwt_buf(0), bwt_d(0), rbwt_d(
				0)
{
#ifdef BLOOM_FPP
	//***bloom filter initialization

	// How many elements roughly do we expect to insert?
	parameters.projected_element_count = VCF_SITES *( WINDOW_SIZE-KMER_SIZE+1+KMER_SIZE) * 2;

	// Maximum tolerable false positive probability? (0,1)
	parameters.false_positive_probability = BLOOM_FPP;// 1 in 10000

	// Simple randomizer (optional)
	parameters.random_seed = 0xA5A5A5A5;

	if (!parameters)
	{
		std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
		exit(1);
	}
	parameters.compute_optimal_parameters();

	filter.set_param(parameters);

	//***bloom filter initialization end
#endif
	InitializeRollHashTable();
}

/******************************/
//hash related functions
bool BwtIndexer::IsKmerInHash(uint64_t kmer)const
{
	int counter(0);
	uint32_t shrinked(0);
	for (int iter = 0; iter != 6; ++iter)
	{
		//printf("the DATUM is : %016llx    masked DATUM:%x    the hash value is :%d as well as mask:%x\n",kmer,KmerShrinkage(kmer, mask[iter]),roll_hash_table[iter][KmerShrinkage(kmer, mask[iter])], mask[iter]);
		shrinked = KmerShrinkage(kmer, mask[iter]);
		if (roll_hash_table[iter][shrinked / 8] & (1 << (shrinked % 8)))
			counter++;
	}
	if (counter > RollParam.thresh)
		return true;
	else
		return false;
}
bool BwtIndexer::IsReadInHash(const ubyte_t * S, int len)const
{
	uint64_t kmer1(0), kmer2(0), kmer3(0);

	for (int i = 0; i != 32; ++i)
	{
		kmer1 = ((kmer1 << 2) | S[0 + i]);
		kmer2 = ((kmer2 << 2) | S[32 + i]);
		kmer3 = ((kmer3 << 2) | S[64 + i]);
	}
	if (IsKmerInHash(kmer1) || IsKmerInHash(kmer2) || IsKmerInHash(kmer3))
		return true;
	else
		return false;

}
bool BwtIndexer::IsReadFiltered(const ubyte_t * S, const ubyte_t * Q, int len)const
{
	if (IsReadHighQ(Q, len)) //pass error
	{
		if (IsReadInHash(S, len))
		{
			return false;
		}
		else
			return true;
	}
	else
	{
		//fprintf(stderr,"low qual\n");
		return false;
	}
}
std::string BwtIndexer::ReverseComplement(const std::string & a)
{

	std::string b(a.rbegin(), a.rend());

	for(uint32_t i = 0; i != b.size(); ++i)
	{
		b[i] = match_table[b[i]];
	}
	return b;
}
uint32_t BwtIndexer::KmerShrinkage(uint64_t a, unsigned int mask)const
{
	uint32_t tmp(0), i(1);
	//printf("the incoming kmer:%016llx\t",a);
	while (mask) //32 circle
	{
		if (mask & 1)
		{
			tmp |= ((a & 3) << (2 * i - 2));
			++i;
		}

		mask >>= 1;
		a >>= 2;
	}
	//printf("the resulted i is :%d\n",i);
	return tmp;
}
void BwtIndexer::InitializeRollHashTable()
{
	RollParam.kmer_size = KMER_SIZE;
	RollParam.read_step_size = READ_STEP_SIZE;
	RollParam.thresh = 3;
	clock_t t = clock();

	hash_table_size = pow(4, 16) / 8;
	for (int i = 0; i != 6; ++i)
		roll_hash_table[i] = (unsigned char *) calloc(hash_table_size,
				sizeof(char));

	notice("Initializing Rolling Hash Table with size: %d bytes in %d sec\n",hash_table_size, (float) (clock() - t) / CLOCKS_PER_SEC);

	assert(roll_hash_table != 0);
}

void BwtIndexer::ReadRollHashTable(const std::string& prefix)
{
	std::string infile = prefix + ".rollhash";
	ifstream fin(infile, std::ofstream::binary);
	if (!fin.is_open())
	{
		error("Open %s failed!\n",infile.c_str());
		exit(1);
	}
	for (int i = 0; i != 6; ++i)
		fin.read(reinterpret_cast<char*>(roll_hash_table[i]), hash_table_size);
	fin.close();
}

void BwtIndexer::DumpRollHashTable(const std::string& prefix)
{
	std::string outfile = prefix + ".rollhash";
	ofstream fout(outfile, std::ofstream::binary);
	if (!fout.is_open())
	{
		error("Open %s failed!\n",outfile.c_str());
		exit(1);
	}
	for (int i = 0; i != 6; ++i)
		fout.write(reinterpret_cast<char*>(roll_hash_table[i]),
				hash_table_size);
	fout.close();
}

void BwtIndexer::DestroyRollHashTable()
{
	for (int i = 0; i != 6; ++i)
		free(roll_hash_table[i]);
}
void BwtIndexer::AddSeq2HashCore(const std::string & Seq, int iter)
{

	uint64_t datum(0);
	unsigned int i = 0;
	uint32_t shrinked(0);
	//std::string tmp=Seq.substr(0,KMER_SIZE);

	for (; i != 32; ++i)
	{
		datum = ((datum << 2) | nst_nt4_table[Seq[i]]);/*
		 & LOMEGA(
		 min_only(RollParam.kmer_size,OVERFLOWED_KMER_SIZE)); */ // N not considered
	}
	shrinked = KmerShrinkage(datum, mask[iter]);
	roll_hash_table[iter][shrinked / 8] |= (1 << (shrinked % 8));

	for (; i != Seq.length() / 2/*last one considered*/; ++i)
	{
		//std::cerr<<"number:"<<i<<"COmparing length:"<<KMER_SIZE<<"and"<<Seq.length()<<"while max:"<<Seq.max_size()<<std::endl;
		datum = ((datum << 2) | nst_nt4_table[Seq[i]]);/*
		 & LOMEGA(
		 min_only(RollParam.kmer_size,OVERFLOWED_KMER_SIZE));*/ // N not considered
		//roll_hash_table[iter][KmerShrinkage(datum, mask[iter])]++;
		shrinked = KmerShrinkage(datum, mask[iter]);
		roll_hash_table[iter][shrinked / 8] |= (1 << (shrinked % 8));
	}
	//printf("the DATUM is : %016llx    masked DATUM:%x    the hash value is :%d as well as mask:%x\n",datum,KmerShrinkage(datum, mask[iter]),roll_hash_table[iter][KmerShrinkage(datum, mask[iter])/8]&(1<<(shrinked%8)), mask[iter]);

	uint64_t tmp = datum;
	for (unsigned int j = i, let = 0; let != 4; ++let, j = i)
	{
		tmp = datum;
		for (; j != Seq.length() / 2 + 32; ++j)
		{
			tmp = ((tmp << 2) | nst_nt4_table["ACGT"[(char) let]]);
			shrinked = KmerShrinkage(tmp, mask[iter]);
			roll_hash_table[iter][shrinked / 8] |= (1 << (shrinked % 8));
		}
	}
	datum = tmp;
	for (i = Seq.length() / 2 + 32; i != Seq.length()/*last one considered*/;
			++i)
	{
		//std::cerr<<"number:"<<i<<"COmparing length:"<<KMER_SIZE<<"and"<<Seq.length()<<"while max:"<<Seq.max_size()<<std::endl;
		datum = ((datum << 2) | nst_nt4_table[Seq[i]]);/*
		 & LOMEGA(
		 min_only(RollParam.kmer_size,OVERFLOWED_KMER_SIZE));*/ // N not considered
		//roll_hash_table[iter][KmerShrinkage(datum, mask[iter])]++;
		shrinked = KmerShrinkage(datum, mask[iter]);
		roll_hash_table[iter][shrinked / 8] |= (1 << (shrinked % 8));
	}
	//printf("the DATUM is : %x    the hash value is :%d as well as:%x\n",datum,roll_hash_table[datum], LOMEGA(min_only(RollParam.kmer_size,OVERFLOWED_KMER_SIZE)));
}
/******************************/
bool BwtIndexer::BuildIndex(RefBuilder & ArtiRef, string & prefix,
		const gap_opt_t * opt)
{
	RefPath=prefix;
	string str;
	Fa2Pac(ArtiRef, prefix.c_str(), opt); //dump .pac
	Fa2RevPac(prefix.c_str()); //dump .rpac

	DumpRollHashTable(prefix);
	//str=prefix+".pac";
	//std::cerr<<"The bwa seq len is:"<<bwa_seq_len(str.c_str())<<std::endl;
	bns_dump(bns, prefix.c_str());
	notice(" Building Pac from Bwt...\n");
	bwt_d = Pac2Bwt(pac_buf);
	//cerr<<"Pac2Bwt..."<<bwt_d->bwt_size<<endl;
	rbwt_d = Pac2Bwt(rpac_buf);
	bwt_gen_cnt_table(bwt_d);
	bwt_gen_cnt_table(rbwt_d);

	bwt_bwtupdate_core(bwt_d);
	//cerr<<"Bwt update..."<<bwt_d->bwt_size<<endl;
	bwt_bwtupdate_core(rbwt_d);
	notice("Dumping Bwt...\n");
	str = prefix + ".bwt";
	bwt_dump_bwt(str.c_str(), bwt_d);
	str = prefix + ".rbwt";
	//cerr<<"The Mapper hs rbwt_d->bwt_size is:"<<rbwt_d->bwt_size<<endl;
	bwt_dump_bwt(str.c_str(), rbwt_d);

	//bwt_d=bwt_gen_cnt_table(bwt_d);
	//rbwt_d=bwt_gen_cnt_table(rbwt_d);
	notice("Calculate SA from Bwt...\n");
	str = prefix + ".sa";
	bwt_cal_sa(bwt_d, 32);
	bwt_dump_sa(str.c_str(), bwt_d);
	str = prefix + ".rsa";
	bwt_cal_sa(rbwt_d, 32);
	bwt_dump_sa(str.c_str(), rbwt_d);

	DBG(fprintf(stderr,"Indexer initialization finished...\n");)
	return true;
}
bool BwtIndexer::LoadIndex(string & prefix)
{
	RefPath=prefix;
	ReadRollHashTable(prefix);
	string str = prefix + ".bwt";
	bwt_d = bwt_restore_bwt(str.c_str());
	str = prefix + ".sa";
	bwt_restore_sa(str.c_str(), bwt_d);
	//bwt_cal_sa(bwt_d, 32);
	str = prefix + ".rbwt";
	if (stat(str.c_str(), &sb) == 0)
		rbwt_d = bwt_restore_bwt(str.c_str());
	else
	{
		str = prefix + ".pac";
		string str2 = prefix + ".rpac";
		bwa_pac_rev_core(str.c_str(), str2.c_str());
		rbwt_d = bwt_pac2bwt(str2.c_str(), 3);
		bwt_gen_cnt_table(rbwt_d);
		bwt_bwtupdate_core(rbwt_d);
		str = prefix + ".rbwt";
		//cerr<<"The Mapper rbwt_d->bwt_size is:"<<rbwt_d->bwt_size<<endl;
		bwt_dump_bwt(str.c_str(), rbwt_d);
	}
	str = prefix + ".rsa";
	bwt_restore_sa(str.c_str(), rbwt_d);
	//bwt_cal_sa(rbwt_d, 32);
	bns = bns_restore(prefix.c_str());
	DBG(fprintf(stderr,"Indexer loading  finished...\n");)
	return true;
}

bool BwtIndexer::Fa2Pac(RefBuilder & ArtiRef, const char *prefix,
		const gap_opt_t* opt)
{
	// initialization
	//seq = kseq_init(fp_fa);
	char name[1024];
	int32_t m_seqs, m_holes;
	uint32_t i;
	bntamb1_t *q;
	FILE *fp;

	bns = (bntseq_t*) calloc(1, sizeof(bntseq_t));
	bns->seed = 11; // fixed seed for random generator
	srand48(bns->seed);
	m_seqs = m_holes = 8;
	bns->anns = (bntann1_t*) calloc(m_seqs, sizeof(bntann1_t));
	bns->ambs = (bntamb1_t*) calloc(m_holes, sizeof(bntamb1_t));
	q = bns->ambs;
	l_buf = 0; //base counter
	m_pac_buf = 0x40000; //limit of number of reads
	int32_t size_pac_buf = 0x10000;
	pac_buf = (unsigned char *) calloc(size_pac_buf, sizeof(unsigned char));
	strcpy(name, prefix);
	strcat(name, ".pac");
	fp = xopen(name, "wb");
	// read sequences
	//while ((l = kseq_read(seq)) >= 0) {
	DBG(fprintf(stderr,"NOTE:Come into Fa2Pac...\n");)
	string RefOutput(prefix);
	RefOutput+=".ref.fa";
	ofstream Fout(RefOutput);
	for (unordered_map<string, uint32_t>::iterator iter =
			ArtiRef.RefTableIndex.begin(); iter != ArtiRef.RefTableIndex.end();
			++iter)
	{

		string CurrentSeqName(iter->first.c_str());
		//if(ArtiRef.longRefTable[CurrentSeqName]==true) continue;// long ref seq would not be indexed
		string CurrentSeq = ArtiRef.SeqVec[iter->second];

		AddSeq2Hash(CurrentSeq, opt);
		AddSeq2Hash(string(CurrentSeq.rbegin(), CurrentSeq.rend()), opt);
		string tmp_rev_cmp = ReverseComplement(CurrentSeq);
		AddSeq2Hash(tmp_rev_cmp, opt);
		AddSeq2Hash(string(tmp_rev_cmp.rbegin(), tmp_rev_cmp.rend()), opt);
		Fout << ">" << CurrentSeqName << "\n" << CurrentSeq << "\n";

		//DBG(fprintf(stderr,"Come into outer loop...%s\n%s\n",CurrentSeqName.c_str(),CurrentSeq.c_str());)
		//DBG(fprintf(stderr,"We got in seqs : %d ...\n%s\n", vec_iter,CurrentSeq.c_str());)
		bntann1_t *p; // tmp ann object
		int lasts;
		if (bns->n_seqs == m_seqs)
		{
			m_seqs <<= 1; // double the space to store number of sequences
			bns->anns = (bntann1_t*) realloc(bns->anns,
					m_seqs * sizeof(bntann1_t));
		}
		p = bns->anns + bns->n_seqs; // next sequence
		//p->name = strdup((char*)seq->name.s);
		p->name = strdup(const_cast<char*>(CurrentSeqName.c_str()));
		//p->anno = seq->comment.s? strdup((char*)seq->comment.s) : strdup("(null)");
		p->anno = strdup("(null)");
		p->gi = 0;
		p->len = CurrentSeq.length(); //length of the sequence
		p->offset = (bns->n_seqs == 0) ? 0 : (p - 1)->offset + (p - 1)->len; // start position of new sequence
		p->n_ambs = 0;
		for (i = 0, lasts = 0; i < CurrentSeq.length(); ++i)
		{
			int c = nst_nt4_table[(int) CurrentSeq[i]];
			if (c >= 4)
			{ // N
				if (lasts == CurrentSeq[i])
				{ // contiguous N
					++q->len;
				}
				else
				{
					if (bns->n_holes == m_holes)
					{
						m_holes <<= 1;
						bns->ambs = (bntamb1_t*) realloc(bns->ambs,
								m_holes * sizeof(bntamb1_t));
					}
					q = bns->ambs + bns->n_holes;
					q->len = 1;
					q->offset = p->offset + i;
					q->amb = CurrentSeq[i];
					++p->n_ambs;
					++bns->n_holes;
				}
			}
			lasts = CurrentSeq[i];
			{ // fill buffer
				if (c >= 4)
					c = lrand48() & 0x3;
				//if (l_buf == 0x40000) {
				if (l_buf == m_pac_buf)
				{
					//fwrite(buf, 1, 0x10000, fp);
					//fwrite(pac_buf, 1, 0x10000, fp);
					//memset(buf, 0, 0x10000);
					//l_buf = 0;
					m_pac_buf += 0x40000;
					int32_t tmp_size_pac_buf = size_pac_buf + 0x10000;
					ubyte_t * tmp_pac_buf = (unsigned char *) calloc(
							tmp_size_pac_buf, sizeof(unsigned char));
					memcpy(tmp_pac_buf, pac_buf, size_pac_buf);
					free(pac_buf);
					pac_buf = tmp_pac_buf;
					size_pac_buf = tmp_size_pac_buf;

					//pac_buf = (unsigned char*)realloc(pac_buf, size_pac_buf * sizeof(char));
				}
				pac_buf[l_buf >> 2] |= c << ((3 - (l_buf & 3)) << 1);
				++l_buf;
			}
		}
		++bns->n_seqs;
		bns->l_pac += CurrentSeq.length();
		//}//inner for-loop
	} //outer for-loop
	  //clean tmp variables
	Fout.close();
	xassert(bns->l_pac, "zero length sequence.");
	{ // finalize .pac file
		ubyte_t ct;
		fwrite(pac_buf, 1, (l_buf >> 2) + ((l_buf & 3) == 0 ? 0 : 1), fp);
		// the following codes make the pac file size always (l_pac/4+1+1)
		if (bns->l_pac % 4 == 0)
		{
			ct = 0;
			fwrite(&ct, 1, 1, fp);
		}
		ct = bns->l_pac % 4;
		fwrite(&ct, 1, 1, fp);
		// close .pac file
		fclose(fp);
	}
	bns->to_pac_buf = pac_buf;
	//DBG(fprintf(stderr,"l_pac: %d\nbns->pac_buf: %x\n",bns->l_pac,bns->to_pac_buf);)
	l_buf = 0;
	return 0;
}
bool BwtIndexer::Fa2RevPac(const char * prefix)
{//TODO: check usage of j
	int64_t seq_len, i;
	bwtint_t pac_len, j;
	ubyte_t *bufin, *bufout, ct;
	FILE *fp;
	char name[1024];
	strcpy(name, prefix);
	strcat(name, ".rpac");
	seq_len = bns->l_pac;
	pac_len = (seq_len >> 2) + 1;
	bufin = pac_buf;
	bufout = rpac_buf = (ubyte_t*) calloc(pac_len, 1);
	//fp = xopen(fn, "rb");
	//fread(bufin, 1, pac_len, fp);
	//fclose(fp);
	for (i = seq_len - 1, j = 0; i >= 0; --i)
	{
		int c = bufin[i >> 2] >> ((~i & 3) << 1) & 3;
		bwtint_t j = seq_len - 1 - i;
		bufout[j >> 2] |= c << ((~j & 3) << 1);
	}
	//free(bufin);
	fp = xopen(name, "wb");
	fwrite(bufout, 1, pac_len, fp);
	ct = seq_len % 4;
	fwrite(&ct, 1, 1, fp);
	fclose(fp);
	//free(bufout);
	return 0;
}

bwt_t* BwtIndexer::Pac2Bwt(unsigned char * pac)
{
	bwt_t *bwt;
	ubyte_t *buf, *buf2;
	uint32_t i; //, pac_size;
	//FILE *fp;
	//DBG(fprintf(stderr,"checkpoint 1....%d\n",bns->l_pac);)
	// initialization
	bwt = (bwt_t*) calloc(1, sizeof(bwt_t));
	bwt->seq_len = bns->l_pac;

	bwt->bwt_size = (bwt->seq_len + 15) >> 4;
	//fp = xopen(fn_pac, "rb");

	// prepare sequence
	//pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
	//buf2 = (ubyte_t*)calloc(pac_size, 1);
	//fread(buf2, 1, pac_size, fp);
	//fclose(fp);
	buf2 = pac;
	memset(bwt->L2, 0, 5 * 4);
	buf = (ubyte_t*) calloc(bwt->seq_len + 1, 1);
	//DBG(fprintf(stderr,"checkpoint 2-1....\n");)
	for (i = 0; i < bwt->seq_len; ++i)
	{
		buf[i] = buf2[i >> 2] >> ((3 - (i & 3)) << 1) & 3;
		//	DBG(fprintf(stderr,"checkpoint 3-1....buf[i]= %d\n",buf[i]);)
		++bwt->L2[1 + buf[i]]; // Number of specific allele plus
	}
	//DBG(fprintf(stderr,"checkpoint 3-1....\n");)
	for (i = 2; i <= 4; ++i)
		bwt->L2[i] += bwt->L2[i - 1]; //cumulative
	//free(buf2);
	//DBG(fprintf(stderr,"checkpoint 3....\n");)
	// Burrows-Wheeler Transform
//		if (use_is) {
	bwt->primary = is_bwt(buf, bwt->seq_len);
	/*		} else {
	 #ifdef _DIVBWT
	 bwt->primary = divbwt(buf, buf, 0, bwt->seq_len);
	 #else
	 err_fatal_simple("libdivsufsort is not compiled in.");
	 #endif
	 }*/

	bwt->bwt = (u_int32_t*) calloc(bwt->bwt_size, 4);
	for (i = 0; i < bwt->seq_len; ++i)
		bwt->bwt[i >> 4] |= buf[i] << ((15 - (i & 15)) << 1);
	free(buf);
	return bwt;
}

#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)
void BwtIndexer::bwt_bwtupdate_core(bwt_t *bwt)
{
	bwtint_t i, k, c[4], n_occ;
	uint32_t *buf;

	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	bwt->bwt_size += n_occ * 4; // the new size
	buf = (uint32_t*) calloc(bwt->bwt_size, 4); // will be the new bwt
	c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = 0; i < bwt->seq_len; ++i)
	{
		if (i % OCC_INTERVAL == 0)
		{
			memcpy(buf + k, c, sizeof(bwtint_t) * 4);
			k += 4;
		}
		if (i % 16 == 0)
			buf[k++] = bwt->bwt[i / 16];
		++c[bwt_B00(bwt, i)];
	}
	// the last element
	memcpy(buf + k, c, sizeof(bwtint_t) * 4);
	xassert(k + 4 == bwt->bwt_size, "inconsistent bwt_size");
	// update bwt
	free(bwt->bwt);
	bwt->bwt = buf;
}

void BwtIndexer::bwt_cal_sa(bwt_t *bwt, int intv)
{
	bwtint_t isa, sa, i; // S(isa) = sa

	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

	if (bwt->sa)
		free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) / intv;
	bwt->sa = (bwtint_t*) calloc(bwt->n_sa, sizeof(bwtint_t));
	// calculate SA value
	isa = 0;
	sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i)
	{
		if (isa % intv == 0)
			bwt->sa[isa / intv] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if (isa % intv == 0)
		bwt->sa[isa / intv] = sa;
	bwt->sa[0] = (bwtint_t) -1; // before this line, bwt->sa[0] = bwt->seq_len
}

BwtIndexer::~BwtIndexer()
{
	// TODO Auto-generated destructor stub
	DBG(fprintf(stderr,"Incoming destroyer...\n");)
	free(bwt_d->sa);
	free(bwt_d->bwt);
	free(bwt_d);
	DBG(fprintf(stderr,"bwt_d delete successfully...\n");)
	free(rbwt_d->sa);
	free(rbwt_d->bwt);
	free(rbwt_d);
	DBG(fprintf(stderr,"rbwt_d delete successfully...\n");)
	if (pac_buf)
		free(pac_buf);
	if (rpac_buf)
		free(rpac_buf);
	if (bwt_buf)
		free(bwt_buf);
	if (rbwt_buf)
		free(rbwt_buf);
	DBG(fprintf(stderr,"buf delete successfully...\n");)
	if (bns == 0)
		return;
	else
	{
		int i;
		if (bns->fp_pac)
			fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i)
		{
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
	DestroyRollHashTable();
}

