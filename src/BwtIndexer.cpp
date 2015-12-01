#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <iostream>
#include <sstream>
#include "BwtIndexer.h"
#include "RefBuilder.h"
#include "stdint.h"
#include "stdio.h"
#include "../libbwa/bwt.h"
#include "../libbwa/utils.h"
#include "../libbwa/bntseq.h"
#include "../libbwa/kseq.h"
#include "../misc/general/Error.h"
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

#define KMER_64_BIT 1
#ifdef KMER_8_BIT
unsigned int mask[6] =
{ 240, //11110000
15,//00001111
195,//11000011
60,//00111100
204,//11001100
51 }; //00110011
#endif
#ifdef KMER_16_BIT
unsigned int mask[6] =
{ 65280, //1111111100000000
255,//0000000011111111
61455,//1111000000001111
4080,//00001111111110000
12528,//1111000011110000
3855 }; //0000111100001111
#endif
#ifdef KMER_32_BIT
unsigned int mask[6] =
{ 0xffff0000, 0xffff, 0xff0000ff, 0xffff00, 0xff00ff00, 0xff00ff };
#endif
#ifdef KMER_64_BIT
const uint64_t mask[6] =
{ 0xffffffff00000000, 0xffffffff, 0xffff00000000ffff, 0xffffffff0000, 0xffff0000ffff0000, 0xffff0000ffff };
#endif

#define NST_NT4_TABLE(X) (nst_nt4_table[(unsigned char)X]<4?nst_nt4_table[(unsigned char)X]:rand()%4)


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

	InitializeRollHashTable(3);

}

BwtIndexer::BwtIndexer(int thresh) :RefPath(),
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

	InitializeRollHashTable(thresh);

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
	InitializeRollHashTable(3);
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
	InitializeRollHashTable(3);
}

/******************************/
//hash related functions
inline bool BwtIndexer::IsKmerInHash(uint64_t kmer)const
{
	int counter(0);
	uint64_t shrinked(0);
	//for (int iter = 0; iter != 6; ++iter)
	//{
	//	//printf("the DATUM is : %016llx    masked DATUM:%x    the hash value is :%d as well as mask:%x\n",kmer,KmerShrinkage(kmer, iter),roll_hash_table[iter][KmerShrinkage(kmer, iter)], iter);
	//	uint32_t shrinked = KmerShrinkage(kmer, iter);
	//	if (roll_hash_table[iter][shrinked / 8] & (1 << (shrinked % 8)))
	//		counter++;
	//}

	/*0xffffffff00000000*/
	shrinked |= (kmer & 0xffffffff00000000) >> 32;
	//fprintf(stderr, "After the mask 0xffffffff00000000 of %016llx \n", shrinked);
	if (roll_hash_table[0][shrinked / 8] & (1 << (shrinked % 8))) { counter++; }
	/*0xffffffff*/
	shrinked = 0;
	shrinked |= kmer & 0xffffffff;
	//fprintf(stderr, "After the mask 0xffffffff of %016llx \n", shrinked);
	if (roll_hash_table[1][shrinked / 8] & (1 << (shrinked % 8))) { counter++;}
	/*0xffff00000000ffff*/
	shrinked = 0;
	shrinked |= (kmer & 0xffff000000000000) >> 32;
	shrinked |= (kmer & 0xffff);
	//fprintf(stderr, "After the mask 0xffff00000000ffff of %016llx \n", shrinked);
	if (roll_hash_table[2][shrinked / 8] & (1 << (shrinked % 8))) { counter++; }
	/*0xffffffff0000*/
	shrinked = 0;
	shrinked |= (kmer & 0xffffffff0000)>>16;
	//fprintf(stderr, "After the mask 0xffffffff0000 of %016llx \n", shrinked);
	if (roll_hash_table[3][shrinked / 8] & (1 << (shrinked % 8))) { counter++;  }
	/*0xffff0000ffff0000*/
	shrinked = 0;
	shrinked |= (kmer & 0xffff000000000000) >> 32;
	shrinked |= (kmer & 0xffff0000)>>16;
	//fprintf(stderr, "After the mask 0xffff0000ffff0000 of %016llx \n", shrinked);
	if (roll_hash_table[4][shrinked / 8] & (1 << (shrinked % 8))) { counter++;  }
	/*0xffff0000ffff*/
	shrinked = 0;
	shrinked |= (kmer & 0xffff00000000) >> 16;
	shrinked |= (kmer & 0xffff);
	//fprintf(stderr, "After the mask 0xffff0000ffff of %016llx \n", shrinked);
	if (roll_hash_table[5][shrinked / 8] & (1 << (shrinked % 8))) { counter++; }
	//fprintf(stderr,"The counter of %016llx is %d\n",kmer,counter);
	if (counter >= RollParam.thresh)
		return true;
	else
		return false;
}
inline int BwtIndexer::CountKmerHitInHash(uint64_t kmer)const
{
	int counter(0);
	uint64_t shrinked(0);
	/*0xffffffff00000000*/
	shrinked |= (kmer & 0xffffffff00000000) >> 32;
	//fprintf(stderr, "After the mask 0xffffffff00000000 of %016llx \n", shrinked);
	if (roll_hash_table[0][shrinked / 8] & (1 << (shrinked % 8))) { counter++; }
	/*0xffffffff*/
	shrinked = 0;
	shrinked |= kmer & 0xffffffff;
	//fprintf(stderr, "After the mask 0xffffffff of %016llx \n", shrinked);
	if (roll_hash_table[1][shrinked / 8] & (1 << (shrinked % 8))) { counter++;}
	/*0xffff00000000ffff*/
	shrinked = 0;
	shrinked |= (kmer & 0xffff000000000000) >> 32;
	shrinked |= (kmer & 0xffff);
	//fprintf(stderr, "After the mask 0xffff00000000ffff of %016llx \n", shrinked);
	if (roll_hash_table[2][shrinked / 8] & (1 << (shrinked % 8))) { counter++; }
	/*0xffffffff0000*/
	shrinked = 0;
	shrinked |= (kmer & 0xffffffff0000)>>16;
	//fprintf(stderr, "After the mask 0xffffffff0000 of %016llx \n", shrinked);
	if (roll_hash_table[3][shrinked / 8] & (1 << (shrinked % 8))) { counter++;  }
	/*0xffff0000ffff0000*/
	shrinked = 0;
	shrinked |= (kmer & 0xffff000000000000) >> 32;
	shrinked |= (kmer & 0xffff0000)>>16;
	//fprintf(stderr, "After the mask 0xffff0000ffff0000 of %016llx \n", shrinked);
	if (roll_hash_table[4][shrinked / 8] & (1 << (shrinked % 8))) { counter++;  }
	/*0xffff0000ffff*/
	shrinked = 0;
	shrinked |= (kmer & 0xffff00000000) >> 16;
	shrinked |= (kmer & 0xffff);
	//fprintf(stderr, "After the mask 0xffff0000ffff of %016llx \n", shrinked);
	if (roll_hash_table[5][shrinked / 8] & (1 << (shrinked % 8))) { counter++; }
	//fprintf(stderr,"The counter of %016llx is %d\n",kmer,counter);
	return counter;
}
uint64_t swap_uint64(uint64_t val)
{
	val = ((val << 8) & 0xFF00FF00FF00FF00ULL) | ((val >> 8) & 0x00FF00FF00FF00FFULL);
	val = ((val << 16) & 0xFFFF0000FFFF0000ULL) | ((val >> 16) & 0x0000FFFF0000FFFFULL);
	return (val << 32) | (val >> 32);
}
typedef uint8_t v16qi __attribute__((vector_size(16)));
//typedef uint64_t v2li __attribute__((vector_size(16)));
bool BwtIndexer::IsReadInHash(const ubyte_t * S, int len, bool more_chunck)const
{
	//	for (int i = 0; i != len;++i)
	//	fprintf(stderr, "%c\t", "ACGT"[S[i]]);
	//	fprintf(stderr, "\n");
	//	double enter_time=realtime();
		//for debug use
		int n_chunk=len/32;
		uint64_t *kmer=new uint64_t [n_chunk];
		for (int i=0;i!=n_chunk;++i)
		{
			kmer[i]=0;
			for (int j = 0; j != 32; ++j)
			{
				kmer[i] = ((kmer[i] << 2) | S[32*i + j]);
			}
			if(IsKmerInHash(kmer[i]))
			{
				return true;
			}
		}

		return false;

}
bool BwtIndexer::IsReadInHash(const ubyte_t * S, int len)const
{	
	uint64_t kmer3[2];
	v16qi mmx1 = { S[0], S[4], S[8], S[12], S[16], S[20], S[24], S[28], S[32], S[36], S[40], S[44], S[48], S[52], S[56], S[60] };
	v16qi mmx2 = { S[1], S[5], S[9], S[13], S[17], S[21], S[25], S[29], S[33], S[37], S[41], S[45], S[49], S[53], S[57], S[61] };
	v16qi mmx3 = { S[2], S[6], S[10], S[14], S[18], S[22], S[26], S[30], S[34], S[38], S[42], S[46], S[50], S[54], S[58], S[62] };
	v16qi mmx4 = { S[3], S[7], S[11], S[15], S[19], S[23], S[27], S[31], S[35], S[39], S[43], S[47], S[51], S[55], S[59], S[63] };
	v16qi mmx5 = { S[64], S[68], S[72], S[76], S[80], S[84], S[88], S[92], 0, 0, 0, 0, 0, 0, 0, 0 };
	v16qi mmx6 = { S[65], S[69], S[73], S[77], S[81], S[85], S[89], S[93], 0, 0, 0, 0, 0, 0, 0, 0 };
	v16qi mmx7 = { S[66], S[70], S[74], S[78], S[82], S[86], S[90], S[94], 0, 0, 0, 0, 0, 0, 0, 0 };
	v16qi mmx8 = { S[67], S[71], S[75], S[79], S[83], S[87], S[91], S[95], 0, 0, 0, 0, 0, 0, 0, 0 };


#define GCC_VERSION (__GNUC__ * 10000 \
               + __GNUC_MINOR__ * 100 \
                + __GNUC_PATCHLEVEL__)	
#if defined(GCC_VERSION)&& GCC_VERSION >= 40700//4.7#endif
//fprintf(stderr,"GCC_VERSION is about 4.7");
	mmx1 = mmx1 << 6;
	mmx1 = mmx1 | (mmx2 << 4);
	mmx1 = mmx1 | (mmx3 << 2);
#else
	v16qi mul64 = {64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64};
	v16qi mul16 = {16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16};
	v16qi mul4 = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
	mmx1 = mmx1 * mul64;
	mmx1 = mmx1 | (mmx2 * mul16);
	mmx1 = mmx1 | (mmx3 * mul4);
#endif

	mmx1 = mmx1 | mmx4;

#if defined(GCC_VERSION)&& GCC_VERSION > 40700//4.7#endif
	mmx5 = mmx5 << 6;
	mmx5 = mmx5 | (mmx6 << 4);
	mmx5 = mmx5 | (mmx7 << 2);
#else
	mmx5 = mmx5 * mul64;
	mmx5 = mmx5 | (mmx6 *mul16);
	mmx5 = mmx5 | (mmx7 *mul4);
#endif

	mmx5 = mmx5 | mmx8;
	memcpy(&kmer3,&mmx1,16);

	//for (int i = 0; i != 32; ++i)
	//	fprintf(stderr, "%c\t", "ACGT"[(swap_uint64(kmer3[0])>> ((31 - i) * 2)) & 3]);
	////fprintf(stderr, "\n");
	//for (int i = 0; i != 32; ++i)
	//	fprintf(stderr, "%c\t", "ACGT"[(swap_uint64(kmer3[1]) >> ((31 - i) * 2)) & 3]);
	//memcpy(&kmer3, &mmx5, 16);
	//for (int i = 0; i != 32; ++i)
	//	fprintf(stderr, "%c\t", "ACGT"[(swap_uint64(kmer3[0]) >> ((31 - i) * 2)) & 3]);
	//fprintf(stderr, "\n");
	//memcpy(&kmer3, &mmx1, 16);
	//fprintf(stderr, "0-31: %016llx\t", swap_uint64(kmer3[0]));
	//fprintf(stderr, "32-63: %016llx\n", swap_uint64(kmer3[1]));
		//fprintf(stderr, "%d\t", (mmx1[0] >> 2) & 3);
	//	//fprintf(stderr, "%d\t", mmx1[0]  & 3);
	//	//fprintf(stderr, "%d\t", (mmx1[1] >> 6) & 3);
	//	//fprintf(stderr, "%d\t", (mmx1[1] >> 4) & 3);
	//	//fprintf(stderr, "%d\t", (mmx1[1] >> 2) & 3);
	//	//fprintf(stderr, "%d\t", mmx1[1] & 3);
	//	//fprintf(stderr, "%d\t", (mmx1[2] >> 6) & 3);
	//	//fprintf(stderr, "%d\t", (mmx1[2] >> 4) & 3);
	//	//fprintf(stderr, "%d\t", (mmx1[2] >> 2) & 3);
	//	//fprintf(stderr, "%d\t", mmx1[2] & 3);
	////for (int i = 0; i != 32; ++i)
	////	fprintf(stderr, "%d\t", (kmer[0][1] >> ((31 - i) * 2)) & 3);
	////for (int i = 0; i != 32; ++i)
	//////	fprintf(stderr, "%d\t", (kmer[1][0] >> ((31 - i) * 2)) & 3);
	////fprintf(stderr, "\n");
	////exit(EXIT_FAILURE);
	if (IsKmerInHash(swap_uint64(kmer3[0])) || IsKmerInHash(swap_uint64(kmer3[1])))
	{

		return true;
	}
	else
	{
		memcpy(&kmer3, &mmx5, 16);
		if (IsKmerInHash(swap_uint64(kmer3[0])))
		{

			return true;
		}
		else
			return false;
	}


}
bool BwtIndexer::IsReadInHashByCount(const ubyte_t *S, int len, bool more_chunck)const
{
	int n_chunk=len/32;
	int count=0;
	uint64_t *kmer=new uint64_t [n_chunk];
	for (int i=0;i!=n_chunk;++i)
	{
		kmer[i]=0;
		for (int j = 0; j != 32; ++j)
		{
			kmer[i] = ((kmer[i] << 2) | S[32*i + j]);
		}
		count+=IsKmerInHash(kmer[i]);
	}
	if (count >= RollParam.thresh)
		return true;
	else
		return false;
}
bool BwtIndexer::IsReadInHashByCount(const ubyte_t * S, int len)const
{
	uint64_t kmer3[2];
	v16qi mmx1 = { S[0], S[4], S[8], S[12], S[16], S[20], S[24], S[28], S[32], S[36], S[40], S[44], S[48], S[52], S[56], S[60] };
	v16qi mmx2 = { S[1], S[5], S[9], S[13], S[17], S[21], S[25], S[29], S[33], S[37], S[41], S[45], S[49], S[53], S[57], S[61] };
	v16qi mmx3 = { S[2], S[6], S[10], S[14], S[18], S[22], S[26], S[30], S[34], S[38], S[42], S[46], S[50], S[54], S[58], S[62] };
	v16qi mmx4 = { S[3], S[7], S[11], S[15], S[19], S[23], S[27], S[31], S[35], S[39], S[43], S[47], S[51], S[55], S[59], S[63] };
	v16qi mmx5 = { S[64], S[68], S[72], S[76], S[80], S[84], S[88], S[92], 0, 0, 0, 0, 0, 0, 0, 0 };
	v16qi mmx6 = { S[65], S[69], S[73], S[77], S[81], S[85], S[89], S[93], 0, 0, 0, 0, 0, 0, 0, 0 };
	v16qi mmx7 = { S[66], S[70], S[74], S[78], S[82], S[86], S[90], S[94], 0, 0, 0, 0, 0, 0, 0, 0 };
	v16qi mmx8 = { S[67], S[71], S[75], S[79], S[83], S[87], S[91], S[95], 0, 0, 0, 0, 0, 0, 0, 0 };


#define GCC_VERSION (__GNUC__ * 10000 \
               + __GNUC_MINOR__ * 100 \
                + __GNUC_PATCHLEVEL__)
#if defined(GCC_VERSION)&& GCC_VERSION >= 40700//4.7#endif
//fprintf(stderr,"GCC_VERSION is about 4.7");
	mmx1 = mmx1 << 6;
	mmx1 = mmx1 | (mmx2 << 4);
	mmx1 = mmx1 | (mmx3 << 2);
#else
	v16qi mul64 = {64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64};
	v16qi mul16 = {16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16};
	v16qi mul4 = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
	mmx1 = mmx1 * mul64;
	mmx1 = mmx1 | (mmx2 * mul16);
	mmx1 = mmx1 | (mmx3 * mul4);
#endif

	mmx1 = mmx1 | mmx4;

#if defined(GCC_VERSION)&& GCC_VERSION > 40700//4.7#endif
	mmx5 = mmx5 << 6;
	mmx5 = mmx5 | (mmx6 << 4);
	mmx5 = mmx5 | (mmx7 << 2);
#else
	mmx5 = mmx5 * mul64;
	mmx5 = mmx5 | (mmx6 *mul16);
	mmx5 = mmx5 | (mmx7 *mul4);
#endif

	mmx5 = mmx5 | mmx8;
	memcpy(&kmer3,&mmx1,16);

	int count=CountKmerHitInHash(swap_uint64(kmer3[0]))+ CountKmerHitInHash(swap_uint64(kmer3[1]));

	memcpy(&kmer3, &mmx5, 16);
	count+=CountKmerHitInHash(swap_uint64(kmer3[0]));

	if (count >= RollParam.thresh)
		return true;
	else
		return false;
}
bool BwtIndexer::IsReadFiltered(const ubyte_t * S, const ubyte_t * Q, int len)const
{
/*	if(len==0)
	{
		warning("Read with empty sequence.");
		return true;
	}*/
	//fprintf(stderr, "enter IsReadFiltered\n");
	//if (IsReadHighQ(Q, len)) //pass error
	{
		//if (IsReadInHash(S, len))
		if ((len>200?IsReadInHashByCount(S, len, true):IsReadInHashByCount(S, len)))
		{
//			if (strncmp(name, "ERR015764.2262",14) == 0)
//				fprintf(stderr, "%s not filtered because of kmer hit\n",name);
			return false;
		}
		else
		{
//			if (strncmp(name, "ERR015764.2262",14) == 0)
//							fprintf(stderr, "%s is indeed filtered because of kmer hit\n",name);
			return true;
		}
	}
	//else
	{
//		if (strncmp(name,"ERR015764.2262",14)==0)
//		fprintf(stderr,"%s not filtered because of low qual\n",name);
	//	return false;
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

void BwtIndexer::InitializeRollHashTable(int thresh=3)
{
	RollParam.kmer_size = KMER_SIZE;
	RollParam.read_step_size = READ_STEP_SIZE;
	RollParam.thresh = thresh;
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

int BwtIndexer::AddSeq2DebugHash(int shrinked)
{
	if(debug_hash.find(shrinked)==debug_hash.end())
	{
		debug_hash[shrinked]=1;
	}
	else
		debug_hash[shrinked]++;
	return 0;
}
void BwtIndexer::AddSeq2HashCore(const std::string & Seq, int iter, const std::vector<char>& alleles )
{

	uint64_t datum(0);
	unsigned int i = 0;
	uint32_t shrinked(0);
	//std::string tmp=Seq.substr(0,KMER_SIZE);
	//uint32_t test = 0x00000000903b7458;
	for (; i != 32; ++i)
	{
		datum = ((datum << 2) | NST_NT4_TABLE((unsigned char)Seq[i]));/*
		 & LOMEGA(
		 min_only(RollParam.kmer_size,OVERFLOWED_KMER_SIZE)); */ // N not considered
	}
	shrinked = KmerShrinkage(datum, iter);
	roll_hash_table[iter][shrinked / 8] |= (1 << (shrinked % 8));
	//AddSeq2DebugHash(shrinked);
	//if (debug_flag)
	//{
	//	fprintf(stderr, "I found it!!!!!\t"); 
	//	for (int k = 0; k != 32;++k)
	//		fprintf(stderr, "%c","ACGT"[(datum>>((31-k)*2))&3]);
	//	fprintf(stderr, "\t");
	//	for (int k = 0; k != 16; ++k)
	//		fprintf(stderr, "%c", "ACGT"[(shrinked >> ((15 - k) * 2)) & 3]);
	//	fprintf(stderr, "\t%llx\n", mask[iter]);
	//}
	for (; i != Seq.length() / 2/*last one considered*/; ++i)
	{
		//std::cerr<<"number:"<<i<<"COmparing length:"<<KMER_SIZE<<"and"<<Seq.length()<<"while max:"<<Seq.max_size()<<std::endl;
		datum = ((datum << 2) | NST_NT4_TABLE((unsigned char)Seq[i]));/*
		 & LOMEGA(
		 min_only(RollParam.kmer_size,OVERFLOWED_KMER_SIZE));*/ // N not considered
		//roll_hash_table[iter][KmerShrinkage(datum, iter)]++;
		shrinked = KmerShrinkage(datum, iter);
		//if (debug_flag)
		//{
		//	fprintf(stderr, "I found it!!!!!\t");
		//	for (int k = 0; k != 32; ++k)
		//		fprintf(stderr, "%c", "ACGT"[(datum >> ((31 - k) * 2)) & 3]);
		//	fprintf(stderr, "\t");
		//	for (int k = 0; k != 16; ++k)
		//		fprintf(stderr, "%c", "ACGT"[(shrinked >> ((15 - k) * 2)) & 3]);
		//	fprintf(stderr, "\t%llx\n", mask[iter]);
		//}
		roll_hash_table[iter][shrinked / 8] |= (1 << (shrinked % 8));
		//AddSeq2DebugHash(shrinked);
	}
	//printf("the DATUM is : %016llx    masked DATUM:%x    the hash value is :%d as well as mask:%x\n",datum,KmerShrinkage(datum, iter),roll_hash_table[iter][KmerShrinkage(datum, iter)/8]&(1<<(shrinked%8)), iter);

	uint64_t tmp = datum;
	for (unsigned int j = i, let = 0; let != alleles.size(); ++let, j = i)
	{
		tmp = datum;
		for (; j != Seq.length() / 2 + 32; ++j)
		{
			if (j==Seq.length()/2)
				tmp = ((tmp << 2) | NST_NT4_TABLE((unsigned char)alleles[let]));
			else
				tmp = ((tmp << 2) | NST_NT4_TABLE((unsigned char)Seq[j]));
			shrinked = KmerShrinkage(tmp, iter);
			//if (debug_flag)
			//{
			//	fprintf(stderr, "I found it!!!!!\t%c\t", nst_nt4_table["ACGT"[(char)let]]);
			//	for (int k = 0; k != 32; ++k)
			//		fprintf(stderr, "%c", "ACGT"[(datum >> ((31 - k) * 2)) & 3]);
			//	fprintf(stderr, "\t");
			//	for (int k = 0; k != 16; ++k)
			//		fprintf(stderr, "%c", "ACGT"[(shrinked >> ((15 - k) * 2)) & 3]);
			//	fprintf(stderr, "\t%llx\n", mask[iter]);
			//}
			roll_hash_table[iter][shrinked / 8] |= (1 << (shrinked % 8));
			//AddSeq2DebugHash(shrinked);
		}
	}
	datum = tmp;
	for (i = Seq.length() / 2 + 32; i != Seq.length()/*last one considered*/;
			++i)
	{
		//std::cerr<<"number:"<<i<<"COmparing length:"<<KMER_SIZE<<"and"<<Seq.length()<<"while max:"<<Seq.max_size()<<std::endl;
		datum = ((datum << 2) | NST_NT4_TABLE((unsigned char)Seq[i]));/*
		 & LOMEGA(
		 min_only(RollParam.kmer_size,OVERFLOWED_KMER_SIZE));*/ // N not considered
		//roll_hash_table[iter][KmerShrinkage(datum, iter)]++;
		shrinked = KmerShrinkage(datum, iter);
		//if (debug_flag)
		//{
		//	fprintf(stderr, "I found it!!!!!\t");
		//	for (int k = 0; k != 32; ++k)
		//		fprintf(stderr, "%c", "ACGT"[(datum >> ((31 - k) * 2)) & 3]);
		//	fprintf(stderr, "\t");
		//	for (int k = 0; k != 16; ++k)
		//		fprintf(stderr, "%c", "ACGT"[(shrinked >> ((15 - k) * 2)) & 3]);
		//	fprintf(stderr, "\t%llx\n", mask[iter]);
		//}
		roll_hash_table[iter][shrinked / 8] |= (1 << (shrinked % 8));
		//AddSeq2DebugHash(shrinked);
	}
	//printf("the DATUM is : %x    the hash value is :%d as well as:%x\n",datum,roll_hash_table[datum], LOMEGA(min_only(RollParam.kmer_size,OVERFLOWED_KMER_SIZE)));
}
/******************************/
bool BwtIndexer::BuildIndex(RefBuilder & ArtiRef, string & OldRef,string & NewRef,
		const gap_opt_t * opt)
{
	//NewRef += ".FASTQuick.fa";
	RefPath = OldRef;
	string str;
	Fa2Pac(ArtiRef, NewRef.c_str(), opt); //dump .pac

	/*debug hash print*/
	//DebugHashPrint();
	/*debug hash print end*/
	Fa2RevPac(NewRef.c_str()); //dump .rpac

	DumpRollHashTable(NewRef);
	//str=NewRef+".pac";
	//std::cerr<<"The bwa seq len is:"<<bwa_seq_len(str.c_str())<<std::endl;
	bns_dump(bns, NewRef.c_str());
	notice("Building Pac from Bwt...\n");
	bwt_d = Pac2Bwt(pac_buf);
	//cerr<<"Pac2Bwt..."<<bwt_d->bwt_size<<endl;
	rbwt_d = Pac2Bwt(rpac_buf);
	bwt_gen_cnt_table(bwt_d);
	bwt_gen_cnt_table(rbwt_d);

	bwt_bwtupdate_core(bwt_d);
	//cerr<<"Bwt update..."<<bwt_d->bwt_size<<endl;
	bwt_bwtupdate_core(rbwt_d);
	notice("Dumping Bwt...\n");
	str = NewRef + ".bwt";
	bwt_dump_bwt(str.c_str(), bwt_d);
	str = NewRef + ".rbwt";
	//cerr<<"The Mapper hs rbwt_d->bwt_size is:"<<rbwt_d->bwt_size<<endl;
	bwt_dump_bwt(str.c_str(), rbwt_d);

	//bwt_d=bwt_gen_cnt_table(bwt_d);
	//rbwt_d=bwt_gen_cnt_table(rbwt_d);
	notice("Calculate SA from Bwt...\n");
	str = NewRef + ".sa";
	bwt_cal_sa(bwt_d, 32);
	bwt_dump_sa(str.c_str(), bwt_d);
	str = NewRef + ".rsa";
	bwt_cal_sa(rbwt_d, 32);
	bwt_dump_sa(str.c_str(), rbwt_d);

	DBG(fprintf(stderr,"Indexer initialization finished...\n");)
	return true;
}
bool BwtIndexer::LoadIndex(string & NewRef)
{
	string str;
	ifstream ParamIN(NewRef + ".param");
	std::getline(ParamIN, str);
	stringstream ss(str);
	ss >> RefPath >> RefPath;
	//RefPath=NewRef;
	ReadRollHashTable(NewRef);
	str = NewRef + ".bwt";
	bwt_d = bwt_restore_bwt(str.c_str());
	str = NewRef + ".sa";
	bwt_restore_sa(str.c_str(), bwt_d);
	//bwt_cal_sa(bwt_d, 32);
	str = NewRef + ".rbwt";
	if (stat(str.c_str(), &sb) == 0)
		rbwt_d = bwt_restore_bwt(str.c_str());
	else
	{
		str = NewRef + ".pac";
		string str2 = NewRef + ".rpac";
		bwa_pac_rev_core(str.c_str(), str2.c_str());
		rbwt_d = bwt_pac2bwt(str2.c_str(), 3);
		bwt_gen_cnt_table(rbwt_d);
		bwt_bwtupdate_core(rbwt_d);
		str = NewRef + ".rbwt";
		//cerr<<"The Mapper rbwt_d->bwt_size is:"<<rbwt_d->bwt_size<<endl;
		bwt_dump_bwt(str.c_str(), rbwt_d);
	}
	str = NewRef + ".rsa";
	bwt_restore_sa(str.c_str(), rbwt_d);
	//bwt_cal_sa(rbwt_d, 32);
	bns = bns_restore(NewRef.c_str());
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
	//RefOutput+=".FASTQuick.fa";
	ofstream Fout(RefOutput);
	for (unordered_map<string, uint32_t>::iterator iter =
			ArtiRef.RefTableIndex.begin(); iter != ArtiRef.RefTableIndex.end();
			++iter)
	{

		string CurrentSeqName(iter->first.c_str());
		std::vector<char> alleles;
		alleles.push_back(CurrentSeqName[CurrentSeqName.find('@')+1]);
		alleles.push_back(CurrentSeqName[CurrentSeqName.find('@')+3]);

		//if (CurrentSeqName == "8:32617714@G/A") debug_flag = true;
		//if(ArtiRef.longRefTable[CurrentSeqName]==true) continue;// long ref seq would not be indexed
		string CurrentSeq = ArtiRef.SeqVec[iter->second];
		//if (debug_flag) std::cerr << CurrentSeq << endl;
		AddSeq2Hash(CurrentSeq, alleles);
		//if (debug_flag) std::cerr << string(CurrentSeq.rbegin(), CurrentSeq.rend()) << endl;
		//AddSeq2Hash(string(CurrentSeq.rbegin(), CurrentSeq.rend()), alleles);
		string tmp_rev_cmp = ReverseComplement(CurrentSeq);
		//if (debug_flag) std::cerr << tmp_rev_cmp << endl;
		AddSeq2Hash(tmp_rev_cmp, alleles);
		//if (debug_flag) std::cerr << string(tmp_rev_cmp.rbegin(), tmp_rev_cmp.rend()) << endl;
		//AddSeq2Hash(string(tmp_rev_cmp.rbegin(), tmp_rev_cmp.rend()), alleles);
		Fout << ">" << CurrentSeqName << "\n" << CurrentSeq << "\n";
		//if (debug_flag) exit(EXIT_FAILURE);
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
		if (bwt_d)
		{
			
		if (bwt_d->sa)
			free(bwt_d->sa);
		if (bwt_d->bwt)
			free(bwt_d->bwt);
			free(bwt_d);
		}

	DBG(fprintf(stderr,"bwt_d delete successfully...\n");)
		if (rbwt_d)
		{
			
		if (rbwt_d->sa)
			free(rbwt_d->sa);
		if (rbwt_d->bwt)
			free(rbwt_d->bwt);
		free(rbwt_d);
		}

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
		if (bns->ambs)
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i)
		{

			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		if (bns->anns)
		free(bns->anns);
		free(bns);
	}
	DestroyRollHashTable();
}

