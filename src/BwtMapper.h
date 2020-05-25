/*The MIT License (MIT)
Copyright (c) 2017 Fan Zhang, Hyun Min Kang
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */
/* Contact: Fan Zhang <fanzhang@umich.edu> */

#ifndef BWTMAPPER_H_
#define BWTMAPPER_H_
#include "../libbwa/bwtaln.h"
#include "../misc/bam/BamInterface.h"
#include "../misc/bam/SamFile.h"
#include "../misc/bam/SamFileHeader.h"
#include "../misc/bam/SamRecord.h"
#include "BwtIndexer.h"
#include "StatCollector.h"
#include <fstream>
#include <iostream>
#ifndef BWAPE_H_
#include "../libbwa/bwape.h"
#endif
#define READ_BUFFER_SIZE 1000000
class BwtMapper {
  StatCollector collector;
  char *bwa_rg_line;
  char *bwa_rg_id;

public:
  BwtMapper();

  BwtMapper(BwtIndexer &BwtIndex, const std::string &FQList,
            const std::string &Fastq_1, const std::string &Fastq_2,
            const std::string &Prefix, const std::string &RefPath,
            const pe_opt_t *popt, gap_opt_t *opt,
            const std::string &targetRegionPath);

  int bwa_cal_pac_pos(BwtIndexer &BwtIndex, int n_seqs, bwa_seq_t *seqs,
                      int max_mm, float fnr);

  // bwa_seq_t* bwa_read_seq_with_hash(BwtIndexer* BwtIndex, bwa_seqio_t *bs,
  // int n_needed, int *n, int mode, int trim_qual, double frac, uint32_t seed);
  // int bwa_read_seq_pair_with_hash(BwtIndexer& BwtIndex, bwa_seqio_t *bs,
  // bwa_seqio_t *bs2, int n_needed, int *n, int mode, int trim_qual, double
  // frac, bwa_seq_t * seqs, bwa_seq_t * seqs2); static void
  // bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt[2], int n_seqs, bwa_seq_t
  // *seqs, const gap_opt_t *opt);
  bool SingleEndMapper(BwtIndexer &BwtIndex, const gap_opt_t *opt,
                       SamFileHeader &SFH, BamInterface &BamIO, IFILE BamFile,
                       StatGenStatus &StatusTracker, std::ofstream &fout,
                       FileStatCollector &FSC);

  bool PairEndMapper(BwtIndexer &BwtIndex, const pe_opt_t *popt, gap_opt_t *opt,
                     SamFileHeader &SFH, BamInterface &BamIO, IFILE BamFile,
                     StatGenStatus &StatusTracker, std::ofstream &fout,
                     FileStatCollector &FSC);

  bool PairEndMapper_without_asyncIO(BwtIndexer &BwtIndex, const pe_opt_t *popt,
                                     gap_opt_t *opt, SamFileHeader &SFH,
                                     BamInterface &BamIO, IFILE BamFile,
                                     StatGenStatus &StatusTracker,
                                     std::ofstream &fout,
                                     FileStatCollector &FSC);

  bool SetSamRecord(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate,
                    SamFileHeader &SFH, SamRecord &SR, const gap_opt_t *opt);

  bool SetSamFileHeader(SamFileHeader &SFH, const BwtIndexer &BwtIndex);

  int bwa_set_rg(const char *s);

  virtual ~BwtMapper();

  inline void revealopt(const gap_opt_t *opt) {
    std::cerr << opt->s_mm << std::endl;
    std::cerr << opt->s_gapo << std::endl;
    std::cerr << opt->s_gape << std::endl;
    std::cerr << opt->mode << std::endl
              << opt->indel_end_skip << std::endl
              << opt->max_del_occ << std::endl
              << opt->max_entries << std::endl
              << opt->fnr << std::endl
              << opt->max_diff << std::endl
              << opt->max_gapo << std::endl
              << opt->max_gape << std::endl
              << opt->max_seed_diff << std::endl
              << opt->seed_len << std::endl
              << opt->n_threads << std::endl
              << opt->max_top2 << std::endl
              << opt->trim_qual << std::endl
              << "opt reveal end" << std::endl;
  }

private:
  BwtMapper(BwtMapper const &);

  void operator=(BwtMapper const &x);
};

#endif /* BWTMAPPER_H_ */
