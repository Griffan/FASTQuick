//
// Created by Fan Zhang on 2/6/20.
//
#include "../src/RegionList.h"
#include "../src/StatCollector.h"
#include "gtest/gtest.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

//gtest configuration test
int add(int a, int b){
  return a+b;
}

TEST(test1, c1){
  EXPECT_EQ(3, add(1,2));
}

//project specific tests
///region tests
TEST(test2,c1){
  RegionList t("/Users/fanzhang/Downloads/FASTQuick_test/EmpQualityInvestigation/20130108.exome.targets.bed.nochr.bed");
  EXPECT_FALSE(t.IsEmpty());
  EXPECT_TRUE(t.IsOverlapped("1",14642));
  EXPECT_TRUE(t.IsOverlapped("1",14742));
  EXPECT_TRUE(t.IsOverlapped("1",14881));
  EXPECT_FALSE(t.IsOverlapped("1",14641));

  EXPECT_TRUE(t.IsOverlapped("1",14944));

  EXPECT_TRUE(t.IsOverlapped("Y",59342899));
}


TEST(test2,c2){
  RegionList a;
  a.AddRegion("chr1", 100, 200);

  RegionList b;
  b.AddRegion("chr1", 150, 160);

  RegionList e;
  e.AddRegion("chr1", 150, 160);

  a.Join(b,false);
//  a.PrintRegion();
  EXPECT_TRUE(a==e);
  EXPECT_TRUE(a.Size() == 11);

  a.Clear();
  a.AddRegion("chr1", 300, 400);

  b.Clear();
  b.AddRegion("chr1", 250, 350);

  e.Clear();
  e.AddRegion("chr1", 300, 350);

  a.Join(b,false);
//  a.PrintRegion();
  EXPECT_TRUE(a==e);
  EXPECT_TRUE(a.Size() == 51);


  a.Clear();
  a.AddRegion("chr1", 100, 155);
  a.AddRegion("chr1", 155, 200);

  a.AddRegion("chr1", 300, 400);
  a.AddRegion("chr1", 500, 600);

  b.Clear();
  b.AddRegion("chr1", 150, 160);
  b.AddRegion("chr1", 250, 350);
  b.AddRegion("chr1", 550, 650);
  b.AddRegion("chr1", 650, 750);

  e.Clear();
  e.AddRegion("chr1", 150, 160);
  e.AddRegion("chr1", 300, 350);
  e.AddRegion("chr1", 550, 600);

  a.Join(b,false);
//  a.PrintRegion();

  EXPECT_TRUE(a==e);
  EXPECT_TRUE(a.Size() == 113);

  a.Clear();
  a.AddRegion("chr1", 990130, 990630);
  a.AddRegion("chr1", 1020346, 1022346);

  // 1       990130  990630
  // 1       1020346 1022346


  b.Clear();
  b.AddRegion("chr1", 989819, 989939);
  b.AddRegion("chr1", 990162, 990402);
  //1       989819  989939
  //1       990162  990402

  e.Clear();
  e.AddRegion("chr1", 990162, 990402);


  a.Join(b,false);
  a.PrintRegion();

  EXPECT_TRUE(a==e);
  EXPECT_TRUE(a.Size() == 241);

}

///Refseq recovering from cigar and MD tag
TEST(test3,c1){
  std::string read = "AAAAAATAAAAAA";

  std::string refseq = StatCollector::RecoverRefseqByMDandCigar(read,"T12","13M");
  EXPECT_TRUE(refseq == "TAAAAATAAAAAA");

  refseq = StatCollector::RecoverRefseqByMDandCigar(read,"4T8","13M");
  EXPECT_TRUE(refseq == "AAAATATAAAAAA");

  refseq = StatCollector::RecoverRefseqByMDandCigar(read,"12T","13M");
  EXPECT_TRUE(refseq == "AAAAAATAAAAAT");

  refseq = StatCollector::RecoverRefseqByMDandCigar(read,"11T0T","13M");
  EXPECT_TRUE(refseq == "AAAAAATAAAATT");

  refseq = StatCollector::RecoverRefseqByMDandCigar(read,"13","4M3I9M");
  EXPECT_TRUE(refseq == "AAAAAATAAAAAA");

  refseq = StatCollector::RecoverRefseqByMDandCigar(read,"9","4S9M");
  EXPECT_TRUE(refseq == "AATAAAAAA");

  refseq = StatCollector::RecoverRefseqByMDandCigar(read,"9","9M4S");
  EXPECT_TRUE(refseq == "AAAAAATAA");

  refseq = StatCollector::RecoverRefseqByMDandCigar(read,"9^TTTT4","9M4D4M");
  EXPECT_TRUE(refseq == "AAAAAATAATTTTAAAA");

  refseq = StatCollector::RecoverRefseqByMDandCigar(read,"6G^GGG5T","7M3D6M");
  EXPECT_TRUE(refseq == "AAAAAAGGGGAAAAAT");

  refseq = StatCollector::RecoverRefseqByMDandCigar(read,"9^TTTT2^TT2","9M4D2M2D2M");
  EXPECT_TRUE(refseq == "AAAAAATAATTTTAATTAA");

}
