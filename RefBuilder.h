/*
 * RefBuilder.h
 *
 *  Created on: 2014Äê7ÔÂ9ÈÕ
 *      Author: Administrator
 */

#ifndef REFBUILDER_H_
#define REFBUILDER_H_
#include "Utility.h"
#include "./libbwa/bwtaln.h"
using namespace std;

class RefBuilder
{
public:
	vector<string> SeqVec;
	unordered_map<string,uint32_t > RefTableIndex;
	RefBuilder();
	RefBuilder(string VcfPath,string RefPath, string MaskPath, const gap_opt_t* opt);
	virtual ~RefBuilder();
};

#endif /* REFBUILDER_H_ */
