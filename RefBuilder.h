/*
 * RefBuilder.h
 *
 *  Created on: 2014��7��9��
 *      Author: Administrator
 */

#ifndef REFBUILDER_H_
#define REFBUILDER_H_
#include "Utility.h"
#include "./libbwa/bwtaln.h"
//using namespace std;

class RefBuilder
{
public:
	std::vector<std::string> SeqVec;
	std::unordered_map<std::string,uint32_t > RefTableIndex;
	//unordered_map<string,bool> longRefTable;
	RefBuilder();
	RefBuilder(std::string VcfPath,std::string RefPath, std::string MaskPath, const gap_opt_t* opt);//, std::unordered_map<std::string,bool>& longRefTable);
	virtual ~RefBuilder();
};

#endif /* REFBUILDER_H_ */
