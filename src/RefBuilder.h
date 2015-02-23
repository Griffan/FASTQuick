/*
 * RefBuilder.h
 *
 *  Created on: 2014Äê7ÔÂ9ÈÕ
 *      Author: Administrator
 */

#ifndef REFBUILDER_H_
#define REFBUILDER_H_
#include "Utility.h"
#include "../libbwa/bwtaln.h"

class RefBuilder
{
public:
	std::vector<std::string> SeqVec;
	std::unordered_map<std::string,uint32_t > RefTableIndex;
	//unordered_map<string,bool> longRefTable;
	RefBuilder();
	RefBuilder(const std::string& VcfPath, const std::string& RefPath, const std::string& DBsnpPath, const std::string& MaskPath, const gap_opt_t* opt, bool reselect);//, std::unordered_map<std::string,bool>& longRefTable);
	virtual ~RefBuilder();
};

#endif /* REFBUILDER_H_ */
