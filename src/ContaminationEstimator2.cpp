/*
 * ContaminationEstimator2.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: fanzhang
 */

#include "ContaminationEstimator2.h"
#include "../libmpu/mpuTool.h"
#include "PopulationIdentifier.h"
#include <iostream>
#include <sstream>

#ifndef MPU_PATH
#define MPU_PATH "mpuTools"
#endif
ContaminationEstimator2::ContaminationEstimator2() {
	// TODO Auto-generated constructor stub

}

ContaminationEstimator2::~ContaminationEstimator2() {
	// TODO Auto-generated destructor stub
}



int ContaminationEstimator2::RunFromVCF(const std::string VcfSiteAFFile,const std::string CurrentMPU, const std::string ReadGroup, const std::string Prefix)
{
	char cmdline[9*1024];
	char** params = new char*[9];
	sprintf(cmdline, MPU_PATH " verify --vcf %s --mpu %s --smID %s --out %s.ctm", VcfSiteAFFile.c_str(), CurrentMPU.c_str(), ReadGroup.c_str(), Prefix.c_str());
	std::stringstream ss(cmdline);
	std::string para;
	ss >> para;
	for (int i = 0; i != 9; ++i)
	{
		ss >> para;
		params[i] = new char[1024];
		strcpy(params[i], para.c_str());
	}
	runVerify(9, params);
	return 0;
}
int ContaminationEstimator2::RunFromSVDMatrix(const std::string UDpath, const std::string PCpath, const std::string Mean, const std::string & MPUpath, const std::string &Bed, const std::string& Prefix,const std::string& ReadGroup)
{
	PopulationIdentifier pop(UDpath, PCpath, Mean, MPUpath, std::string("Empty"),Bed);
	pop.OptimizeLLK();
	pop.writeVcfFile(MPUpath+".vcf");
	/*preparations done, then call RunFromVCF*/
	RunFromVCF(std::string(MPUpath+".vcf"),MPUpath,ReadGroup,Prefix);
	return 0;
}



