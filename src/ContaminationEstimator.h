/*
 * ContaminationEstimator.h
 *
 *  Created on: Oct 30, 2015
 *      Author: fanzhang
 */

#ifndef CONTAMINATIONESTIMATOR_H_
#define CONTAMINATIONESTIMATOR_H_
#include <string>
class ContaminationEstimator {
public:

	/*Initialize from existed UD*/
	/*This assumes the markers are the same as the selected vcf*/
	int RunFromSVDMatrix(const std::string UDpath, const std::string PCpath, const std::string Mean, const std::string & GLpath, const std::string &Bed,const std::string& Prefix,const std::string& ReadGroup);
	/*Directly obtain AFs from specified VCF files*/
	int RunFromVCF(const std::string VcfSiteAFFile,const std::string CurrentMPU, const std::string ReadGroup, const std::string Prefix);

public:
	ContaminationEstimator();
	virtual ~ContaminationEstimator();
};

#endif /* CONTAMINATIONESTIMATOR_H_ */
