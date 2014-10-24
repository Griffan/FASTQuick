#ifndef __F_PED_H
#define __F_PED_H

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <set>
#include <map>

#include "pFile.h"
#include "wFile.h"
#include "Error.h"
#include "boolParser.h"
#include "PhredHelper.h"

class fPed {
 public:
  std::vector<int> dads;
  std::vector<int> moms;
  std::vector<int> sexes;
  std::map<std::string,int> icols;
  std::map<std::string,int> iinds;
  std::vector<std::string> data;
  int ninds;
  int ncols;

  std::string& getPheno(const char* pheno, const char* ind) {
    if ( icols.find(pheno) == icols.end() ) error("Cannot find phenotype name '%s'",pheno);
    if ( iinds.find(ind) == iinds.end() ) error("Cannot find individual ID '%s'",ind);
    return data[icols[pheno] + ncols * iinds[ind]];
  }

  std::string& getPheno(const char* pheno, int iind) {
    //notice("foo");
    //notice("icols[pheno] = %d",icols[pheno]);
    if ( icols.find(pheno) == icols.end() ) error("Cannot find phenotype name %s",pheno);
    return data[icols[pheno] + ncols * iind];
  }

  fPed(const char* ped, const char* dat = NULL) {
    pFile tfped;
    const char* line;
    std::vector<std::string> colnames, tokens;
    int i;

    tfped.load(ped);
    if ( ( dat != NULL ) && ( strlen(dat) > 0 ) ) {
      pFile tfdat(dat);  
      colnames.push_back("FAM_ID");
      colnames.push_back("IND_ID");
      colnames.push_back("DAD_ID");
      colnames.push_back("MOM_ID");
      colnames.push_back("SEX");
      while( (line = (char*)tfdat.getLine()) != NULL ) {
	pFile::tokenizeLine(line," \t\n",tokens);
	if ( tokens.size() != 2 ) error("Cannot parse a line in dat file %s\n",dat);
	colnames.push_back(tokens[0]);
      }
      tfdat.close();
    }
    else {
      line = (char*)tfped.getLine();
      if ( line[0] == '#' ) {
	pFile::tokenizeLine(line+1," \t\n",tokens);
	for(i=0; i < (int)tokens.size(); ++i) {
	  switch(i) {
	  case 0:
	    if ( tokens[i] != "FAM_ID" ) 
	      error("The first column must be named as FAM_ID");
	    break;
	  case 1:
	    if ( tokens[i] != "IND_ID" ) 
	      error("The second column must be named as IND_ID");
	    break;
	  case 2:
	    if ( ( tokens[i] != "FAT_ID" ) && ( tokens[i] != "DAD_ID" ) )
	      error("The third column must be named as DAD_ID");
	    break;
	  case 3:
	    if ( ( tokens[i] != "MOT_ID" ) && ( tokens[i] != "MOM_ID" ) )
	      error("The fourth column must be named as MOM_ID");
	    break;
	  case 4:
	    if ( tokens[i] != "SEX" ) 
	      error("The fifth column must be named as SEX");
	    break;
	  default:
	    break;
	    // do nothing
	  }
	  colnames.push_back(tokens[i]);
	}
      }
      else {
	error("Without a separate .dat file specified, we expect a header '#' in the first line of ped file %s",ped);
      }
    }

    for(i=0; i < (int)colnames.size(); ++i) {
      if ( icols.find(colnames[i]) != icols.end() ) error("Duplicate column name %s",colnames[i].c_str());
      icols[colnames[i]] = i;
    }

    tokens.clear();
    ncols = (int)colnames.size();

    ninds = 0;
    while( (line = (char*)tfped.getLine()) != NULL ) {
      pFile::tokenizeLine(line," \t\n",tokens);
      if ( tokens.size() != icols.size() ) error("Number of tokens differ between header (%d) at %s",icols.size(),line);
      for(i=0; i < (int)tokens.size() ; ++i) {
	data.push_back(tokens[i]);
      }

      if ( iinds.find(data[ninds*ncols+1]) != iinds.end() ) error("Duplicate individual ID %d",data[ninds*ncols+1].c_str());
      iinds[data[ninds*ncols+1]] = ninds;
      ++ninds;
    }

    tfped.close();

    for(int i=0; i < ninds; ++i) {
      if ( data[i*ncols+2].compare("0") != 0 ) {
	if ( iinds.find(data[i*ncols+2]) == iinds.end() ) {
	  warning("Father's ID %s do not exist in current PED. Ignoring..",data[i*ncols+2].c_str());
	  dads.push_back(-1);
	}
	else {
	  dads.push_back(iinds[data[i*ncols+2]]);
	}
      }
      else {
	dads.push_back(-1);
      }

      if ( data[i*ncols+3].compare("0") != 0 ) {
	if ( iinds.find(data[i*ncols+3]) == iinds.end() ) {
	  warning("Mother's ID %s do not exist in current PED. Ignoring..",data[i*ncols+3].c_str());
	  moms.push_back(-1);
	}
	else {
	  moms.push_back(iinds[data[i*ncols+3]]);
	}
      }
      else {
	moms.push_back(-1);
      }

      if ( data[i*ncols+4].compare("0") == 0 ) sexes.push_back(0);
      else if ( data[i*ncols+4].compare("1") == 0 ) sexes.push_back(1);
      else if ( data[i*ncols+4].compare("2") == 0 ) sexes.push_back(2);
      else {
	warning("Cannot recognize sex %s. Setting to zero..",data[i*ncols+4].c_str());
	sexes.push_back(0);
      }
    }

    //error("ncols=%d, ninds=%d, data[5] = %s, data[11] = %s",ncols,ninds,data[5].c_str(),data[5].c_str());
    //error("getPheno(MPU,1) = %s",getPheno("MPU",1).c_str());
  }
};

#endif // __FPED_H
