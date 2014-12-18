#include "fVcf.h"

const float fVcf::NAN_FLT = sqrtf(-1.);     // assign float  NAN value
const double fVcf::NAN_DBL = sqrt(-1.);  // assign double NAN value
std::string fVcf::xLabel = "X";
std::string fVcf::yLabel = "Y";
std::string fVcf::mitoLabel = "MT";
int fVcf::xStart = 2699520;
int fVcf::xStop = 154931044;
