#ifndef EOS_H
#define EOS_H
#include "AMReX_REAL.H"

class Eos {
public:

  std::string eos_name;
  static Real smallt, smalld; // these should be deprecated for mintemp, mindens below ...

  static Real mintemp;
  static Real maxtemp;
  static Real mindens;
  static Real maxdens;
  static Real minx;
  static Real maxx;
  static Real minye;
  static Real maxye;
  static Real mine;
  static Real maxe;
  static Real minp;
  static Real maxp;
  static Real mins;
  static Real maxs;
  static Real minh;
  static Real maxh;
};
#endif
