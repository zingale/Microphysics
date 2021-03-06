#ifndef _conductivity_H_
#define _conductivity_H_

#include <eos_type.H>
#include <actual_conductivity.H>

using namespace amrex;

AMREX_FORCE_INLINE
void conductivity_init() {
  actual_conductivity_init();
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void conductivity(eos_t& state) {
  actual_conductivity(state);
}
#endif
