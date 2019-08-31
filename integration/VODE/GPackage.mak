ifdef SDC
  F90sources += vode_integrator_sdc.F90
  F90sources += vode_rhs_sdc.F90
  F90sources += vode_type_sdc.F90
else
  f90sources += vode_integrator.f90
  F90sources += vode_rhs.F90
  F90sources += vode_type.F90
endif
