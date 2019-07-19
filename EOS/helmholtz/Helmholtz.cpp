#include "Helmholtz.H"

/* Define Helmholtz static class data */

Array1D<amrex::Real, Helmholtz::imax> Helmholtz::d;
Array1D<amrex::Real, Helmholtz::jmax> Helmholtz::t;

// ..for the helmholtz free energy tables (2D)
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::f;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::fd;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::ft;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::fdd;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::ftt;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::fdt;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::fddt;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::fdtt;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::fddtt;

// ..for the pressure derivative with density ables
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::dpdf;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::dpdfd;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::dpdft;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::dpdfdt;

// ..for chemical potential tables
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::ef;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::efd;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::eft;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::efdt;

// ..for the number density tables
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::xf;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::xfd;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::xft;
Array2D<amrex::Real, Helmholtz::imax, Helmholtz::jmax> Helmholtz::xfdt;

// ..for storing the differences
Array1D<amrex::Real, Helmholtz::jmax> Helmholtz::dt_sav;
Array1D<amrex::Real, Helmholtz::jmax> Helmholtz::dt2_sav;
Array1D<amrex::Real, Helmholtz::jmax> Helmholtz::dti_sav;
Array1D<amrex::Real, Helmholtz::jmax> Helmholtz::dt2i_sav;

Array1D<amrex::Real, Helmholtz::imax> Helmholtz::dd_sav;
Array1D<amrex::Real, Helmholtz::imax> Helmholtz::dd2_sav;
Array1D<amrex::Real, Helmholtz::imax> Helmholtz::ddi_sav;
Array1D<amrex::Real, Helmholtz::imax> Helmholtz::dd2i_sav;
