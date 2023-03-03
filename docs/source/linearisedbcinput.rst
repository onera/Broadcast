.. _linearisedbcinput:

Inputs for the linearised boundary conditions
==============


* :func:`bc_general_2d_d`: (*w, wd, loc, interf, field, im, jm*).

* :func:`bc_no_reflexion_2d_d`: (*w, wd, wbd, loc, interf, nx, ny, gam, gh, im, jm*).

* :func:`bc_supandsubinlet_2d_d`: (*w, wd0, loc, interf, field, nx, ny, gam, gh, im, jm*).

* :func:`bc_extrapolate_o2_2d_d`: (*w, wd, loc, interf, im, jm, gh*).

* :func:`bc_symmetry_2d_d`: (*w, wd, loc, interf, nx, ny, gh, im, jm*).

* :func:`bc_antisymmetry_2d_d`: (*w, wd, loc, interf, nx, ny, gh, im, jm*).

* :func:`jn_match_2d_d`: (*wr, wrd, prr, gh1r, gh2r, gh3r, gh4r, imr, jmr, wd, wdd, prd, gh1d, gh2d, gh3d, gh4d, imd, jmd, tr*).

* :func:`bc_wall_viscous_adia_2d_d`: (*w, wd, loc, gam, interf, gh, im, jm*).

* :func:`bc_wall_viscous_iso_2d_d`: (*w, wd, twall, loc, gam, rgaz, interf, gh, im, jm*).

* :func:`bc_wall_viscous_iso_profile_2d_d`: (*w, wd, twallprof, twallprofd, loc, gam, gamd, rgaz, rgazd, interf, gh, im, jm*).

* :func:`bc_wall_blow_profile_2d_d`: (*w, wd, velprof, velprofd, loc, gam, gamd, rgaz, interf, gh, im, jm*).
