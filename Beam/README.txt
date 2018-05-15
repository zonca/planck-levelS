Polarised Beam Analysis Package
===============================

Mark Ashdown, CPAC

This package contains tools for reading beam data from Grasp files (in
a limited set of "cut" and "grid" sub-formats), converting polarised
beam data into Stokes parameters, calculating effective beams
including cross-polar leakage in the detectors and extracting the
spherical harmonic coefficients from the Stokes parameters.

The sources of the modules are commented, as are the F90 modules on
which they are based.  For descriptions of the parameters required by
the modules, see the corresponding text files, for example,
grasp2stokes.par.txt.

The five modules that make up this package are -

grasp2stokes
------------

Reads beam amplitudes from a Grasp file and converts them to Stokes
parameters.  This module can read four particular sub-formats of Grasp
files which are denoted as follows -

 grd_polar_hfi: grid format file containing beam data on a spherical
                polar theta-phi grid.  Used by HFI for main beam
                simulations.  This type of input will be converted
                into Stokes parameters on a polar grid.

grd_square_hfi: grid format file containing beam data on a square u-v
                grid.  Used by HFI for main beam simulations.  This
                type of input will be converted into Stokes parameters
                on a square grid.

grd_square_lfi: grid format file containing beam data on a square u-v
                grid.  Used by LFI for main beam simulations.  This
                type of input will be converted into Stokes parameters
                on a square grid.

           cut: cut format file containing beam data on constant-phi
                cuts.  Usually contains full-sky data.  This type of
                input will be converted into Stokes parameters on a
                polar grid.
                NOTE: The cuts must have an odd number of samples, i.e.
                there must be samples exactly at the North Pole!


crosspol
--------

Combines two perfectly polarised beams (as produced by grasp2stokes)
to create an effective beam with detector cross-polar leakage
included.  The "co-polar" beam which measures the optical response of
the detector along its direction of polarisation is combined with the
"cross-polar" beam which measures the optical response in the
orthogonal direction.  It is possible to use only the co-polar beam to
simulate the effective beam: in this case it is assumed that the
cross-polar beam has the same shape as the co-polar beam and only
differs in the direction of polarisation.

In order to be able to combine co- and cross-polar beam patterns, they
must be simulated in coordinate systems with the same (nominal)
pointing centres, that is, the z-axes of the coordinate systems must
be aligned.  The x- and y- axes need not be aligned, but there are
some restrictions; see below.  The "angle" parameter gives the angle
through which the cross-polar beam must be rotated (in a right-handed
sense about the z-axis) in order to align it in the co-polar
coordinate system.  In order to be able to make the effective beam
without interpolating the cross-polar beam, permissible values of
"angle" are restricted to those which align the grids on which the
beam patterns are sampled.  For beams on a square grid, this restricts
"angle" to be a multiple of 90 degrees.  For beams on a polar grid, it
restricts "angle" to be a multiple of the phi spacing of the grid.


beam2alm
--------

Reads beam Stokes parameters and then calculates their spherical
harmonic coefficients.

This program can take its input from one or two objects:

1) Main beam at high resolution (the "north polar cap");
2) Full-sky beam at lower resolution.

Both are optional, but there must be at least one input object!

The main beam may be a polar or square beam object.  If it is a square
beam object, it will be interpolated onto a polar grid before
extracting the multipoles.  The full-sky beam must be in a polar beam
object.

If there is both a main beam and a full-sky beam then the full-sky
beam will be interpolated in the theta direction so that it has the
same theta resolution as the full-sky beam. This is to give better
extraction of the multipoles.


gaussbeampol
------------

Calculates the multipoles of an linearly-polarised elliptical Gaussian
beam with following properties -

  - beam points in z-direction (towards "north pole").

  - beam shape is described by three parameters:
    1) mean FWHM, defined as fwhm = sqrt(fwhm_max*fwhm_min);
    2) ellipticity, defined as fwhm_max/fwhm_min;
    3) orientation angle. Major axis is oriented at angle psi_ell to
       the x-axis;

  - beam is polarised along direction at angle psi_pol to the
     x-axis;

  - beam has cross-polar leakage of degree epsilon.


stokes_extract
--------------

Extracts the Stokes parameters from the DMC (or a FITS file) and
writes them to a text file.  This is used for visualisation and
checking purposes only.
