# Fortran implementation of the geodesic routines in geodesic_module

This is a library to solve geodesic problems on an ellipsoid model of
the earth.

The two tools ngsforward and ngsinverse are replacements for the tools
FORWARD and INVERSE available from the
[NGS](http://www.ngs.noaa.gov/PC_PROD/Inv_Fwd/)

Licensed under the MIT/X11 License; see
[LICENSE.txt](https://geodesic_module.sourceforge.io/LICENSE.txt).

The algorithms are documented in

* C. F. F. Karney,
  [Algorithms for geodesics](https://doi.org/10.1007/s00190-012-0578-z),
  J. Geodesy **87**(1), 43–55 (2013);
  [Addenda](https://geodesic_module.sourceforge.io/geod-addenda.html).

## Other links:

* Library documentation: https://geodesic_module.sourceforge.io/Fortran/doc
* Change log: https://geodesic_module.sourceforge.io/Fortran/doc/changes.html
* GIT repository: https://github.com/geodesic_module/geodesic_module-fortran
  Releases are tagged in git as, e.g., [`v1.52`](../../tree/v1.52),
  [`v2.0`](../../tree/v2.0), etc.
* Source distribution:
  https://sourceforge.net/projects/geodesic_module/files/distrib-Fortran
* geodesic_module: https://geodesic_module.sourceforge.io
* Author: Charles Karney, <charles@karney.com>
