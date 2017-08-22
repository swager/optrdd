When using our method for regression discontinuity designs with a multivariate
running variable, we recommend using the optimizer MOSEK. This is a commerical
optimizer that needs to be installed separately; however, MOSEK is free for academics.
MOSEK is called from R via the RMosek interface.

MOSEK can be installed as follows

1. [Download](https://www.mosek.com/resources/downloads) MOSEK for the relevant platform.
2. Place the downloaded files in the home directory (e.g., /home/<userid>/mosek/)
3. Request a licence. If you're an [academic](https://www.mosek.com/resources/academic-license),
you can get one for free.
4. Place the license in the mosek directory (e.g., /home/<userid>/mosek/mosek.lic)
5. Install the RMosek interface. The following (Mac OSX specific) command worked for me:

```R
install.packages("Rmosek", type="source", repos="http://download.mosek.com/R/8", configure.vars=c("PKG_MOSEKHOME=~/mosek/8/tools/platform/osx64x86", "PKG_MOSEKLIB=mosek64"))
```