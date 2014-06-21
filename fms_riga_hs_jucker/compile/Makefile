# Makefile for Experiment 'HSt42'
include intel.mk

fms_HSt42.x: libfms_spectral_solo.a
	$(LD) $^ $(LDFLAGS) -o $@ 

libfms_spectral_solo.a:  FORCE
	make  NETCDF=3  -f Makefile.fms_spectral_solo $@

FORCE:

clean:
	make  NETCDF=3  -f Makefile.fms_spectral_solo clean

localize:
	make  NETCDF=3  -f Makefile.fms_spectral_solo localize

