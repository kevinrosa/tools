PROGRAM = glabel_matlab nsfces_matlab

FC=gfortran

FFLAGS = -O4 -ffixed-line-length-132

.KEEP_STATE:

.PRECIOUS :

LIBS       =  ../gamma.a -lnetcdff -lnetcdf

all: $(PROGRAM)

$(PROGRAM) : $@
	$(LINK.f) $@.f -o $@ $(LIBS)

clean :
	rm -f *% core $@
