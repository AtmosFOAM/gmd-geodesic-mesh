# $1 -- NGRIDS value
define ngrids

gengrid_hex.$1.f: gengrid_hex.f
	NGRIDS=$1 envsubst < $$< > $$@

gengrid_hex.$1: gengrid_hex.$1.f
	gfortran -mcmodel=medium $$< -o $$@

endef

refinements := 3 4 5 6 7 8 9

.PHONY: all

all: $(foreach n,$(refinements),gengrid_hex.$n)

$(foreach n,$(refinements),$(eval $(call ngrids,$n)))
