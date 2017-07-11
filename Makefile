# $1 -- NGRIDS value
define ngrids

gengrid_hex.$1.f: gengrid_hex.f
	NGRIDS=$1 envsubst < $$< > $$@

gengrid_hex.$1: gengrid_hex.$1.f
	gfortran -mcmodel=medium $$< -o $$@

endef

refinements := 3 4 5 6 7 8 9

.PHONY: clean all

all: $(foreach n,$(refinements),gengrid_hex.$n)

clean:
	$(RM) gengrid_hex.[0-9]
	$(RM) gengrid_hex.?.f

$(foreach n,$(refinements),$(eval $(call ngrids,$n)))
