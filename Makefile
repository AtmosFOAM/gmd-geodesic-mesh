gengrid_hex: gengrid_hex.f
	gfortran -mcmodel=medium $< -o $@
