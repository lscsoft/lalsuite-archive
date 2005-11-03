
UPLOAD_FILES=make_run.tcl split.tcl get_param.tcl make_1800_sfts.tcl \
	get_lines.tcl make_power_injection_run.tcl split_pi.tcl reshape.tcl

all:


upload:
	scp -S globus-ssh -C $(UPLOAD_FILES) hydra.phys.uwm.edu:PowerFlux/scripts/

upload.h:
	scp -S globus-ssh -C $(UPLOAD_FILES) ldas-grid.ligo-wa.caltech.edu:PowerFlux/scripts/

upload.l:
	scp -S globus-ssh -C $(UPLOAD_FILES) ldas-grid.ligo-la.caltech.edu:PowerFlux/scripts/
