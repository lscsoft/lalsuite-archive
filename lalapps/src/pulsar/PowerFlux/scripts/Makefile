
UPLOAD_FILES=make_run.tcl split.tcl get_param.tcl make_1800_sfts.tcl get_lines.tcl

all:


upload:
	scp -S globus-ssh -C $(UPLOAD_FILES) hydra.phys.uwm.edu:PowerFlux/scripts/

upload.h:
	scp -S globus-ssh -C $(UPLOAD_FILES) ldas-grid.ligo-wa.caltech.edu:PowerFlux/scripts/

