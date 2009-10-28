
UPLOAD_FILES=make_run.tcl split.tcl get_param.tcl make_1800_sfts.tcl \
	get_lines.tcl make_power_injection_run.tcl split_pi.tcl reshape.tcl \
	generate_injection_batch.tcl make_sft_injection_run.tcl split_dataset.tcl \
	split_sfti.tcl make_second_sft_injection_run.tcl normalize_sft_name.tcl \
	make_native_sft_injection_run.tcl reshape_sfti.tcl archive.tcl keyed_split.tcl \
	make_split_run.tcl find_outlier_matches.tcl make_outlier_histogram.tcl

all:


upload:
	globus-scp -C $(UPLOAD_FILES) hydra.phys.uwm.edu:PowerFlux/scripts/

upload.h:
	globus-scp -C $(UPLOAD_FILES) ldas-grid.ligo-wa.caltech.edu:PowerFlux/scripts/

upload.l:
	globus-scp -C $(UPLOAD_FILES) ldas-grid.ligo-la.caltech.edu:PowerFlux/scripts/

upload.c:
	globus-scp -C $(UPLOAD_FILES) ldas-grid.ligo.caltech.edu:PowerFlux/scripts/

upload.a:
	globus-scp -C $(UPLOAD_FILES) h2.atlas.aei.uni-hannover.de:PowerFlux/scripts/
