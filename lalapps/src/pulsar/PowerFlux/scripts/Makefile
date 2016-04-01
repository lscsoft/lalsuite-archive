
UPLOAD_FILES=make_run.tcl split.tcl get_param.tcl make_1800_sfts.tcl \
	get_lines.tcl make_power_injection_run.tcl split_pi.tcl reshape.tcl \
	generate_injection_batch.tcl make_sft_injection_run.tcl split_dataset.tcl \
	split_sfti.tcl make_second_sft_injection_run.tcl normalize_sft_name.tcl \
	make_native_sft_injection_run.tcl make_native_sft_injection_run2.tcl reshape_sfti.tcl archive.tcl keyed_split.tcl \
	make_split_run.tcl find_outlier_matches.tcl find_outlier_matches2.tcl find_outlier_matches3.tcl make_outlier_histogram.tcl analyze_coinc.R \
	pick_job.tcl shuffle.tcl analyze_coinc.S6.R

all:


upload:
	rsync -avz --rsh=gsissh $(UPLOAD_FILES) pcdev2.phys.uwm.edu:PowerFlux/scripts/

upload.h:
	gsiscp -C $(UPLOAD_FILES) ldas-grid.ligo-wa.caltech.edu:PowerFlux/scripts/

upload.l:
	gsiscp -C $(UPLOAD_FILES) ldas-grid.ligo-la.caltech.edu:PowerFlux/scripts/

upload.c:
	rsync -avz --rsh=gsissh $(UPLOAD_FILES) ldas-grid.ligo.caltech.edu:PowerFlux/scripts/

upload.a:
	rsync -avz --rsh=gsissh $(UPLOAD_FILES) atlas8.atlas.aei.uni-hannover.de:PowerFlux/scripts/
