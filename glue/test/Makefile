all : test01 test02 lal_verify segmentsUtils_verify segments_verify iterutils_verify offsetvector_verify ligolw_sqlite_test

offsetvector_verify :
	@echo "=== start offsetvector_verify ===>"
	python offsetvector_verify.py
	@echo "<=== end offsetvector_verify ==="

iterutils_verify :
	@echo "=== start iterutils_verify ===>"
	python iterutils_verify.py
	@echo "<=== end iterutils_verify ==="

lal_verify :
	@echo "=== start lal_verify ===>"
	python lal_verify.py
	@echo "<=== end lal_verify ==="

segmentsUtils_verify :
	@echo "=== start segmentsUtils_verify ===>"
	python segmentsUtils_verify.py
	@echo "<=== end segmentsUtils_verify ==="

segments_verify :
	@echo "=== start segments_verify ===>"
	python segments_verify.py
	@echo "<=== end segments_verify ==="

test01 :
	@echo "=== start test01 ===>"
	python ligo_lw_test_01.py
	@echo "<=== end test01 ==="

test02 :
	@echo "=== start test02 ===>"
	bash ligo_lw_test_02.sh
	@echo "<=== end test02 ==="

ligo_lw_test_03 :
	# test reading compressed and non-compressed files from stdin
	@echo "=== start ligo_lw_test_03 ===>"
	cat inspiral_event_id_test_in1.xml.gz | ligolw_print >/dev/null && echo "Pass" || { echo "Fail"; false; }
	cat inspiral_event_id_test_in2.xml | ligolw_print >/dev/null && echo "Pass" || { echo "Fail"; false; }
	@echo "<=== end ligo_lw_test_03 ==="

ligolw_sqlite_test:
	@echo "=== start ligolw_sqlite_test ===>"
	bash ligolw_sqlite_test.sh
	@echo "<=== end ligolw_sqlite_test ==="

glue_ligolw_ilwd_verify :
	@echo "=== start glue_ligolw_ilwd_verify ===>"
	python glue_ligolw_ilwd_verify.py
	@echo "<=== end glue_ligolw_ilwd_verify ==="

clean :
	rm -f ligo_lw_test_01*png
