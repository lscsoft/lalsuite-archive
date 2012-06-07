all : glue_ligolw_ilwd_verify iterutils_verify lal_verify ligo_lw_test01 ligo_lw_test02 ligo_lw_test_03 ligolw_sqlite_test offsetvector_verify segmentsUtils_verify segments_verify

define printpassfail
	{ echo "Pass" ; true ; } || { echo "Fail" ; false ; }
endef

glue_ligolw_ilwd_verify :
	@echo "=== start glue_ligolw_ilwd_verify ===>"
	python glue_ligolw_ilwd_verify.py && $(printpassfail)
	@echo "<=== end glue_ligolw_ilwd_verify ==="

iterutils_verify :
	@echo "=== start iterutils_verify ===>"
	python iterutils_verify.py && $(printpassfail)
	@echo "<=== end iterutils_verify ==="

lal_verify :
	@echo "=== start lal_verify ===>"
	python lal_verify.py && $(printpassfail)
	@echo "<=== end lal_verify ==="

ligo_lw_test01 :
	@echo "=== start test01 ===>"
	python ligo_lw_test_01.py && $(printpassfail)
	@echo "<=== end test01 ==="

ligo_lw_test02 :
	@echo "=== start test02 ===>"
	sh ligo_lw_test_02.sh && $(printpassfail)
	@echo "<=== end test02 ==="

ligo_lw_test_03 :
	# test reading compressed and non-compressed files from stdin
	@echo "=== start ligo_lw_test_03 ===>"
	{ cat inspiral_event_id_test_in1.xml.gz | ligolw_print >/dev/null ; } && $(printpassfail)
	{ cat inspiral_event_id_test_in2.xml | ligolw_print >/dev/null ; } && $(printpassfail)
	@echo "<=== end ligo_lw_test_03 ==="

ligolw_sqlite_test :
	@echo "=== start ligolw_sqlite_test ===>"
	sh ligolw_sqlite_test.sh && $(printpassfail)
	@echo "<=== end ligolw_sqlite_test ==="

offsetvector_verify :
	@echo "=== start offsetvector_verify ===>"
	python offsetvector_verify.py && $(printpassfail)
	@echo "<=== end offsetvector_verify ==="

segmentsUtils_verify :
	@echo "=== start segmentsUtils_verify ===>"
	python segmentsUtils_verify.py && $(printpassfail)
	@echo "<=== end segmentsUtils_verify ==="

segments_verify :
	@echo "=== start segments_verify ===>"
	python segments_verify.py && $(printpassfail)
	@echo "<=== end segments_verify ==="

clean :
	rm -f ligo_lw_test_01*png
