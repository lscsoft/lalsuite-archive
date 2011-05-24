all : test01 test02 lal_verify segmentsUtils_verify segments_verify iterutils_verify

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

clean :
	rm -f ligo_lw_test_01*png
