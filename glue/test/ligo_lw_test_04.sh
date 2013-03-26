# ligolw_cut must produce a byte-identical copy
ligolw_cut --delete-table sngl_inspiral <inspiral_event_id_test_in1.xml.gz | cmp ligolw_cut_proof.xml
