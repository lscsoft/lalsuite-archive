--Beginning of script 2--   DatabaseDB2LUOW (SEG_LHO) [WARNING***Please do not alter this line]--
-- CONNECT TO SEG_LHO USER XXXX using XXXX;
INSERT INTO ASN.IBMQREP_RECVQUEUES
 (repqmapname, recvq, sendq, adminq, capture_alias, capture_schema,
 num_apply_agents, memory_limit, state, description, capture_server)
 VALUES
 ('SEG_CIT_ASN_TO_SEG_LHO_ASN', 'ASN.QM3_TO_QM1.DATAQ',
 'ASN.QM3_TO_QM1.DATAQ', 'ASN.QM3.ADMINQ', 'SEG_CIT', 'ASN', 16, 2,
 'A', '', 'SEG_CIT');
-- COMMIT;