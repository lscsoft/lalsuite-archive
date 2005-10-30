--Beginning of script 1--   DatabaseDB2LUOW (SEG_LLO) [WARNING***Please do not alter this line]--
-- CONNECT TO SEG_LLO USER XXXX using XXXX;
INSERT INTO ASN.IBMQREP_SENDQUEUES
 (pubqmapname, sendq, message_format, msg_content_type, state,
 error_action, heartbeat_interval, max_message_size, description,
 apply_alias, apply_schema, recvq, apply_server)
 VALUES
 ('SEG_LLO_ASN_TO_SEG_LHO_ASN', 'ASN.QM2_TO_QM1.DATAQ', 'C', 'T', 'A',
 'S', 60, 64, '', 'SEG_LHO', 'ASN', 'ASN.QM2_TO_QM1.DATAQ', 'SEG_LHO');
-- COMMIT;