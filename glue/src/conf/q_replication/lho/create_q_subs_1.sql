-- This file is part of the Grid LSC User Environment (GLUE)
-- 
-- GLUE is free software: you can redistribute it and/or modify it under the
-- terms of the GNU General Public License as published by the Free Software
-- Foundation, either version 3 of the License, or (at your option) any later
-- version.
-- 
-- This program is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
-- FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
-- details.
-- 
-- You should have received a copy of the GNU General Public License along
-- with this program.  If not, see <http://www.gnu.org/licenses/>.

--Beginning of script 1--   DatabaseDB2LUOW (SEG_LHO) [WARNING***Please do not alter this line]--
-- CONNECT TO SEG_LHO USER XXXX using XXXX#
CREATE BUFFERPOOL BP8K SIZE 125 PAGESIZE 8192#
CREATE TABLESPACE QPASN PAGESIZE 8192 MANAGED BY SYSTEM USING 
('QPASN_TSC') BUFFERPOOL BP8K#
CREATE TABLE ASN.IBMQREP_DELTOMB
(
 TARGET_OWNER VARCHAR(30) NOT NULL,
 TARGET_NAME VARCHAR(128) NOT NULL,
 VERSION_TIME TIMESTAMP NOT NULL,
 VERSION_NODE SMALLINT NOT NULL,
 KEY_HASH INTEGER NOT NULL,
 PACKED_KEY VARCHAR(4096) FOR BIT DATA NOT NULL
)
 IN QPASN#
ALTER TABLE ASN.IBMQREP_DELTOMB
 VOLATILE CARDINALITY#
CREATE INDEX ASN.IX1DELTOMB ON ASN.IBMQREP_DELTOMB
(
 VERSION_TIME ASC,
 TARGET_NAME ASC,
 TARGET_OWNER ASC,
 KEY_HASH ASC
)#
ALTER TABLE LDBD.PROCESS
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'PROCESS', 'LDBD',
 'PROCESS', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F', 'Q', '000001'
, 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'VERSION', 'VERSION', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CVS_REPOSITORY',
 'CVS_REPOSITORY', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CVS_ENTRY_TIME',
 'CVS_ENTRY_TIME', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'IS_ONLINE', 'IS_ONLINE', 'N'
, 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'NODE', 'NODE', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'USERNAME', 'USERNAME', 'N',
 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'UNIX_PROCID', 'UNIX_PROCID',
 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N',
 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'JOBID', 'JOBID', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'DOMAIN', 'DOMAIN', 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'Y', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PARAM_SET', 'PARAM_SET', 'N'
, 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'IFOS', 'IFOS', 'N', 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 18)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 19)#
CREATE TRIGGER LDBD.APROCESSESS NO CASCADE BEFORE INSERT ON
 LDBD.PROCESS
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BPROCESSESS NO CASCADE BEFORE UPDATE ON
 LDBD.PROCESS
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.PROCESS_PARAMS
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD',
 'PROCESS_PARAMS', 'LDBD', 'PROCESS_PARAMS', 'SEG_LLO', 'SEG_LLO', 1,
 'I', 'P', 'V', 'F', 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PARAM', 'PARAM', 'Y',
 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'TYPE', 'TYPE', 'N', 4
)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'VALUE', 'VALUE', 'N',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 8)#
CREATE TRIGGER LDBD.APROCESS_PARAMAMS NO CASCADE BEFORE INSERT ON
 LDBD.PROCESS_PARAMS
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BPROCESS_PARAMAMS NO CASCADE BEFORE UPDATE ON
 LDBD.PROCESS_PARAMS
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.LFN
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'LFN', 'LDBD', 'LFN',
 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F', 'Q', '000001', 1, 2, 0,
 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB', 'Y',
 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID', 'N',
 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'LFN_ID', 'LFN_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'NAME', 'NAME', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME', 'START_TIME', 'N',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 9)#
CREATE TRIGGER LDBD.ALFNLFN NO CASCADE BEFORE INSERT ON LDBD.LFN
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BLFNLFN NO CASCADE BEFORE UPDATE ON LDBD.LFN
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.GRIDCERT
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('GRIDCERT0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'GRIDCERT', 'LDBD',
 'GRIDCERT', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F', 'Q',
 '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'Y', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0001', 'ASN.QM2_TO_QM1.DATAQ', 'DN', 'DN', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 5)#
CREATE TRIGGER LDBD.AGRIDCERTERT NO CASCADE BEFORE INSERT ON
 LDBD.GRIDCERT
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BGRIDCERTERT NO CASCADE BEFORE UPDATE ON
 LDBD.GRIDCERT
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.FILTER
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'FILTER', 'LDBD',
 'FILTER', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F', 'Q', '000001',
 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'FILTER_NAME', 'FILTER_NAME',
 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'FILTER_ID', 'FILTER_ID', 'Y',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'PARAM_SET', 'PARAM_SET', 'N',
 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 10)#
CREATE TRIGGER LDBD.AFILTERTER NO CASCADE BEFORE INSERT ON LDBD.FILTER
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BFILTERTER NO CASCADE BEFORE UPDATE ON LDBD.FILTER
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.FILTER_PARAMS
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'FILTER_PARAMS'
, 'LDBD', 'FILTER_PARAMS', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F'
, 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'FILTER_NAME',
 'FILTER_NAME', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'FILTER_ID',
 'FILTER_ID', 'Y', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PARAM', 'PARAM', 'Y',
 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'TYPE', 'TYPE', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'VALUE', 'VALUE', 'N',
 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 9)#
CREATE TRIGGER LDBD.AFILTER_PARAMSAMS NO CASCADE BEFORE INSERT ON
 LDBD.FILTER_PARAMS
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BFILTER_PARAMSAMS NO CASCADE BEFORE UPDATE ON
 LDBD.FILTER_PARAMS
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.FRAMESET_CHANLIST
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD',
 'FRAMESET_CHANLIST', 'LDBD', 'FRAMESET_CHANLIST', 'SEG_LLO',
 'SEG_LLO', 1, 'I', 'P', 'V', 'F', 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME',
 'START_TIME', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME',
 'END_TIME', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'CHANLIST_ID',
 'CHANLIST_ID', 'Y', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'CHANLIST',
 'CHANLIST', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'CHANLIST_LENGTH',
 'CHANLIST_LENGTH', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 10)#
CREATE TRIGGER LDBD.AFRAMESET_CHANIST NO CASCADE BEFORE INSERT ON
 LDBD.FRAMESET_CHANLIST
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BFRAMESET_CHANIST NO CASCADE BEFORE UPDATE ON
 LDBD.FRAMESET_CHANLIST
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.FRAMESET_WRITER
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD',
 'FRAMESET_WRITER', 'LDBD', 'FRAMESET_WRITER', 'SEG_LLO', 'SEG_LLO', 1
, 'I', 'P', 'V', 'F', 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'Y', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'DATA_SOURCE',
 'DATA_SOURCE', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'IFOS', 'IFOS', 'N',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'COMMENT', 'COMMENT',
 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 9)#
CREATE TRIGGER LDBD.AFRAMESET_WRITTER NO CASCADE BEFORE INSERT ON
 LDBD.FRAMESET_WRITER
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BFRAMESET_WRITTER NO CASCADE BEFORE UPDATE ON
 LDBD.FRAMESET_WRITER
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.FRAMESET
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'FRAMESET', 'LDBD',
 'FRAMESET', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F', 'Q',
 '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'N', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'CHANLIST_CDB',
 'CHANLIST_CDB', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'CHANLIST_ID', 'CHANLIST_ID'
, 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'Y', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'Y',
 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME_NS', 'END_TIME_NS'
, 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'N_FRAMES', 'N_FRAMES', 'N',
 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'MISSING_FRAMES',
 'MISSING_FRAMES', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'N_BYTES', 'N_BYTES', 'N',
 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'NAME', 'NAME', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 15)#
CREATE TRIGGER LDBD.AFRAMESETSET NO CASCADE BEFORE INSERT ON
 LDBD.FRAMESET
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BFRAMESETSET NO CASCADE BEFORE UPDATE ON
 LDBD.FRAMESET
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.FRAMESET_LOC
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'FRAMESET_LOC',
 'LDBD', 'FRAMESET_LOC', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F',
 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'N', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'NAME', 'NAME', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'MEDIA_TYPE',
 'MEDIA_TYPE', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'NODE', 'NODE', 'Y', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'MEDIA_LOC', 'MEDIA_LOC'
, 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'MEDIA_STATUS',
 'MEDIA_STATUS', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'FULLNAME', 'FULLNAME',
 'Y', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'FILENUM', 'FILENUM',
 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'DECOMPRESS_CMD',
 'DECOMPRESS_CMD', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 11)#
CREATE TRIGGER LDBD.AFRAMESET_LOCLOC NO CASCADE BEFORE INSERT ON
 LDBD.FRAMESET_LOC
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BFRAMESET_LOCLOC NO CASCADE BEFORE UPDATE ON
 LDBD.FRAMESET_LOC
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SEGMENT_DEFINER
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD',
 'SEGMENT_DEFINER', 'LDBD', 'SEGMENT_DEFINER', 'SEG_LLO', 'SEG_LLO', 1
, 'I', 'P', 'V', 'F', 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'RUN', 'RUN', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'IFOS', 'IFOS', 'N',
 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'NAME', 'NAME', 'N',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'VERSION', 'VERSION',
 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'COMMENT', 'COMMENT',
 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'STATE_VEC_MAJOR',
 'STATE_VEC_MAJOR', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'STATE_VEC_MINOR',
 'STATE_VEC_MINOR', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 12)#
CREATE TRIGGER LDBD.ASEGMENT_DEFINNER NO CASCADE BEFORE INSERT ON
 LDBD.SEGMENT_DEFINER
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSEGMENT_DEFINNER NO CASCADE BEFORE UPDATE ON
 LDBD.SEGMENT_DEFINER
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SEGMENT
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'SEGMENT', 'LDBD',
 'SEGMENT', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F', 'Q', '000001'
, 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_ID', 'SEGMENT_ID',
 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME_NS', 'END_TIME_NS',
 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'ACTIVE', 'ACTIVE', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGNUM', 'SEGNUM', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 11)#
CREATE TRIGGER LDBD.ASEGMENTENT NO CASCADE BEFORE INSERT ON
 LDBD.SEGMENT
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSEGMENTENT NO CASCADE BEFORE UPDATE ON
 LDBD.SEGMENT
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SEGMENT_DEF_MAP
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD',
 'SEGMENT_DEF_MAP', 'LDBD', 'SEGMENT_DEF_MAP', 'SEG_LLO', 'SEG_LLO', 1
, 'I', 'P', 'V', 'F', 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEG_DEF_MAP_ID',
 'SEG_DEF_MAP_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_CDB',
 'SEGMENT_CDB', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_ID',
 'SEGMENT_ID', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_CDB',
 'SEGMENT_DEF_CDB', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'STATE_VEC_MAP',
 'STATE_VEC_MAP', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 10)#
CREATE TRIGGER LDBD.ASEGMENT_DEF_MMAP NO CASCADE BEFORE INSERT ON
 LDBD.SEGMENT_DEF_MAP
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSEGMENT_DEF_MMAP NO CASCADE BEFORE UPDATE ON
 LDBD.SEGMENT_DEF_MAP
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SEGMENT_LFN_MAP
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD',
 'SEGMENT_LFN_MAP', 'LDBD', 'SEGMENT_LFN_MAP', 'SEG_LLO', 'SEG_LLO', 1
, 'I', 'P', 'V', 'F', 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEG_LFN_MAP_ID',
 'SEG_LFN_MAP_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_CDB',
 'SEGMENT_CDB', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_ID',
 'SEGMENT_ID', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'LFN_CDB', 'LFN_CDB',
 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'LFN_ID', 'LFN_ID',
 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 9)#
CREATE TRIGGER LDBD.ASEGMENT_LFN_MMAP NO CASCADE BEFORE INSERT ON
 LDBD.SEGMENT_LFN_MAP
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSEGMENT_LFN_MMAP NO CASCADE BEFORE UPDATE ON
 LDBD.SEGMENT_LFN_MAP
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SUMM_VALUE
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'SUMM_VALUE',
 'LDBD', 'SUMM_VALUE', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F',
 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB'
, 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'SUMM_VALUE_ID',
 'SUMM_VALUE_ID', 'Y', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM', 'N',
 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID'
, 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME', 'START_TIME'
, 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME',
 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME_NS',
 'END_TIME_NS', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'IFO', 'IFO', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'NAME', 'NAME', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'VALUE', 'VALUE', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'ERROR', 'ERROR', 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'INTVALUE', 'INTVALUE',
 'N', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N',
 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 18)#
CREATE TRIGGER LDBD.ASUMM_VALUELUE NO CASCADE BEFORE INSERT ON
 LDBD.SUMM_VALUE
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSUMM_VALUELUE NO CASCADE BEFORE UPDATE ON
 LDBD.SUMM_VALUE
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SUMM_STATISTICS
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD',
 'SUMM_STATISTICS', 'LDBD', 'SUMM_STATISTICS', 'SEG_LLO', 'SEG_LLO', 1
, 'I', 'P', 'V', 'F', 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'N', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME',
 'START_TIME', 'Y', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME',
 'END_TIME', 'Y', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME_NS',
 'END_TIME_NS', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMES_USED',
 'FRAMES_USED', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'SAMPLES', 'SAMPLES',
 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'CHANNEL', 'CHANNEL',
 'Y', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'MIN_VALUE',
 'MIN_VALUE', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'MAX_VALUE',
 'MAX_VALUE', 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'MIN_DELTA',
 'MIN_DELTA', 'N', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'MAX_DELTA',
 'MAX_DELTA', 'N', 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'MIN_DELTADELTA',
 'MIN_DELTADELTA', 'N', 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'MAX_DELTADELTA',
 'MAX_DELTADELTA', 'N', 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'MEAN', 'MEAN', 'N',
 18)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'VARIANCE',
 'VARIANCE', 'N', 19)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'RMS', 'RMS', 'N', 20
)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'SKEWNESS',
 'SKEWNESS', 'N', 21)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'KURTOSIS',
 'KURTOSIS', 'N', 22)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 23)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 24)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 25)#
CREATE TRIGGER LDBD.ASUMM_STATISTIICS NO CASCADE BEFORE INSERT ON
 LDBD.SUMM_STATISTICS
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSUMM_STATISTIICS NO CASCADE BEFORE UPDATE ON
 LDBD.SUMM_STATISTICS
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SUMM_SPECTRUM
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'SUMM_SPECTRUM'
, 'LDBD', 'SUMM_SPECTRUM', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F'
, 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'SUMM_SPECTRUM_ID',
 'SUMM_SPECTRUM_ID', 'Y', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME',
 'START_TIME', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME',
 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME_NS',
 'END_TIME_NS', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMES_USED',
 'FRAMES_USED', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_FREQUENCY',
 'START_FREQUENCY', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'DELTA_FREQUENCY',
 'DELTA_FREQUENCY', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'MIMETYPE', 'MIMETYPE',
 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'CHANNEL', 'CHANNEL',
 'N', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'SPECTRUM_TYPE',
 'SPECTRUM_TYPE', 'N', 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'SPECTRUM', 'SPECTRUM',
 'N', 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'SPECTRUM_LENGTH',
 'SPECTRUM_LENGTH', 'N', 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 18)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 19)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 20)#
CREATE TRIGGER LDBD.ASUMM_SPECTRUMRUM NO CASCADE BEFORE INSERT ON
 LDBD.SUMM_SPECTRUM
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSUMM_SPECTRUMRUM NO CASCADE BEFORE UPDATE ON
 LDBD.SUMM_SPECTRUM
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SUMM_CSD
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'SUMM_CSD', 'LDBD',
 'SUMM_CSD', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F', 'Q',
 '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'SUMM_CSD_ID', 'SUMM_CSD_ID'
, 'Y', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM', 'N', 2
)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N',
 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME_NS', 'END_TIME_NS'
, 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMES_USED', 'FRAMES_USED'
, 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_FREQUENCY',
 'START_FREQUENCY', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'DELTA_FREQUENCY',
 'DELTA_FREQUENCY', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'MIMETYPE', 'MIMETYPE', 'N',
 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'CHANNEL1', 'CHANNEL1', 'N',
 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'CHANNEL2', 'CHANNEL2', 'N',
 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'SPECTRUM_TYPE',
 'SPECTRUM_TYPE', 'N', 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'SPECTRUM', 'SPECTRUM', 'N',
 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'SPECTRUM_LENGTH',
 'SPECTRUM_LENGTH', 'N', 18)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 19)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 20)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 21)#
CREATE TRIGGER LDBD.ASUMM_CSDCSD NO CASCADE BEFORE INSERT ON
 LDBD.SUMM_CSD
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSUMM_CSDCSD NO CASCADE BEFORE UPDATE ON
 LDBD.SUMM_CSD
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SUMM_COMMENT
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'SUMM_COMMENT',
 'LDBD', 'SUMM_COMMENT', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F',
 'Q', '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'SUBMITTER', 'SUBMITTER'
, 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME',
 'START_TIME', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME',
 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME_NS',
 'END_TIME_NS', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'IFO', 'IFO', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'TEXT', 'TEXT', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'SUMM_COMMENT_ID',
 'SUMM_COMMENT_ID', 'Y', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 15)#
CREATE TRIGGER LDBD.ASUMM_COMMENTENT NO CASCADE BEFORE INSERT ON
 LDBD.SUMM_COMMENT
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSUMM_COMMENTENT NO CASCADE BEFORE UPDATE ON
 LDBD.SUMM_COMMENT
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
ALTER TABLE LDBD.SUMM_MIME
 ADD "ibmqrepVERTIME" TIMESTAMP NOT NULL WITH DEFAULT
 '0001-01-01-00.00.00'
 ADD "ibmqrepVERNODE" SMALLINT NOT NULL WITH DEFAULT 0
 DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'LDBD', 'SUMM_MIME', 'LDBD'
, 'SUMM_MIME', 'SEG_LLO', 'SEG_LLO', 1, 'I', 'P', 'V', 'F', 'Q',
 '000001', 1, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'ORIGIN', 'ORIGIN', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'FILENAME', 'FILENAME', 'N'
, 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'SUBMITTER', 'SUBMITTER',
 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'SUBMIT_TIME',
 'SUBMIT_TIME', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N'
, 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'END_TIME_NS',
 'END_TIME_NS', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'CHANNEL', 'CHANNEL', 'N',
 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'DESCRIP', 'DESCRIP', 'N',
 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'MIMEDATA', 'MIMEDATA', 'N'
, 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'MIMEDATA_LENGTH',
 'MIMEDATA_LENGTH', 'N', 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'MIMETYPE', 'MIMETYPE', 'N'
, 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N',
 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'SUMM_MIME_ID',
 'SUMM_MIME_ID', 'Y', 18)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 19)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 20)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0001', 'ASN.QM2_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 21)#
CREATE TRIGGER LDBD.ASUMM_MIMEIME NO CASCADE BEFORE INSERT ON
 LDBD.SUMM_MIME
 REFERENCING NEW AS new FOR EACH ROW MODE DB2SQL 
 WHEN (new."ibmqrepVERNODE" = 0)
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = (CURRENT TIMESTAMP - CURRENT TIMEZONE), 
new."ibmqrepVERNODE" = 2*4; 
END#
CREATE TRIGGER LDBD.BSUMM_MIMEIME NO CASCADE BEFORE UPDATE ON
 LDBD.SUMM_MIME
 REFERENCING NEW AS new OLD AS old FOR EACH ROW MODE DB2SQL 
 WHEN ((new."ibmqrepVERTIME" = old."ibmqrepVERTIME")
AND (((new."ibmqrepVERNODE")/4) = ((old."ibmqrepVERNODE")/4))) 
BEGIN ATOMIC 
SET new."ibmqrepVERTIME" = 
CASE 
WHEN (CURRENT TIMESTAMP - CURRENT TIMEZONE) < old."ibmqrepVERTIME"
THEN old."ibmqrepVERTIME" + 00000000000000.000001 
ELSE (CURRENT TIMESTAMP - CURRENT TIMEZONE) 
END, 
new."ibmqrepVERNODE" = 2*4; 
END#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS0002', 'LDBD', 'PROCESS', 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N',
 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO'
, 'LDBD', 'PROCESS', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'VERSION', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'CVS_REPOSITORY', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'CVS_ENTRY_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'IS_ONLINE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'NODE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'USERNAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'UNIX_PROCID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'JOBID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'DOMAIN', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'PROCESS_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'PARAM_SET', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'IFOS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS_PARAMS0002', 'LDBD', 'PROCESS_PARAMS',
 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1,
 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO', 'LDBD', 'PROCESS_PARAMS', 1, 'ASN'
)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0002', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0002', 'PROCESS_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0002', 'PARAM', 3)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0002', 'TYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0002', 'VALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('LFN0002', 'LDBD', 'LFN', 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N'
, 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO', 'LDBD',
 'LFN', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0002', 'LFN_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0002', 'NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0002', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0002', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0002', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('GRIDCERT0002', 'LDBD', 'GRIDCERT', 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N'
, 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO',
 'SEG_LLO', 'LDBD', 'GRIDCERT', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0002', 'CREATOR_DB', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0002', 'PROCESS_ID', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0002', 'DN', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FILTER0002', 'LDBD', 'FILTER', 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N',
 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO'
, 'LDBD', 'FILTER', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'FILTER_NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'FILTER_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'PARAM_SET', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FILTER_PARAMS0002', 'LDBD', 'FILTER_PARAMS', 'ASN.QM1_TO_QM2.DATAQ'
, 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO'
, 'SEG_LLO', 'LDBD', 'FILTER_PARAMS', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0002', 'FILTER_NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0002', 'FILTER_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0002', 'PARAM', 3)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0002', 'TYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0002', 'VALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FRAMESET_CHANLIST0002', 'LDBD', 'FRAMESET_CHANLIST',
 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1,
 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO', 'LDBD', 'FRAMESET_CHANLIST', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'CREATOR_DB', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'CHANLIST_ID', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'CHANLIST', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'CHANLIST_LENGTH', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FRAMESET_WRITER0002', 'LDBD', 'FRAMESET_WRITER',
 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1,
 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO', 'LDBD', 'FRAMESET_WRITER', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0002', 'CREATOR_DB', 3)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0002', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0002', 'PROCESS_ID', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0002', 'FRAMESET_GROUP', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0002', 'DATA_SOURCE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0002', 'IFOS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0002', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FRAMESET0002', 'LDBD', 'FRAMESET', 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N'
, 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO',
 'SEG_LLO', 'LDBD', 'FRAMESET', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'CREATOR_DB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'FRAMESET_GROUP', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'CHANLIST_CDB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'CHANLIST_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'START_TIME', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'END_TIME', 3)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'N_FRAMES', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'MISSING_FRAMES', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'N_BYTES', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FRAMESET_LOC0002', 'LDBD', 'FRAMESET_LOC', 'ASN.QM1_TO_QM2.DATAQ',
 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO',
 'SEG_LLO', 'LDBD', 'FRAMESET_LOC', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'CREATOR_DB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'MEDIA_TYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'NODE', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'MEDIA_LOC', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'MEDIA_STATUS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'FULLNAME', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'FILENUM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'DECOMPRESS_CMD', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_DEFINER0002', 'LDBD', 'SEGMENT_DEFINER',
 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1,
 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO', 'LDBD', 'SEGMENT_DEFINER', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'SEGMENT_DEF_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'RUN', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'IFOS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'VERSION', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'STATE_VEC_MAJOR', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'STATE_VEC_MINOR', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT0002', 'LDBD', 'SEGMENT', 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N',
 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO'
, 'LDBD', 'SEGMENT', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'SEGMENT_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'ACTIVE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'SEGNUM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'LDBD', 'SEGMENT_DEF_MAP',
 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1,
 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO', 'LDBD', 'SEGMENT_DEF_MAP', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'SEG_DEF_MAP_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'SEGMENT_CDB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'SEGMENT_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'SEGMENT_DEF_CDB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'STATE_VEC_MAP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'LDBD', 'SEGMENT_LFN_MAP',
 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1,
 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO', 'LDBD', 'SEGMENT_LFN_MAP', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'SEG_LFN_MAP_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'SEGMENT_CDB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'SEGMENT_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'LFN_CDB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'LFN_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_VALUE0002', 'LDBD', 'SUMM_VALUE', 'ASN.QM1_TO_QM2.DATAQ', 'P',
 'N', 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO',
 'SEG_LLO', 'LDBD', 'SUMM_VALUE', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'SUMM_VALUE_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'IFO', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'VALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'ERROR', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'INTVALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_STATISTICS0002', 'LDBD', 'SUMM_STATISTICS',
 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1,
 'NNNN', 'N', 'SEG_LLO', 'SEG_LLO', 'LDBD', 'SUMM_STATISTICS', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'CREATOR_DB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'START_TIME', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'END_TIME', 3)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'FRAMES_USED', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'SAMPLES', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'CHANNEL', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'MIN_VALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'MAX_VALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'MIN_DELTA', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'MAX_DELTA', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'MIN_DELTADELTA', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'MAX_DELTADELTA', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'MEAN', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'VARIANCE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'RMS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'SKEWNESS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'KURTOSIS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_SPECTRUM0002', 'LDBD', 'SUMM_SPECTRUM', 'ASN.QM1_TO_QM2.DATAQ'
, 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO'
, 'SEG_LLO', 'LDBD', 'SUMM_SPECTRUM', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'SUMM_SPECTRUM_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'FRAMES_USED', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'START_FREQUENCY', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'DELTA_FREQUENCY', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'MIMETYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'CHANNEL', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'SPECTRUM_TYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'SPECTRUM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'SPECTRUM_LENGTH', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_CSD0002', 'LDBD', 'SUMM_CSD', 'ASN.QM1_TO_QM2.DATAQ', 'P', 'N'
, 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO',
 'SEG_LLO', 'LDBD', 'SUMM_CSD', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'SUMM_CSD_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'FRAMES_USED', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'START_FREQUENCY', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'DELTA_FREQUENCY', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'MIMETYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'CHANNEL1', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'CHANNEL2', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'SPECTRUM_TYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'SPECTRUM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'SPECTRUM_LENGTH', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_COMMENT0002', 'LDBD', 'SUMM_COMMENT', 'ASN.QM1_TO_QM2.DATAQ',
 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO',
 'SEG_LLO', 'LDBD', 'SUMM_COMMENT', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'SUBMITTER', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'IFO', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'TEXT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'SUMM_COMMENT_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_MIME0002', 'LDBD', 'SUMM_MIME', 'ASN.QM1_TO_QM2.DATAQ', 'P',
 'N', 'N', 'N', 'E', 'I', '000001', 2, 1, 'NNNN', 'N', 'SEG_LLO',
 'SEG_LLO', 'LDBD', 'SUMM_MIME', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'ORIGIN', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'FILENAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'SUBMITTER', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'SUBMIT_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'CHANNEL', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'DESCRIP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'MIMEDATA', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'MIMEDATA_LENGTH', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'MIMETYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'SUMM_MIME_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0002', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'PROCESS', 'LDBD',
 'PROCESS', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F', 'Q', '000001'
, 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'VERSION', 'VERSION', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CVS_REPOSITORY',
 'CVS_REPOSITORY', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CVS_ENTRY_TIME',
 'CVS_ENTRY_TIME', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'IS_ONLINE', 'IS_ONLINE', 'N'
, 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'NODE', 'NODE', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'USERNAME', 'USERNAME', 'N',
 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'UNIX_PROCID', 'UNIX_PROCID',
 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N',
 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'JOBID', 'JOBID', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'DOMAIN', 'DOMAIN', 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'Y', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PARAM_SET', 'PARAM_SET', 'N'
, 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'IFOS', 'IFOS', 'N', 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 18)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 19)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD',
 'PROCESS_PARAMS', 'LDBD', 'PROCESS_PARAMS', 'SEG_CIT', 'SEG_CIT', 1,
 'I', 'P', 'V', 'F', 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PARAM', 'PARAM', 'Y',
 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'TYPE', 'TYPE', 'N', 4
)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'VALUE', 'VALUE', 'N',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('PROCESS_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 8)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'LFN', 'LDBD', 'LFN',
 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F', 'Q', '000001', 3, 2, 0,
 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB', 'Y',
 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID', 'N',
 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'LFN_ID', 'LFN_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'NAME', 'NAME', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME', 'START_TIME', 'N',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('LFN0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 9)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('GRIDCERT0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'GRIDCERT', 'LDBD',
 'GRIDCERT', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F', 'Q',
 '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'Y', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0003', 'ASN.QM3_TO_QM1.DATAQ', 'DN', 'DN', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('GRIDCERT0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 5)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'FILTER', 'LDBD',
 'FILTER', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F', 'Q', '000001',
 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'FILTER_NAME', 'FILTER_NAME',
 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'FILTER_ID', 'FILTER_ID', 'Y',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'PARAM_SET', 'PARAM_SET', 'N',
 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 10)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'FILTER_PARAMS'
, 'LDBD', 'FILTER_PARAMS', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F'
, 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'FILTER_NAME',
 'FILTER_NAME', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'FILTER_ID',
 'FILTER_ID', 'Y', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PARAM', 'PARAM', 'Y',
 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'TYPE', 'TYPE', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'VALUE', 'VALUE', 'N',
 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FILTER_PARAMS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 9)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD',
 'FRAMESET_CHANLIST', 'LDBD', 'FRAMESET_CHANLIST', 'SEG_CIT',
 'SEG_CIT', 1, 'I', 'P', 'V', 'F', 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME',
 'START_TIME', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME',
 'END_TIME', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'CHANLIST_ID',
 'CHANLIST_ID', 'Y', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'CHANLIST',
 'CHANLIST', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'CHANLIST_LENGTH',
 'CHANLIST_LENGTH', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_CHANLIST0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 10)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD',
 'FRAMESET_WRITER', 'LDBD', 'FRAMESET_WRITER', 'SEG_CIT', 'SEG_CIT', 1
, 'I', 'P', 'V', 'F', 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'Y', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'DATA_SOURCE',
 'DATA_SOURCE', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'IFOS', 'IFOS', 'N',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'COMMENT', 'COMMENT',
 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_WRITER0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 9)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'FRAMESET', 'LDBD',
 'FRAMESET', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F', 'Q',
 '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'N', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'CHANLIST_CDB',
 'CHANLIST_CDB', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'CHANLIST_ID', 'CHANLIST_ID'
, 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'Y', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'Y',
 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME_NS', 'END_TIME_NS'
, 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'N_FRAMES', 'N_FRAMES', 'N',
 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'MISSING_FRAMES',
 'MISSING_FRAMES', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'N_BYTES', 'N_BYTES', 'N',
 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'NAME', 'NAME', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 15)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'FRAMESET_LOC',
 'LDBD', 'FRAMESET_LOC', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F',
 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'N', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'NAME', 'NAME', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'MEDIA_TYPE',
 'MEDIA_TYPE', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'NODE', 'NODE', 'Y', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'MEDIA_LOC', 'MEDIA_LOC'
, 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'MEDIA_STATUS',
 'MEDIA_STATUS', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'FULLNAME', 'FULLNAME',
 'Y', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'FILENUM', 'FILENUM',
 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'DECOMPRESS_CMD',
 'DECOMPRESS_CMD', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('FRAMESET_LOC0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 11)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD',
 'SEGMENT_DEFINER', 'LDBD', 'SEGMENT_DEFINER', 'SEG_CIT', 'SEG_CIT', 1
, 'I', 'P', 'V', 'F', 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'RUN', 'RUN', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'IFOS', 'IFOS', 'N',
 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'NAME', 'NAME', 'N',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'VERSION', 'VERSION',
 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'COMMENT', 'COMMENT',
 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'STATE_VEC_MAJOR',
 'STATE_VEC_MAJOR', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'STATE_VEC_MINOR',
 'STATE_VEC_MINOR', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEFINER0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 12)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'SEGMENT', 'LDBD',
 'SEGMENT', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F', 'Q', '000001'
, 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_ID', 'SEGMENT_ID',
 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N',
 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME_NS', 'END_TIME_NS',
 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'ACTIVE', 'ACTIVE', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGNUM', 'SEGNUM', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 11)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD',
 'SEGMENT_DEF_MAP', 'LDBD', 'SEGMENT_DEF_MAP', 'SEG_CIT', 'SEG_CIT', 1
, 'I', 'P', 'V', 'F', 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEG_DEF_MAP_ID',
 'SEG_DEF_MAP_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_CDB',
 'SEGMENT_CDB', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_ID',
 'SEGMENT_ID', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_CDB',
 'SEGMENT_DEF_CDB', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'STATE_VEC_MAP',
 'STATE_VEC_MAP', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_DEF_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 10)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD',
 'SEGMENT_LFN_MAP', 'LDBD', 'SEGMENT_LFN_MAP', 'SEG_CIT', 'SEG_CIT', 1
, 'I', 'P', 'V', 'F', 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEG_LFN_MAP_ID',
 'SEG_LFN_MAP_ID', 'Y', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_CDB',
 'SEGMENT_CDB', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_ID',
 'SEGMENT_ID', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'LFN_CDB', 'LFN_CDB',
 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'LFN_ID', 'LFN_ID',
 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SEGMENT_LFN_MAP0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 9)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'SUMM_VALUE',
 'LDBD', 'SUMM_VALUE', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F',
 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB'
, 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'SUMM_VALUE_ID',
 'SUMM_VALUE_ID', 'Y', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM', 'N',
 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID'
, 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME', 'START_TIME'
, 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME',
 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME_NS',
 'END_TIME_NS', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'IFO', 'IFO', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'NAME', 'NAME', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'VALUE', 'VALUE', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'ERROR', 'ERROR', 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'INTVALUE', 'INTVALUE',
 'N', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N',
 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_VALUE0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 18)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD',
 'SUMM_STATISTICS', 'LDBD', 'SUMM_STATISTICS', 'SEG_CIT', 'SEG_CIT', 1
, 'I', 'P', 'V', 'F', 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'N', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME',
 'START_TIME', 'Y', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME',
 'END_TIME', 'Y', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME_NS',
 'END_TIME_NS', 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMES_USED',
 'FRAMES_USED', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'SAMPLES', 'SAMPLES',
 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'CHANNEL', 'CHANNEL',
 'Y', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'MIN_VALUE',
 'MIN_VALUE', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'MAX_VALUE',
 'MAX_VALUE', 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'MIN_DELTA',
 'MIN_DELTA', 'N', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'MAX_DELTA',
 'MAX_DELTA', 'N', 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'MIN_DELTADELTA',
 'MIN_DELTADELTA', 'N', 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'MAX_DELTADELTA',
 'MAX_DELTADELTA', 'N', 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'MEAN', 'MEAN', 'N',
 18)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'VARIANCE',
 'VARIANCE', 'N', 19)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'RMS', 'RMS', 'N', 20
)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'SKEWNESS',
 'SKEWNESS', 'N', 21)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'KURTOSIS',
 'KURTOSIS', 'N', 22)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 23)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 24)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_STATISTICS0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 25)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'SUMM_SPECTRUM'
, 'LDBD', 'SUMM_SPECTRUM', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F'
, 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'SUMM_SPECTRUM_ID',
 'SUMM_SPECTRUM_ID', 'Y', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME',
 'START_TIME', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME',
 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME_NS',
 'END_TIME_NS', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMES_USED',
 'FRAMES_USED', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_FREQUENCY',
 'START_FREQUENCY', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'DELTA_FREQUENCY',
 'DELTA_FREQUENCY', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'MIMETYPE', 'MIMETYPE',
 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'CHANNEL', 'CHANNEL',
 'N', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'SPECTRUM_TYPE',
 'SPECTRUM_TYPE', 'N', 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'SPECTRUM', 'SPECTRUM',
 'N', 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'SPECTRUM_LENGTH',
 'SPECTRUM_LENGTH', 'N', 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 18)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 19)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_SPECTRUM0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 20)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'SUMM_CSD', 'LDBD',
 'SUMM_CSD', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F', 'Q',
 '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'SUMM_CSD_ID', 'SUMM_CSD_ID'
, 'Y', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM', 'N', 2
)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N',
 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME_NS', 'END_TIME_NS'
, 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMES_USED', 'FRAMES_USED'
, 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_FREQUENCY',
 'START_FREQUENCY', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'DELTA_FREQUENCY',
 'DELTA_FREQUENCY', 'N', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'MIMETYPE', 'MIMETYPE', 'N',
 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'CHANNEL1', 'CHANNEL1', 'N',
 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'CHANNEL2', 'CHANNEL2', 'N',
 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'SPECTRUM_TYPE',
 'SPECTRUM_TYPE', 'N', 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'SPECTRUM', 'SPECTRUM', 'N',
 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'SPECTRUM_LENGTH',
 'SPECTRUM_LENGTH', 'N', 18)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 19)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 20)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_CSD0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 21)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'SUMM_COMMENT',
 'LDBD', 'SUMM_COMMENT', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F',
 'Q', '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB',
 'CREATOR_DB', 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROGRAM', 'PROGRAM',
 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID',
 'PROCESS_ID', 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'SUBMITTER', 'SUBMITTER'
, 'N', 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME',
 'START_TIME', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME',
 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME_NS',
 'END_TIME_NS', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'IFO', 'IFO', 'N', 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'TEXT', 'TEXT', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'SUMM_COMMENT_ID',
 'SUMM_COMMENT_ID', 'Y', 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_COMMENT0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 15)#
INSERT INTO ASN.IBMQREP_TARGETS
 (subname, recvq, source_owner, source_name, target_owner, target_name
, source_server, source_alias, target_type, state, subtype,
 conflict_rule, conflict_action, error_action, subgroup, source_node,
 target_node, load_type, has_loadphase)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'LDBD', 'SUMM_MIME', 'LDBD'
, 'SUMM_MIME', 'SEG_CIT', 'SEG_CIT', 1, 'I', 'P', 'V', 'F', 'Q',
 '000001', 3, 2, 0, 'E')#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'CREATOR_DB', 'CREATOR_DB',
 'Y', 0)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'ORIGIN', 'ORIGIN', 'N', 1)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'PROCESS_ID', 'PROCESS_ID',
 'N', 2)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'FILENAME', 'FILENAME', 'N'
, 3)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'SUBMITTER', 'SUBMITTER',
 'N', 4)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'SUBMIT_TIME',
 'SUBMIT_TIME', 'N', 5)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'FRAMESET_GROUP',
 'FRAMESET_GROUP', 'N', 6)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'SEGMENT_DEF_ID',
 'SEGMENT_DEF_ID', 'N', 7)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME', 'START_TIME',
 'N', 8)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'START_TIME_NS',
 'START_TIME_NS', 'N', 9)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME', 'END_TIME', 'N'
, 10)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'END_TIME_NS',
 'END_TIME_NS', 'N', 11)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'CHANNEL', 'CHANNEL', 'N',
 12)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'DESCRIP', 'DESCRIP', 'N',
 13)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'MIMEDATA', 'MIMEDATA', 'N'
, 14)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'MIMEDATA_LENGTH',
 'MIMEDATA_LENGTH', 'N', 15)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'MIMETYPE', 'MIMETYPE', 'N'
, 16)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'COMMENT', 'COMMENT', 'N',
 17)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'SUMM_MIME_ID',
 'SUMM_MIME_ID', 'Y', 18)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'INSERTION_TIME',
 'INSERTION_TIME', 'N', 19)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERTIME',
 'ibmqrepVERTIME', 'N', 20)#
INSERT INTO ASN.IBMQREP_TRG_COLS
 (subname, recvq, target_colname, source_colname, is_key, target_colNo)
 VALUES
 ('SUMM_MIME0003', 'ASN.QM3_TO_QM1.DATAQ', 'ibmqrepVERNODE',
 'ibmqrepVERNODE', 'N', 21)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS0004', 'LDBD', 'PROCESS', 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N',
 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT'
, 'LDBD', 'PROCESS', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'VERSION', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'CVS_REPOSITORY', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'CVS_ENTRY_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'IS_ONLINE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'NODE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'USERNAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'UNIX_PROCID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'JOBID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'DOMAIN', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'PROCESS_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'PARAM_SET', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'IFOS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('PROCESS_PARAMS0004', 'LDBD', 'PROCESS_PARAMS',
 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3,
 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT', 'LDBD', 'PROCESS_PARAMS', 1, 'ASN'
)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0004', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0004', 'PROCESS_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0004', 'PARAM', 3)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0004', 'TYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0004', 'VALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('PROCESS_PARAMS0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('LFN0004', 'LDBD', 'LFN', 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N'
, 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT', 'LDBD',
 'LFN', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0004', 'LFN_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0004', 'NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0004', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0004', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0004', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('LFN0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('GRIDCERT0004', 'LDBD', 'GRIDCERT', 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N'
, 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT',
 'SEG_CIT', 'LDBD', 'GRIDCERT', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0004', 'CREATOR_DB', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0004', 'PROCESS_ID', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0004', 'DN', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('GRIDCERT0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FILTER0004', 'LDBD', 'FILTER', 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N',
 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT'
, 'LDBD', 'FILTER', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'FILTER_NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'FILTER_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'PARAM_SET', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FILTER_PARAMS0004', 'LDBD', 'FILTER_PARAMS', 'ASN.QM1_TO_QM3.DATAQ'
, 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT'
, 'SEG_CIT', 'LDBD', 'FILTER_PARAMS', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0004', 'FILTER_NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0004', 'FILTER_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0004', 'PARAM', 3)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0004', 'TYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0004', 'VALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FILTER_PARAMS0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FRAMESET_CHANLIST0004', 'LDBD', 'FRAMESET_CHANLIST',
 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3,
 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT', 'LDBD', 'FRAMESET_CHANLIST', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'CREATOR_DB', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'CHANLIST_ID', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'CHANLIST', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'CHANLIST_LENGTH', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_CHANLIST0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FRAMESET_WRITER0004', 'LDBD', 'FRAMESET_WRITER',
 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3,
 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT', 'LDBD', 'FRAMESET_WRITER', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0004', 'CREATOR_DB', 3)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0004', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0004', 'PROCESS_ID', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0004', 'FRAMESET_GROUP', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0004', 'DATA_SOURCE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0004', 'IFOS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0004', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_WRITER0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FRAMESET0004', 'LDBD', 'FRAMESET', 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N'
, 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT',
 'SEG_CIT', 'LDBD', 'FRAMESET', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'CREATOR_DB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'FRAMESET_GROUP', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'CHANLIST_CDB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'CHANLIST_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'START_TIME', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'END_TIME', 3)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'N_FRAMES', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'MISSING_FRAMES', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'N_BYTES', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('FRAMESET_LOC0004', 'LDBD', 'FRAMESET_LOC', 'ASN.QM1_TO_QM3.DATAQ',
 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT',
 'SEG_CIT', 'LDBD', 'FRAMESET_LOC', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'CREATOR_DB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'MEDIA_TYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'NODE', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'MEDIA_LOC', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'MEDIA_STATUS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'FULLNAME', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'FILENUM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'DECOMPRESS_CMD', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('FRAMESET_LOC0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_DEFINER0004', 'LDBD', 'SEGMENT_DEFINER',
 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3,
 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT', 'LDBD', 'SEGMENT_DEFINER', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'SEGMENT_DEF_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'RUN', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'IFOS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'VERSION', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'STATE_VEC_MAJOR', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'STATE_VEC_MINOR', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEFINER0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT0004', 'LDBD', 'SEGMENT', 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N',
 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT'
, 'LDBD', 'SEGMENT', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'SEGMENT_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'ACTIVE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'SEGNUM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'LDBD', 'SEGMENT_DEF_MAP',
 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3,
 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT', 'LDBD', 'SEGMENT_DEF_MAP', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'SEG_DEF_MAP_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'SEGMENT_CDB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'SEGMENT_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'SEGMENT_DEF_CDB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'STATE_VEC_MAP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_DEF_MAP0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'LDBD', 'SEGMENT_LFN_MAP',
 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3,
 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT', 'LDBD', 'SEGMENT_LFN_MAP', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'SEG_LFN_MAP_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'SEGMENT_CDB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'SEGMENT_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'LFN_CDB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'LFN_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SEGMENT_LFN_MAP0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_VALUE0004', 'LDBD', 'SUMM_VALUE', 'ASN.QM1_TO_QM3.DATAQ', 'P',
 'N', 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT',
 'SEG_CIT', 'LDBD', 'SUMM_VALUE', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'SUMM_VALUE_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'IFO', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'NAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'VALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'ERROR', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'INTVALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_VALUE0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_STATISTICS0004', 'LDBD', 'SUMM_STATISTICS',
 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3,
 'NNNN', 'N', 'SEG_CIT', 'SEG_CIT', 'LDBD', 'SUMM_STATISTICS', 1,
 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'CREATOR_DB', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'START_TIME', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'END_TIME', 3)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'FRAMES_USED', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'SAMPLES', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'CHANNEL', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'MIN_VALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'MAX_VALUE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'MIN_DELTA', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'MAX_DELTA', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'MIN_DELTADELTA', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'MAX_DELTADELTA', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'MEAN', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'VARIANCE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'RMS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'SKEWNESS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'KURTOSIS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_STATISTICS0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_SPECTRUM0004', 'LDBD', 'SUMM_SPECTRUM', 'ASN.QM1_TO_QM3.DATAQ'
, 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT'
, 'SEG_CIT', 'LDBD', 'SUMM_SPECTRUM', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'SUMM_SPECTRUM_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'FRAMES_USED', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'START_FREQUENCY', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'DELTA_FREQUENCY', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'MIMETYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'CHANNEL', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'SPECTRUM_TYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'SPECTRUM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'SPECTRUM_LENGTH', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_SPECTRUM0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_CSD0004', 'LDBD', 'SUMM_CSD', 'ASN.QM1_TO_QM3.DATAQ', 'P', 'N'
, 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT',
 'SEG_CIT', 'LDBD', 'SUMM_CSD', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'SUMM_CSD_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'FRAMES_USED', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'START_FREQUENCY', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'DELTA_FREQUENCY', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'MIMETYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'CHANNEL1', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'CHANNEL2', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'SPECTRUM_TYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'SPECTRUM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'SPECTRUM_LENGTH', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_CSD0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_COMMENT0004', 'LDBD', 'SUMM_COMMENT', 'ASN.QM1_TO_QM3.DATAQ',
 'P', 'N', 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT',
 'SEG_CIT', 'LDBD', 'SUMM_COMMENT', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'PROGRAM', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'SUBMITTER', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'IFO', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'TEXT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'SUMM_COMMENT_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_COMMENT0004', 'ibmqrepVERNODE', 0)#
INSERT INTO ASN.IBMQREP_SUBS
 (subname, source_owner, source_name, sendq, subtype, all_changed_rows
, before_values, changed_cols_only, has_loadphase, state, subgroup,
 source_node, target_node, options_flag, suppress_deletes,
 target_server, target_alias, target_owner, target_name, target_type,
 apply_schema)
 VALUES
 ('SUMM_MIME0004', 'LDBD', 'SUMM_MIME', 'ASN.QM1_TO_QM3.DATAQ', 'P',
 'N', 'N', 'N', 'E', 'I', '000001', 2, 3, 'NNNN', 'N', 'SEG_CIT',
 'SEG_CIT', 'LDBD', 'SUMM_MIME', 1, 'ASN')#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'CREATOR_DB', 2)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'ORIGIN', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'PROCESS_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'FILENAME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'SUBMITTER', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'SUBMIT_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'FRAMESET_GROUP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'SEGMENT_DEF_ID', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'START_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'START_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'END_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'END_TIME_NS', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'CHANNEL', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'DESCRIP', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'MIMEDATA', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'MIMEDATA_LENGTH', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'MIMETYPE', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'COMMENT', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'SUMM_MIME_ID', 1)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'INSERTION_TIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'ibmqrepVERTIME', 0)#
INSERT INTO ASN.IBMQREP_SRC_COLS
 (subname, src_colname, is_key)
 VALUES
 ('SUMM_MIME0004', 'ibmqrepVERNODE', 0)#
-- COMMIT#
