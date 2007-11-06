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
