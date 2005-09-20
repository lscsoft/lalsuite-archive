-- 
-- Give only "select" privilege to the public
--
-- to run, first connect to db2 via command line connect 
-- e.g. db2 connect to cit_1 user db2sol7s using <password>
--
echo granting priveledges tables;
-- ----------------------------------------------------------------------
-- misc tables
-- ----------------------------------------------------------------------
echo table calib_info;
revoke all on table calib_info from public;
grant select on table calib_info to public;
-- ----------------------------------------------------------------------
echo table runlist;
revoke all on table runlist from public;
grant select on table runlist to public;
-- ----------------------------------------------------------------------
echo table search_summvars;
revoke all on table search_summvars from public;
grant select on table search_summvars to public;
-- ----------------------------------------------------------------------
echo table search_summary;
revoke all on table search_summary from public;
grant select on table search_summary to public;
-- ----------------------------------------------------------------------
-- multi events
-- ----------------------------------------------------------------------
echo table exttrig_search;
revoke all on table exttrig_search from public;
grant select on table exttrig_search to public;
-- ----------------------------------------------------------------------
echo table coinc_sngl;
revoke all on table coinc_sngl from public;
grant select on table coinc_sngl to public;
-- ----------------------------------------------------------------------
echo table multi_burst;
revoke all on table multi_burst from public;
grant select on table multi_burst to public;
-- ----------------------------------------------------------------------
echo table multi_inspiral;
revoke all on table multi_inspiral from public;
grant select on table multi_inspiral to public;
-- ----------------------------------------------------------------------
-- single events 
-- ----------------------------------------------------------------------
echo table waveburst_mime;
revoke all on table waveburst_mime from public;
grant select on table waveburst_mime to public;
-- ----------------------------------------------------------------------
echo table waveburst;
revoke all on table waveburst from public;
grant select on table waveburst to public;
-- ----------------------------------------------------------------------
echo table gds_trigger;
revoke all on table gds_trigger from public;
grant select on table gds_trigger to public;
-- ----------------------------------------------------------------------
echo table sngl_burst;
revoke all on table sngl_burst from public;
grant select on table sngl_burst to public;
-- ----------------------------------------------------------------------
echo table sngl_block;
revoke all on table sngl_block from public;
grant select on table sngl_block to public;
-- ----------------------------------------------------------------------
echo table sngl_dperiodic;
revoke all on table sngl_dperiodic from public;
grant select on table sngl_dperiodic to public;
-- ----------------------------------------------------------------------
echo table sngl_inspiral;
revoke all on table sngl_inspiral from public;
grant select on table sngl_inspiral to public;
-- ----------------------------------------------------------------------
echo table sngl_ringdown;
revoke all on table sngl_ringdown from public;
grant select on table sngl_ringdown to public;
-- ----------------------------------------------------------------------
echo table sngl_unmodeled;
revoke all on table sngl_unmodeled from public;
grant select on table sngl_unmodeled to public;
-- ----------------------------------------------------------------------
echo table sngl_unmodeled_v;
revoke all on table sngl_unmodeled_v from public;
grant select on table sngl_unmodeled_v to public;
-- ----------------------------------------------------------------------
echo table sngl_mime;
revoke all on table sngl_mime from public;
grant select on table sngl_mime to public;
-- ----------------------------------------------------------------------
echo table sngl_transdata;
revoke all on table sngl_transdata from public;
grant select on table sngl_transdata to public;
-- ----------------------------------------------------------------------
echo table sngl_datasource;
revoke all on table sngl_datasource from public;
grant select on table sngl_datasource to public;
-- ----------------------------------------------------------------------
-- summary info
-- ----------------------------------------------------------------------
echo table summ_mime;
revoke all on table summ_mime from public;
grant select on table summ_mime to public;
-- ----------------------------------------------------------------------
echo table summ_comment;
revoke all on table summ_comment from public;
grant select on table summ_comment to public;
-- ----------------------------------------------------------------------
echo table summ_spectrum;
revoke all on table summ_spectrum from public;
grant select on table summ_spectrum to public;
-- ----------------------------------------------------------------------
echo table summ_csd;
revoke all on table summ_csd from public;
grant select on table summ_csd to public;
-- ----------------------------------------------------------------------
echo table summ_statistics;
revoke all on table summ_statistics from public;
grant select on table summ_statistics to public;
-- ----------------------------------------------------------------------
echo table summ_value;
revoke all on table summ_value from public;
grant select on table summ_value to public;
-- ----------------------------------------------------------------------
echo table summ_csd;
revoke all on table summ_csd from public;
grant select on table summ_csd to public;
-- ----------------------------------------------------------------------
-- ----------------------------------------------------------------------
echo table frameset_loc;
revoke all on table frameset_loc from public;
grant select on table frameset_loc to public;
-- ----------------------------------------------------------------------
echo table frameset;
revoke all on table frameset from public;
grant select on table frameset to public;
-- ----------------------------------------------------------------------
echo table frameset_chanlist;
revoke all on table frameset_chanlist from public;
grant select on table frameset_chanlist to public;
-- ----------------------------------------------------------------------
echo table frameset_writer;
revoke all on table frameset_writer from public;
grant select on table frameset_writer to public;
-- ----------------------------------------------------------------------
echo table segment;
revoke all on table segment from public;
grant select on table segment to public;
-- ----------------------------------------------------------------------
echo table segment_definer;
revoke all on table segment_definer from public;
grant select on table segment_definer to public;
-- ----------------------------------------------------------------------
-- simulation tables
-- ----------------------------------------------------------------------
echo table sim_inst_params;
revoke all on table sim_inst_params from public;
grant select on table sim_inst_params to public;
-- ----------------------------------------------------------------------
echo table sim_inst;
revoke all on table sim_inst from public;
grant select on table sim_inst to public;
-- ----------------------------------------------------------------------
echo table sim_type_params;
revoke all on table sim_type_params from public;
grant select on table sim_type_params to public;
-- ----------------------------------------------------------------------
echo table sim_type;
revoke all on table sim_type from public;
grant select on table sim_type to public;
-- ----------------------------------------------------------------------
-- filter and params
-- ----------------------------------------------------------------------
echo table filter;
revoke all on table filter from public;
grant select on table filter to public;
-- ----------------------------------------------------------------------
echo table filter_params;
revoke all on table filter_params from public;
grant select on table filter_params to public;
-- ----------------------------------------------------------------------
echo table process_params;
revoke all on table process_params from public;
grant select on table process_params to public;
-- ----------------------------------------------------------------------
echo table gridcert;
revoke all on table gridcert from public;
grant select on table gridcert to public;
-- ----------------------------------------------------------------------
echo table lfn;
revoke all on table lfn from public;
grant select on table lfn to public;
-- ----------------------------------------------------------------------
-- handle process table last
-- ----------------------------------------------------------------------
echo table process;
revoke all on table process from public;
grant select on table process to public;
-- ----------------------------------------------------------------------
