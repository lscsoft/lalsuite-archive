-- 
-- select process and process_params from test2.ilwd
--
-- to run, first connect to db2 via command line connect 
-- e.g. db2 connect to cit_1 user db2sol7s using <password>
--

echo table process;
select * from process where unix_procid in (949096406,949096407);
echo table process_params;
select * from process_params where param in ('Param_949096406','Param_949096407') ;




