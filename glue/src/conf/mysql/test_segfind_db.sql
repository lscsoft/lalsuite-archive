-- Issue some SQL to test a new LSCSegFindServer database
-- Duncan Brown <dbrown@ligo.caltech.edu>
-- $Id$

use segfind;

-- create some human readable state vector information (using our fake SV bitmask values)
insert into state_vec (version,value,state) values(0,0,'Unlocked');
insert into state_vec (version,value,state) values(0,1024,'Up');
insert into state_vec (version,value,state) values(0,2048,'Badgers');
insert into state_vec (version,value,state) values(0,65535,'Science');
insert into state_vec (version,value,state) values(0,65534,'Injection');

-- create some fake information. it's a pain doing this in sql, the python
-- code to do this is much simpler as it can use LAST_INSERT_ID() to get
-- the frame_file_id; here we just select it by querying on the frame name
-- this is the price we pay for wanting to do a pure database test

-- the instrument starts unlocked
insert into frame_file (lfn,start_time,end_time) values('L-R-700000000-32.gwf',700000000,700000032);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000000, 700000032, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000000-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 0;
  
insert into frame_file (lfn,start_time,end_time) values('L-R-700000032-32.gwf',700000032,700000064);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000032, 700000050, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000032-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 0;

-- then they get it up
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000050, 700000064, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000032-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 1024;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000064-32.gwf',700000064,700000096);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000064, 700000075, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000064-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 1024;

-- and into science mode
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000075, 700000096, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000064-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65535;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000096-32.gwf',700000096,700000128);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000096, 700000128, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000096-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65535;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000128-32.gwf',700000128,700000160);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000128, 700000160, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000128-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65535;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000160-32.gwf',700000160,700000192);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000160, 700000192, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000160-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65535;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000192-32.gwf',700000192,700000223);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000192, 700000223, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000192-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65535;

-- now lets do some injections
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000223, 700000224, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000192-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65534;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000224-32.gwf',700000224,700000256);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000224, 700000256, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000224-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65534;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000256-32.gwf',700000256,700000288);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000256, 700000288, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000256-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65534;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000288-32.gwf',700000288,700000320);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000288, 700000320, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000288-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65534;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000320-32.gwf',700000320,700000352);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000320, 700000324, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000320-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65534;

-- and back into science mode
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000324, 700000352, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000320-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65535;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000352-32.gwf',700000352,700000384);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000352, 700000384, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000352-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65535;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000384-32.gwf',700000384,700000416);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000384, 700000416, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000384-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65535;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000416-32.gwf',700000416,700000448);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000416, 700000448, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000416-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65535;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000448-32.gwf',700000448,700000480);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000448, 700000452, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000448-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 65535;

-- here come the badgers!
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000452, 700000453, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000448-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 2048;

-- curse them, lost lock...
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000453, 700000480, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000448-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 0;

insert into frame_file (lfn,start_time,end_time) values('L-R-700000480-32.gwf',700000480,700000512);
insert into state_segment (ifo, start_time, end_time, state_vec_id, frame_file_id)
  SELECT 'L1', 700000480, 700000512, state_vec.state_vec_id, frame_file.frame_file_id
  FROM state_vec, frame_file
  WHERE frame_file.lfn = 'L-R-700000480-32.gwf'
  AND state_vec.version = 0
  AND state_vec.value = 0;

-- now get some information back from the database

-- greg wants the science mode segments from l1
select state_segment.ifo, state_vec.state, state_segment.start_time, state_segment.end_time from state_segment, state_vec where state_segment.ifo = 'L1' and state_vec.state = 'Science' and state_segment.state_vec_id = state_vec.state_vec_id order by state_segment.start_time;

-- inspiral running at livingston wants all science and injection segments
select state_segment.start_time, state_segment.end_time from state_segment, state_vec where state_segment.ifo = 'L1' and ( state_vec.state = 'Science' or state_vec.state = 'Injection' ) and state_segment.state_vec_id = state_vec.state_vec_id order by state_segment.start_time;

-- stuart wants all the frame file lfns for science mode data
select distinct frame_file.lfn from frame_file,state_segment,state_vec where state_vec.state = 'Science' and state_segment.state_vec_id = state_vec.state_vec_id and frame_file.frame_file_id = state_segment.frame_file_id;

-- dan wants all the distinct lfns
select distinct lfn from frame_file;

-- i want to know what interferometers are available 
-- (cf the observatories query in LSCdataFind)
select distinct ifo from state_segment;

-- and i want the differnt states available 
-- (cf the types query in LSCdataFind)
select distinct state from state_vec;

-- what was our duty cycle?

