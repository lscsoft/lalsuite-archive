-- Issue some SQL to test a new LSCSegFindServer database
-- Duncan Brown <dbrown@ligo.caltech.edu>
-- $Id$

use segfind;

-- create some frame files
insert into frame_file (lfn,start_time,end_time) values('L-R-700000000-32.gwf',700000000,700000032); -- 1
insert into frame_file (lfn,start_time,end_time) values('L-R-700000032-32.gwf',700000032,700000064); -- 2
insert into frame_file (lfn,start_time,end_time) values('L-R-700000064-32.gwf',700000064,700000096); -- 3
insert into frame_file (lfn,start_time,end_time) values('L-R-700000096-32.gwf',700000096,700000128); -- 4
insert into frame_file (lfn,start_time,end_time) values('L-R-700000128-32.gwf',700000128,700000160); -- 5
insert into frame_file (lfn,start_time,end_time) values('L-R-700000160-32.gwf',700000160,700000192); -- 6
insert into frame_file (lfn,start_time,end_time) values('L-R-700000192-32.gwf',700000192,700000224); -- 7
insert into frame_file (lfn,start_time,end_time) values('L-R-700000224-32.gwf',700000224,700000256); -- 8
insert into frame_file (lfn,start_time,end_time) values('L-R-700000256-32.gwf',700000256,700000288); -- 9
insert into frame_file (lfn,start_time,end_time) values('L-R-700000288-32.gwf',700000288,700000320); -- 10
insert into frame_file (lfn,start_time,end_time) values('L-R-700000320-32.gwf',700000320,700000352); -- 11
insert into frame_file (lfn,start_time,end_time) values('L-R-700000352-32.gwf',700000352,700000384); -- 12
insert into frame_file (lfn,start_time,end_time) values('L-R-700000384-32.gwf',700000384,700000416); -- 13
insert into frame_file (lfn,start_time,end_time) values('L-R-700000416-32.gwf',700000416,700000448); -- 14
insert into frame_file (lfn,start_time,end_time) values('L-R-700000448-32.gwf',700000448,700000480); -- 15
insert into frame_file (lfn,start_time,end_time) values('L-R-700000480-32.gwf',700000480,700000512); -- 16

-- create some state segments
--     for this test we manually insert the frame_file_id, but in the real
--     code this will be fetched by LAST_INSERT_ID()
  -- the instrument starts unlocked
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000000,700000032,0,1);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000032,700000050,0,2);
  -- then they get it up
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000050,700000064,1024,2);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000064,700000075,1024,3);
  -- and into science mode
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000075,700000096,65535,3);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000096,700000128,65535,4);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000128,700000160,65535,5);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000160,700000192,65535,6);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000192,700000224,65535,7);
  -- now lets do some injections
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000224,700000250,65534,8);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000250,700000256,65534,8);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000256,700000288,65534,9);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000288,700000320,65534,10);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000320,700000323,65534,11);
  -- and back into science mode
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000323,700000352,65534,11);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000352,700000384,65535,12);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000384,700000416,65535,13);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000416,700000448,65535,14);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000448,700000480,65535,15);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000480,700000487,65535,16);
  -- oops, lost lock
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('L1',700000487,700000512,0,16);

-- some frames from hanford
insert into frame_file (lfn,start_time,end_time) values('H-R-700000000-32.gwf',700000000,700000032); -- 16
insert into frame_file (lfn,start_time,end_time) values('H-R-700000032-32.gwf',700000032,700000064); -- 17
insert into frame_file (lfn,start_time,end_time) values('H-R-700000064-32.gwf',700000064,700000096); -- 18
insert into frame_file (lfn,start_time,end_time) values('H-R-700000096-32.gwf',700000096,700000128); -- 19
insert into frame_file (lfn,start_time,end_time) values('H-R-700000128-32.gwf',700000128,700000160); -- 20

  -- the instrument starts unlocked
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000000,700000032,0,16);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000032,700000047,0,17);
  -- then they get it up
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000047,700000064,1024,17);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000064,700000075,1024,18);
  -- and into science mode
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000075,700000096,65535,18);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000096,700000128,65535,19);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000128,700000160,65535,20);

  -- H2 is not working 'cos it's crap
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000000,700000032,0,16);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000032,700000064,0,17);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000064,700000096,0,18);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000096,700000096,0,19);
insert into state_segment (ifo,start_time,end_time,state_vec_id,frame_file_id) values('H1',700000128,700000160,0,20);

-- create some human readable state vector information (using our fake SV bitmask values)
insert into state_vec (state_vec_id,state) values(0,'Unlocked');
insert into state_vec (state_vec_id,state) values(1024,'Up');
insert into state_vec (state_vec_id,state) values(65535,'Science');
insert into state_vec (state_vec_id,state) values(65534,'Injection');


-- now get some information back from the database

-- greg wants the science mode segments from l1
select state_segment.ifo, state_vec.state, state_segment.start_time, state_segment.end_time from state_segment, state_vec where state_segment.ifo = 'L1' and state_vec.state = 'Science' and state_segment.state_vec_id = state_vec.state_vec_id order by state_segment.start_time;

-- inspiral running at livingston wants all science and injection segments
select state_segment.start_time, state_segment.end_time from state_segment, state_vec where state_segment.ifo = 'L1' and ( state_vec.state = 'Science' or state_vec.state = 'Injection' ) and state_segment.state_vec_id = state_vec.state_vec_id order by state_segment.start_time;

-- power running at hanford wants the h1 science segments
select state_segment.start_time, state_segment.end_time from state_segment, state_vec where state_segment.ifo = 'H1' and state_vec.state = 'Science' and state_segment.state_vec_id = state_vec.state_vec_id order by state_segment.start_time;

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
