-- Create the LSCsegFindServer database
-- Duncan Brown <dbrown@ligo.caltech.edu>
-- $Id$

DROP DATABASE /*! IF EXISTS */ segfind;
CREATE DATABASE segfind;

USE segfind;

create table state_segment (
-- A state_segment is a time interval for which a given interferometer is in
-- a particular state. For example: science mode, injection mode, etc.

-- Unique id of the state segment.
    state_segment_id integer(11) not null auto_increment,

-- The IFO which this state segment describes.
    ifo char(2) not null,

-- Segment start and end times, in GPS seconds and nanoseconds.
    start_time integer not null,
    start_time_ns integer not null default 0,
    end_time integer not null,
    end_time_ns integer not null default 0,

-- The value of the state vector for this segment.
    state_vec_id integer(11) not null references state_vec(state_vec_id),

-- The frame file from which this state segment was derived.
    frame_file_id integer(11) not null references frame_file(frame_file_id),

-- Insertion time (automatically assigned by the database)
    insertion_time timestamp default null,

    primary key( state_segment_id )
  ) type=myISAM;

create table frame_file (
-- A frame file used to construct state information
    frame_file_id integer(11) not null auto_increment,

-- The logial file name of the frame file
    lfn varchar(255) not null,

-- The start and end times of the frame file, in GPS seconds
    start_time integer not null,
    end_time integer not null,

-- Insertion time (automatically assigned by the database)
    insertion_time timestamp default null,

    primary key( frame_file_id )
  ) type=myISAM;

create table state_vec (
-- Human readable representations of the state vector bits.

-- The value of the state vector
    state_vec_id integer(11) not null,

-- Description of the IFO state
    state varchar(255),

-- Insertion time (automatically assigned by the database)
    insertion_time timestamp default null,

    primary key (state_vec_id)
  ) type=myISAM;
