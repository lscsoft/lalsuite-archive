-- Create the LSCsegFindServer database
-- Duncan Brown <dbrown@ligo.caltech.edu>
-- $Id$
--
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

-- The file from which this state segment was derived.
    lfn_id integer(11) not null references lfn(lfn_id),

-- Insertion time (automatically assigned by the database)
    insertion_time timestamp default null,

    primary key( ifo, start_time, start_time_ns, end_time, end_time_ns, state_vec_id ),
    unique index( state_segment_id ),
    index( start_time, end_time, state_vec_id )
  ) type=myISAM;

create table lfn (
-- A LFN from which some state information is derived

-- A unique identifier of this LFN
    lfn_id integer(11) not null auto_increment,

-- The logial file name of the file
    lfn varchar(255) not null,

-- The start and end times of the file, in GPS seconds
    start_time integer not null,
    end_time integer not null,

-- Insertion time (automatically assigned by the database)
    insertion_time timestamp default null,

    primary key( lfn, start_time, end_time ),
    index( lfn ),
    unique index( lfn_id )

  ) type=myISAM;

create table state_vec (
-- Human readable representations of the state vector bits.

-- Unique identifier for the state vector version and state
    state_vec_id integer(11) not null auto_increment,

-- Version number of the state vector
    version integer(11) not null,
    
-- The value of the state vector
    value integer(11) not null,

-- Description of the IFO state
    state varchar(255),

-- Insertion time (automatically assigned by the database)
    insertion_time timestamp default null,

    primary key( version, value ),
    unique index( state_vec_id ),
    unique index( version, value, state ),
    index( state )

  ) type=myISAM;
