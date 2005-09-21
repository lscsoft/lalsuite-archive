CREATE TABLE summ_value
(
-- Table to record a value about a particular time interval

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH RECORDED THE VALUE
-- Program name
      program            CHAR(16) NOT NULL,
-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- TIME INTERVAL FROM WHICH THIS VALUE WAS CALCULATED
-- Group name for frameset which determined this time interval, if any
      frameset_group     VARCHAR(48),
-- segment type which determined this time interval, if any
      segment_def_id     CHAR(13) FOR BIT DATA,
-- Start and end times (in GPS seconds and nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,

-- THE SUMMARY VALUE
-- Site or interferometer to which this applies (H0, H1, H2, L0, L1)
      ifo                CHAR(2) NOT NULL,
-- Descriptive name
      name               VARCHAR(128) NOT NULL,
-- The value itself (must be a real number)
      value              REAL,
-- Optional uncertainty on the value
      error              REAL,
-- An optional 4-byte integer value or bitmask
      intvalue           INTEGER,
-- Optional comment
      comment            VARCHAR(80),
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,
  
      CONSTRAINT summval_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db),

-- If segment_group or frameset_group is non-null, make sure there is a
-- corresponding entry in the appropriate table.  If null, then no
-- foreign-key check is performed.
      CONSTRAINT summval_fk_seg
      FOREIGN KEY (segment_def_id)
          REFERENCES segment_definer(segment_def_id),

      CONSTRAINT summval_fk_fs
      FOREIGN KEY (frameset_group, start_time, end_time)
          REFERENCES frameset(frameset_group, start_time, end_time)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time
CREATE INDEX summval_ind_time ON summ_value(start_time, end_time)
;
-- Create an index based on program name
CREATE INDEX summval_ind_prog ON summ_value(program, start_time, name)
;
-- Create an index based on process_id
CREATE INDEX summval_ind_pid ON summ_value(process_id)
;
-- Create an index based on frameset_group
CREATE INDEX summval_ind_fsg ON summ_value(frameset_group)
;
-- Create an index based on segment_group
CREATE INDEX summval_ind_sgrp ON summ_value(segment_def_id)
;
