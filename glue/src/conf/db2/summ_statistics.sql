CREATE TABLE summ_statistics
(
-- Table to contain minimum, maximum, mean, rms, etc. for a single channel
-- (or pseudo-channel) for a specific time interval.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH PRODUCED THESE STATISTICS
-- Program name
      program            CHAR(16) NOT NULL,
-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- TIME INTERVAL FOR WHICH THESE STATISTICS WERE CALCULATED
-- Group name for frameset which determined this time interval, if any
      frameset_group     VARCHAR(48),
-- Group and version of segment which determined this time interval, if any
      segment_def_id     CHAR(13) FOR BIT DATA,
-- Start and end times (in GPS seconds and nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,
-- Number of frames actually used to calculate statistics
      frames_used        INTEGER,
-- Number of samples from which sums were accumulated (needed so that one
-- can convert between raw sums and mean, rms, etc.)
      samples            INTEGER NOT NULL,

-- CHANNEL (OR PSEUDO-CHANNEL) NAME
-- The channel name should indicate the interferometer or site
      channel            VARCHAR(240) NOT NULL,

-- STATISTICS INFO
-- Minimum and maximum value of the channel during this time interval
      min_value          DOUBLE,
      max_value          DOUBLE,
-- Minimum and maximum CHANGE in the value of the channel
      min_delta          DOUBLE,
      max_delta          DOUBLE,
-- Minimum and maximum second-order finite difference
      min_deltadelta     DOUBLE,
      max_deltadelta     DOUBLE,
-- Mean, rms, etc.
      mean               DOUBLE,
      variance           DOUBLE,
      rms                DOUBLE,
      skewness           DOUBLE,
      kurtosis           DOUBLE,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT summstat_pk
      PRIMARY KEY (channel, start_time, end_time),

      CONSTRAINT summstat_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db),

-- If segment_group or frameset_group is non-null, make sure there is a
-- corresponding entry in the appropriate table.  If null, then no
-- foreign-key check is performed.
      CONSTRAINT summstat_fk_seg
      FOREIGN KEY (segment_def_id,creator_db)
          REFERENCES segment_definer(segment_def_id,creator_db),

      CONSTRAINT summstat_fk_fs
      FOREIGN KEY (frameset_group, start_time, end_time)
          REFERENCES frameset(frameset_group, start_time, end_time)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time
CREATE INDEX summstat_ind_time ON summ_statistics(start_time, end_time)
;
-- Create an index based on channel
CREATE INDEX summstat_ind_chan ON summ_statistics(channel, start_time)
;
-- Create an index based on process_id
CREATE INDEX summstat_ind_pid ON summ_statistics(process_id)
;
-- Create an index based on frameset_group
CREATE INDEX summstat_ind_fsg ON summ_statistics(frameset_group)
;
-- Create an index based on segment_group
CREATE INDEX summstat_ind_sgrp ON summ_statistics(segment_def_id)
;
