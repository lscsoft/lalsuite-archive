CREATE TABLE summ_csd
(
-- Table to contain a cross-spectral density derived from the data in
-- a specific time interval.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH PRODUCED THIS CSD
-- Program name
      program            CHAR(16) NOT NULL,
-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- TIME INTERVAL FROM WHICH THIS CSD WAS DERIVED
-- Group name for frameset which determined this time interval, if any
      frameset_group     VARCHAR(48),
-- Group and version of segment which determined this time interval, if any
      segment_group      VARCHAR(64),
      version            INTEGER,
-- Start and end times (in GPS seconds and nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,
-- Number of frames actually used to create CSD
      frames_used        INTEGER,
-- Start and end frequencies (stop frequency is derived from these)
      start_frequency    DOUBLE NOT NULL,
      delta_frequency    DOUBLE NOT NULL,
      mimetype        	 VARCHAR(64) NOT NULL,
	  
-- CHANNEL (OR PSEUDO-CHANNEL) NAMES
-- The channel/pseudo-channel names should indicate the interferometer or site
      channel1           VARCHAR(240) NOT NULL,
      channel2           VARCHAR(240) NOT NULL,

-- CSD DATA AND ASSOCIATED INFO
-- Spectrum type (descriptive name)
      spectrum_type      VARCHAR(128) NOT NULL,
-- The spectrum itself is stored in a Binary Large OBject (BLOB).
-- We specify COMPACT since we do not expect this ever to be updated.
      spectrum           BLOB(1M) COMPACT NOT NULL,
-- Length of the spectrum (in bytes)
      spectrum_length    INTEGER NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT summcsd_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db),

-- If segment_group or frameset_group is non-null, make sure there is a
-- corresponding entry in the appropriate table.  If null, then no
-- foreign-key check is performed.
      CONSTRAINT summcsd_fk_seg
      FOREIGN KEY (segment_group, version, start_time, end_time)
          REFERENCES segment(segment_group, version, start_time, end_time),

      CONSTRAINT summcsd_fk_fs
      FOREIGN KEY (frameset_group, start_time, end_time)
          REFERENCES frameset(frameset_group, start_time, end_time)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create a clustering index based on spectrum_type
CREATE INDEX summcsd_cind ON summ_csd(spectrum_type) CLUSTER
;
-- Create an index based on channel names
CREATE INDEX summcsd_ind_chan ON summ_csd(channel1, channel2)
;
-- Create an index based on time
CREATE INDEX summcsd_ind_time ON summ_csd(start_time, end_time)
;
-- Create an index based on process_id
CREATE INDEX summcsd_ind_pid ON summ_csd(process_id)
;
-- Create an index based on frameset_group
CREATE INDEX summcsd_ind_fsg ON summ_csd(frameset_group)
;
-- Create an index based on segment_group
CREATE INDEX summcsd_ind_sgrp ON summ_csd(segment_group, version)
;
