CREATE TABLE summ_comment
(
-- Table to attach a comment to a particular time interval

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH ADDED THE COMMENT (may be null)
-- Program name
      program            CHAR(16),
-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA,

-- ORIGIN OF THE COMMENT
-- Name of person who made comment
      submitter          VARCHAR(48) NOT NULL,

-- TIME INTERVAL TO WHICH THIS COMMENT APPLIES
-- Group name for frameset which determined this time interval, if any
      frameset_group     VARCHAR(48),
-- Group and version of segment which determined this time interval, if any
       segment_def_id    CHAR(13) FOR BIT DATA,
-- Start and end times (in GPS seconds and nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,

-- COMMENT AND ASSOCIATED INFO
-- Interferometer or site to which the comment applies (if appropriate)
      ifo                CHAR(2),
-- The comment itself
      text               VARCHAR(1000) NOT NULL,
-- Unique identifier for this comment (needed for primary key)
      summ_comment_id    CHAR(13) FOR BIT DATA NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT summcomm_pk
      PRIMARY KEY (summ_comment_id, creator_db),

-- Note that process_id is allowed to be null, in which case no check is made.
      CONSTRAINT summcomm_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db),

-- If segment_group or frameset_group is non-null, make sure there is a
-- corresponding entry in the appropriate table.  If null, then no
-- foreign-key check is performed.
      CONSTRAINT summcomm_fk_seg
      FOREIGN KEY (segment_def_id)
          REFERENCES segment_definer(segment_def_id),

      CONSTRAINT summcomm_fk_fs
      FOREIGN KEY (frameset_group, start_time, end_time)
          REFERENCES frameset(frameset_group, start_time, end_time)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time the comment refers to
CREATE INDEX summcomm_ind_time ON summ_comment(start_time, end_time)
;
-- Create an index based on frameset_group
CREATE INDEX summcomm_ind_fsg ON summ_comment(frameset_group)
;
-- Create an index based on segment_group
CREATE INDEX summcomm_ind_sgrp ON summ_comment(segment_def_id)
;
