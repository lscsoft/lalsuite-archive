CREATE TABLE segment
(
-- A "segment" is a time interval which is meaningful for some reason.  For
-- example, it may indicate a period during which an interferometer is locked.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Unique process ID of the process which defined this segment
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Unique segment ID for this segment
      segment_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Segment group (e.g. 'H2-locked') and version to which this segment belongs
      segment_group      VARCHAR(64) NOT NULL,
      version            INTEGER NOT NULL,

-- INFORMATION ABOUT THIS SEGMENT
-- Segment start and end times, in GPS seconds.
      start_time         INTEGER NOT NULL,
      end_time           INTEGER NOT NULL,

-- Activity bit for segment. If zero then this time has been analyzed but
-- the segment should not be applied. If positive then this segment should
-- be added to the list of segments (by a union), if negative then this
-- segment should be removed from the list of segments (by a negation and an
-- intersection).
      active             INTEGER NOT NULL,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

-- Note that (segment_group,version,start_time) should be sufficient to
-- uniquely identify a segment, but we include end_time in the primary key
-- to facilitate faster queries.
      CONSTRAINT segment_pk
      PRIMARY KEY (segment_id, start_time, end_time),

      CONSTRAINT segment_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create a clustering index for quicker scanning of a given segment_id
CREATE INDEX segment_cind ON segment(segment_id) CLUSTER
;
-- Create an index based on time
CREATE INDEX segment_ind_time ON segment(start_time, end_time)
;
-- Create an index based on process_id
CREATE INDEX segment_ind_pid ON segment(process_id)
;
