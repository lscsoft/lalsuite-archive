CREATE TABLE segment_definer
(
-- List of processes which define segments.  Note that multiple processes
-- can define segments in the same segment_group.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROGRAM WHICH IS DEFINING SEGMENTS
-- Program name
      program            CHAR(16) NOT NULL,
-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- INFORMATION ABOUT THE SEGMENTS BEING DEFINED
-- Descriptive name for this group of segments (e.g. 'H2-locked')
      segment_group      VARCHAR(64) NOT NULL,
-- Version number for segment group (to allow re-evaluation)
      version            INTEGER NOT NULL,
-- Interferometer(s) for which these segments are meaningful
      ifos               CHAR(12),

-- Optional user comment about this segment_group
      comment            VARCHAR(240),
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT segdef_pk
      PRIMARY KEY (process_id, creator_db, segment_group, version),

      CONSTRAINT segdef_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create a clustering index for quicker scanning of a given segment_group
CREATE INDEX segdef_cind ON segment_definer(segment_group, version) CLUSTER
;
-- Create an index based on program name
CREATE INDEX segdef_ind_program ON segment_definer(program)
;
