CREATE TABLE segment_map
(
-- Create a map between segments and segment definitions

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Segment ID
      segment_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Segment definition ID
      segment_def_id     CHAR(13) FOR BIT DATA NOT NULL,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT segmap_pk
      PRIMARY KEY (process_id, segment_id, segment_def_id),

      CONSTRAINT segdef_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on process ID
CREATE INDEX segmap_pid on segment_map(process_id);
;
-- Create an index based on segment ID
CREATE INDEX segmap_sid on segment_map(segment_id);
;
-- Create an index based on segment definintion ID
CREATE INDEX segmap_sdid on segment_map(segment_def_id);
;
