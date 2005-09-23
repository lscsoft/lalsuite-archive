CREATE TABLE segment_def_map
(
-- Create a map between segments and segment definitions

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Process ID of the program which created this map
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- ID of the segment
      segment_id         CHAR(13) FOR BIT DATA NOT NULL,

-- ID of the segment definer
      segment_def_id     CHAR(13) FOR BIT DATA NOT NULL,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT segdefmap_pk
      PRIMARY KEY (process_id, segment_id, segment_def_id),

      CONSTRAINT segdefmap_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)

      CONSTRAINT segdefmap_fk_sid
      FOREIGN KEY (segment_id)
          REFERENCES segment(segment_id),

      CONSTRAINT segdefmap_fk_sdid
      FOREIGN KEY (segment_def_id)
          REFERENCES segment_definer(segment_def_id)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on process ID
CREATE INDEX segdefmap_pid on segment_def_map(process_id);
;
-- Create an index based on segment ID
CREATE INDEX segdefmap_sid on segment_def_map(segment_id);
;
-- Create an index based on segment definintion ID
CREATE INDEX segdefmap_sdid on segment_def_map(segment_def_id);
;
