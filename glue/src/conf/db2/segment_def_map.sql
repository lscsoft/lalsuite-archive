CREATE TABLE segment_def_map
(
-- Create a map between segments and segment definitions

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Process ID of the program which created this map
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Unique ID for this map
      seg_def_map_id     CHAR(13) FOR BIT DATA NOT NULL,

-- ID of the segment
      segment_cdb        INTEGER NOT NULL WITH DEFAULT 1,
      segment_id         CHAR(13) FOR BIT DATA NOT NULL,

-- ID of the segment definer
      segment_def_cdb    INTEGER NOT NULL WITH DEFAULT 1,
      segment_def_id     CHAR(13) FOR BIT DATA NOT NULL,

-- Flag to indicate if this is a state vector map
      state_vec_map      INTEGER,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT segdefmap_pk
      PRIMARY KEY (seg_def_map_id, creator_db),

-- This map must be created by some process
      CONSTRAINT segdefmap_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE,

-- It must point to a segment
      CONSTRAINT segdefmap_fk_sid
      FOREIGN KEY (segment_id, segment_cdb)
          REFERENCES segment(segment_id, creator_db),

-- It must point to a definition
      CONSTRAINT segdefmap_fk_sdid
      FOREIGN KEY (segment_def_id, segment_def_cdb)
          REFERENCES segment_definer(segment_def_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on segment ID
CREATE INDEX segdefmap_sid on segment_def_map(segment_id, segment_cdb)
;
-- Create an index based on segment definintion ID
CREATE INDEX segdefmap_sdid on segment_def_map(segment_def_id, segment_def_cdb)
;
