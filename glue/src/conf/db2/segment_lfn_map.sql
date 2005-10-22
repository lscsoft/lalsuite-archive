CREATE TABLE segment_lfn_map
(
-- Create a map between segments and logical file names (LFNs)

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Process ID of the program which created this map
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- ID of the segment
      segment_id         CHAR(13) FOR BIT DATA NOT NULL,

-- ID of the logiacal file name
      lfn_id             CHAR(13) FOR BIT DATA NOT NULL,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT seglfnmap_pk
      PRIMARY KEY (creator_db, process_id, segment_id, lfn_id),

      CONSTRAINT seglfnmap_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db),

      CONSTRAINT seglfnmap_fk_sid
      FOREIGN KEY (segment_id, creator_db)
          REFERENCES segment(segment_id, creator_db),

      CONSTRAINT seglfnmap_fk_sdid
      FOREIGN KEY (lfn_id)
          REFERENCES lfn(lfn_id)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on process ID
CREATE INDEX seglfnmap_pid on segment_lfn_map(process_id)
;
-- Create an index based on segment ID
CREATE INDEX seglfnmap_sid on segment_lfn_map(segment_id)
;
-- Create an index based on segment lfninintion ID
CREATE INDEX seglfnmap_lfnid on segment_lfn_map(lfn_id)
;
