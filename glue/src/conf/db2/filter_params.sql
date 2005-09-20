CREATE TABLE filter_params
(
-- This table contains parameters for filters.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Filter name
      filter_name        CHAR(32) NOT NULL,
-- Unique process ID (not unix process ID)
      filter_id          CHAR(13) FOR BIT DATA NOT NULL,

-- The triplet will store a single parameter.
-- One example might be name = "pole1", type="REAL"
-- and value = "12345.6789"
      param              VARCHAR(32) NOT NULL,
      type               VARCHAR(16) NOT NULL,
      value              VARCHAR(1024) NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT filtpar_pk
      PRIMARY KEY (filter_id, creator_db, param),

-- Foreign key relationship to filter table.  The 'ON DELETE CASCADE'
-- modifier means that if a row in the filter table is deleted, then
-- all its associated parameters are deleted too.
      CONSTRAINT filtpar_fk_fid
      FOREIGN KEY (filter_id, creator_db)
          REFERENCES filter(filter_id, creator_db)
          ON DELETE CASCADE,

-- Foreign key relationship to process table.
      CONSTRAINT filtpar_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on filter name
CREATE INDEX filtpar_ind_name ON filter_params(filter_name)
;
-- Create an index based on process_id
CREATE INDEX filtpar_ind_pid ON filter_params(process_id)
;
