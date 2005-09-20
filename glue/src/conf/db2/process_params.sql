CREATE TABLE process_params
(
-- This table contains input parameters for programs.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Program name
      program            CHAR(16) NOT NULL,
-- Unique process ID (not unix process ID)
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- The triplet will store the value of a single parameter.
-- One example might be param = "mass", type="REAL"
-- and value = "12345.6789"
      param              VARCHAR(32) NOT NULL,
      type               VARCHAR(16) NOT NULL,
      value              VARCHAR(1024) NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT procpar_pk
      PRIMARY KEY (process_id, creator_db, param),

-- Foreign key relationship to process table.  The 'ON DELETE CASCADE'
-- modifier means that if a row in the process table is deleted, then
-- all its associated parameters are deleted too.
      CONSTRAINT procpar_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on program name
CREATE INDEX procpar_ind_name ON process_params(program)
;
