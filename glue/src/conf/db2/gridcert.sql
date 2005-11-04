CREATE TABLE gridcert
(
-- An entry in the process table can be associated with a grid certificate to
-- to track the distinguished name (DN) of the person or program who created
-- the entry.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Unique process ID of the process which wrote this lfn
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Distinguished name of the user or program
      dn                 VARCHAR(255) NOT NULL,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT gridcert_pk
      PRIMARY KEY (creator_db,process_id),

-- This entry must map to an entry in the process table
      CONSTRAINT gridcert_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on distinguished name
CREATE INDEX gridcert_ind_dn ON gridcert(dn)
;
