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

--  Each entry in the process table can by associated with only one DN.
      CONSTRAINT gridcert_pk
      PRIMARY KEY (process_id,dn)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on process_id
CREATE INDEX gridcert_ind_pid ON gridcert(process_id)
;
-- Create an index based on distinguished name
CREATE INDEX gridcert_ind_dn ON gridcert(dn)
;
