CREATE TABLE lfn
(
-- An entry in the process table may have one or more logial file names (LFNs)
-- associated with it. These LFNs can be used in a query to a replica location
-- service (RLS) to obtain the physical file names (PFNs) associated with the
-- entry in the process table. A process to LFN mapping is typically used to 
-- allow storage of information that is too large to be kept in the database
-- (e.g. raw trigger files)

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Unique process ID of the process which wrote this lfn
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- logical file name
      lfn                VARCHAR(255) NOT NULL,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

-- LFNs in this table must be unique
      CONSTRAINT lfn_pk
      PRIMARY KEY (lfn)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on process_id
CREATE INDEX lfn_ind_pid ON lfn(process_id)
;
