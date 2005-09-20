CREATE TABLE runlist
(
-- This table contains the list of data runs.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Site (H for Hanford, L for Livingston, etc.)
-- If more than one DAQS system operates at a site, use different codes.
      site               CHAR(4) NOT NULL,

-- Run number
      run                INTEGER NOT NULL,

-- Time range (GPS seconds)
      start_time         INTEGER NOT NULL,
      end_time           INTEGER,

-- Optional additional information about this run
      run_type           VARCHAR(64),
      comment            VARCHAR(1024),
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT runlist_pk
      PRIMARY KEY (site, run)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on site and time
CREATE INDEX runlist_ind_time ON runlist(site, start_time)
;
