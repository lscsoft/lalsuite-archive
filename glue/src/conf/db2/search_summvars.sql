CREATE TABLE search_summvars
(
-- This table contains search-specific summary variables in the form of
-- name/value pairs.  Any given search can create an arbitrary number of
-- entries in this table.
-- Created by Peter, 21 Feb 2002

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH RAN THIS SEARCH
-- Process which generated this event
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Name/value pairs.  The value can be either a string or a number (expressed
-- as a double-precision real number, even if the value is an integer).
-- To do this, two columns are provided; just fill the appropriate one and
-- leave the other one blank.
      name               VARCHAR(64) NOT NULL,
      string             VARCHAR(256),
      value              DOUBLE,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

-- Require this to correspond to an entry in the search_summary table
      CONSTRAINT s_summvar_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES search_summary(process_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on name
CREATE INDEX s_summvar_ind_name ON search_summvars(name)
;
