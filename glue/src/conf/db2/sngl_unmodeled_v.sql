CREATE TABLE sngl_unmodeled_v
(
-- Generic table to store additional values describing events found by a
-- search for "unmodeled" sources.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Unique identifier for the source-independent event
      event_id           CHAR(13) FOR BIT DATA NOT NULL,
-- Site or interferometer to which this applies (H0, H1, H2, L0, L1)
      ifo                CHAR(2) NOT NULL,

-- Descriptive name of the result variable
      name               VARCHAR(32) NOT NULL,
-- The value of the result (must be a real number)
      value              REAL NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT s_unmodv_pk
      PRIMARY KEY (event_id, creator_db, name),

      CONSTRAINT s_unmodv_fk_unmod
      FOREIGN KEY (event_id, creator_db)
          REFERENCES sngl_unmodeled(event_id, creator_db)
          ON DELETE CASCADE,

-- Foreign key relationship to process table.
      CONSTRAINT s_unmodv_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on name
CREATE INDEX s_unmodv_ind_name ON sngl_unmodeled_v(name, event_id)
;
-- Create an index based on process_id
CREATE INDEX s_unmodv_ind_pid ON sngl_unmodeled_v(process_id)
;
