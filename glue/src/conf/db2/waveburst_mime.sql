CREATE TABLE waveburst_mime  
(
-- Table with smallest rectangle of pixels containing the cluster

-- Process which generated this entry
      process_id         CHAR(13) FOR BIT DATA NOT NULL, 
-- LIGO/LSC site that created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1, 
-- Unique identifier for this event
      event_id           CHAR(13) FOR BIT DATA NOT NULL, 
-- Description of BLOB format if different from t_size by f_size array
-- of REAL4 entries
      mimetype            VARCHAR(64), 
-- Dimensions of the rectangle
-- Number of time steps
      t_size             INTEGER, 
-- Number of frequency steps
      f_size             INTEGER, 
-- Rectangle contains original amplitudes of each pixel
      cluster_o          BLOB(102400) NOT LOGGED COMPACT, 
-- Rectangle contains percentile amplitudes of each pixel
      cluster_p          BLOB(102400) NOT LOGGED COMPACT, 
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,
      
      CONSTRAINT wbm_fk_wb
      FOREIGN KEY (event_id, creator_db)
          REFERENCES waveburst(event_id,creator_db)
          ON DELETE CASCADE

)   
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;

