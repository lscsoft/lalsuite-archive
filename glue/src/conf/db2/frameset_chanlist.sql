CREATE TABLE frameset_chanlist
(
-- List of channels included in a frameset.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Unique process ID of the process which submits the chanlist info
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Frameset group for which this channel list applies
-- (But note that one frameset group could have several different channel
-- lists for data taken at different times)
      frameset_group     CHAR(48) NOT NULL,
-- Validity range (in GPS seconds) for this channel list
      start_time         INTEGER,
      end_time           INTEGER,

-- Unique identifier for this channel list
      chanlist_id        CHAR(13) FOR BIT DATA NOT NULL,

-- Channel list, with modifiers and sampling rates.  Separated by spaces.
-- Examples of items in the list:
--      'H2:PEM-LVEA_SEISX 256'                 Raw data stream
--      'H2:PEM-LVEA_SEISX 16'                  Decimated
--      'H2:PEM-LVEA_SEISX.LOWPASS(10.0) 16'    Filtered and decimated
-- List is stored in a Character Large OBject (CLOB).
      chanlist           CLOB(512K),
-- Length of channel list in bytes (0 if the CLOB is empty)
      chanlist_length    INTEGER NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT fschanlist_pk
      PRIMARY KEY (creator_db, chanlist_id),

      CONSTRAINT fschanlist_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- The following line ensures that replication will work properly on any
-- LONG VARCHAR columns that we might add to this table in the future.
ALTER TABLE frameset_chanlist
    DATA CAPTURE CHANGES INCLUDE LONGVAR COLUMNS
;
