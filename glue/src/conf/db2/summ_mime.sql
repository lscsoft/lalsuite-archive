CREATE TABLE summ_mime
(
-- Table to store arbitary binary data, identified by a MIME type,
-- relevant to a particular time interval

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- ORIGIN OF THE DATA
-- Generator of the data (program name, 'video camera', etc.)
      origin             VARCHAR(64),
-- Unique process ID.  May be null, if entry was added manually
      process_id         CHAR(13) FOR BIT DATA,
-- Original filename (optional)
      filename           VARCHAR(64),
-- Name of person who entered data into database
      submitter          VARCHAR(48),
-- Timestamp at submission (automatically set by DB2)
      submit_time        TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

-- TIME INTERVAL TO WHICH THIS DATA APPLIES
-- Group name for frameset which determined this time interval, if any
      frameset_group     VARCHAR(48),
-- Group and version of segment which determined this time interval, if any
      segment_def_id     CHAR(13) FOR BIT DATA,
-- Start and end times (in GPS seconds and nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,

-- CHANNEL (OR PSEUDO-CHANNEL) NAME (if appropriate)
-- The channel/pseudo-channel name should indicate the interferometer or site
      channel            VARCHAR(240),
-- Brief description of the contents, e.g. "power spectrum"
      descrip            VARCHAR(64),

-- BINARY DATA AND ASSOCIATED INFO
-- The data itself is stored in a Binary Large OBject (BLOB).
-- We specify COMPACT since we do not expect this ever to be updated.
      mimedata           BLOB(1M) COMPACT NOT NULL,
-- Length of the binary data (in bytes)
      mimedata_length    INTEGER NOT NULL,
-- MIME content-type (e.g. 'image/jpeg', 'application/ligo-spectrum')
      mimetype           VARCHAR(64) NOT NULL,
-- Optional comment about this data
      comment            VARCHAR(240),
-- Unique identifier for this data (needed for primary key)
      summ_mime_id       CHAR(13) FOR BIT DATA NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,
      
      CONSTRAINT summmime_pk
      PRIMARY KEY (summ_mime_id, creator_db),

-- Note that process_id is allowed to be null, in which case no check is made.
      CONSTRAINT summmime_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db),

-- If segment_group or frameset_group is non-null, make sure there is a
-- corresponding entry in the appropriate table.  If null, then no
-- foreign-key check is performed.
      CONSTRAINT summmime_fk_seg
      FOREIGN KEY (segment_def_id,creator_db)
          REFERENCES segment_definer(segment_def_id,creator_db),

      CONSTRAINT summmime_fk_fs
      FOREIGN KEY (frameset_group, start_time, end_time)
          REFERENCES frameset(frameset_group, start_time, end_time)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time
CREATE INDEX summmime_ind_time ON summ_mime(start_time, end_time)
;
-- Create an index based on mimetype
CREATE INDEX summmime_ind_type ON summ_mime(mimetype)
;
-- Create an index based on description
CREATE INDEX summmime_ind_desc ON summ_mime(descrip, start_time)
;
-- Create an index based on process_id
CREATE INDEX summmime_ind_pid ON summ_mime(process_id)
;
-- Create an index based on frameset_group
CREATE INDEX summmime_ind_fsg ON summ_mime(frameset_group)
;
-- Create an index based on segment_group
CREATE INDEX summmime_ind_sgrp ON summ_mime(segment_def_id)
;
