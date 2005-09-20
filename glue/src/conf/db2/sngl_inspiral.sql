CREATE TABLE sngl_inspiral
(
-- Event table for single-interferometer binary-inspiral search.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH GENERATED THIS EVENT
-- Process which generated this event
      process_id         CHAR(13) FOR BIT DATA NOT NULL,
-- Filter identifier (indicates type of filter, plus parameters).  May be null
      filter_id          CHAR(13) FOR BIT DATA,
-- Interferometer
      ifo                CHAR(2) NOT NULL,
-- Brief keyword to identify the search technique, e.g. "template" or "FCT"
      search             VARCHAR(24) NOT NULL,
-- Channel that was analyzed
      channel            VARCHAR(64),

-- TIME OF THE EVENT
-- The coalescence time of this inspiral event (GPS seconds and nanoseconds)
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,
      end_time_gmst      DOUBLE,
-- Time of impulse that would generate a peak filter output at this time
-- (end_time-impulse_time is a fixed quantity for a given template)
      impulse_time       INTEGER,
      impulse_time_ns    INTEGER,
-- Duration of the template (in seconds) used to identify this event
      template_duration  DOUBLE,
-- Duration of the event (amount of time over threshold, in seconds)
      event_duration     DOUBLE,

-- PROPERTIES OF THE EVENT
-- Absolute strain amplitude of fitted template at 100 Hz
      amplitude          REAL NOT NULL,
-- Effective distance, assuming optimally oriented binary (in megaparsecs)
      eff_distance       REAL,
-- Coalescence phase angle (radians)
      coa_phase          REAL,

-- Mass of the larger object (in solar mass units)
      mass1              REAL,
-- Mass of the smaller (or equal) object
      mass2              REAL,
-- Chirp mass [ (m1+m2)*eta^(3/5) ] and eta [ m1*m2/(m1+m2)^2 ]
      mchirp             REAL,
      mtotal             REAL,
      eta                REAL,
-- Chirp-time parameters
      tau0               REAL,
      tau2               REAL,
      tau3               REAL,
      tau4               REAL,
      tau5               REAL,
-- ttotal is defined as tau0+tau2-tau3+tau4-tau5
      ttotal             REAL,

-- BCV detection template family waveform parameters
      psi0               REAL,
      psi3               REAL,
      f_final            REAL,
      
-- Amplitude signal-to-noise ratio
      snr                REAL,
-- Chi-squared statistic, and number of degrees of freedom
      chisq              REAL,
      chisq_dof          INTEGER,
-- Variance of filter output (nominally normalized snr) near time of event
      sigmasq            DOUBLE,
-- Measured value of BCV shape parameter
      alpha              REAL,
      alpha1             REAL,
      alpha2             REAL,
      alpha3             REAL,
      alpha4             REAL,
      alpha5             REAL,
      alpha6             REAL,
      beta		 REAL,

-- Unique identifier for this event
      event_id           CHAR(13) FOR BIT DATA NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT s_inspiral_pk
      PRIMARY KEY (event_id, creator_db),

      CONSTRAINT s_inspiral_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db),

-- Note that filter_id is allowed to be null, in which case no check is made.
      CONSTRAINT s_inspiral_fk_filt
      FOREIGN KEY (filter_id, creator_db)
          REFERENCES filter(filter_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create a clustering index based on search algorithm
CREATE INDEX s_inspiral_cind ON sngl_inspiral(search) CLUSTER
;
-- Create an index based on search algorithm plus time and snr
CREATE INDEX s_inspiral_ind_sea ON sngl_inspiral(search, end_time, snr)
;
-- Create an index based on time and snr
CREATE INDEX s_inspiral_ind_tim ON sngl_inspiral(end_time, snr)
;
-- Create an index based on process_id
CREATE INDEX s_inspiral_ind_pid ON sngl_inspiral(process_id, end_time)
;
-- Create an index based on filter_id
CREATE INDEX s_inspiral_ind_fid ON sngl_inspiral(filter_id, end_time)
;
-- Create an SQL trigger so that if a sngl_inspiral entry is deleted, any
-- associated sngl_datasource, sngl_transdata, and/or sngl_mime entries
-- are deleted too.
-- Must be done this way because there is no foreign-key relationship.
-- Run script sngl_inspiral_tr_del.sql to create delete trigger for 
-- sngl_datasource, sngl_transdata, and sngl_mime records.
