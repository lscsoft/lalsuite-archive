CREATE TABLE multi_burst
(
-- Event table for multi-interferometer burst-event search.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH GENERATED THIS EVENT
-- Process which generated this event
      process_id         CHAR(13) FOR BIT DATA NOT NULL,
-- Filter identifier (indicates type of filter, plus parameters).  May be null
      filter_id          CHAR(13) FOR BIT DATA,
-- Interferometers used for this search
      ifos               CHAR(12) NOT NULL,

-- TIME OF THE EVENT
-- The start time of this burst event (in GPS seconds and nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
-- The time duration of this burst event (seconds)
      duration           REAL NOT NULL,

-- PROPERTIES OF THE EVENT
-- Center of frequency band in which observation is made (Hz)
      central_freq       REAL,
-- Range of frequency observed (Hz)
      bandwidth          REAL,
-- Absolute signal amplitude (fractional strain)
      amplitude          REAL NOT NULL,
-- Signal to noise ratio
      snr                REAL,
-- Confidence variable
      confidence         REAL,

-- Direction of Hanford-to-Livingston ray at time of event
-- (i.e. the central axis of the cone on which the source lies)
      ligo_axis_ra       REAL,
      ligo_axis_dec      REAL,
-- Wave arrival angle with respect to Hanford-to-Livingston ray, and error
      ligo_angle         REAL,
      ligo_angle_sig     REAL,

-- Unique identifier for this event
      event_id           CHAR(13) FOR BIT DATA NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT m_burst_pk
      PRIMARY KEY (event_id, creator_db),

      CONSTRAINT m_burst_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db),

-- Note that filter_id is allowed to be null, in which case no check is made.
      CONSTRAINT m_burst_fk_filt
      FOREIGN KEY (filter_id, creator_db)
          REFERENCES filter(filter_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time alone
CREATE INDEX m_burst_ind_tim ON multi_burst(start_time, snr)
;
-- Create an index based on process_id
CREATE INDEX m_burst_ind_pid ON multi_burst(process_id, start_time)
;
-- Create an index based on filter_id
CREATE INDEX m_burst_ind_fid ON multi_burst(filter_id, start_time)
;
