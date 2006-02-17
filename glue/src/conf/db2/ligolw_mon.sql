CREATE TABLE ligolw_mon
(
-- Table used to communicate with the DMT LIGOLwMon program.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH GENERATED THIS DATA
-- Process which generated this trigger
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- TIME OF THE TRIGGER
-- The time at which to display the trigger (GPS seconds and nanoseconds)
      time               INTEGER NOT NULL,
      time_ns            INTEGER NOT NULL,

-- PROPERTIES OF THE TRIGGER
-- amplitude is mandatory
      amplitude          REAL NOT NULL,
-- optional extra information about the event
      confidence         REAL,
      frequency          REAL,

-- Unique identifier for this event
      event_id           CHAR(13) FOR BIT DATA NOT NULL,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT ligolw_mon_pk
      PRIMARY KEY (event_id, creator_db),

      CONSTRAINT ligolw_mon_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE,
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time
CREATE INDEX ligolw_mon_ind_t ON ligolw_mon(time)
;
-- Create an index based on time and time_ns
CREATE INDEX ligolw_mon_ind_tns ON ligolw_mon(time, time_ns)
;
