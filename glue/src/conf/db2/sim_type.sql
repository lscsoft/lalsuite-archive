CREATE TABLE sim_type
(
-- Classification of simulations. This table is not changed by DSO but by hand.
-- This table should migrate from one database to another together with data.

-- Simulations types are enumerated here automatically by DB2
      sim_type          INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY ( START WITH 0 , INCREMENT BY 1 , NO CACHE ), 
-- Simulation name
      name               VARCHAR(24) NOT NULL, 
-- Simulation description.
      comment            VARCHAR(64),
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,
      
      CONSTRAINT s_type_pk
      PRIMARY KEY (sim_type)
)   

-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on name
CREATE INDEX s_type_ind_name ON sim_type(name)
;

