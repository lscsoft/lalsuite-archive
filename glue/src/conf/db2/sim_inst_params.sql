CREATE TABLE sim_inst_params
(
-- Parameters of the simulation instance
	
-- Simulation instance id
      simulation_id               CHAR(13) FOR BIT DATA NOT NULL, 
-- LIGO/LSC site that created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1, 
-- Parameter name
      name               VARCHAR(24) NOT NULL, 
-- Parameter descrition
      comment            VARCHAR(64),
-- Parameter value
      value              DOUBLE,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT sim_in_par_fk_si
      FOREIGN KEY (simulation_id, creator_db)
          REFERENCES sim_inst(simulation_id, creator_db)
          ON DELETE CASCADE
)   
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;


