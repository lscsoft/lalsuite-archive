CREATE TABLE sim_type_params
(
-- Parameters for the simulation types from sim_type table

-- Foreign key refering to sim_type.type
      sim_type          INTEGER NOT NULL, 
-- Parameter name
      name               VARCHAR(24) NOT NULL, 
-- Parameter description
      comment            VARCHAR(64), 
-- Parameter value
      value              DOUBLE,

      CONSTRAINT sim_type_par_pk
      PRIMARY KEY (sim_type, name),

      CONSTRAINT sim_type_par_fk_st 
      FOREIGN KEY (sim_type)
          REFERENCES sim_type(sim_type)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;

