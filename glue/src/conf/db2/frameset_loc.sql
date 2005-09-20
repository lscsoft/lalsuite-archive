CREATE TABLE frameset_loc
(
-- Table to keep track of frameset locations.  There can be more than one
-- location for a given frameset.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Frameset name
      name               VARCHAR(80) NOT NULL,

-- Media type (disk, hpss, 8mm, dvd,...)
      media_type         CHAR(16) NOT NULL,
-- The node with the data (i.e. computer on which disk is mounted, or hpss
-- server); OR, in the case of removable media, this should be the media label
      node               VARCHAR(48) NOT NULL,
-- Physical location of removable media (e.g. LHO, LLO, transit, CACR, ...)
      media_loc          VARCHAR(48),
-- Status of the media (e.g. OK, broken, ...)
      media_status       CHAR(8) WITH DEFAULT 'OK',
-- Full path and actual file name.  OR, in the case of a tape without named
-- files, this should just be the frameset name again.
      fullname           VARCHAR(128) NOT NULL,
-- File number (for tapes without named files)
      filenum            INTEGER,
-- Decompression command (e.g. 'gunzip *.tar.gz; tar xf *.tar', where '*'
-- is replaced by the frameset name)
      decompress_cmd     VARCHAR(128),
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT fsloc_pk
      PRIMARY KEY (node, fullname),

      CONSTRAINT fsloc_fk_frameset
      FOREIGN KEY (name) 
	  	REFERENCES frameset(name)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Also create a clustering index for quicker name lookup
CREATE INDEX fsloc_ind_name ON frameset_loc(name) CLUSTER;
