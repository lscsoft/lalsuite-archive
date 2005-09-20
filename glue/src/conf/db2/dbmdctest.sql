CREATE TABLE dbmdctest
(
  -- This table is used by Peter's database MDC scripts.  It should not be
  -- created automatically, but should be created manually when needed.

  seqid      INTEGER NOT NULL,
  tshort     SMALLINT,
  tint       INTEGER,
  tlong      BIGINT,
  treal      REAL,
  tdouble    DOUBLE,
  tchar      CHAR(20),
  tvarchar   VARCHAR(20),
  tcharbin   CHAR(13) FOR BIT DATA,
  tblob      BLOB(32),
  tbigblob   BLOB(1048576),

      CONSTRAINT dbmdctest_pk
      PRIMARY KEY (seqid)
);
