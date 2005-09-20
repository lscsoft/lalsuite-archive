CREATE TRIGGER gdstrig_tr_del \
    AFTER DELETE ON gds_trigger \
    REFERENCING OLD_TABLE AS oldtable \
    FOR EACH STATEMENT MODE DB2SQL \
    BEGIN ATOMIC \
        DELETE FROM sngl_datasource \
           WHERE (event_id, creator_db) in \
           ( select event_id, creator_db from oldtable ) ; \
        DELETE FROM sngl_transdata \
           WHERE (event_id, creator_db) in \
           ( select event_id, creator_db from oldtable ) ; \
        DELETE FROM sngl_mime \
           WHERE (event_id, creator_db) in \
           ( select event_id, creator_db from oldtable ) ; \
    END
