CREATE TRIGGER s_block_tr_del \
    AFTER DELETE ON sngl_block \
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
