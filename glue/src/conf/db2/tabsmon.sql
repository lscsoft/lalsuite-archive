set event monitor tabsmon state = 0;
drop event monitor tabsmon;
create event monitor tabsmon 
for tables, statements 
write to FILE '/tmp/db2' 
replace;
set event monitor tabsmon state = 1 ;
