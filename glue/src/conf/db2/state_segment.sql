create view state_segment (process_id,segment_id,
        start_time,start_time_ns,end_time,end_time_ns,segment_def_id,lfn_id) as
    select segment.process_id,segment.segment_id,
        segment.start_time,segment.start_time_ns,
        segment.end_time,segment.end_time_ns,
        segment_definer.segment_def_id,
        segment_lfn_map.lfn_id from
    segment,segment_def_map,segment_definer,segment_lfn_map where
        segment.segment_id = segment_def_map.segment_id 
        and segment_def_map.segment_def_id = segment_definer.segment_def_id 
        and segment.segment_id = segment_lfn_map.segment_id
        and segment_definer.state_vec_major is not null 
        and segment_definer.state_vec_minor is not null
@
create trigger state_segment_i instead of insert on state_segment
    referencing new as n 
    for each row mode db2sql
        begin atomic
        values ( 
        case when ( (select count(s.segment_id) from 
            segment as s, segment_def_map as m where
            m.segment_def_id = n.segment_def_id 
            and s.start_time = n.start_time 
            and s.start_time_ns = n.start_time_ns 
            and s.end_time = n.end_time 
            and s.end_time_ns = n.end_time_ns 
            and s.segment_id = m.segment_id) > 0 ) 
        then raise_error ( '70001', 'state_segment rows must be unique' ) 
        else 0 end );
        insert into segment (process_id,segment_id,active,
            start_time,start_time_ns,end_time,end_time_ns)
        values (n.process_id,n.segment_id,1,
            n.start_time,n.start_time_ns,n.end_time,n.end_time_ns);
        insert into segment_def_map (process_id,segment_id,segment_def_id)
        values (n.process_id,n.segment_id,n.segment_def_id);
        insert into segment_lfn_map (process_id,segment_id,lfn_id)
        values (n.process_id,n.segment_id,n.lfn_id);
    end
@
