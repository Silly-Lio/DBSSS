WITH traj_cp AS (
	    SELECT pt.*, cpt.geom AS cp
	    FROM {dataset} pt
	    INNER JOIN cp_trajs cpt
	    ON pt.tazid = cpt.tazid
     ), 
     taz_std AS (
	    SELECT tazid, st_scale(st_translate(geom, st_x(cp), st_y(cp)), std_x, std_y) AS geom
	    FROM traj_cp
     ), 
     target_std AS ( -- 0. 表格生成与存储
	    SELECT *
	    FROM taz_std pt
	    WHERE pt.tazid = {target_tazid}
     ), 
     
     table_area AS (
	    SELECT taz_std.tazid, taz_std.geom, ABS(st_area(st_envelope(taz_std.geom))- st_area(st_envelope(target_std.geom))) AS diff_area
	    FROM taz_std
	    CROSS JOIN target_std
	    ORDER BY diff_area ASC
	    LIMIT {k*2}
     ),
     
     table_f_dist AS (
	    SELECT table_area.tazid, table_area.geom
	    FROM table_area
	    CROSS JOIN target_std
	    ORDER BY st_frechetdistance(table_area.geom, target_std.geom) ASC
	    LIMIT {k} 
     )

    SELECT tazid
    FROM table_f_dist; 
