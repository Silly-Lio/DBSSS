-- Parameters
-- {dataset} - target dataset
-- {target_tazid} - of the query series
-- {k} - the retrieval number K

WITH traj_cp AS ( -- 1. Add columns of TABLE <traj_cp> to the target dataset. TABLE <traj_cp> contains the mean value of time series in X and Y axis, and the standard deviation of time series in X and Y axis
	    SELECT pt.*, cpt.geom AS cp
	    FROM {dataset} pt
	    INNER JOIN cp_trajs cpt
	    ON pt.tazid = cpt.tazid
     ), 
     taz_std AS ( -- 2. Translate and rescale the time series in original dataset, which is equivalent to the z-score normalization
	    SELECT tazid, st_scale(st_translate(geom, st_x(cp), st_y(cp)), std_x, std_y) AS geom
	    FROM traj_cp
     ), 
     target_std AS ( -- 3. Get the query series
	    SELECT *
	    FROM taz_std pt
	    WHERE pt.tazid = {target_tazid}
     ), 
     
     table_area AS ( -- Filter1: Areal difference
	    SELECT taz_std.tazid, taz_std.geom, ABS(st_area(st_envelope(taz_std.geom))- st_area(st_envelope(target_std.geom))) AS diff_area
	    FROM taz_std
	    CROSS JOIN target_std
	    ORDER BY diff_area ASC
	    LIMIT {k*2}
     ),
     
     table_f_dist AS ( -- Filter2: Frechet Distance
	    SELECT table_area.tazid, table_area.geom
	    FROM table_area
	    CROSS JOIN target_std
	    ORDER BY st_frechetdistance(table_area.geom, target_std.geom) ASC
	    LIMIT {k} 
     )

SELECT tazid
FROM table_f_dist; 
