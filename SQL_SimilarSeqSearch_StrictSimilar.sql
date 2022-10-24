WITH target_taz AS ( -- 1. Get the query series
	SELECT *
	FROM {dataset} trajs
	WHERE trajs.tazid = {target_tazid}
     ), 
     
     -- Filter1: Distance between centroids 
     cp_dist AS ( -- Filter1.1: Distance filter with spatial index. Operator "<->" and keyword "LATERAL" can do (rough) NN-query with spatial index
	SELECT matchtrajs.tazid
	FROM (  SELECT * 
		FROM cp_trajs
		WHERE tazid = {target_tazid}) target, 
	LATERAL (SELECT *
		 FROM cp_trajs cp_trajs
		 ORDER BY target.geom <-> cp_trajs.geom ASC 
		 LIMIT {k*5}) as matchtrajs
     ), 
     table_cp_dist AS ( -- Filter1.2: Accurate distance between centroids
	SELECT trajs.*
	FROM cp_dist 
	INNER JOIN trajs trajs 
	ON cp_dist.tazid = trajs.tazid 
     ), 
     
     table_area AS ( -- Filter2: Areal difference
	SELECT trajs.tazid, trajs.geom
	FROM table_cp_dist trajs, target_taz target
	ORDER BY ABS(st_area(st_envelope(target.geom)) - st_area(st_envelope(trajs.geom))) ASC 
	LIMIT {k*2}
     ), 
     
     table_f_dist AS ( -- Filter3: Frechet Distance
	SELECT trajs.tazid, trajs.geom
	FROM table_area trajs, target_taz target
	ORDER BY st_frechetdistance(target.geom, trajs.geom) ASC 
	LIMIT {k}
     )

SELECT tazid
FROM table_f_dist;
