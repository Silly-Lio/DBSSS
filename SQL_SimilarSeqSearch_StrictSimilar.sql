WITH target_taz AS ( -- 0. 表格生成与存储
	SELECT *
	FROM {dataset} trajs
	WHERE trajs.tazid = {target_tazid}
     ), 
     cp_dist AS (
	SELECT matchtrajs.tazid
	FROM (  SELECT * 
		FROM cp_trajs
		WHERE tazid = {target_tazid}) target, 
	LATERAL (SELECT *
		 FROM cp_trajs cp_trajs
		 ORDER BY target.geom <-> cp_trajs.geom ASC 
		 LIMIT {k*5}) as matchtrajs
     ), table_cp_dist AS (
	SELECT trajs.*
	FROM cp_dist 
	INNER JOIN trajs trajs 
	ON cp_dist.tazid = trajs.tazid -- ORDER BY st_centroid(target.geom) <-> st_centroid(trajs.geom) ASC 
     ), table_area AS ( -- 2. 多边形面积
	SELECT trajs.tazid, trajs.geom
	FROM table_cp_dist trajs, target_taz target
	ORDER BY ABS(st_area(st_envelope(target.geom)) - st_area(st_envelope(trajs.geom))) ASC 
	LIMIT {k*2}
     ), 
	table_f_dist AS ( -- 3. Frechet 距离
	SELECT trajs.tazid, trajs.geom
	FROM table_area trajs, target_taz target
	ORDER BY st_frechetdistance(target.geom, trajs.geom) ASC 
	LIMIT {k}
     )

SELECT tazid
FROM table_f_dist;
