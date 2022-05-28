-- 757 北京西站
-- 60 527 525 526 三里屯
-- 84 天安门
-- 1565 1566 永定河水库
-- 657 660 望京
-- 669 朝阳公园
-- 924 天通苑

-- 自定义参数
-- {dataset} 被查找的集合名称
-- {target_tazid} 查询序列的ID
-- {k} 查找数量


WITH target_taz AS ( -- 1. 表格生成与存储
		SELECT *
		FROM {dataset} trajs
		WHERE trajs.tazid = {target_tazid}
	), 
	-- Filter1 中心点距离约束
	bc_dist AS ( -- Filter1.1 使用索引预先进行判断：对应外包矩形的中心点距离约束
		SELECT matchtrajs.tazid, matchtrajs.geom
		FROM target_taz target, 
		LATERAL (SELECT *, target.geom <-> trajs.geom as dist
					FROM {dataset} trajs
					ORDER BY dist ASC 
					LIMIT {k*10}) as matchtrajs
	), table_cp_dist AS ( -- Filter1.2 真实轨迹的中心点距离约束
		SELECT trajs.tazid, trajs.geom
		FROM target_taz target, bc_dist trajs
		ORDER BY st_centroid(target.geom) <-> st_centroid(trajs.geom) ASC 
		LIMIT {k*5}
	),

	table_area AS ( -- Filter2 外包矩形的面积约束
		SELECT trajs.tazid, trajs.geom
		FROM table_cp_dist trajs, target_taz target
		ORDER BY ABS(st_area(st_envelope(target.geom)) - st_area(st_envelope(trajs.geom))) ASC 
		LIMIT {k*2}
	), 

	table_f_dist AS ( -- Filter3 Frechet 距离约束
		SELECT trajs.tazid, trajs.geom
		FROM table_area trajs, target_taz target
		ORDER BY st_frechetdistance(target.geom, trajs.geom) ASC 
		LIMIT {k}
	)

SELECT tazid
FROM table_f_dist;