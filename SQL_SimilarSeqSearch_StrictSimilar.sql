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


WITH target_taz AS ( -- 1. 获取查询序列的对应轨迹
		SELECT *
		FROM {dataset} pt
		WHERE pt.tazid = {target_tazid}
	), 

	traj_cp AS ( -- 2. 获取查询轨迹中心点
		SELECT *, st_centroid(geom) AS cp
		FROM {dataset} pt
	), 
	taz_std AS ( -- 3. 对查询轨迹进行平移与缩放，相当于对原始的时间序列进行z-score标准化
		SELECT tazid, st_scale(st_translate(geom, st_x(cp), st_y(cp)), std_x, std_y) AS geom
		FROM traj_cp
	), 
	target_cp AS ( -- 4. 获取全数据集所包含轨迹的中心点
		SELECT *, st_centroid(geom) AS cp
		FROM target_taz
	),
	target_std AS ( -- 5. 对全数据集的轨迹进行平移与缩放，相当于对原始的时间序列进行z-score标准化
		SELECT tazid, st_scale(st_translate(geom, st_x(cp), st_y(cp)), std_x, std_y) AS geom
		FROM target_cp
	),

	table_area AS ( -- Filter1. 外包矩形的面积约束
		SELECT taz_std.tazid, taz_std.geom, ABS(st_area(st_envelope(taz_std.geom))- st_area(st_envelope(target_std.geom))) AS diff_area
		FROM taz_std
		CROSS JOIN target_std
		ORDER BY diff_area ASC
		LIMIT {k*2}
	), 

	table_f_dist AS ( -- Filter2. Frechet 距离约束
		SELECT table_area.tazid, table_area.geom
		FROM table_area
		CROSS JOIN target_std
		ORDER BY st_frechetdistance(table_area.geom, target_std.geom) ASC
		LIMIT {k} 
	)

SELECT tazid
FROM table_f_dist; 
