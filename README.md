## Update (2021.06.27)

- [x] 采用更简洁方式改写convertKernel;
- [x] 重构接口，将Viewer与库隔离;
- [x] 添加“需定义`PolygonCenterLine_Implementation`宏”的提示;
- [x] 将`PolygonCenterlineSolver`中较长的函数定义放在模板类外;
- [ ] 继续Debug

## Progress

![output 0](README.assets/output 0.png)

![output 2](README.assets/output 2.png)



## 目前已知的bug

- [ ] 个别点的end_loc未被赋值；

- [ ] 个别边的两个端点end_loc被放在了同一位置

  ![output 1](README.assets/output 1.png)



------------------------------------------

# Old

- 上次的问题：没有考虑到可能在

在interval增加一个属性is_endpoint，1表示src沿着中线走，2表示target沿着中线走。
求edge与line的交点。
若端点location均小于等于origin的location，continue。
//（若线段运动方向背对直线，如果依然能相交，说明运动方向偏转超过中线，说明相邻边偏转超过对位边，所以相邻边不会是对位边的相邻边，那么显然这个相邻边更能代表这一点的碰撞，所以不用考虑？）（若两边在中线同侧，则正对中线的一方一定是靠近另一边的一侧先接触中线，而对应的相邻边一定更接近背对中线的边，且角度更接近）
若交点为空：

- 若inner小于0，取target的path；
- 否则取src的path，与bisecter求交点。
- 若交点不存在，则不会相交；
- 否则交点就是相交的起始点。

若交点非空：

- 则交点就是相交区间的起始点。

- 两个端点ray分别与bisector求交点。
  若交点均不存在，则延伸到无穷远；
  若一个存在另一个不存在：若存在的那个是point，则point为结束点；若是ray，则
  若交点落在起始点之前(<=)，则相交区间不存在；
  否则，交点就是相交区间的结束点。

初始化区间。

## TODO

- [ ] 支持Polygon_with_holes
- [ ] 