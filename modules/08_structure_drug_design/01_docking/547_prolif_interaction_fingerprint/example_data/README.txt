本模块的 turnkey demo 直接使用 prolif 官方打包的真实蛋白-配体 MD 复合物
(plf.datafiles.TOP / plf.datafiles.TRAJ:一个 250 帧、含盐桥/氢键/疏水/π 相互作用
的真实复合物),无需在此放置示例数据,故本目录默认为空。

若要换自己的体系,把拓扑与轨迹放这里并用 --top/--traj 指向即可,例如:
  python 547_prolif_interaction_fingerprint.py --top example_data/my.pdb --traj example_data/my.xtc
