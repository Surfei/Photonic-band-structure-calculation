# photonic band structure calculation
## 光子晶体能带计算
### 程序介绍
- 使用Matlab程序  
- 主程序 `photonic_crystal.m` 计算本征频率并画出能带图  
- 调用外部函数 `ecrcepsilon.m` 计算介电常数在到空间中的展开系数  
### 使用方法：  
- 将`photonic_crystal.m`和`ecrcepsilon.m`两个文件放在同一文件夹中，用Matlab运行主程序
- 设置计算模式`mode=0`，可以得到二维正方晶格TE模式的能带图  
- 设置计算模式`mode=1`，可以得到二维正方晶格TM模式的能带图  
- 经测试，程序可在Matlab 2016a及Matlab 2017a版本上成功运行
### 运行结果
- TE模式，如图`TE modes.png`
- TM模式，如图`TM modes.png`
- 与文献 `Meade R D V , Johnson S G , Winn J N . Photonic Crystals: Molding the Flow of Light - Second Edition[M]. 2008.`相比较,如图`reference.png`，可以看出计算结果与其相符。
