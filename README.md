# Heat-exchanger
Simulation of a CHX using Matlab </br>
划分换热器微元，分别使用焓值（Untitled）、温度（test）迭代计算换热器出口温度</br>
未考虑压力变化
最终计算结果与EES进行对比，CFHX计算较为准确

2020/05/15更新
考虑辐射与压降，test2采用经验关联式进行逆流式套管换热器仿真的弯管修正。
对于套管式换热器进行适当简化——通道内流体一维流动，壁面一维导热。
