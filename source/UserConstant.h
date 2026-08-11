//
// 文件名:     UserConstant.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2017-12-25
// 描述　:     自定义宏及数据常量
// 类别　:     %Internal File% ( Don't delete this line )
//
#ifndef _USERCONSTANT_H_
#define _USERCONSTANT_H_
#define AREA(a, b, c)                      \
  ((((b)[1] - (a)[1]) * ((c)[2] - (a)[2]) -((b)[2] - (a)[2]) * ((c)[1] - (a)[1]))*(((b)[1] - (a)[1]) * ((c)[2] - (a)[2]) -((b)[2] - (a)[2]) * ((c)[1] - (a)[1]))+\
  (((b)[2] - (a)[2]) * ((c)[0] - (a)[0]) -((b)[0] - (a)[0]) * ((c)[2] - (a)[2]))*(((b)[2] - (a)[2]) * ((c)[0] - (a)[0]) -((b)[0] - (a)[0]) * ((c)[2] - (a)[2]))+\
  (((b)[0] - (a)[0]) * ((c)[1] - (a)[1]) -((b)[1] - (a)[1]) * ((c)[0] - (a)[0]))*(((b)[0] - (a)[0]) * ((c)[1] - (a)[1]) -((b)[1] - (a)[1]) * ((c)[0] - (a)[0])))
/*************************************************************************
  *                  常量定义
  ************************************************************************/
const double hop_distantce=0.1e-9;//跃迁距离
const double escape_freq=1e12;//逃逸频率
const double Diff0=2e-7;//扩散系数（m2/s）
const double K_b=1.380650400000000E-23;//波尔兹曼常数
const double EA=1.0*1.602176620800000E-19;//激活能
const double Echarge=1.602176620800000E-19;//电子电荷量
const double KHfO0=0.5;//(w/(m*K)) 氧化铪的热导率
const double KHf = 23;//(w/(m*K)) 铪热导率
const double nd_int[2]={1.0e20,12e26};//初始掺杂浓度
const double eps1=1e-5;//有无漂移流的判据，当两个结点的电势差小于该值的时候，判断无漂移流
const double Kte=71.6;//上电极的热导率 pt
const double Sigmate=8.9e6;//上电极的电导率 pt
const double Kbe=29.31;//下电极的热导率 w
const double Sigmabe=5e6;//下电极的导率
const double Kdiodeon=22;//二极管热导率 sio2
const double Sdiodeon=30700.;//二极管电导率
const double Kdiodeoff=22;//二极管热导率 fengzhuang
const double Sdiodeoff=0.005;//二极管电导率
const double Kme=300;//二极管热导率 tin
const double Sigmame=1e-5;//二极管电导率
const double EPS= 1e-13;///(2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2*2);//很小的正数




/*************************************************************************
  *                  线性函数
  ************************************************************************/
//导体激活能线性函数
double Eac_var(double nD);
//电导率线性函数
double Sigma_0(double nD);
//热导率线性函数
double Kth_var(double nD);
/*************************************************************************
  *                  电压脉冲函数
  ************************************************************************/
//过原点的直线
double voltagePulse(double time);
#endif
