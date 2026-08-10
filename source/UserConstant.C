/*
 *  UserContant.C
 */

#include "UserConstant.h"
/*************************************************************************
  *                  线性函数
  ************************************************************************/
//导体激活能线性函数（nD表示浓度）
double Eac_var(double nD)
{
    double Eac=0;
    if(nD<0) {Eac=0.05;}
    if(nD>=0&&nD<2.0e26) { Eac=((nD/(2.0e26-0))*(0-0.05)+0.05);}
    if(nD>=2.0e26) {Eac=0;}
    return Eac*1.602176620800000E-19;//换算为焦耳
}
//电导率线性函数（nD表示浓度）
double Sigma_0(double nD)
{
    double sigma=1000;
    if(nD<0) {sigma=1000;}
    if(nD>=0&&nD<1.2e27) {sigma=((nD/(1.2e27-0))*(330000-1000)+1000);}
    if(nD>=1.2e27) {sigma=330000;}
    return sigma;
}
//热导率线性函数（nD表示浓度）
double Kth_var(double nD)
{
    double Kth=0;
    if(nD<0) {Kth=0.5;}
    if(nD>=0&&nD<1.2e27) { Kth=((nD/(1.2e27-0))*(23-0.5)+0.5);}
    if(nD>=1.2e27) {Kth=23;}
    return Kth;
}
/*************************************************************************
  *                  电压脉冲函数
  ************************************************************************/
//过原点的直线
double voltagePulse(double time){
    double voltage=0;
    voltage = 1*time;
    if(time<0.6) {voltage=1.0*time;}
    if(time>=0.6 && time<=1.0) { voltage = 1.0*time;}

    return voltage;
}
