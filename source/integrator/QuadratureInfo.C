//
// 文件名:     QuadratureInfo.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:44:19 2014 $
// 描述　:     有限元积分信息类
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "DoubleVector.h"
#include "QuadratureInfo.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
template <int DIM>
QuadratureInfo<DIM>::QuadratureInfo(int algebric_accuracy) {}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
template <int DIM>
QuadratureInfo<DIM>::~QuadratureInfo() {}
template class QuadratureInfo<NDIM>;
