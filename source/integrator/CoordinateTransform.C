//
// 文件名:     CoordinateTransform.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:44:51 2014 $
// 描述　:     坐标变换基类
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "CoordinateTransform.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
template <int DIM>
CoordinateTransform<DIM>::CoordinateTransform() {}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
template <int DIM>
CoordinateTransform<DIM>::~CoordinateTransform() {}
template class CoordinateTransform<NDIM>;
