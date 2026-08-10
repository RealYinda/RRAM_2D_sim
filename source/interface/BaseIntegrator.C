//
// 文件名:     BaseIntegrator.C
// 软件包:     JAUMIN fem
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-06-17 16:43:13 +0800 (日, 2012-06-17) $
// 描述　:     积分器基类实现.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "DoubleVector.h"
#include "BaseIntegrator.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
template <int DIM>
BaseIntegrator<DIM>::BaseIntegrator(const std::string& name) {
  d_object_name = name;
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
template <int DIM>
BaseIntegrator<DIM>::~BaseIntegrator() {}

/*****************************************************************************
 * 获取对象名字.
 *****************************************************************************/
template <int DIM>
const std::string& BaseIntegrator<DIM>::getName() {
  return d_object_name;
}

/*****************************************************************************
 * 判断给定名字是否以本对象名字匹配.
 *****************************************************************************/
template <int DIM>
bool BaseIntegrator<DIM>::isMatching(const std::string& name) {
  return (d_object_name == name);
}
template class BaseIntegrator<NDIM>;
