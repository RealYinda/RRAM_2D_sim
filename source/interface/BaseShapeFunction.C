//
// 文件名:     BaseShapeFunction.C
// 软件包:     JAUMIN fem
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-02-16 12:55:07 +0800 (四, 2012-02-16) $
// 描述　:     有限元基函数类.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "DoubleVector.h"
#include "BaseShapeFunction.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
template <int DIM>
BaseShapeFunction<DIM>::BaseShapeFunction(const string& name) {
  d_object_name = name;
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
template <int DIM>
BaseShapeFunction<DIM>::~BaseShapeFunction() {}

/*****************************************************************************
 * 判断给定名字是否以本对象名字匹配.
 *****************************************************************************/
template <int DIM>
bool BaseShapeFunction<DIM>::isMatching(const string& name) {
  return (d_object_name == name);
}

/*****************************************************************************
 * 获取对象名字.
 *****************************************************************************/
template <int DIM>
const string& BaseShapeFunction<DIM>::getName() {
  return d_object_name;
}

template <int DIM>
tbox::Array<int> BaseShapeFunction<DIM>::getNumberOfDofOnEntity() {
  TBOX_WARNING(
      "Function getNumberOfDofOnEntity is not implemented, please make sure "
      "your code is right! ");
  tbox::Array<int> rt;
  return rt;
}
template class BaseShapeFunction<NDIM>;
