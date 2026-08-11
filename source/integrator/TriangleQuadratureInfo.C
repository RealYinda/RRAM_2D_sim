//
// 文件名:     TriangleQuadratureInfo.C
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-02-16 12:55:07 +0800 (四, 2012-02-16) $
// 描述　:     有限元积分信息类实现.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "DoubleVector.h"
#include "TriangleQuadratureInfo.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
TriangleQuadratureInfo::TriangleQuadratureInfo(int algebric_accuracy) {
  d_volume = 0.5;
  if (algebric_accuracy == 1) {
    /// 若积分精度为 1 则设置积分信息为:
    d_number_quad_points = 1;
    d_algebric_accuracy = algebric_accuracy;
    d_quad_points.resizeArray(1);
    d_quad_weights.resizeArray(1);
    for (int i = 0; i < d_number_quad_points; ++i) {
      d_quad_weights[i] = 1.0;
    }
    d_quad_points[0][0] = 0.33333333333333333333;
    d_quad_points[0][1] = 0.33333333333333333333;
  }

  if (algebric_accuracy == 2) {
    /// 若积分精度为 2 则设置积分信息为:
    d_number_quad_points = 3;
    d_algebric_accuracy = algebric_accuracy;
    d_quad_points.resizeArray(3);
    d_quad_weights.resizeArray(3);
    for (int i = 0; i < d_number_quad_points; ++i) {
      d_quad_weights[i] = 0.33333333333333333333;
    }
    d_quad_points[0][0] = 0.16666666666666666667;
    d_quad_points[0][1] = 0.16666666666666666667;
    d_quad_points[1][0] = 0.16666666666666666667;
    d_quad_points[1][1] = 0.66666666666666666667;
    d_quad_points[2][0] = 0.66666666666666666667;
    d_quad_points[2][1] = 0.16666666666666666667;
  }
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
TriangleQuadratureInfo::~TriangleQuadratureInfo() {}

/*****************************************************************************
 * 获取积分点.
 *****************************************************************************/
const tbox::Array<hier::DoubleVector<NDIM> >&
TriangleQuadratureInfo::getQuadraturePoints() {
  return d_quad_points;
}

/*****************************************************************************
 * 获取积分点权重.
 *****************************************************************************/
const tbox::Array<double>& TriangleQuadratureInfo::getQuadratureWeights() {
  return d_quad_weights;
}

/*****************************************************************************
 * 获取积分点数目.
 *****************************************************************************/
int TriangleQuadratureInfo::getNumberOfQuadraturePoints() {
  return d_number_quad_points;
}

/*****************************************************************************
 * 获取积分单元体积.
 *****************************************************************************/
double TriangleQuadratureInfo::getElementVolume() { return d_volume; }
