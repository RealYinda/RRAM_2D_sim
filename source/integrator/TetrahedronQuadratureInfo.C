//
// 文件名:     TetrahedronQuadratureInfo.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:42:35 2014 $
// 描述　:     四面体有限元积分信息类实现.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "TetrahedronQuadratureInfo.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
TetrahedronQuadratureInfo::TetrahedronQuadratureInfo(int algebric_accuracy) {
  d_volume = 0.1666666666666667;
  if (algebric_accuracy == 1) {
    /// 若积分精度为 1 则设置积分信息为:
    d_number_quad_points = 1;
    d_algebric_accuracy = algebric_accuracy;
    d_quad_points.resizeArray(1);
    d_quad_weights.resizeArray(1);
    for (int i = 0; i < d_number_quad_points; ++i) {
      d_quad_weights[i] = 1.0;
    }
    d_quad_points[0][0] = 0.25;
    d_quad_points[0][1] = 0.25;
    d_quad_points[0][2] = 0.25;
  }

  if (algebric_accuracy == 2) {
    /// 若积分精度为 2 则设置积分信息为:
    d_number_quad_points = 4;
    d_algebric_accuracy = algebric_accuracy;
    d_quad_points.resizeArray(4);
    d_quad_weights.resizeArray(4);
    for (int i = 0; i < d_number_quad_points; ++i) {
      d_quad_weights[i] = 0.25;
    }
    d_quad_points[0][0] = 0.13819660;
    d_quad_points[0][1] = 0.13819660;
    d_quad_points[0][2] = 0.13819660;
    d_quad_points[1][0] = 0.58541020;
    d_quad_points[1][1] = 0.13819660;
    d_quad_points[1][2] = 0.13819660;
    d_quad_points[2][0] = 0.13819660;
    d_quad_points[2][1] = 0.58541020;
    d_quad_points[2][2] = 0.13819660;
    d_quad_points[3][0] = 0.13819660;
    d_quad_points[3][1] = 0.13819660;
    d_quad_points[3][2] = 0.58541020;
  }
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
TetrahedronQuadratureInfo::~TetrahedronQuadratureInfo() {}

/*****************************************************************************
 * 获取积分点.
 *****************************************************************************/
const tbox::Array<hier::DoubleVector<NDIM> >&
TetrahedronQuadratureInfo::getQuadraturePoints() {
  return d_quad_points;
}

/*****************************************************************************
 * 获取积分点权重.
 *****************************************************************************/
const tbox::Array<double>& TetrahedronQuadratureInfo::getQuadratureWeights() {
  return d_quad_weights;
}

/*****************************************************************************
 * 获取积分点数目.
 *****************************************************************************/
int TetrahedronQuadratureInfo::getNumberOfQuadraturePoints() {
  return d_number_quad_points;
}

/*****************************************************************************
 * 获取积分单元体积.
 *****************************************************************************/
double TetrahedronQuadratureInfo::getElementVolume() { return d_volume; }
