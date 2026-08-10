//
// 文件名:     TriangleCoordTran.C
// 软件包:     JAUMIN fem
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-02-16 12:55:07 +0800 (四, 2012-02-16) $
// 描述　:
// 类别　:     %Internal File% ( Don't delete this line )
//

#define AREA(a, b, c) \
  ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))

#include "DoubleVector.h"
#include "TriangleCoordTran.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
TriangleCoordTran::TriangleCoordTran(
    tbox::Array<hier::DoubleVector<NDIM> > &template_vertex) {
  d_template_vertex = template_vertex;
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
TriangleCoordTran::~TriangleCoordTran() {}

/*****************************************************************************
 * 把一个模板单元上的局部积分点映射到实际网格上的积分点.
 *****************************************************************************/
void TriangleCoordTran::local2Global(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex,
    const hier::DoubleVector<NDIM> &lp, hier::DoubleVector<NDIM> &gp) {
  double lambda[3];
  double area =
      AREA(d_template_vertex[0], d_template_vertex[1], d_template_vertex[2]);
  lambda[0] = AREA(lp, d_template_vertex[1], d_template_vertex[2]) / area;
  lambda[1] = AREA(lp, d_template_vertex[2], d_template_vertex[0]) / area;
  lambda[2] = AREA(lp, d_template_vertex[0], d_template_vertex[1]) / area;
  gp[0] = lambda[0] * real_vertex[0][0] + lambda[1] * real_vertex[1][0] +
          lambda[2] * real_vertex[2][0];
  gp[1] = lambda[0] * real_vertex[0][1] + lambda[1] * real_vertex[1][1] +
          lambda[2] * real_vertex[2][1];
}

/*****************************************************************************
 * 把一个实际单元上的积分点映射到模板单元上的局部积分点.
 *****************************************************************************/
void TriangleCoordTran::global2Local(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex,
    const hier::DoubleVector<NDIM> &gp, hier::DoubleVector<NDIM> &lp) {
  double lambda[3];
  double area =
      AREA(d_template_vertex[0], d_template_vertex[1], d_template_vertex[2]);
  lambda[0] = AREA(gp, d_template_vertex[1], d_template_vertex[2]) / area;
  lambda[1] = AREA(gp, d_template_vertex[2], d_template_vertex[0]) / area;
  lambda[2] = AREA(gp, d_template_vertex[0], d_template_vertex[1]) / area;
  lp[0] = lambda[0] * d_template_vertex[0][0] +
          lambda[1] * d_template_vertex[1][0] +
          lambda[2] * d_template_vertex[2][0];
  lp[1] = lambda[0] * d_template_vertex[0][1] +
          lambda[1] * d_template_vertex[1][1] +
          lambda[2] * d_template_vertex[2][1];
}

/*****************************************************************************
 * 从模板单元到实际网格Jacibian变换行列式. .
 *****************************************************************************/
double TriangleCoordTran::Local2GlobalJacobian(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {
  double larea =
      AREA(d_template_vertex[0], d_template_vertex[1], d_template_vertex[2]);
  double garea = AREA(real_vertex[0], real_vertex[1], real_vertex[2]);
  return garea / larea;
}

/*****************************************************************************
 * 从实际网格到模板单元Jacibian变换行列式. .
 *****************************************************************************/
double TriangleCoordTran::Global2LocalJacobian(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {
  double larea =
      AREA(d_template_vertex[0], d_template_vertex[1], d_template_vertex[2]);
  double garea = AREA(real_vertex[0], real_vertex[1], real_vertex[2]);
  return larea / garea;
}
