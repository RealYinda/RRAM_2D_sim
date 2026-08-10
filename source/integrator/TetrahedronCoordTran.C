//
// 文件名:     TetrahedronCoordTran.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:43:41 2014 $
// 描述　:     四面体仿射变换类
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "TetrahedronCoordTran.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
TetrahedronCoordTran::TetrahedronCoordTran(
    tbox::Array<hier::DoubleVector<NDIM> > &template_vertex) {
  d_template_vertex = template_vertex;
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
TetrahedronCoordTran::~TetrahedronCoordTran() {}

/*****************************************************************************
 * 把一个模板单元上的局部积分点映射到实际网格上的积分点.
 *****************************************************************************/
void TetrahedronCoordTran::local2Global(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex,
    const hier::DoubleVector<NDIM> &lp, hier::DoubleVector<NDIM> &gp) {
  double lambda[4];
  double volume = get_volume(d_template_vertex[0], d_template_vertex[1],
                             d_template_vertex[2], d_template_vertex[3]);
  lambda[0] = get_volume(lp, d_template_vertex[1], d_template_vertex[2],
                         d_template_vertex[3]) /
              volume;
  lambda[1] = get_volume(d_template_vertex[0], lp, d_template_vertex[2],
                         d_template_vertex[3]) /
              volume;
  lambda[2] = get_volume(d_template_vertex[0], d_template_vertex[1], lp,
                         d_template_vertex[3]) /
              volume;
  lambda[3] = get_volume(d_template_vertex[0], d_template_vertex[1],
                         d_template_vertex[2], lp) /
              volume;
  gp[0] = lambda[0] * real_vertex[0][0] + lambda[1] * real_vertex[1][0] +
          lambda[2] * real_vertex[2][0] + lambda[3] * real_vertex[3][0];
  gp[1] = lambda[0] * real_vertex[0][1] + lambda[1] * real_vertex[1][1] +
          lambda[2] * real_vertex[2][1] + lambda[3] * real_vertex[3][1];
  gp[2] = lambda[0] * real_vertex[0][2] + lambda[1] * real_vertex[1][2] +
          lambda[2] * real_vertex[2][2] + lambda[3] * real_vertex[3][2];
}

/*****************************************************************************
 * 把一个实际单元上的积分点映射到模板单元上的局部积分点.
 *****************************************************************************/
void TetrahedronCoordTran::global2Local(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex,
    const hier::DoubleVector<NDIM> &gp, hier::DoubleVector<NDIM> &lp) {
  double lambda[4];
  double volume = get_volume(real_vertex[0], real_vertex[1], real_vertex[2],
                             real_vertex[3]);
  lambda[0] =
      get_volume(gp, real_vertex[1], real_vertex[2], real_vertex[3]) / volume;
  lambda[1] =
      get_volume(real_vertex[0], gp, real_vertex[2], real_vertex[3]) / volume;
  lambda[2] =
      get_volume(real_vertex[0], real_vertex[1], gp, real_vertex[3]) / volume;
  lambda[3] =
      get_volume(real_vertex[0], real_vertex[1], real_vertex[2], gp) / volume;
  lp[0] = lambda[0] * d_template_vertex[0][0] +
          lambda[1] * d_template_vertex[1][0] +
          lambda[2] * d_template_vertex[2][0] +
          lambda[3] * d_template_vertex[3][0];
  lp[1] = lambda[0] * d_template_vertex[0][1] +
          lambda[1] * d_template_vertex[1][1] +
          lambda[2] * d_template_vertex[2][1] +
          lambda[3] * d_template_vertex[3][1];
  lp[2] = lambda[0] * d_template_vertex[0][2] +
          lambda[1] * d_template_vertex[1][2] +
          lambda[2] * d_template_vertex[2][2] +
          lambda[3] * d_template_vertex[3][2];
}

/*****************************************************************************
 * 从模板单元到实际网格Jacibian变换行列式.
 *****************************************************************************/
double TetrahedronCoordTran::Local2GlobalJacobian(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {
  double lvolume = get_volume(d_template_vertex[0], d_template_vertex[1],
                              d_template_vertex[2], d_template_vertex[3]);
  double gvolume = get_volume(real_vertex[0], real_vertex[1], real_vertex[2],
                              real_vertex[3]);
  return gvolume / lvolume;
}

/*****************************************************************************
 * 从实际网格到模板单元Jacibian变换行列式.
 *****************************************************************************/
double TetrahedronCoordTran::Global2LocalJacobian(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {
  double lvolume = get_volume(d_template_vertex[0], d_template_vertex[1],
                              d_template_vertex[2], d_template_vertex[3]);
  double gvolume = get_volume(real_vertex[0], real_vertex[1], real_vertex[2],
                              real_vertex[3]);
  return lvolume / gvolume;
}

/*****************************************************************************
 * 计算由四个点构成的四面体体积.
 *****************************************************************************/
double TetrahedronCoordTran::get_volume(const hier::DoubleVector<NDIM> &v0,
                                        const hier::DoubleVector<NDIM> &v1,
                                        const hier::DoubleVector<NDIM> &v2,
                                        const hier::DoubleVector<NDIM> &v3) {
  return ((v1[0] - v0[0]) * (v2[1] - v0[1]) * (v3[2] - v0[2]) +
          (v1[1] - v0[1]) * (v2[2] - v0[2]) * (v3[0] - v0[0]) +
          (v1[2] - v0[2]) * (v2[0] - v0[0]) * (v3[1] - v0[1]) -
          (v1[0] - v0[0]) * (v2[2] - v0[2]) * (v3[1] - v0[1]) -
          (v1[1] - v0[1]) * (v2[0] - v0[0]) * (v3[2] - v0[2]) -
          (v1[2] - v0[2]) * (v2[1] - v0[1]) * (v3[0] - v0[0]));
}
