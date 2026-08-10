//
// 文件名:     TetrahedronShapeFunction.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:41:56 2014 $
// 描述　:     四面体有限元形函数类实现.
// 类别　:     %Internal File% ( Don't delete this line )
//

#define co_det(v, m, n)                                               \
  (((m % 2 == 0) ? -1. : 1.) *                                        \
   ((v[(m + 2) % 4][(n + 1) % 3] - v[(m + 1) % 4][(n + 1) % 3]) *     \
        (v[(m + 3) % 4][(n + 2) % 3] - v[(m + 1) % 4][(n + 2) % 3]) - \
    (v[(m + 2) % 4][(n + 2) % 3] - v[(m + 1) % 4][(n + 2) % 3]) *     \
        (v[(m + 3) % 4][(n + 1) % 3] - v[(m + 1) % 4][(n + 1) % 3])))

#include "TetrahedronShapeFunction.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
TetrahedronShapeFunction::TetrahedronShapeFunction(const string& name)
    : BaseShapeFunction<NDIM>(name) {
  d_num_dof = 4;
  d_num_dof_on_entity.resizeArray(4);
  d_num_dof_on_entity[0] = 1;
  d_num_dof_on_entity[1] = 0;
  d_num_dof_on_entity[2] = 0;
  d_num_dof_on_entity[3] = 0;
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
TetrahedronShapeFunction::~TetrahedronShapeFunction() {}

/*****************************************************************************
 * 获取某一个基函数在某一点的值.
 *****************************************************************************/
double TetrahedronShapeFunction::value(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  double tmp_val;
  double volume = get_volume(real_vertex[0], real_vertex[1], real_vertex[2],
                             real_vertex[3]);
  tmp_val = ((j % 2 == 0) ? 1. : -1.) *
            get_volume(pnt, real_vertex[(j + 1) % 4], real_vertex[(j + 2) % 4],
                       real_vertex[(j + 3) % 4]);
  tmp_val /= volume;
  return tmp_val;
}

/*****************************************************************************
 * 获取某一个基函数在某些点的值.
 *****************************************************************************/
tbox::Array<double> TetrahedronShapeFunction::value(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<double> tmp_val(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tmp_val[i] = value(j, real_vertex, pnts[i]);
  }
  return tmp_val;
}

/*****************************************************************************
 * 获取基函数在某一点的值.
 *****************************************************************************/
tbox::Array<double> TetrahedronShapeFunction::value(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<double> tmp_val(4);
  double volume = get_volume(real_vertex[0], real_vertex[1], real_vertex[2],
                             real_vertex[3]);
  tmp_val[0] = get_volume(pnt, real_vertex[1], real_vertex[2], real_vertex[3]);
  tmp_val[0] /= volume;

  tmp_val[1] = get_volume(real_vertex[0], pnt, real_vertex[2], real_vertex[3]);
  tmp_val[1] /= volume;

  tmp_val[2] = get_volume(real_vertex[0], real_vertex[1], pnt, real_vertex[3]);
  tmp_val[2] /= volume;

  tmp_val[3] = get_volume(real_vertex[0], real_vertex[1], real_vertex[2], pnt);
  tmp_val[3] /= volume;
  return tmp_val;
}

/*****************************************************************************
 * 获取基函数在某些点的值.
 *****************************************************************************/
tbox::Array<tbox::Array<double> > TetrahedronShapeFunction::value(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<tbox::Array<double> > tmp_val(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tbox::Array<double> val_i = value(real_vertex, pnts[i]);
    tmp_val[i] = val_i;
  }
  return tmp_val;
}

/*****************************************************************************
 * 获取某一个基函数在某一点的梯度值.
 *****************************************************************************/
tbox::Array<double> TetrahedronShapeFunction::gradient(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<double> tmp_grad(3);

  double volume = get_volume(real_vertex[0], real_vertex[1], real_vertex[2],
                             real_vertex[3]);

  tmp_grad[0] = co_det(real_vertex, j, 0) / volume;
  tmp_grad[1] = co_det(real_vertex, j, 1) / volume;
  tmp_grad[2] = co_det(real_vertex, j, 2) / volume;

  return tmp_grad;
}

/*****************************************************************************
 * 获取某一个基函数在某些点的梯度值.
 *****************************************************************************/
tbox::Array<tbox::Array<double> > TetrahedronShapeFunction::gradient(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<tbox::Array<double> > tmp_grad(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tbox::Array<double> grad = gradient(j, real_vertex, pnts[i]);
    tmp_grad[i] = grad;
  }
  return tmp_grad;
}

/*****************************************************************************
 * 获取基函数在某一点的梯度值.
 *****************************************************************************/
tbox::Array<tbox::Array<double> > TetrahedronShapeFunction::gradient(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<tbox::Array<double> > tmp_grad(4);

  double volume = get_volume(real_vertex[0], real_vertex[1], real_vertex[2],
                             real_vertex[3]);

  for (int i = 0; i < 4; ++i) {
    tmp_grad[i].resizeArray(3);
    tmp_grad[i][0] = co_det(real_vertex, i, 0) / volume;
    tmp_grad[i][1] = co_det(real_vertex, i, 1) / volume;
    tmp_grad[i][2] = co_det(real_vertex, i, 2) / volume;
  }
  return tmp_grad;
}

/*****************************************************************************
 * 获取基函数在某些点的梯度值.
 *****************************************************************************/
tbox::Array<tbox::Array<tbox::Array<double> > >
TetrahedronShapeFunction::gradient(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<tbox::Array<tbox::Array<double> > > tmp_grad(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tbox::Array<tbox::Array<double> > grad = gradient(real_vertex, pnts[i]);
    tmp_grad[i] = grad;
  }
  return tmp_grad;
}

/*****************************************************************************
 * 计算由四个点构成的四面体体积.
 *****************************************************************************/
double TetrahedronShapeFunction::get_volume(
    const hier::DoubleVector<NDIM>& v0, const hier::DoubleVector<NDIM>& v1,
    const hier::DoubleVector<NDIM>& v2, const hier::DoubleVector<NDIM>& v3) {
  return ((v1[0] - v0[0]) * (v2[1] - v0[1]) * (v3[2] - v0[2]) +
          (v1[1] - v0[1]) * (v2[2] - v0[2]) * (v3[0] - v0[0]) +
          (v1[2] - v0[2]) * (v2[0] - v0[0]) * (v3[1] - v0[1]) -
          (v1[0] - v0[0]) * (v2[2] - v0[2]) * (v3[1] - v0[1]) -
          (v1[1] - v0[1]) * (v2[0] - v0[0]) * (v3[2] - v0[2]) -
          (v1[2] - v0[2]) * (v2[1] - v0[1]) * (v3[0] - v0[0]));
}

/****************************************************************************
 * 获取自由度数目.
 ****************************************************************************/
int TetrahedronShapeFunction::getNumberOfDof() { return d_num_dof; }

tbox::Array<int> TetrahedronShapeFunction::getNumberOfDofOnEntity() {
  return d_num_dof_on_entity;
}
