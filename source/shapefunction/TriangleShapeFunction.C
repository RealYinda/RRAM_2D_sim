//
// 文件名:     TriangleShapeFunction.C
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-02-16 12:55:07 +0800 (四, 2012-02-16) $
// 描述　:     有限元基函数类.
// 类别　:     %Internal File% ( Don't delete this line )
//

#define AREA(a, b, c)                      \
  (((b)[0] - (a)[0]) * ((c)[1] - (a)[1]) - \
   ((b)[1] - (a)[1]) * ((c)[0] - (a)[0]))

#define GRADIENT(gradient, p0, p1, vol)          \
  {                                              \
    (gradient)[0] = ((p0)[1] - (p1)[1]) / (vol); \
    (gradient)[1] = ((p1)[0] - (p0)[0]) / (vol); \
  }

#include "DoubleVector.h"
#include "TriangleShapeFunction.h"

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
TriangleShapeFunction::TriangleShapeFunction(const string& name)
    : BaseShapeFunction<NDIM>(name) {
  d_num_dof = 3;
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
TriangleShapeFunction::~TriangleShapeFunction() {}

/*****************************************************************************
 * 获取某一个基函数在某一点的值.
 *****************************************************************************/
double TriangleShapeFunction::value(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  double tmp_val;
  double area = AREA(real_vertex[0], real_vertex[1], real_vertex[2]);
  tmp_val = AREA(pnt, d_real_vertex[(j + 1) % 3], d_real_vertex[(j + 2) % 3]);
  tmp_val /= area;
  return tmp_val;
}

/*****************************************************************************
 * 获取某一个基函数在某些点的值.
 *****************************************************************************/
tbox::Array<double> TriangleShapeFunction::value(
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
tbox::Array<double> TriangleShapeFunction::value(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<double> tmp_val(3);
  double area = AREA(real_vertex[0], real_vertex[1], real_vertex[2]);
  tmp_val[0] = AREA(pnt, real_vertex[1], real_vertex[2]);
  tmp_val[0] /= area;

  tmp_val[1] = AREA(real_vertex[0], pnt, real_vertex[2]);
  tmp_val[1] /= area;

  tmp_val[2] = AREA(real_vertex[0], real_vertex[1], pnt);
  tmp_val[2] /= area;
  return tmp_val;
}

/*****************************************************************************
 * 获取基函数在某些点的值.
 *****************************************************************************/
tbox::Array<tbox::Array<double> > TriangleShapeFunction::value(
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
tbox::Array<double> TriangleShapeFunction::gradient(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<double> tmp_grad(2);

  double area = AREA(real_vertex[0], real_vertex[1], real_vertex[2]);
  GRADIENT(tmp_grad, real_vertex[(j + 1) % 3], real_vertex[(j + 2) % 3], area);

  return tmp_grad;
}

/*****************************************************************************
 * 获取某一个基函数在某些点的梯度值.
 *****************************************************************************/
tbox::Array<tbox::Array<double> > TriangleShapeFunction::gradient(
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
tbox::Array<tbox::Array<double> > TriangleShapeFunction::gradient(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<tbox::Array<double> > tmp_grad(3);

  double area = AREA(real_vertex[0], real_vertex[1], real_vertex[2]);

  tbox::Array<double> val_0(2);
  GRADIENT(val_0, real_vertex[1], real_vertex[2], area);
  tmp_grad[0] = val_0;

  tbox::Array<double> val_1(2);
  GRADIENT(val_1, real_vertex[2], real_vertex[0], area);
  tmp_grad[1] = val_1;

  tbox::Array<double> val_2(2);
  GRADIENT(val_2, real_vertex[0], real_vertex[1], area);
  tmp_grad[2] = val_2;

  return tmp_grad;
}

/*****************************************************************************
 * 获取基函数在某些点的梯度值.
 *****************************************************************************/
tbox::Array<tbox::Array<tbox::Array<double> > > TriangleShapeFunction::gradient(
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

/****************************************************************************
 * 获取自由度数目.
 ****************************************************************************/
int TriangleShapeFunction::getNumberOfDof() { return d_num_dof; }
