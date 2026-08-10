//
// 文件名:     TetrahedronIntegrator.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:42:50 2014 $
// 描述　:     四面体线性连续有限元积分器实现.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "DoubleVector.h"
#include "TetrahedronShapeFunction.h"
#include "TetrahedronIntegrator.h"

/****************************************************************************
 * 构造函数.
 ****************************************************************************/
TetrahedronIntegrator::TetrahedronIntegrator(int algebric_accuracy,
                                             const string &name)
    : BaseIntegrator<NDIM>(name) {
  d_num_vertex = 4;
  d_template_vertex.resizeArray(d_num_vertex);

  d_template_vertex[0][0] = 0.0;
  d_template_vertex[0][1] = 0.0;
  d_template_vertex[0][2] = 0.0;
  d_template_vertex[1][0] = 1.0;
  d_template_vertex[1][1] = 0.0;
  d_template_vertex[1][2] = 0.0;
  d_template_vertex[2][0] = 0.0;
  d_template_vertex[2][1] = 1.0;
  d_template_vertex[2][2] = 0.0;
  d_template_vertex[3][0] = 0.0;
  d_template_vertex[3][1] = 0.0;
  d_template_vertex[3][2] = 1.0;

  /// 坐标变换.
  d_coord_transform = new TetrahedronCoordTran(d_template_vertex);
  /// 积分信息.
  d_quad_info = new TetrahedronQuadratureInfo(algebric_accuracy);
}

/****************************************************************************
 * 析构函数.
 ****************************************************************************/
TetrahedronIntegrator::~TetrahedronIntegrator() {}

/****************************************************************************
 * 获取积分点.
 ****************************************************************************/
tbox::Array<hier::DoubleVector<NDIM> >
TetrahedronIntegrator::getQuadraturePoints(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {
  tbox::Array<hier::DoubleVector<NDIM> > l_pnt =
      d_quad_info->getQuadraturePoints();
  int num_pnt = l_pnt.getSize();
  tbox::Array<hier::DoubleVector<NDIM> > g_pnt(num_pnt);
  for (int i = 0; i < num_pnt; ++i) {
    d_coord_transform->local2Global(real_vertex, l_pnt[i], g_pnt[i]);
  }
  return g_pnt;
}

/****************************************************************************
 * 获取积分点权重.
 ****************************************************************************/
const tbox::Array<double> &TetrahedronIntegrator::getQuadratureWeights() {
  return d_quad_info->getQuadratureWeights();
}

/****************************************************************************
 * 获取积分点数目.
 ****************************************************************************/
int TetrahedronIntegrator::getNumberOfQuadraturePoints() {
  return d_quad_info->getNumberOfQuadraturePoints();
}

/****************************************************************************
 * 获取单元的面积(体积).
 ****************************************************************************/
double TetrahedronIntegrator::getElementVolume() {
  return d_quad_info->getElementVolume();
}

/****************************************************************************
 * 获取仿射变换的矩阵行列式.
 ****************************************************************************/
double TetrahedronIntegrator::getLocal2GlobalJacobian(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {
  return d_coord_transform->Local2GlobalJacobian(real_vertex);
}

/****************************************************************************
 * 获取仿射变换的矩阵行列式.
 ****************************************************************************/
double TetrahedronIntegrator::getGlobal2LocalJacobian(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {
  return d_coord_transform->Global2LocalJacobian(real_vertex);
}

/****************************************************************************
 * 获取结点数目.
 ****************************************************************************/
int TetrahedronIntegrator::getNumberOfVertex() { return d_num_vertex; }
