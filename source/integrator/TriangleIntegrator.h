//
// 文件名:     TriangleIntegrator.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     三角形积分器.
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_TriangleIntegrator
#define included_fem_TriangleIntegrator

#include "Array.h"
#include "Pointer.h"
#include "BaseIntegrator.h"

namespace JAUMIN {
namespace hier {
template <int DIM>
class DoubleVector;
}
}
class TriangleQuadratureInfo;
class TriangleCoordTran;

/*!
 * @brief 积分器类，提供单元上有限元积分计算所需要的信息。该类根据用户提
 * 供的网格单元的结点信息和积分精度初始化有限元模板单元并获取有限元积分计算需要的
 * 信息，其中用到了TriangleQuadratureInfo，TriangleCoordTran两个成员类。
 *
 * @see fem::TriangleQuadratureInfo, fem::TriangleCoordTran
 */

class TriangleIntegrator : public BaseIntegrator<NDIM> {
public:
  /**
   * @brief 构造函数。
   *
   * @param algebric_accuracy 输入参数，整型，积分精度。
   *
   * @param name              输入参数, 字符串, 积分器的名字.
   *
   * @return
   */
  TriangleIntegrator(int algebric_accuracy, const string &name);

  /**
   * @brief 析构函数。
   */
  ~TriangleIntegrator();

  /**
   * @brief 获取单元的积分点。
   *
   * @return 坐标数组，单元的积分点数组。
   */
  virtual tbox::Array<hier::DoubleVector<NDIM> > getQuadraturePoints(
      tbox::Array<hier::DoubleVector<NDIM> > &real_vertex);

  /**
   * @brief 获取积分点权重。
   *
   * @return 双精度浮点型数组，积分点权重数组。
   */
  virtual const tbox::Array<double> &getQuadratureWeights();

  /**
   * @brief 获取积分点数目。
   *
   * @return 整型，积分点数目。
   */
  virtual int getNumberOfQuadraturePoints();

  /**
   * @brief 获取单元的面积(体积)。
   *
   * @return 双精度浮点型，单元面积(体积)。
   */
  virtual double getElementVolume();

  /**
   * @brief 获取由局部到全局仿射变换的矩阵行列式。
   *
   * @param real_vertex 输入参数, 坐标数组, 单元的结点坐标.
   *
   * @return 双精度浮点型，仿射变换的矩阵行列式。
   */
  virtual double getLocal2GlobalJacobian(
      tbox::Array<hier::DoubleVector<NDIM> > &real_vertex);

  /**
   * @brief 获取由全局到局部仿射变换的矩阵行列式。
   *
   * @param real_vertex 输入参数, 坐标数组, 单元的结点坐标.
   *
   * @return 双精度浮点型，仿射变换的矩阵行列式。
   */
  virtual double getGlobal2LocalJacobian(
      tbox::Array<hier::DoubleVector<NDIM> > &real_vertex);

  /**
   * @brief 获取单元结点数目。
   *
   * @return 整型，单元结点数目。
   */
  virtual int getNumberOfVertex();

  /**
   * @brief 获取局部编号为i的结点的实际坐标.
   *
   * @param i 输入参数, 整型, 结点的局部编号.
   *
   * @return 局部编号为i的结点的坐标.
   */
  hier::DoubleVector<NDIM> getCoordinate(int i);

private:
  tbox::Pointer<TriangleCoordTran> d_coord_transform; /**< 坐标变换 */
  tbox::Pointer<TriangleQuadratureInfo> d_quad_info;  /**< 积分信息 */
  tbox::Array<hier::DoubleVector<NDIM> >
      d_template_vertex; /**< 有限元模板单元的结点 */
  int d_num_vertex;      /**< 结点数目 */
};
#endif
