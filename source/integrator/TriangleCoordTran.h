//
// 文件名:     TriangleCoordTran.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-02-16 12:55:07 +0800 (四, 2012-02-16) $
// 描述　:     有限元仿射变换类。
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_TriangleCoordTran
#define included_fem_TriangleCoordTran

#include "Array.h"
#include "Pointer.h"
#include "CoordinateTransform.h"

namespace JAUMIN {
namespace hier {
template <int DIM>
class DoubleVector;
}
}
using namespace JAUMIN;
/*!
 * @brief 有限元模板单元坐标变换类，把模板单元上的积分点转换成实际网格单元上的积
 * 分点。
 *
 */

class TriangleCoordTran : public CoordinateTransform<NDIM> {
public:
  /**
   * @brief 构造函数。
   *
   * @param template_vertex 输入参数，坐标数组，模板单元的结点坐标。
   */
  TriangleCoordTran(tbox::Array<hier::DoubleVector<NDIM> > &template_vertex);

  /**
   * @brief 析构函数。
   */
  ~TriangleCoordTran();

  /**
   * @brief 把一个局部积分点映射到一个全局积分点。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   * @param lp              输入参数，坐标，局部积分点。
   * @param gp              输出参数，坐标，全局积分点。
   *
   */
  void local2Global(tbox::Array<hier::DoubleVector<NDIM> > &real_vertex,
                    const hier::DoubleVector<NDIM> &lp,
                    hier::DoubleVector<NDIM> &gp);

  /**
   * @brief 把一个全局积分点映射到一个局部积分点。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   * @param gp              输出参数，坐标，全局积分点。
   * @param lp              输入参数，坐标，局部积分点。
   *
   */
  void global2Local(tbox::Array<hier::DoubleVector<NDIM> > &real_vertex,
                    const hier::DoubleVector<NDIM> &gp,
                    hier::DoubleVector<NDIM> &lp);

  /**
   * @brief 局部到全局仿射变换的Jacibian行列式。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   *
   * @return 双精度浮点型，局部到全局仿射变换的Jacibian行列式。
   */
  double Local2GlobalJacobian(
      tbox::Array<hier::DoubleVector<NDIM> > &real_vertex);

  /**
   * @brief 全局到局部仿射变换的Jacibian行列式。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   *
   * @return 双精度浮点型，全局到局部仿射变换的Jacibian行列式。
   */
  double Global2LocalJacobian(
      tbox::Array<hier::DoubleVector<NDIM> > &real_vertex);

private:
  tbox::Array<hier::DoubleVector<NDIM> >
      d_template_vertex; /**< 模板单元的结点坐标 */
};
#endif
