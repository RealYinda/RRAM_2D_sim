//
// 文件名:     CoordinateTransform.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:44:30 2014 $
// 描述　:     坐标变换基类
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_CoordinateTransform
#define included_fem_CoordinateTransform

#include "Array.h"
#include "Pointer.h"

namespace JAUMIN {
namespace hier {
template <int DIM>
class DoubleVector;
}
}
using namespace JAUMIN;
/*!
 * @brief 有限元单元坐标变换类，提供坐标变换的接口, 各个形状单元的坐标变换继承该
 *        基类, 并实现接口函数。
 *
 */

template <int DIM>
class CoordinateTransform {
public:
  /**
   * @brief 构造函数。
   *
   */
  CoordinateTransform();

  /**
   * @brief 析构函数。
   */
  virtual ~CoordinateTransform();

  /**
   * @brief 把一个局部积分点映射到一个全局积分点。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   * @param lp              输入参数，坐标，局部积分点。
   * @param gp              输出参数，坐标，全局积分点。
   *
   */
  virtual void local2Global(tbox::Array<hier::DoubleVector<DIM> >& real_vertex,
                            const hier::DoubleVector<DIM>& lp,
                            hier::DoubleVector<DIM>& gp) = 0;

  /**
   * @brief 把一个全局积分点映射到一个局部积分点。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   * @param gp              输出参数，坐标，全局积分点。
   * @param lp              输入参数，坐标，局部积分点。
   *
   */
  virtual void global2Local(tbox::Array<hier::DoubleVector<DIM> >& real_vertex,
                            const hier::DoubleVector<DIM>& gp,
                            hier::DoubleVector<DIM>& lp) = 0;

  /**
   * @brief 局部到全局仿射变换的Jacibian行列式。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   *
   * @return 双精度浮点型，局部到全局仿射变换的Jacibian行列式。
   */
  virtual double Local2GlobalJacobian(
      tbox::Array<hier::DoubleVector<DIM> >& real_vertex) = 0;

  /**
   * @brief 全局到局部仿射变换的Jacibian行列式。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   *
   * @return 双精度浮点型，全局到局部仿射变换的Jacibian行列式。
   */
  virtual double Global2LocalJacobian(
      tbox::Array<hier::DoubleVector<DIM> >& real_vertex) = 0;

private:
};
#endif
