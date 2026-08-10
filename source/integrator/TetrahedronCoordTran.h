//
// 文件名:     TetrahedronCoordTran.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:43:18 2014 $
// 描述　:     四面体仿射变换类
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_TetrahedronCoordTran
#define included_fem_TetrahedronCoordTran

#include "Array.h"
#include "Pointer.h"
#include "CoordinateTransform.h"
#include "DoubleVector.h"
using namespace JAUMIN;

/*!
 * @brief 有限元坐标变换类，支持模板单元上的积分点和实际网格单元上的积分点转换。
 *
 */

class TetrahedronCoordTran : public CoordinateTransform<NDIM> {
public:
  /**
   * @brief 构造函数。
   *
   * @param template_vertex 输入参数，坐标数组，模板单元的结点坐标。
   */
  TetrahedronCoordTran(tbox::Array<hier::DoubleVector<NDIM> >& template_vertex);

  /**
   * @brief 析构函数。
   */
  ~TetrahedronCoordTran();

  /**
   * @brief 把一个局部积分点映射到一个全局积分点。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   * @param lp              输入参数，坐标，局部积分点。
   * @param gp              输出参数，坐标，全局积分点。
   *
   */
  void local2Global(tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
                    const hier::DoubleVector<NDIM>& lp,
                    hier::DoubleVector<NDIM>& gp);

  /**
   * @brief 把一个全局积分点映射到一个局部积分点。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   * @param gp              输出参数，坐标，全局积分点。
   * @param lp              输入参数，坐标，局部积分点。
   *
   */
  void global2Local(tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
                    const hier::DoubleVector<NDIM>& gp,
                    hier::DoubleVector<NDIM>& lp);

  /**
   * @brief 局部到全局仿射变换的Jacibian行列式。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   *
   * @return 双精度浮点型，局部到全局仿射变换的Jacibian行列式。
   */
  double Local2GlobalJacobian(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex);

  /**
   * @brief 全局到局部仿射变换的Jacibian行列式。
   *
   * @param real_vertex     输入参数，坐标数组，实际单元的结点坐标.
   *
   * @return 双精度浮点型，全局到局部仿射变换的Jacibian行列式。
   */
  double Global2LocalJacobian(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex);

private:
  /**
   * @brief 计算四个点构成的四面体的体积.
   *
   * @param v0 输入参数, 双精度数组, 坐标.
   * @param v1 输入参数, 双精度数组, 坐标.
   * @param v2 输入参数, 双精度数组, 坐标.
   * @param v3 输入参数, 双精度数组, 坐标.
   *
   * @return 双精度浮点型, 四面体的体积.
   */
  double get_volume(const hier::DoubleVector<NDIM>& v0,
                    const hier::DoubleVector<NDIM>& v1,
                    const hier::DoubleVector<NDIM>& v2,
                    const hier::DoubleVector<NDIM>& v3);

  tbox::Array<hier::DoubleVector<NDIM> >
      d_template_vertex; /**< 模板单元的结点坐标 */
};
#endif
