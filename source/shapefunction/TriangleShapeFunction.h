//
// 文件名:     TriangleShapeFunction.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-02-16 12:55:07 +0800 (四, 2012-02-16) $
// 描述　:     有限元基函数类。
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_TriangleShapeFunction
#define included_fem_TriangleShapeFunction

#include "Array.h"
#include "Pointer.h"
#include "BaseShapeFunction.h"

namespace JAUMIN {
namespace hier {
template <int DIM>
class DoubleVector;
}
}
/*!
 * @brief 有限元模板单元基函数类，为有限元积分计算提供基函数在积分点上的函数值和
 * 梯度值。
 *
 */

class TriangleShapeFunction : public BaseShapeFunction<NDIM> {
public:
  /**
   * @brief 构造函数。
   *
   */
  TriangleShapeFunction(const string& name);

  /**
   * @brief 析构函数。
   */
  ~TriangleShapeFunction();

  /**
   * @brief 取出某一个积分点上的基函数值。
   *
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数，坐标，积分点坐标。
   *
   * @return 双精度型数组，积分点的基函数值。
   */
  tbox::Array<double> value(tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
                            const hier::DoubleVector<NDIM>& pnt);

  /**
   * @brief 取出某一个积分点上的某一个基函数值。
   *
   * @param j            输入参数, 整型, 积分点编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数，坐标，积分点坐标。
   *
   * @return 双精度型，积分点的基函数值。
   */
  double value(int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
               const hier::DoubleVector<NDIM>& pnt);

  /**
   * @brief 取出若干个积分点上的某一个基函数值。
   *
   * @param j            输入参数, 整型, 基函数编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts         输入参数，坐标数组，积分点坐标数组。
   *
   * @return 双精度型数组，积分点数组的基函数值。
   */
  tbox::Array<double> value(int j,
                            tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
                            const tbox::Array<hier::DoubleVector<NDIM> >& pnts);

  /**
   * @brief 取出若干个积分点上的基函数值。
   *
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts 输入参数，坐标数组，积分点坐标数组。
   *
   * @return 双精度型，积分点数组的基函数值。
   */
  tbox::Array<tbox::Array<double> > value(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<NDIM> >& pnts);

  /**
   * @brief 取出某一个积分点上的基函数梯度。
   *
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt 输入参数，坐标，积分点坐标。
   *
   * @return 双精度型数组，积分点的基函数梯度。
   */
  tbox::Array<tbox::Array<double> > gradient(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const hier::DoubleVector<NDIM>& pnt);

  /**
   * @brief 取出某一个积分点上的某一个基函数梯度。
   *
   * @param j            输入参数, 整型, 基函数编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数，坐标，积分点坐标.
   *
   * @return 双精度型数组，积分点的基函数梯度。
   */
  tbox::Array<double> gradient(
      int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const hier::DoubleVector<NDIM>& pnt);

  /**
   * @brief 取出某若干积分点上的基函数梯度。
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts 输入参数，坐标，积分点坐标数组。
   *
   * @return 二维双精度型数组，积分点数组的基函数梯度。
   */
  tbox::Array<tbox::Array<tbox::Array<double> > > gradient(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<NDIM> >& pnts);

  /**
   * @brief 取出某若干积分点上的某一个基函数梯度。
   *
   * @param j 输入参数, 整型, 基函数编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts 输入参数，坐标，积分点坐标数组。
   *
   * @return 二维双精度型数组，积分点数组的基函数梯度。
   */
  virtual tbox::Array<tbox::Array<double> > gradient(
      int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<NDIM> >& pnts);

  /**
   * @brief 获取单元上自由度的数目.
   *
   * @return 整型, 单元自由度的数目.
   */
  virtual int getNumberOfDof();

private:
  tbox::Array<hier::DoubleVector<NDIM> >
      d_real_vertex; /**< 实际网格单元的结点坐标 */
  int d_num_dof;
  string d_object_name;
};
#endif
