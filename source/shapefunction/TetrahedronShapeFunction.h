//
// 文件名:     TetrahedronShapeFunction.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:41:38 2014 $
// 描述　:     四面体有限元形函数类.
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_TetrahedronShapeFunction
#define included_TetrahedronShapeFunction

#include "Array.h"
#include "Pointer.h"
#include "BaseShapeFunction.h"
#include "DoubleVector.h"

using namespace JAUMIN;

/*!
 * @brief 有限元形四面体线性形函数类，为有限元计算中的单元积分计算提供形函数在积
 *        分点上的函数值和梯度值。
 *
 */

class TetrahedronShapeFunction : public BaseShapeFunction<NDIM> {
public:
  /**
   * @brief 构造函数。
   *
   * @param name 输入参数, 字符串, 形函数的名字.
   */
  TetrahedronShapeFunction(const string& name);

  /**
   * @brief 析构函数。
   */
  ~TetrahedronShapeFunction();

  /**
   * @brief 取出某一个积分点上的形函数值。
   *
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数，坐标，积分点坐标。
   *
   * @return 双精度型数组，积分点的形函数值。
   */
  tbox::Array<double> value(tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
                            const hier::DoubleVector<NDIM>& pnt);

  /**
   * @brief 取出某一个积分点上的某一个形函数值。
   *
   * @param j            输入参数, 整型, 积分点编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数，坐标，积分点坐标。
   *
   * @return 双精度型，积分点的形函数值。
   */
  double value(int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
               const hier::DoubleVector<NDIM>& pnt);

  /**
   * @brief 取出若干个积分点上的某一个形函数值。
   *
   * @param j            输入参数, 整型, 形函数编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts         输入参数，坐标数组，积分点坐标数组。
   *
   * @return 双精度型数组，积分点数组的形函数值。
   */
  tbox::Array<double> value(int j,
                            tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
                            const tbox::Array<hier::DoubleVector<NDIM> >& pnts);

  /**
   * @brief 取出若干个积分点上的形函数值。
   *
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts         输入参数，坐标数组，积分点坐标数组。
   *
   * @return 双精度型，积分点数组的形函数值。
   */
  tbox::Array<tbox::Array<double> > value(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<NDIM> >& pnts);

  /**
   * @brief 取出某一个积分点上的形函数梯度。
   *
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数，坐标，积分点坐标。
   *
   * @return 双精度型数组，积分点的形函数梯度。
   */
  tbox::Array<tbox::Array<double> > gradient(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const hier::DoubleVector<NDIM>& pnt);

  /**
   * @brief 取出某一个积分点上的某一个形函数梯度。
   *
   * @param j            输入参数, 整型,     形函数编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数, 坐标,     积分点坐标.
   *
   * @return 双精度型数组，积分点的形函数梯度。
   */
  tbox::Array<double> gradient(
      int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const hier::DoubleVector<NDIM>& pnt);

  /**
   * @brief 取出某若干积分点上的形函数梯度。
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts         输入参数，坐标，    积分点坐标数组。
   *
   * @return 二维双精度型数组，积分点数组的形函数梯度。
   */
  tbox::Array<tbox::Array<tbox::Array<double> > > gradient(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<NDIM> >& pnts);

  /**
   * @brief 取出某若干积分点上的某一个形函数梯度。
   *
   * @param j            输入参数, 整型,     形函数编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts         输入参数，坐标，    积分点坐标数组。
   *
   * @return 二维双精度型数组，积分点数组的形函数梯度。
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

  tbox::Array<int> getNumberOfDofOnEntity();

private:
  /**
   * @brief 计算四个点构成的四面体的体积.
   *
   * @param v0 输入参数, 双精度数组, 坐标.
   * @param v1 输入参数, 双精度数组, 坐标.
   * @param v2 输入参数, 双精度数组, 坐标.
   * @param v3 输入参数, 双精度数组, 坐标.
   *
   * @return 双精度浮点型, 四个点构成的四面体的体积.
   */
  double get_volume(const hier::DoubleVector<NDIM>& v0,
                    const hier::DoubleVector<NDIM>& v1,
                    const hier::DoubleVector<NDIM>& v2,
                    const hier::DoubleVector<NDIM>& v3);

  tbox::Array<int> d_num_dof_on_entity;
  int d_num_dof; /**< 该形函数上自由度的数目 */
  tbox::Array<hier::DoubleVector<NDIM> >
      d_real_vertex; /**< 实际网格单元的结点坐标 */
};
#endif
