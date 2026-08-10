//
// 文件名:     BaseShapeFunction.h
// 软件包:     JAUMIN fem
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-02-16 12:55:07 +0800 (四, 2012-02-16) $
// 描述　:     有限元基函数基类。
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_BaseShapeFunction
#define included_fem_BaseShapeFunction

#include "Array.h"
#include "Pointer.h"
namespace JAUMIN {
namespace hier {
template <int DIM>
class DoubleVector;
}
}

using namespace std;
using namespace JAUMIN;

/*!
 * @brief 有限元模板单元基函数类，为有限元积分计算提供基函数在积分点上的函数值和
 * 梯度值。
 *
 */

template <int DIM>
class BaseShapeFunction {
public:
  /**
   * @brief 构造函数。
   *
   * @param name 输入参数, 字符串, 对象名称.
   */
  BaseShapeFunction(const string& name);

  /**
   * @brief 析构函数。
   */
  virtual ~BaseShapeFunction();

  /**
   * @brief 取出某一个积分点上的基函数值。
   *
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数，坐标，积分点坐标。
   *
   * @return 双精度型数组，积分点的基函数值。
   */
  virtual tbox::Array<double> value(
      tbox::Array<hier::DoubleVector<DIM> >& real_vertex,
      const hier::DoubleVector<DIM>& pnt) = 0;

  /**
   * @brief 取出某一个积分点上的某一个基函数值。
   *
   * @param j            输入参数, 整型, 积分点编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数，坐标，积分点坐标。
   *
   * @return 双精度型，积分点的基函数值。
   */
  virtual double value(int j,
                       tbox::Array<hier::DoubleVector<DIM> >& real_vertex,
                       const hier::DoubleVector<DIM>& pnt) = 0;

  /**
   * @brief 取出若干个积分点上的某一个基函数值。
   *
   * @param j            输入参数, 整型, 基函数编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts         输入参数，坐标数组，积分点坐标数组。
   *
   * @return 双精度型数组，积分点数组的基函数值。
   */
  virtual tbox::Array<double> value(
      int j, tbox::Array<hier::DoubleVector<DIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<DIM> >& pnts) = 0;

  /**
   * @brief 取出若干个积分点上的基函数值。
   *
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts         输入参数，坐标数组，积分点坐标数组。
   *
   * @return 双精度型，积分点数组的基函数值。
   */
  virtual tbox::Array<tbox::Array<double> > value(
      tbox::Array<hier::DoubleVector<DIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<DIM> >& pnts) = 0;

  /**
   * @brief 取出某一个积分点上的基函数梯度。
   *
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数，坐标，积分点坐标。
   *
   * @return 双精度型数组，积分点的基函数梯度。
   */
  virtual tbox::Array<tbox::Array<double> > gradient(
      tbox::Array<hier::DoubleVector<DIM> >& real_vertex,
      const hier::DoubleVector<DIM>& pnt) = 0;

  /**
   * @brief 取出某一个积分点上的某一个基函数梯度。
   *
   * @param j            输入参数, 整型, 基函数编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnt          输入参数，坐标，积分点坐标.
   *
   * @return 双精度型数组，积分点的基函数梯度。
   */
  virtual tbox::Array<double> gradient(
      int j, tbox::Array<hier::DoubleVector<DIM> >& real_vertex,
      const hier::DoubleVector<DIM>& pnt) = 0;

  /**
   * @brief 取出某若干积分点上的基函数梯度。
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts         输入参数，坐标，积分点坐标数组。
   *
   * @return 二维双精度型数组，积分点数组的基函数梯度。
   */
  virtual tbox::Array<tbox::Array<tbox::Array<double> > > gradient(
      tbox::Array<hier::DoubleVector<DIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<DIM> >& pnts) = 0;

  /**
   * @brief 取出某若干积分点上的某一个基函数梯度。
   *
   * @param j            输入参数, 整型, 基函数编号.
   * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
   * @param pnts         输入参数，坐标，积分点坐标数组。
   *
   * @return 二维双精度型数组，积分点数组的基函数梯度。
   */
  virtual tbox::Array<tbox::Array<double> > gradient(
      int j, tbox::Array<hier::DoubleVector<DIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<DIM> >& pnts) = 0;

  /**
   * @brief 判断给定的名字是否为该对象的名字.
   *
   * @param name 输入参数, 字符串, 约束名字.
   *
   * @return 布尔型, 是否一致.
   */
  virtual bool isMatching(const string& name);

  /**
   * @brief 获取对象的名字.
   *
   * @return 字符串, 对象的名字.
   */
  virtual const string& getName();

  /**
   * @brief 获取单元上自由度数目.
   *
   * @return 整型, 单元上自由度数目.
   */
  virtual int getNumberOfDof() = 0;

  /**
   * @brief 获取形函数在各个网格实体上自由度数目.
   *
   * @return 整型数组, 形函数在各个网格实体上自由度数目.
   */
  virtual tbox::Array<int> getNumberOfDofOnEntity();

private:
  string d_object_name;
};
#endif
