//
// 文件名:     BaseIntegrator.h
// 软件包:     JAUMIN fem
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-06-17 16:43:13 +0800 (日, 2012-06-17) $
// 描述　:     积分器基类
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_BaseIntegrator
#define included_fem_BaseIntegrator

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
 * @brief 积分器基类，提供单元上积分功能. 各种形状单元积分器继承该基类.
 *
 */

template <int DIM>
class BaseIntegrator {
public:
  /**
   * @brief 构造函数。
   *
   * @param name 输入参数, 字符串, 对象名称.
   */
  BaseIntegrator(const std::string &name);

  /**
   * @brief 析构函数。
   */
  virtual ~BaseIntegrator();

  /**
   * @brief 获取积分点。
   *
   * @param real_vertex 输入参数, 坐标数组, 单元的结点坐标.
   *
   * @return 坐标数组，积分点数组。
   */
  virtual tbox::Array<hier::DoubleVector<DIM> > getQuadraturePoints(
      tbox::Array<hier::DoubleVector<DIM> > &real_vertex) = 0;

  /**
   * @brief 获取积分点权重。
   *
   * @return 双精度浮点型数组，积分点权重数组。
   */
  virtual const tbox::Array<double> &getQuadratureWeights() = 0;

  /**
   * @brief 获取积分点数目。
   *
   * @return 整型，积分点数目。
   */
  virtual int getNumberOfQuadraturePoints() = 0;

  /**
   * @brief 获取单元结点数目。
   *
   * @return 整型，单元结点数目。
   */
  virtual int getNumberOfVertex() = 0;

  /**
   * @brief 获取单元的面积(体积)。
   *
   * @return 双精度浮点型，单元面积(体积)。
   */
  virtual double getElementVolume() = 0;

  /**
   * @brief 获取由局部到全局仿射变换的矩阵行列式。
   *
   * @param real_vertex 输入参数, 坐标数组, 单元的结点坐标.
   *
   * @return 双精度浮点型，仿射变换的矩阵行列式。
   */
  virtual double getLocal2GlobalJacobian(
      tbox::Array<hier::DoubleVector<DIM> > &real_vertex) = 0;

  /**
   * @brief 获取由全局到局部仿射变换的矩阵行列式。
   *
   * @param real_vertex 输入参数, 坐标数组, 单元的结点坐标.
   *
   * @return 双精度浮点型，仿射变换的矩阵行列式。
   */
  virtual double getGlobal2LocalJacobian(
      tbox::Array<hier::DoubleVector<DIM> > &real_vertex) = 0;

  /**
   * @brief 判断给定的名字是否与本对象的匹配.
   *
   * @param name  输入参数, 字符串, 对象名字.
   *
   * @return 布尔型, 是否匹配.
   */
  virtual bool isMatching(const std::string &ele_name);

  /**
   * @brief 获取对象的名字.
   *
   * @return 字符串, 对象的名字.
   */
  virtual const std::string &getName();

private:
  std::string d_object_name;
};
#endif
