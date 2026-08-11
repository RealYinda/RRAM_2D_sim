//
// 文件名:     TriangleQuadratureInfo.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2012-02-16 12:55:07 +0800 (四, 2012-02-16) $
// 描述　:     有限元积分信息类。
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_TriangleQuadratureInfo
#define included_fem_TriangleQuadratureInfo

#include "Array.h"
#include "Pointer.h"
#include "QuadratureInfo.h"

namespace JAUMIN {
namespace hier {
template <int DIM>
class DoubleVector;
}
}

/*!
 * @brief 有限元模板单元积分信息类，辅助模板单元类提供有限元积分计算需要的信息。
 * 该类根据用户提供的积分精度，在有限元计算中提供积分点以及积分点权重信息。
 *
 */

class TriangleQuadratureInfo : public QuadratureInfo<NDIM> {
public:
  /**
   * @brief 构造函数。
   *
   * @param algebric_accuracy 输入参数，整型，积分精度。
   *
   * @note 参数缺省值为1，表示默认积分精度为一阶。
   */
  TriangleQuadratureInfo(int algebric_accuracy = 1);

  /**
   * @brief 析构函数。
   */
  ~TriangleQuadratureInfo();

  /**
   * @brief 获取积分点。
   *
   *
   * @return 坐标数组，积分点数组。
   */
  virtual const tbox::Array<hier::DoubleVector<NDIM> >& getQuadraturePoints();

  /**
   * @brief 获取积分点权重。
   *
   *
   * @return 双精度浮点数组，积分点权重数组。
   */
  virtual const tbox::Array<double>& getQuadratureWeights();

  /**
   * @brief 获取积分点数目。
   *
   *
   * @return 整型，积分点数目。
   */
  virtual int getNumberOfQuadraturePoints();

  /**
   * @brief 获取模板单元的面积(体积)。
   *
   *
   * @return 双精度浮点型，模板单元面积(体积)。
   */
  virtual double getElementVolume();

private:
  int d_algebric_accuracy;  /**< 积分精度 */
  int d_number_quad_points; /**< 积分点个数 */
  double d_volume;          /**< 模板单元的面积(体积) */
  tbox::Array<hier::DoubleVector<NDIM> > d_quad_points; /**< 积分点数组 */
  tbox::Array<double> d_quad_weights; /**< 积分点权重数组 */
};
#endif
