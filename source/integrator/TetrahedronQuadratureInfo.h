//
// 文件名:     TetrahedronQuadratureInfo.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:42:22 2014 $
// 描述　:     四面体有限元积分信息类.
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_TetrahedronQuadratureInfo
#define included_fem_TetrahedronQuadratureInfo

#include "Array.h"
#include "Pointer.h"
#include "QuadratureInfo.h"
#include "DoubleVector.h"

using namespace JAUMIN;

/*!
 * @brief 有限元积分信息类，辅助积分器类提供有限元积分积分计算需要的信息。该类根
 *        据用户提供的积分精度，在有限元计算中提供积分点以及积分点权重信息。
 *
 */

class TetrahedronQuadratureInfo : public QuadratureInfo<NDIM> {
public:
  /**
   * @brief 构造函数。
   *
   * @param algebric_accuracy 输入参数，整型，积分精度。
   *
   * @note 参数缺省值为1，表示默认积分精度为一阶。
   */
  TetrahedronQuadratureInfo(int algebric_accuracy = 1);

  /**
   * @brief 析构函数。
   */
  ~TetrahedronQuadratureInfo();

  /**
   * @brief 获取积分点。
   *
   *
   * @return 坐标数组，积分点数组。
   */
  const tbox::Array<hier::DoubleVector<NDIM> >& getQuadraturePoints();

  /**
   * @brief 获取积分点权重。
   *
   *
   * @return 双精度浮点数组，积分点权重数组。
   */
  const tbox::Array<double>& getQuadratureWeights();

  /**
   * @brief 获取积分点数目。
   *
   *
   * @return 整型，积分点数目。
   */
  int getNumberOfQuadraturePoints();

  /**
   * @brief 获取模板单元的面积(体积)。
   *
   *
   * @return 双精度浮点型，模板单元面积(体积)。
   */
  double getElementVolume();

private:
  int d_algebric_accuracy;  /**< 积分精度 */
  int d_number_quad_points; /**< 积分点个数 */
  double d_volume;          /**< 模板单元的面积(体积) */
  tbox::Array<hier::DoubleVector<NDIM> > d_quad_points; /**< 积分点数组 */
  tbox::Array<double> d_quad_weights; /**< 积分点权重数组 */
};
#endif
