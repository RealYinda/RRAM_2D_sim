//
// 文件名:     QuadratureInfo.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:44:03 2014 $
// 描述　:     有限元积分信息类
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_QuadratureInfo
#define included_fem_QuadratureInfo

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
 * @brief 单元积分信息类，辅助积分器类提供有限元积分计算。
 *
 */

template <int DIM>
class QuadratureInfo {
public:
  /**
   * @brief 构造函数。
   *
   * @param algebric_accuracy 输入参数，整型，积分精度。
   *
   * @note 参数缺省值为1，表示默认积分精度为一阶。
   */
  QuadratureInfo(int algebric_accuracy = 1);

  /**
   * @brief 析构函数。
   */
  virtual ~QuadratureInfo() = 0;

  /**
   * @brief 获取积分点。
   *
   *
   * @return 坐标数组，积分点数组。
   */
  virtual const tbox::Array<hier::DoubleVector<DIM> >&
  getQuadraturePoints() = 0;

  /**
   * @brief 获取积分点权重。
   *
   *
   * @return 双精度浮点数组，积分点权重数组。
   */
  virtual const tbox::Array<double>& getQuadratureWeights() = 0;

  /**
   * @brief 获取积分点数目。
   *
   *
   * @return 整型，积分点数目。
   */
  virtual int getNumberOfQuadraturePoints() = 0;

  /**
   * @brief 获取模板单元的面积(体积)。
   *
   *
   * @return 双精度浮点型，模板单元面积(体积)。
   */
  virtual double getElementVolume() = 0;

private:
};
#endif
