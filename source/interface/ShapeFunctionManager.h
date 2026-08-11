//
// 文件名:     ShapeFunctionManager.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元形函数管理类
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_ShapeFunctionManager
#define included_fem_ShapeFunctionManager

#include "Pointer.h"
#include "Array.h"
#include <string>
#include <vector>
#include <list>
using namespace std;

using namespace JAUMIN;
template <int DIM>
class BaseShapeFunction;

/**
 * @brief 有限元形函数管理类 ShapeFunctionManager,
 *为JAUMIN框架管理所有的形函数。
 *
 * 该类维护一个数组, 存储用户注册到JAUMIN框架的形函数。
 *
 * 该类提供如下函数, 支持用户从数组中获取指定名字的形函数。 - getShapeFunction()
 * : 查询数据片
 *
 * @see adpt::BaseShapeFunction
 */

template <int DIM>
class ShapeFunctionManager {
public:
  /**
   * @brief 构造函数。
   */
  ShapeFunctionManager();

  /**
   * @brief 虚析构函数。
   */
  ~ShapeFunctionManager();

  /**
   * @brief  获取形函数理器.
   *
   * @return 指针, 形函数管理器.
   */
  static tbox::Pointer<ShapeFunctionManager<DIM> > getManager();

  /**
   * @brief 将有形函数注册到链表中。
   *
   * @param shape_function 输入参数, 指针, 指向形函数.
   */
  void addShapeFunction(tbox::Pointer<BaseShapeFunction<DIM> > Shape_function);

  /**
   * @brief 获取匹配于指定名称的形函数.
   * @param name 输入参数, 字符串, 形函数名字.
   * @return 指针, 指向匹配的形函数.
   */
  tbox::Pointer<BaseShapeFunction<DIM> > getShapeFunction(
      const string& name) const;

  /**
   * @brief 将该对象的信息输出到指定的输出流.
   * @param os 输入参数, 输出流类型, 指定输出流.
   */
  void printClassData(ostream& os) const;

private:
  static tbox::Pointer<ShapeFunctionManager<DIM> >
      s_shape_func_manager_instance;

  /**
   * 存储积分器的数组。
   */
  std::vector<tbox::Pointer<BaseShapeFunction<DIM> > > d_vec_shape_function;
};

#endif
