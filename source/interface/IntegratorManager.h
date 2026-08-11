//
// 文件名:     IntegratorManager.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元积分器管理类
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_IntegratorManager
#define included_fem_IntegratorManager

#include "Pointer.h"
#include <string>
#include <vector>
#include <list>
using namespace std;

template <int DIM>
class BaseIntegrator;
/**
 * @brief 有限元数据片管理类 IntegratorManager, 为JAUMIN框架管理所有的积分器。
 *
 * 该类维护一个数组, 存储用户注册到JAUMIN框架的积分器。
 *
 * 该类提供如下函数, 支持用户从数组中获取指定名字的积分器。 - getIntegrator()
 * : 查询数据片
 *
 * @see adpt::Integrator
 */

template <int DIM>
class IntegratorManager {
public:
  /**
   * @brief 构造函数。
   */
  IntegratorManager();

  /**
   * @brief 虚析构函数。
   */
  virtual ~IntegratorManager();

  /**
   * @brief  获取积分管理器.
   *
   * @return 指针, 积分管理器.
   */
  static tbox::Pointer<IntegratorManager<DIM> > getManager();

  /**
   * @brief 将有积分器添加到链表中.
   *
   * @param integrator 输入参数, 指针, 指向积分器.
   */
  void addIntegrator(tbox::Pointer<BaseIntegrator<DIM> > integrator);

  /**
   * @brief 获取匹配于指定名称的积分器.
   * @param integrator_name 输入参数, 字符串, 积分器名字.
   * @return 指针, 指向匹配的积分器.
   */
  tbox::Pointer<BaseIntegrator<DIM> > getIntegrator(
      const string& integrator_name) const;

  /**
   * @brief 将该对象的信息输出到指定的输出流.
   * @param os 输入参数, 输出流类型, 指定输出流.
   */
  virtual void printClassData(ostream& os) const;

private:
  static tbox::Pointer<IntegratorManager<DIM> > s_integrator_manager_instance;
  /**
   * 存储积分器的数组.
   */
  std::vector<tbox::Pointer<BaseIntegrator<DIM> > > d_vec_integrator;
};

#endif
