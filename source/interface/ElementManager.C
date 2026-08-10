//
// 文件名:     ElementManager.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元单元管理器实现.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "JAUMIN_config.h"
#include "Pointer.h"
#include "Utilities.h"
#include <string>
#include <list>
#include "ElementManager.h"
#include "BaseElement.h"
#ifndef NULL
#define NULL (0)
#endif

template <int DIM>
tbox::Pointer<ElementManager<DIM> >
    ElementManager<DIM>::s_element_manager_instance = 0;

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for ElementManager objects.      *
*                                                                       *
*************************************************************************
*/

template <int DIM>
ElementManager<DIM>::ElementManager() {}

template <int DIM>
ElementManager<DIM>::~ElementManager() {}

template <int DIM>
tbox::Pointer<ElementManager<DIM> > ElementManager<DIM>::getManager() {
  if (s_element_manager_instance.isNull()) {
    s_element_manager_instance = new ElementManager<DIM>();
  }
  return (s_element_manager_instance);
}

/************************************************************************
*                                                                       *
* Add element to list.                                                     *
*                                                                       *
*************************************************************************/
template <int DIM>
void ElementManager<DIM>::addElement(tbox::Pointer<BaseElement<DIM> > ele) {
  const string &name = ele->getName();
  bool exist = false;
  unsigned int size = d_vec_element.size();
  for (unsigned int i = 0; i < size; ++i) {
    exist = d_vec_element[i]->isMatching(name);
    if (exist) break;
  }
  if (!exist)
    d_vec_element.push_back(ele);
  else
    TBOX_ERROR("element \"" << name
                            << "\" is existed. please choose another name!!"
                            << endl);
}

/************************************************************************
*                                                                       *
* Search element lists for element matching request.                          *
*                                                                       *
*************************************************************************/

template <int DIM>
tbox::Pointer<BaseElement<DIM> > ElementManager<DIM>::getElement(
    const string &ele_name) const {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!ele_name.empty());
#endif

  tbox::Pointer<BaseElement<DIM> > element = NULL;
  bool found = false;
  int vec_size = d_vec_element.size();
  for (int i = 0; i < vec_size; ++i) {
    tbox::Pointer<BaseElement<DIM> > tmp_element = d_vec_element[i];
    if (tmp_element->isMatching(ele_name)) {
      element = tmp_element;
      found = true;
    }
  }

  if (!found) {
    TBOX_ERROR("fem element \"" << ele_name << "\" is not found. " << endl);
  }

  if (element.isNull()) {
    TBOX_ERROR("fem element \"" << ele_name << "\" is NULL. " << endl);
  }

  return (element);
}

template <int DIM>
void ElementManager<DIM>::setElement(tbox::Array<string> &element_names,
                                     tbox::Array<int> &element_marks) {
  d_element_names = element_names;
  d_element_marks = element_marks;
}

/************************************************************************
 * 获取形函数名字数组.
*************************************************************************/
template <int DIM>
tbox::Array<string> ElementManager<DIM>::getElementNames() {
  return d_element_names;
}

/************************************************************************
 * 获取形函数的单元标识数组.
*************************************************************************/
template <int DIM>
tbox::Array<int> ElementManager<DIM>::getElementMarks() {
  return d_element_marks;
}

template class ElementManager<NDIM>;
