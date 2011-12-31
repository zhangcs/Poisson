/**
 *       @file  rcm.hpp
 *      @brief  Reverse Cuthill Mckee for csr matrix
 *
 *     @author  Feiteng Huang
 *
 *   @internal
 *     created  2011年12月23日 08时46分38秒
 *    compiler  gcc/g++
 *       email  hftenger@126.com
 * institution  Sichuan University
 *
 * =====================================================================================
 */
#ifndef _RCM_HPP_
#define _RCM_HPP_

#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

int RCM( fsls_CSRMatrix *A, fsls_XVector *f , fsls_XVector *u, int nt);

#endif
