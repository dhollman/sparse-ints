/*

 * binfiles.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#include "sparse_ints.h"

#ifndef BINFILES_H_
#define BINFILES_H_

void write_header(std::ofstream& o, const sc::Ref<sc::GaussianBasisSet>& basis, sparse_ints::SparseIntOptions opts);


#endif /* BINFILES_H_ */
