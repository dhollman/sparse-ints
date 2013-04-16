/*

 * binfiles.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#include "sparse_ints.h"

#ifndef BINFILES_H_
#define BINFILES_H_



void write_header(
	std::ofstream& o,
	const sc::Ref<sc::GaussianBasisSet>& basis,
	sparse_ints::SparseIntOptions opts
);

void write_header(
	std::ofstream& o,
	const sc::Ref<sc::GaussianBasisSet>& bs1,
	const sc::Ref<sc::GaussianBasisSet>& bs2,
	const sc::Ref<sc::GaussianBasisSet>& bs3,
	const sc::Ref<sc::GaussianBasisSet>& bs4,
	const sparse_ints::binfile_type ty
);

void write_density_binfile(
	sc::RefSymmSCMatrix P,
	std::string filename,
	sc::Ref<sc::GaussianBasisSet> basis
);


#endif /* BINFILES_H_ */
