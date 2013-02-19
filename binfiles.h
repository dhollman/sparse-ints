/*

 * binfiles.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#include "sparse_ints.h"

#ifndef BINFILES_H_
#define BINFILES_H_


typedef enum binfile_type {
	MaxAbs = 2,
	Average = 3,
	StdDev = 4,
	Median = 5
} binfile_type;

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
	const binfile_type ty
);


#endif /* BINFILES_H_ */
