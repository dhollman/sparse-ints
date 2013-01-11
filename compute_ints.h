/*
 * compute_ints.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#ifndef COMPUTE_INTS_H_
#define COMPUTE_INTS_H_

#include "sparse_ints.h"

namespace sparse_ints{


void
compute_half_trans_ints(
        const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
        sc::TwoBodyOper::type* otypes, std::string* descs, int num_types,
        std::string prefix, std::string tmpdir,
        const sc::Ref<sc::GaussianBasisSet>& bs13,
        const sc::Ref<sc::GaussianBasisSet>& bs24,
    	DensityMap dens_pairs,
    	const sc::Ref<sc::SCMatrixKit>& inkit = (sc::Ref<sc::SCMatrixKit>)0
);

void
compute_full_trans_ints(
        const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
        sc::TwoBodyOper::type otype,
        std::string prefix, std::string tmpdir,
        const sc::Ref<sc::GaussianBasisSet>& bs1,
        const sc::Ref<sc::GaussianBasisSet>& bs2,
        const sc::Ref<sc::GaussianBasisSet>& bs3,
        const sc::Ref<sc::GaussianBasisSet>& bs4,
        sc::RefSymmSCMatrix& P1,
        sc::RefSymmSCMatrix& P2,
        sc::RefSymmSCMatrix& P3,
        sc::RefSymmSCMatrix& P4,
    	const sc::Ref<sc::SCMatrixKit>& inkit = (sc::Ref<sc::SCMatrixKit>)0
);

} // end namespace sparse_ints


#endif /* COMPUTE_INTS_H_ */
