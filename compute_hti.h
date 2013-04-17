/*
 * compute_hti.h
 *
 *  Created on: Apr 11, 2013
 *      Author: dhollman
 */

#ifndef COMPUTE_HTI_H_
#define COMPUTE_HTI_H_

#include "sparse_ints.h"
#include <math/scmat/local.h>

namespace sparse_ints{

void
compute_hti_threaded(
		sc::Ref<sc::MessageGrp> msg,
		sc::Ref<sc::ThreadGrp> thr,
        const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
        sc::TwoBodyOper::type* otypes, std::string* descs, int num_types,
        std::string prefix, std::string tmpdir,
        const sc::Ref<sc::GaussianBasisSet>& bs13,
        const sc::Ref<sc::GaussianBasisSet>& bs24,
    	DensityMap dens_pairs,
    	sc::Ref<sc::LocalSCMatrixKit>& kit
);

} // end namespace sparse_ints


#endif /* COMPUTE_HTI_H_ */
