/*
 * compute_untrans.h
 *
 *  Created on: Apr 11, 2013
 *      Author: dhollman
 */

#ifndef COMPUTE_UNTRANS_H_
#define COMPUTE_UNTRANS_H_

#include "sparse_ints.h"
#include <math/scmat/local.h>

namespace sparse_ints{

void
compute_untrans_threaded(
		sc::Ref<sc::MessageGrp> msg,
		sc::Ref<sc::ThreadGrp> thr,
        const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
        sc::TwoBodyOper::type* otypes, std::string* descs, int num_types,
        std::string prefix, std::string tmpdir,
        const sc::Ref<sc::GaussianBasisSet>& bs1,
        const sc::Ref<sc::GaussianBasisSet>& bs2,
        const sc::Ref<sc::GaussianBasisSet>& bs3,
        const sc::Ref<sc::GaussianBasisSet>& bs4,
    	sc::Ref<sc::LocalSCMatrixKit>& kit
);

} // end namespace sparse_ints


#endif /* COMPUTE_UNTRANS_H_ */
