/*
 * compute_ints.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#ifndef COMPUTE_INTS_H_
#define COMPUTE_INTS_H_

#include "sparse_ints.h"
#include <math/scmat/local.h>

namespace sparse_ints{


void
compute_full_trans_ints(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
        const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
        sc::TwoBodyOper::type otype,
        std::string prefix, std::string matname, std::string tmpdir,
        const sc::Ref<sc::GaussianBasisSet>& bs1,
        const sc::Ref<sc::GaussianBasisSet>& bs2,
        const sc::Ref<sc::GaussianBasisSet>& bs3,
        const sc::Ref<sc::GaussianBasisSet>& bs4,
        sc::RefSymmSCMatrix& P1,
        sc::RefSymmSCMatrix& P2,
        sc::RefSymmSCMatrix& P3,
        sc::RefSymmSCMatrix& P4,
    	const sc::Ref<sc::LocalSCMatrixKit>& inkit = (sc::Ref<sc::LocalSCMatrixKit>)0
);

} // end namespace sparse_ints


#endif /* COMPUTE_INTS_H_ */
