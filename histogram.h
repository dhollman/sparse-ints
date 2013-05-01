/*
 * histogram.h
 *
 *  Created on: Apr 29, 2013
 *      Author: dhollman
 */

#include "sparse_ints.h"
#include <math.h>

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

namespace sparse_ints {

class LogHistogram {

	double slot_spacing_;

	std::map<int, long> slots_;
	std::set<int> keys_;

public:

	LogHistogram(double slot_spacing = 1.0) :
		slot_spacing_(slot_spacing)
	{ };

	template <class T>
	void insert(T value){
		double lval, slot;
		if((double)value == 0.0) {
			slot = -INFINITY;
		}
		else {
			lval = log10(fabs((double)value));
			slot = floor(lval * 1.0 / slot_spacing_) * slot_spacing_;
		}
		keys_.insert(slot);
		++slots_[slot];
	}

	void insert_matrix(sc::RefSCMatrix mat){
		for(int r = 0; r < mat.nrow(); ++r){
			for(int c = 0; c < mat.ncol(); ++c){
				insert(mat.get_element(r, c));
			}
		}
	}

	void write(ofstream& o){
		int nslots = keys_.size();
		o.write((char*)&nslots, sizeof(int));
		for(std::set<int>::iterator it = keys_.begin(); it != keys_.end(); ++it){
			int bin = *it;
			o.write((char*)&bin, sizeof(int));
			long count = slots_[bin];
			assert(count!=0);
			o.write((char*)&count, sizeof(long));
		}
	}
};


} // end namespace sparse_ints

#endif /* HISTOGRAM_H_ */
