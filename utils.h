/*
 * utils.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#include "sparse_ints.h"
#include <stdlib.h>

#ifndef UTILS_H_
#define UTILS_H_

using namespace std;

//template <class T=value_t>
//inline T
//max_abs(sc::RefSCMatrix& mat){
//	int nrows = mat->nrow();
//	int ncols = mat->ncol();
//	T maxval = 0.0;
//	for_each(row,nrows, col,ncols){
//		T tmpval = (T)fabs(mat->get_element(row, col));
//		if(tmpval > maxval){
//			maxval = tmpval;
//		}
//	}
//	return maxval;
//}
//
//template <class T=value_t>
//inline T
//mean(sc::RefSCMatrix& mat){
//	int nrows = mat->nrow();
//	int ncols = mat->ncol();
//	T total = 0.0;
//	for_each(row,nrows, col,ncols){
//		T tmpval = (T)fabs(mat->get_element(row, col));
//		total += tmpval;
//	}
//	return total/((T)(nrows*ncols));
//}
//
//inline int compare(const void * a, const void * b)
//{
//  return ( *(int*)a - *(int*)b );
//}
//
//
//template <class T=value_t>
//inline T
//median(sc::RefSCMatrix& mat){
//	int nrows = mat->nrow();
//	int ncols = mat->ncol();
//	double *data = new double[nrows*ncols];
//	for_each(idata, nrows*ncols){
//		data[idata] = fabs(data[idata]);
//	}
//	mat->convert(data);
//	qsort(data, nrows*ncols, sizeof(double), compare);
//	T rv = (T)data[nrows*ncols/2];
//	delete[] data;
//	return rv;
//}
//
//template <class T=value_t, class index_type=idx_t>
//T
//max_abs(sc::RefSCMatrix& mat, index_type &maxrow, index_type &maxcol){
//	int nrows = mat->nrow();
//	int ncols = mat->ncol();
//	T maxval = 0.0;
//	for_each(row,nrows, col,ncols){
//		T tmpval = (T)fabs(mat->get_element(row, col));
//		if(tmpval > maxval){
//			maxrow = (index_type)row;
//			maxcol = (index_type)col;
//			maxval = tmpval;
//		}
//	}
//	return maxval;
//}

inline off_t
file_size(string file){
	ifstream i(file.c_str(), ios::in | ios::binary | ios::ate);
	off_t length = (off_t)i.tellg();
	i.close();
	return length;
}

inline string
memory_size_string(double length){
	stringstream sstr;
	if(length < 8.e2)
		sstr << length << " B";
	else if(length < 8.e5)
		sstr << setprecision(2) << fixed << length/1024.0 << " KB";
	else if(length < 8.e8)
		sstr << setprecision(2) << fixed << length/1024.0/1024.0 << " MB";
	else
		sstr << setprecision(2) << fixed << length/1024.0/1024.0/1024.0 << " GB";
	string str = sstr.str();
	return str;
}

inline string
memory_size_string(int length){
	return memory_size_string(double(length));
}

inline string
file_size_string(string file){
	double length = (double)file_size(file);
	return memory_size_string(length);
}

inline double
fake_int(int fake_num, int idx1=0, int idx2=0, int idx3=0, int idx4=0){
	switch(fake_num){
	case 1:
		return 1.0;
	default:
		return 0.0;
	}
}

inline void
copy_buffer(ifstream& i, ofstream& o){

	assert(i.is_open());
	int timeout = 10;
	i.seekg(0, ios::end);
	long length = i.tellg();
	int count = 0;
	while(length < 0 && count < timeout) {
		count++;
		cout << "  Failed to get size of file...retrying..." << endl;
		sleep(1);
		i.seekg(0, ios::end);
		length = i.tellg();
	}
	assert(length >= 0);
	i.seekg(0, ios::beg);
	char* buffer = new char[BUFF_SIZE];
	while(length > BUFF_SIZE){
		i.read(buffer, BUFF_SIZE);
		o.write(buffer, BUFF_SIZE);
		length -= BUFF_SIZE;
	}
	i.read(buffer, length);
	o.write(buffer, length);
	delete[] buffer;

}


#endif /* UTILS_H_ */
