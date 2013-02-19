/*
 * binfiles.cc
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#include "binfiles.h"

using namespace std;
using namespace sc;
using namespace sparse_ints;

void
write_header_common(
	ofstream& o
){
	// Write the int_size, the double size, endianness, and number of atoms
    int8_t big_endian;
    if ( htonl(47) == 47 ) {
      big_endian = (int8_t)1;
    } else {
      big_endian = (int8_t)0;
    }
    o.write((char*)&big_endian, sizeof(int8_t));
    int8_t int_size = (int8_t)sizeof(int);
    o.write((char*)&int_size, sizeof(int8_t));
    int8_t value_size;
    value_size = (int8_t)sizeof(value_t);
	o.write((char*)&value_size, sizeof(int8_t));
}

void
write_header(ofstream& o, const Ref<GaussianBasisSet>& basis, SparseIntOptions opts)
{
	write_header_common(o);
    int8_t max8 = (int8_t)opts.max_only;
    o.write((char*)&max8, sizeof(int8_t));

    // Write a description of the basis
    int natoms = basis->ncenter();
    o.write((char*)&natoms, sizeof(int));
    int nshell = basis->nshell();
    o.write((char*)&nshell, sizeof(int));
    // note that at this point, we are 64-bit aligned

    // Write the info about each shell
    for_each(ish, nshell){
    	GaussianShell sh = basis->shell(ish);
    	// Write the center number
    	int tmp = basis->shell_to_center(ish);
    	o.write((char*)&tmp, sizeof(int));
    	// Write the number of functions
    	tmp = sh.nfunction();
    	o.write((char*)&tmp, sizeof(int));
    	// Write the maximum and minimum angular momentum
    	tmp = sh.max_am();
    	o.write((char*)&tmp, sizeof(int));
    	tmp = sh.min_am();
    	o.write((char*)&tmp, sizeof(int));
    	// TODO Write some sort of "average contraction exponent" here...
    }
}


void
write_header(
	ofstream& o,
	const Ref<GaussianBasisSet>& bs1,
	const Ref<GaussianBasisSet>& bs2,
	const Ref<GaussianBasisSet>& bs3,
	const Ref<GaussianBasisSet>& bs4,
	const binfile_type ty
)
{
	write_header_common(o);
    int8_t ty8 = (int8_t)ty;
    o.write((char*)&ty8, sizeof(int8_t));

    vector<Ref<GaussianBasisSet> > sets;
    sets.push_back(bs1); sets.push_back(bs2);
    sets.push_back(bs3); sets.push_back(bs4);
    for_each(ibas, 4){
    	Ref<GaussianBasisSet> basis = sets[ibas];
		// Write a description of the basis
		int natoms = basis->ncenter();
		o.write((char*)&natoms, sizeof(int));
		int nshell = basis->nshell();
		o.write((char*)&nshell, sizeof(int));
		// note that at this point, we are 64-bit aligned

		// Write the info about each shell
		for_each(ish, nshell){
			GaussianShell sh = bs1->shell(ish);
			// Write the center number
			int tmp = basis->shell_to_center(ish);
			o.write((char*)&tmp, sizeof(int));
			// Write the number of functions
			tmp = sh.nfunction();
			o.write((char*)&tmp, sizeof(int));
			// Write the maximum and minimum angular momentum
			tmp = sh.max_am();
			o.write((char*)&tmp, sizeof(int));
			tmp = sh.min_am();
			o.write((char*)&tmp, sizeof(int));
			// TODO Write some sort of "average contraction exponent" here...
		}
    }
}

