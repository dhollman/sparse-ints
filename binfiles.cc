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

#define BINFILE_VERSION 1

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
write_basis(
		ofstream& o,
		const Ref<GaussianBasisSet> basis
)
{
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

			// Write the number of primitives
			tmp = basis->shell(ish).nprimitive();
			o.write((char*)&tmp, sizeof(int));

			// Write the number of contractions in the shell
			//   just a placeholder since generally contracted basis sets
			//   are not supported yet
			tmp = basis->shell(ish).ncontraction();
			o.write((char*)&tmp, sizeof(int));

			// TODO write coefficients and exponents
		}

}

void
write_density_binfile(
		RefSymmSCMatrix P,
		string filename,
		Ref<GaussianBasisSet> basis
)
{
	ofstream o(filename.c_str(), ios::out | ios::binary);
	write_header_common(o);
    int8_t ty8 = (int8_t)DensityMatrix;
    o.write((char*)&ty8, sizeof(int8_t));

    // write a version number
    // First indicate that we're not writing the number of atoms
    int16_t new_file_marker = -1;
    o.write((char*)&new_file_marker, sizeof(int16_t));
    // Now write the version number
    int16_t version_number = BINFILE_VERSION;
    o.write((char*)&version_number, sizeof(int16_t));

    write_basis(o, basis);

    int nrc = P.n();
    value_t vals[nrc*(nrc+1)/2];
    int ispot = 0;
    for_each(row,nrc){
		for_each(col,row+1){
			vals[ispot] = (value_t)P.get_element(row, col);
			ispot++;
		}
    }
    assert(ispot == nrc*(nrc+1)/2);
    o.write((char*)vals, ispot*sizeof(value_t));
    o.close();
}

void
write_header(
	ofstream& o,
	const Ref<GaussianBasisSet>& bs1,
	const Ref<GaussianBasisSet>& bs2,
	const Ref<GaussianBasisSet>& bs3,
	const Ref<GaussianBasisSet>& bs4,
	const int ty
)
{
	write_header_common(o);
    int8_t ty8 = (int8_t)ty;
    o.write((char*)&ty8, sizeof(int8_t));

    // write a version number
    // First indicate that we're not writing the number of atoms
    int16_t new_file_marker = -1;
    o.write((char*)&new_file_marker, sizeof(int16_t));
    // Now write the version number
    int16_t version_number = BINFILE_VERSION;
    o.write((char*)&version_number, sizeof(int16_t));

    // 64-bit aligned at this point

    write_basis(o, bs1);
    write_basis(o, bs2);
    write_basis(o, bs3);
    write_basis(o, bs4);

    // TODO write molecule information so that results can be duplicated in the future
}

