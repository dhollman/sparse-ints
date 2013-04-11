from __future__ import division
from operator import attrgetter
from struct import unpack
import os
import sys
from os.path import exists, getmtime, getsize
import numpy as np
from math import log10

def is_newer(filea, fileb):
    return getmtime(filea) - getmtime(fileb) > 0

def make_matrix(shape, mm_filename, bin_filename=None):
    shape = tuple(shape)
    if exists(mm_filename) and (bin_filename is None or is_newer(mm_filename, bin_filename)):
        rv = np.memmap(mm_filename, mode="r", shape=shape, dtype='float32')
        return rv, True
    else:
        rv = np.memmap(mm_filename, mode="w+", shape=shape, dtype='float32')
        return rv, False

def make_mag_matrix(shape, mags_filename, mm_filename=None, force_new = False):
    shape = tuple(shape)
    if exists(mags_filename) and is_newer(mags_filename, mm_filename) and not force_new:
        rv = np.memmap(mags_filename, mode="r", shape=shape, dtype='float32')
        return rv, True
    else:
        if exists(mags_filename):
            print "    Removing outdated memmap mags file {}".format(mags_filename)
            os.unlink(mags_filename)
        rv = np.memmap(mags_filename, mode="w+", shape=shape, dtype='float32')
        return rv, False


class ProgressBar(object):
    def __init__(self, max_val, width=80, out=sys.stdout, fill_char='=', indent='', arrow = True):
        self.max_val = max_val
        self.current_percent = 0
        self.width = width
        self.out = out
        self.fill_char = fill_char
        self.indent = indent
        self.arrow = arrow
    def update(self, new_val):
        new_percent = int(float(new_val/self.max_val)*100)
        if new_percent != self.current_percent:
            self.current_percent = new_percent
            self.print_bar()
    def print_bar(self):
        bar_width = self.width - 7
        filled = int(round(float(self.current_percent)/100.0 * bar_width))
        filler = self.fill_char * filled
        if self.arrow and filled > 0 and filled < bar_width:
            filler = filler[:-1] + '>'
        space = " " * (bar_width - filled)
        self.out.write("\r{3}[{0}{1}] {2:>3d}%".format(filler, space, self.current_percent, self.indent))
        self.out.flush()

class BasisSet(object):

    def __init__(self, shells):
        # Sort is guarenteed to be stable
        self.shells = shells
        if not shells == sorted(shells, key=attrgetter('center')):
            raise NotImplementedError()
        self.center_to_function = []
        self.center_to_shell = []
        self.shell_to_function = []
        prev_center = None
        curr_function = 0
        for ish, sh in enumerate(shells):
            if sh.center != prev_center:
                self.center_to_function.append(curr_function)
                self.center_to_shell.append(ish)
                prev_center = sh.center
            self.shell_to_function.append(curr_function)
            sh.first_bf = curr_function
            curr_function += sh.nfunction

    @property
    def nbf(self):
        return sum(s.nfunction for s in self.shells)

    @property
    def nshell(self):
        return len(self.shells)

    @property
    def ncenter(self):
        return len(self.center_to_function)

    def __eq__(self, other):
        return self.shells == other.shells

class Shell(object):

    def __init__(self, idx, center, nbf, maxam, minam):
        self.index = idx
        self.center = center
        self.nfunction = nbf
        self.min_am = minam
        self.max_am = maxam

    def __str__(self):
        return "Shell #{} on center {}:  nfunction={}, min_am={}, max_am={}".format(
            self.index, self.center, self.nfunction, self.min_am, self.max_am
        )

    @property
    def slice(self):
        return slice(self.first_bf, self.first_bf + self.nfunction)

    def __eq__(self, other):
        def signature(s):
            return s.index, s.center, s.nfunction, s.min_am, s.max_am
        return signature(self) == signature(other)


class BasisFunctionTensor(object):

    def __init__(self, basis_sets, int_file_obj):
        self.basis_sets = basis_sets
        self.shape = tuple(b.nbf for b in self.basis_sets)
        self.int_file_obj = int_file_obj
        self.array, self.matrix_filled = make_matrix(
            self.shape,
            int_file_obj.memmap_filename,
            bin_filename=int_file_obj.filename
        )
        self.mag_matrix, self.mags_filled = None, False

    def init_mag_matrix(self):
        self.mag_array, self.mags_filled = make_mag_matrix(
            self.shape,
            self.int_file_obj.mags_memmap_filename,
            self.int_file_obj.memmap_filename,
            force_new=not self.matrix_filled
        )

    @property
    def matrix(self):
        return self.array
    @matrix.setter
    def matrix(self,newval):
        self.array = newval

    @property
    def mag_matrix(self):
        return self.mag_array
    @mag_matrix.setter
    def mag_matrix(self,newval):
        self.mag_array = newval

class BasisFunctionMatrix(BasisFunctionTensor):

    def __init__(self, bs1, bs3, int_file_obj):
        super(BasisFunctionMatrix, self).__init__((bs1,bs3), int_file_obj)

    @property
    def nrow(self):
        return self.shape[0]

    @property
    def ncol(self):
        return self.shape[1]

    def draw_lines_string(self, image_width):
        cum_shell_rows = [
            sum(s.nfunction for s in self.basis_sets[0].shells if s.center < i)
                for i in xrange(self.basis_sets[0].ncenter)
        ]
        cum_shell_cols = [
            sum(s.nfunction for s in self.basis_sets[1].shells if s.center < i)
                for i in xrange(self.basis_sets[1].ncenter)
        ]
        image_height = image_width * self.ncol//self.nrow
        return " ".join("M 0,{ncol} L {th},{ncol} M {nrow},0 L {nrow},{tw}".format(
                nrow=int(c*float(image_height)/self.nrow),
                ncol=int(d*float(image_width)/self.ncol),
                tw=image_width, th=image_height
            ) for c,d in zip(cum_shell_rows[1:],cum_shell_cols[1:])
        )

class LoadedIntFile(object):

    def __init__(self, filename, force_reload=False):
        self.filename = filename
        self.force_reload = force_reload
        self._load()

    @property
    def memmap_filename(self):
        if self.filename[-4:] == ".bin":
            return self.filename[:-4] + ".npmm"
        else:
            return self.filename + ".npmm"

    @property
    def mags_memmap_filename(self):
        return self.memmap_filename + "_mags"

    @property
    def loaded_array(self):
        return self.loaded_matrix
    @loaded_array.setter
    def loaded_array(self, val):
        self.loaded_matrix = val

    def delete_mags(self):
        del self.loaded_matrix.mag_array
        self.loaded_matrix.mag_array = None
        self.loaded_matrix.mags_filled = False
        os.unlink(self.mags_memmap_filename)

    def load_mags(self, mag_function_or_min_val, max_val=None, force_reload=False):
        self.loaded_array.init_mag_matrix()
        if self.loaded_array.mags_filled and not force_reload:
            return
        else:
            if callable(mag_function_or_min_val):
                mag = mag_function_or_min_val
            else:
                min_val = mag_function_or_min_val
                def mag(val):
                    if min_val is not None and abs(val) < pow(10.0, -min_val):
                        return 0
                    elif max_val is not None and abs(val) > pow(10.0, -max_val):
                        return -min_val-max_val
                    return log10(abs(val)) + min_val
            mag_vect = np.vectorize(mag)
            self.loaded_matrix.mag_array = mag_vect(self.loaded_matrix.array)
            self.loaded_matrix.mags_filled = True



    def _load(self):
        print "  Loading matrix from {}".format(self.filename)
        with open(self.filename, 'rb') as f:
            #========================================#
            # first load all of the metadata
            # Get sizes and stuff
            big_endian = unpack("@b", f.read(1))[0]
            endian_char = ">" if big_endian else "<"
            # a convenience function factory
            def getter(n, ch, size):
                val = f.read(n*size)
                if len(val) == 0:
                    raise EOFError
                elif len(val) != n*size:
                    raise IOError("File {} not aligned".format(self.filename))
                return unpack(endian_char+str(n)+ch, val)
                # byte and unsigned short functions
            self.byte_ = byte_ = lambda x: getter(x, 'b', 1)
            self.short_ = short_ = lambda x: getter(x, 'H', 2)
            # Read the int_size, float_size, and the type of data file this is
            int_size, float_size, ints_type = byte_(3)
            self.ints_type = ints_type
            # create convenience functions
            if int_size not in (4, 8): raise ValueError("int size of {} not acceptable".format(int_size))
            int_char = 'i' if int_size == 4 else 'q'
            self.int_ = int_ = lambda x: getter(x, int_char, int_size)
            if float_size not in (4, 8): raise ValueError("float size of {} not acceptable".format(float_size))
            float_char = "f" if float_size == 4 else "d"
            self.float_ = float_ = lambda x: getter(x, float_char, float_size)
            #========================================#
            if ints_type == 0:
                # File contains all integrals
                self._load_basis(f)
                self._load_all(f)
            #========================================#
            elif ints_type == 1:
                # Legacy support for bs1==bs2==bs3==bs4 max only files
                self._load_basis_legacy(f)
                self._load_maxes(f)
            #========================================#
            else:
                # File only contains maxima/average/std_dev/whatever
                self._load_basis(f)
                self._load_maxes(f)

    def _load_basis_legacy(self, f):
        # Legacy support for bs1==bs2==bs3==bs4 max only files
        # read in the shell descriptions
        int_ = self.int_
        natoms, nshell = int_(2)
        shells = []
        tot_bf = 0
        for ish in xrange(nshell):
            center, nbf, max_am, min_am = int_(4)
            sh = Shell(ish, center, nbf, max_am, min_am)
            sh.first_bf = tot_bf
            tot_bf += nbf
            shells.append(sh)
        self.basis_sets = [BasisSet(shells)] * 4

    def _load_basis(self, f):
        int_ = self.int_
        self.basis_sets = []
        for __ in xrange(4):
            natoms, nshell = int_(2)
            shells = []
            tot_bf = 0
            for ish in xrange(nshell):
                center, nbf, max_am, min_am = int_(4)
                sh = Shell(ish, center, nbf, max_am, min_am)
                sh.first_bf = tot_bf
                tot_bf += nbf
                sh.nprimitive = int_(1)
                sh.ncontraction = int_(1)
                shells.append(sh)
            self.basis_sets.append(BasisSet(shells))

    def _load_maxes(self, f):
        metadata_size = f.tell()
        tot_size = getsize(self.filename)
        data_size = tot_size - metadata_size
        float_, short_ = self.float_, self.short_
        pbar = ProgressBar(data_size, indent='    ')
        self.loaded_matrix = BasisFunctionMatrix(self.basis_sets[0], self.basis_sets[2], self)
        if not self.loaded_matrix.matrix_filled and not self.force_reload:
            while True:
                try:
                    if not self.loaded_matrix.matrix_filled:
                        si1, _, si3, __ = short_(4)
                        s1, s3 = self.basis_sets[0].shells[si1], self.basis_sets[2].shells[si3]
                        nfunc = s1.nfunction * s3.nfunction
                        buff = np.array(float_(nfunc))
                        buff = buff.reshape((s1.nfunction, s3.nfunction))
                        self.loaded_matrix.array[s1.slice, s3.slice] = buff
                        pbar.update(f.tell()-metadata_size)
                except EOFError:
                    # We've reached the end of the file, so break
                    break
                except Exception:
                    # First delete the mem_map_file, since the whole thing was not successfully loaded
                    os.unlink(self.memmap_filename)
                    # then reraise the exception
                    raise
        self.loaded_matrix.matrix_filled = True

    def _load_all(self, f):
        metadata_size = f.tell()
        tot_size = getsize(self.filename)
        data_size = tot_size - metadata_size
        float_, short_ = self.float_, self.short_
        pbar = ProgressBar(data_size, indent='    ')
        self.loaded_matrix = BasisFunctionTensor(self.basis_sets, self)
        if not self.loaded_matrix.matrix_filled and not self.force_reload:
            while True:
                try:
                    si1, si2, si3, si4 = short_(4)
                    s1 = self.basis_sets[0].shells[si1]
                    s2 = self.basis_sets[1].shells[si2]
                    s3 = self.basis_sets[2].shells[si3]
                    s4 = self.basis_sets[3].shells[si4]
                    nfunc = s1.nfunction * s2.nfunction * s3.nfunction * s4.nfunction
                    buff = np.array(float_(nfunc))
                    buff = buff.reshape((s1.nfunction, s2.nfunction, s3.nfunction, s4.nfunction))
                    self.loaded_matrix.array[s1.slice, s2.slice, s3.slice, s4.slice] = buff
                    pbar.update(f.tell()-metadata_size)
                except EOFError:
                    # We've reached the end of the file, so break
                    break
                except Exception:
                    # First delete the mem_map_file, since the whole thing was not successfully loaded
                    os.unlink(self.memmap_filename)
                    # then reraise the exception
                    raise
        self.loaded_matrix.matrix_filled = True


