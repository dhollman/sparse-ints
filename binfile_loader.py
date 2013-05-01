from __future__ import division
from itertools import product
from operator import attrgetter
from struct import unpack
import os
import sys
from os.path import exists, getmtime, getsize
import numpy as np
from math import log10


def tensor_from_binfile(filename, force_reload=False, verbose=False):
    int_file = LoadedIntFile(filename, force_reload=force_reload, verbose=verbose)
    return int_file.loaded_array.array

def tensor_and_basis_sets_from_binfile(filename, force_reload=False, verbose=False):
    int_file = LoadedIntFile(filename, force_reload=force_reload, verbose=verbose)
    return int_file.loaded_array.array, int_file.loaded_array.basis_sets

def bsarray_from_binfile(filename, force_reload=False, verbose=False):
    int_file = LoadedIntFile(filename, force_reload=force_reload, verbose=verbose)
    return int_file.loaded_array

def tensor_and_mags_from_binfile(filename, min_value_or_mag_function, max_value=None, force_reload=False, force_reload_mags=False):
    int_file = LoadedIntFile(filename, force_reload=force_reload)
    int_file.load_mags(min_value_or_mag_function, max_value, force_reload=force_reload_mags)
    return int_file.loaded_array.array, int_file.loaded_array.mag_array

#--------------------------------------------------------------------------------#

def is_newer(filea, fileb):
    return getmtime(filea) - getmtime(fileb) > 0

def make_matrix(shape, mm_filename, bin_filename=None, force_reload=False):
    shape = tuple(shape)
    if not force_reload and (exists(mm_filename) and (bin_filename is None or is_newer(mm_filename, bin_filename))):
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

#--------------------------------------------------------------------------------#

class ProgressBar(object):
    def __init__(self, max_val, width=80, out=sys.stdout, fill_char='=', indent='', arrow = True):
        self.max_val = max_val
        self.current_percent = 0
        self.current_value = 0
        self.width = width
        self.out = out
        self.fill_char = fill_char
        self.indent = indent
        self.arrow = arrow

    def update(self, new_val=None):
        if new_val is not None:
            self.current_value = new_val
        new_percent = int(float(self.current_value/self.max_val)*100)
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

    def __iadd__(self, other):
        self.current_value += other
        self.update()

#--------------------------------------------------------------------------------#

class BasisSet(object):

    def __init__(self, shells):
        # Sort is guarenteed to be stable
        self.shells = shells
        if not shells == sorted(shells, key=attrgetter('center')):
            raise NotImplementedError("Shell center ordering was: {}".format(
                ", ".join(str(s.center) for s in shells)
            ))
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

    def __iter__(self):
        for sh in self.shells:
            yield sh

    def __eq__(self, other):
        return self.shells == other.shells

#--------------------------------------------------------------------------------#

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

#--------------------------------------------------------------------------------#

class IntStats(object):

    def __init__(self):
        self.values = None
        self.max = None
        self.average = None
        self.median = None
        self.stddev = None
        self.histogram = None


#--------------------------------------------------------------------------------#

class BasisFunctionTensor(object):

    def __init__(self, basis_sets, int_file_obj, force_reload=False, stat_type="max"):
        self.basis_sets = basis_sets
        self.shape = tuple(b.nbf for b in self.basis_sets)
        self.int_file_obj = int_file_obj
        self.mag_arrays = {}
        self.stat_type = stat_type
        self.array, self.matrix_filled = make_matrix(
            self.shape,
            int_file_obj.memmap_filename(stat_type),
            bin_filename=int_file_obj.filename,
            force_reload=force_reload
        )
        self.mag_matrix, self.mags_filled = None, False

    def init_mag_matrix(self, min_val_or_name, max_val=None):
        if self.mags_filled:
            self.mag_arrays[self.mag_filename] = self.mag_array
            self.mags_filled = False
        self.mag_filename = self.int_file_obj.mags_memmap_filename(min_val_or_name, max_val, stat_type=self.stat_type)
        self.mag_array, self.mags_filled = make_mag_matrix(
            self.shape,
            self.mag_filename,
            self.int_file_obj.memmap_filename(stat_type=self.stat_type),
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

#--------------------------------------------------------------------------------#

class BasisFunctionMatrix(BasisFunctionTensor):

    def __init__(self, bs1, bs3, int_file_obj, **kwargs):
        super(BasisFunctionMatrix, self).__init__((bs1,bs3), int_file_obj, **kwargs)

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

#--------------------------------------------------------------------------------#

class LoadedIntFile(object):

    stat_types = {
        "max" : 2,
        "average": 4,
        "stddev": 8,
        "median": 16,
        "histogram": 128
    }

    def __init__(self, filename, force_reload=False, verbose=False):
        self.filename = filename
        self.force_reload = force_reload
        self.mag_filenames = []
        self.verbose = verbose
        self._load()

    def memmap_filename(self, stat_type="max"):
        parts = self.filename.split('/')
        dot_filename = "/".join(parts[:-1]) + "/." + parts[-1]
        if dot_filename[-4:] == ".bin":
            dot_filename = dot_filename[:-4]
        if stat_type != "max":
            dot_filename += "_" + stat_type
        return dot_filename + ".npmm"

    def mags_memmap_filename(self, minval_or_name, max_val=None, stat_type="max"):
        rv = self.memmap_filename(stat_type) + "_mags_" + str(minval_or_name)
        if max_val is not None:
            rv += "_" + str(max_val)
        self.mag_filenames.append(rv)
        return rv

    @property
    def loaded_array(self):
        return self.loaded_matrix
    @loaded_array.setter
    def loaded_array(self, val):
        self.loaded_matrix = val

    def delete_mags(self):
        """
        Don't use this.  It's unsafe.  It would be better to
        tell the memmap object to delete its associated file
        when the object is deleted.
        """
        raise NotImplementedError()
        #del self.loaded_matrix.mag_array
        #for key in self.loaded_array.mag_arrays:
        #    self.loaded_array.mag_arrays.pop(key)
        #self.loaded_array.mag_array = None
        #self.loaded_array.mags_filled = False
        #for mmfn in self.mag_filenames:
        #    os.unlink(mmfn)

    def load_mags(self, mag_function_or_min_val, max_val=None, force_reload=False):
        bsmats = [getattr(self.loaded_stats, ty) for ty in self.available_stat_types]
        for bsmat, stat_type in zip(bsmats, self.available_stat_types):
            bsmat.init_mag_matrix(mag_function_or_min_val, max_val)
            if not force_reload and (bsmat.mags_filled
                    or self.mags_memmap_filename(mag_function_or_min_val, max_val, stat_type=stat_type) in bsmat.mag_arrays):
                if self.mags_memmap_filename(mag_function_or_min_val, max_val, stat_type=stat_type) in bsmat.mag_arrays:
                    # Swap the active mag_array
                    bsmat.mag_arrays[bsmat.mag_filename] = bsmat.mag_array
                    bsmat.mag_array = bsmat.mag_arrays[self.mags_memmap_filename(mag_function_or_min_val, max_val, stat_type=stat_type)]
                    bsmat.mag_filename = self.mags_memmap_filename(mag_function_or_min_val, max_val, stat_type=stat_type)
                #if self.verbose: print "    Magnitudes already loaded."
                return
            else:
                if callable(mag_function_or_min_val):
                    mag = mag_function_or_min_val
                    raise NotImplementedError("Naming scheme for mag function not yet implemented")
                else:
                    min_val = mag_function_or_min_val
                    def mag(val):
                        if min_val is not None and abs(val) < pow(10.0, -min_val):
                            return 0
                        elif max_val is not None and abs(val) > pow(10.0, -max_val):
                            return min_val-max_val
                        return log10(abs(val)) + min_val
                mag_vect = np.vectorize(mag)
                #print "    Loading magnitudes."
                bsmat.mag_array = mag_vect(self.loaded_matrix.array)
                bsmat.mags_filled = True

    def _load(self):
        with open(self.filename, 'rb') as f:
            self.did_load = False
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
            self.short_ = short_ = lambda x: getter(x, 'h', 2)
            self.ushort_ = ushort_ = lambda x: getter(x, 'H', 2)
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
            # Check to see if this is a "versioned" file
            #   or from a time before I used the version
            #   metadata entry
            chkval, self.version = short_(2)
            if chkval != -1:
                # seek backwards 4 bytes from the current
                #   position (which is what the 1 means)
                f.seek(-4, 1)
                self.version = 0
            #========================================#
            if ints_type == 0:
                # File contains all integrals
                self._load_basis(f)
                self._load_all(f)
            #========================================#
            elif ints_type == 1:
                # Legacy support for bs1==bs2==bs3==bs4 max only files
                self._load_basis_legacy(f)
                self._load_stats(f, 2)
            #========================================#
            elif any(ints_type & ty for ty in self.stat_types.values()):
                # File only contains maxima/average/std_dev/whatever
                self._load_basis(f)
                self._load_stats(f, ints_type)
            #========================================#
            elif ints_type == 32:
                # Special untransformed integrals format
                self._load_basis(f)
                self._load_untrans_all(f)
            #========================================#
            elif ints_type == 64:
                # Density matrix format
                self._load_dens_basis(f)
                self._load_dens(f)
            #========================================#
            else:
                raise NotImplementedError("Unrecognized ints_type {}".format(ints_type))
            #========================================#
            if self.verbose and self.did_load: print "    Done loading bin file into memeory mapped structure."

    def _load_dens_basis(self, f):
        int_ = self.int_
        natoms, nshell = int_(2)
        shells = []
        tot_bf = 0
        for ish in xrange(nshell):
            center, nbf, max_am, min_am = int_(4)
            sh = Shell(ish, center, nbf, max_am, min_am)
            sh.first_bf = tot_bf
            tot_bf += nbf
            sh.nprimitive = int_(1)[0]
            sh.ncontraction = int_(1)[0]
            shells.append(sh)
        self.basis_sets = [BasisSet(shells)]

    def _load_dens(self, f):
        metadata_size = f.tell()
        tot_size = getsize(self.filename)
        data_size = tot_size - metadata_size
        pbar = ProgressBar(data_size, indent='    ')
        float_ = self.float_
        # ugly hack
        self.loaded_stats = IntStats()
        self.loaded_matrix = BasisFunctionMatrix(self.basis_sets[0], self.basis_sets[0], self, force_reload=self.force_reload)
        self.loaded_stats.max = self.loaded_matrix
        if not self.loaded_matrix.matrix_filled or self.force_reload:
            for row in xrange(self.basis_sets[0].nbf):
                for col in xrange(row + 1):
                    val = float_(1)[0]
                    self.loaded_matrix.array[row, col] = val
                    self.loaded_matrix.array[col, row] = val
                    if self.verbose: pbar.update(f.tell()-metadata_size)
            # New line at the end of the progress bar
            print
        else:
            pass
            #print "    Warning:  Matrix not reloaded"
        self.loaded_matrix.matrix_filled = True
        # Ugly hack again...
        self.available_stat_types = ['max']

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
                sh.nprimitive = int_(1)[0]
                sh.ncontraction = int_(1)[0]
                shells.append(sh)
            self.basis_sets.append(BasisSet(shells))

    def _load_stats(self, f, ints_type):
        metadata_size = f.tell()
        tot_size = getsize(self.filename)
        data_size = tot_size - metadata_size
        float_, short_ = self.float_, self.short_
        pbar = ProgressBar(data_size, indent='    ')
        self.loaded_stats = IntStats()
        for stat_type in self.stat_types:
            if ints_type & self.stat_types[stat_type]:
                tmp = BasisFunctionMatrix(
                    self.basis_sets[0],
                    self.basis_sets[2],
                    self,
                    force_reload=self.force_reload,
                    stat_type=stat_type
                )
                setattr(self.loaded_stats, stat_type, tmp)
        self.loaded_matrix = self.loaded_stats.max
        stat_types_sorted = sorted(self.stat_types.keys(), key=lambda k: self.stat_types[k])
        needed_stats = [ty for ty in stat_types_sorted if self.stat_types[ty] & ints_type]
        if "histogram" in needed_stats:
            raise NotImplementedError()
        bsmats = [getattr(self.loaded_stats, ty) for ty in needed_stats]
        if not all(bsmat.matrix_filled for bsmat in bsmats) or self.force_reload:
            self.did_load = True
            if self.verbose: print "  Loading matrix from {}".format(self.filename)
            while True:
                try:
                    # Read the identifier
                    si1, _, si3, __ = short_(4)
                    s1, s3 = self.basis_sets[0].shells[si1], self.basis_sets[2].shells[si3]
                    nfunc = s1.nfunction * s3.nfunction

                    # loop over the possible types of statistics
                    for bsmat in bsmats:
                        # read the data regardless
                        buff = np.array(float_(nfunc))
                        # only write to the matrix if it's writable
                        if not bsmat.matrix_filled:
                            if any(abs(b) > 1e10 for b in buff):
                                raise ValueError("Strange value {} for indices {}, {}".format(
                                    np.amax(abs(buff)), si1, si3
                                ))
                            buff = buff.reshape((s1.nfunction, s3.nfunction))
                            bsmat.matrix[s1.slice, s3.slice] = buff
                        if self.verbose: pbar.update(f.tell()-metadata_size)
                except EOFError:
                    # We've reached the end of the file, so break
                    break
                except Exception:
                    # First delete the mem_map_file, since the whole thing was not successfully loaded
                    for ty in needed_stats:
                        os.unlink(self.memmap_filename(stat_type=ty))
                    # then reraise the exception
                    raise
            # New line at the end of the progress bar
            print
        else:
            pass
            #print "    Warning:  Matrix not reloaded"
        for bsmat in bsmats:
            bsmat.matrix_filled = True
        self.available_stat_types = needed_stats

    def _load_all(self, f):
        metadata_size = f.tell()
        tot_size = getsize(self.filename)
        data_size = tot_size - metadata_size
        float_, short_ = self.float_, self.short_
        pbar = ProgressBar(data_size, indent='    ')
        self.loaded_matrix = BasisFunctionTensor(self.basis_sets, self, force_reload=self.force_reload)
        if not self.loaded_matrix.matrix_filled or self.force_reload:
            if self.verbose: print "  Loading matrix from {}".format(self.filename)
            while True:
                try:
                    si1, _, si3, __ = short_(4)
                    if si1 > self.basis_sets[0].nshell:
                        raise IndexError("got index si1 = {} which is out of range".format(si1))
                    if si3 > self.basis_sets[2].nshell:
                        raise IndexError("got index si3 = {} which is out of range".format(si3))
                    s1 = self.basis_sets[0].shells[si1]
                    s3 = self.basis_sets[2].shells[si3]
                    nfunc = s1.nfunction * s3.nfunction * self.basis_sets[1].nbf * self.basis_sets[3].nbf
                    buff = np.array(float_(nfunc))
                    buff = buff.reshape((s1.nfunction, s3.nfunction, self.basis_sets[1].nbf, self.basis_sets[3].nbf))
                    self.loaded_matrix.array[s1.slice, :, s3.slice, :] = buff.transpose([0,2,1,3])
                    if self.verbose: pbar.update(f.tell()-metadata_size)
                except EOFError:
                    # We've reached the end of the file, so break
                    break
                except Exception:
                    # First delete the mem_map_file, since the whole thing was not successfully loaded
                    os.unlink(self.memmap_filename())
                    # then reraise the exception
                    raise
            # New line at the end of the progress bar
            print
        else:
            print "    Warning:  Matrix not reloaded"
        self.loaded_matrix.matrix_filled = True

    def _load_untrans_all(self, f):
        metadata_size = f.tell()
        tot_size = getsize(self.filename)
        data_size = tot_size - metadata_size
        float_, short_ = self.float_, self.short_
        pbar = ProgressBar(data_size, indent='    ')
        self.loaded_matrix = BasisFunctionTensor(self.basis_sets, self, force_reload=self.force_reload)
        loaded_quartets = np.ndarray(
            shape=(
                self.basis_sets[0].nshell,
                self.basis_sets[1].nshell,
                self.basis_sets[2].nshell,
                self.basis_sets[3].nshell
            ),
            dtype=bool
        )
        loaded_quartets[...] = False
        if not self.loaded_matrix.matrix_filled or self.force_reload:
            if self.verbose: print "  Loading matrix from {}".format(self.filename)
            while True:
                try:
                    si1, si2, si3, si4 = short_(4)
                    loaded_quartets[si1,si2,si3,si4] = True
                    s1 = self.basis_sets[0].shells[si1]
                    s2 = self.basis_sets[1].shells[si2]
                    s3 = self.basis_sets[2].shells[si3]
                    s4 = self.basis_sets[3].shells[si4]
                    nfunc = s1.nfunction * s2.nfunction * s3.nfunction * s4.nfunction
                    buff = np.array(float_(nfunc))
                    buff = buff.reshape((s1.nfunction, s2.nfunction, s3.nfunction, s4.nfunction))
                    self.loaded_matrix.array[s1.slice, s2.slice, s3.slice, s4.slice] = buff
                    if self.verbose: pbar.update(f.tell()-metadata_size)
                except EOFError:
                    # We've reached the end of the file, so break
                    break
                except Exception:
                    # First delete the mem_map_file, since the whole thing was not successfully loaded
                    os.unlink(self.memmap_filename())
                    # then reraise the exception
                    raise
            # New line at the end of the progress bar
            print
            if not np.all(loaded_quartets.ravel()):
                print "  Warning: not all quartets were found."
                numnl = 0
                for i,j,k,l in product(*[xrange(s) for s in loaded_quartets.shape]):
                    if not loaded_quartets[i,j,k,l]:
                        print "    Not loaded: {}".format((i,j,k,l))
                        numnl+= 1
                print "    Total shell quartets not loaded: {}".format(numnl)
        else:
            print "    Warning:  Matrix not reloaded"

        self.loaded_matrix.matrix_filled = True
