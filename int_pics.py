#!/usr/bin/env python
from __future__ import division
from numbers import Number
from os.path import getmtime, exists, join as join_path
import math
import pylab as pl
from matplotlib.colors import Normalize
from itertools import product
from math import log10
from struct import unpack
import os
import numpy as np
import sys
from tempfile import mkdtemp
from glob import glob
import traceback
from binfile_loader import LoadedIntFile, BasisFunctionMatrix
import atexit

convert="/usr/bin/convert"

matrix_folder = "matrices/"

# Figure out what stuff to do
if len(sys.argv) == 2:
    sets = [(sys.argv[1], "cc-pVDZ-F12", "cc-pVDZ-F12-CABS")]
elif len(sys.argv) == 3:
    sets = [(sys.argv[1], sys.argv[2], "cc-pVDZ-F12-CABS")]
elif len(sys.argv) == 4:
    sets = [(sys.argv[1], sys.argv[2], sys.argv[3])]
else:
    # Get everything we can
    # TODO we can do better than this.  It's becoming a mess
    sets = []
    for filen in glob(matrix_folder + "*"):
        parts = filen.split('_')
        if len(parts) >= 4:
            tmp = list(parts[:3])
            qualifiers = ["allmax", "allints", "average", "median", "stddev"]
            qual = [q for q in qualifiers if q in parts[3:-1]]
            if len(qual) > 1:
                raise ValueError("Huh?")
            if len(parts) >= 5 and len(qual) == 1:
                tmp[2] += "_" + qual[0]
            if len(parts) >= 5 and "ri" in parts[3:-1]:
                tmp[2] += "_ri"
            tmp[0] = tmp[0].split('/')[1]
            tmp = tuple(tmp)
            if tmp not in sets:
                sets.append(tuple(tmp))

#TODO find all file names available



min_val = 10
max_val = -1
image_width=1000
do_grid=True
do_scale = False
scale_denom = 30
cmap_name = "hot_r"
png_folder = matrix_folder + "png/"
border=5
bc="blue"
grid_color='blue'
debug = False
grid_thickness=5
pair_load_print = 500000
shell_load_print = 500000
global load_pending
load_pending = None
tmpdir = mkdtemp(dir=os.getenv("SCRATCH", default='/tmp'))
int_mats = ['F', 'g', 'Fg', 'Fsq', 'DC']
dens_mats = ['P', 'Q', 'O', 'I', 'R', '0']
#max_mm_size = float("inf")
max_mm_size = 10 * 1024 * 1024
#max_mm_size = 3000


#--------------------------------------------------------------------------------#

def data_size_str(bytes):
    v = float(bytes)
    if bytes < 800:
        return "{} B".format(bytes)
    elif bytes < 8e5:
        return "{:.2f} KB".format(v/1024)
    elif bytes < 8e8:
        return "{:.2f} MB".format(v/1024/1024)
    elif bytes < 8e11:
        return "{:.2f} GB".format(v/1024/1024/1024)
    else:
        return "{:.2f} TB".format(v/1024/1024/1024/1024)

def load_binfile(mat):
    filename = bin_name(mat)
    #print "  Loading matrix {} from {}".format(mat, filename)
    ints = LoadedIntFile(filename, verbose=True)
    ints.load_mags(min_val, max_val)
    #print "    Loading complete!  Memmap file size is {}".format(
    #    data_size_str(os.path.getsize(mem_map_name(mat)))
    #)
    if ints.ints_type == 0 or ints.ints_type == 32:
        shp = ints.loaded_array.mag_array.shape
        ints.loaded_array.mag_array = ints.loaded_array.mag_array.reshape((shp[0]*shp[1], shp[2]*shp[3]))
    return ints


def add_scale(mat):
    minv = None
    maxv = None
    if min_val is not None:
        minv = min_val
    if max_val is not None:
        maxv = max_val
    if minv is None or maxv is None:
        raise NotImplementedError
    width = mat.shape[1] // scale_denom
    rv = np.ndarray(shape=(mat.shape[0], mat.shape[1] + width))
    rv[:,:mat.shape[1]] = mat
    for col in xrange(mat.shape[1], mat.shape[1]+width):
        for row in xrange(mat.shape[0]):
            rv[row, col] = (maxv-minv) * row/mat.shape[0] - (maxv-minv)
    return rv

global prefix
prefix = None

def png_dir(mat, stuff):
    rv = join_path(matrix_folder, "png", stuff[0])
    extra_stuff = stuff[2].split('_')
    ribas = extra_stuff[0]
    basis = stuff[1]
    qualifiers = ["allmax", "allints", "average", "median", "stddev"]
    rimats = ["R", "O"]
    if "ri" in extra_stuff or any(m in mat[:-1] for m in rimats) or mat in rimats:
        second_dir = basis + "+" + ribas
    else:
        second_dir = basis
    if mat in int_mats:
        if any(q in extra_stuff for q in qualifiers):
            rv = join_path(rv, second_dir, "untrans")
        else:
            rv = join_path(rv, second_dir, "pair_ints")
            raise NotImplementedError("Parsing of pair_ints is not longer supported (encountered in {})".format((stuff,mat)))
    elif mat in ['P', 'Q']:
        rv = join_path(rv, second_dir, "density_mats")
    elif mat == 'O':
        second_dir = basis + "+" + ribas
        rv = join_path(rv, second_dir, "density_mats")
    #elif any(mat == a+b+c for a,b,c in product(dens_mats, dens_mats, int_mats)):
    elif all(a in dens_mats for a in mat[:2]) and mat[2:] in int_mats:
        rv = join_path(rv, second_dir, "half_trans")
    #elif any(mat == a+b+c+d+e for a,b,c,d,e in product(dens_mats, dens_mats, dens_mats, dens_mats, int_mats)):
    elif all(a in dens_mats for a in mat[:4]) and mat[4:] in int_mats:
        rv = join_path(rv, second_dir, "full_trans")
    else:
        raise NotImplementedError
    qual = [q for q in qualifiers[1:] if q in extra_stuff]
    if len(qual) > 1:
        raise ValueError("Huh?")
    elif len(qual) == 1:
        rv = join_path(rv, qual[0])
    return rv

def png_name(mat, stuff, grid=False, resized=False, stat_type="max"):
    if grid and resized:
        raise ValueError
    rv = join_path(png_dir(mat, stuff), mat + "_")
    if stat_type != "max":
        rv += stat_type + "_"
    rv += str(min_val)
    rv += "_" + str(max_val)
    if grid:
        rv += "_grid"
    if do_scale:
        rv += "_scale"
    rv += "_" + cmap_name
    if resized:
        rv += "_small"
    rv += ".png"
    return rv

def bin_name(mat):
    return matrix_folder + prefix + mat + ".bin"

def mem_map_name(mat, mags=False):
    parts = prefix.split("/")
    dot_prefix = "/".join(parts[:-1]) + "/." + parts[-1]
    if mags:
        return matrix_folder + dot_prefix + "_" + str(min_val) + "_" + str(max_val) + "_" + mat + ".npmm_mags"
    else:
        return matrix_folder + dot_prefix + mat + ".npmm"

def is_newer(filea, fileb):
    return getmtime(filea) - getmtime(fileb) > 0

doing_printed = False
def print_doing_if_needed(stuff):
    global doing_printed
    if not doing_printed:
        print "Doing {}...".format(stuff)
        doing_printed = True

#--------------------------------------------------------------------------------#

num_mats = 0
num_processed = 0

for stuff in sets:
    prefix = "_".join(stuff) + "_"
    all_stuff = glob(matrix_folder + "_".join(stuff) + "*.bin")
    mats = [".".join(a.split(".")[:-1]).split("_")[-1] for a in all_stuff]
    mats = list(set(mats))

    doing_printed = False

    num_mats += len([mat for mat in mats if exists(bin_name(mat))])

    for mat in mats:

        if debug: print_doing_if_needed(stuff)

        # Skip if the matrix is not available
        if not exists(bin_name(mat)):
            if debug: print "  Skipping {} because {} does not exist.".format(mat, bin_name(mat))
            continue

        bin = bin_name(mat)
        main_type="max"
        grd = png_name(mat, stuff, grid=True, stat_type=main_type)
        rsz = png_name(mat, stuff, resized=True, stat_type=main_type)
        png = png_name(mat, stuff, stat_type=main_type)

        if not exists(png):
            print_doing_if_needed(stuff)
            print "  Doing {} because full-res png {} does not exist".format(mat, png)
        elif not exists(rsz):
            print_doing_if_needed(stuff)
            print "  Doing {} because resized png {} does not exist".format(mat, rsz)
        elif not "allints" in bin and not exists(grd):
            print_doing_if_needed(stuff)
            print "  Doing {} because gridded png {} does not exist".format(mat, grd)
        elif not "allints" in bin and is_newer(rsz, grd):
            print_doing_if_needed(stuff)
            print "  Doing {} because resized png {} is newer than gridded png {}".format(mat, rsz, grd)
        elif is_newer(png, rsz):
            print_doing_if_needed(stuff)
            print "  Doing {} because full-res png {} is newer than resized png {}".format(mat, png, rsz)
        elif is_newer(bin, png):
            print_doing_if_needed(stuff)
            print "  Doing {} because bin file {} is newer than png file {}".format(mat, bin, png)
        else:
            if debug: print "  Skipping {} because nothing needs to be done.".format(mat)
            continue

        # load the matrix
        try:
            if "allints" in bin_name(mat):
                continue
            lif = load_binfile(mat)
            if lif.ints_type == 0 or lif.ints_type == 32:
                raise NotImplementedError()
        except Exception as e:
            print "  Could not load binary file {}.  The following error was raised:\n    {}".format(
                bin_name(mat), type(e).__name__ + ": " + e.message
            )
            print "    Backtrace:"
            traceback.print_exc(file=sys.stdout)
            continue

        for ty in lif.available_stat_types:
            m = getattr(lif.loaded_stats, ty)
            grd = png_name(mat, stuff, grid=True, stat_type=ty)
            rsz = png_name(mat, stuff, resized=True, stat_type=ty)
            png = png_name(mat, stuff, stat_type=ty)

            num_processed += 1

            # Create a png directory if we need to
            if not exists(png_dir(mat, stuff)):
                os.makedirs(png_dir(mat, stuff))

            # Only write a new png using pylab if we have to
            if not exists(png) or not is_newer(png, bin_name(mat)):
                if not m.mags_filled:
                    raise NotImplementedError("Something went wrong.  Shouldn't get here")
                cmap = pl.get_cmap(cmap_name)
                norm = None
                if min_val is not None and max_val is not None:
                    norm = Normalize(0, min_val-max_val)
                if do_scale:
                    raise NotImplementedError
                mmsize = os.path.getsize(mem_map_name(mat))
                mat_to_use = m.mag_matrix
                if mmsize > max_mm_size:
                    max_chunk = int(math.floor(math.sqrt(max_mm_size / 4.0)))
                    # For now, assume square matrix
                    if mat_to_use.shape[0] != mat_to_use.shape[1]:
                        raise NotImplementedError
                    num_chunks = int(math.ceil(mat_to_use.shape[0]/float(max_chunk)))
                    print "  Running matshow from pylab using {} different {}x{} blocks...".format(
                        num_chunks*num_chunks, max_chunk, max_chunk)
                    tmp_rows = []
                    for ichnk in xrange(num_chunks):
                        tmpfiles = []
                        istart, iend = ichnk*max_chunk, (ichnk+1)*max_chunk
                        for jchnk in xrange(num_chunks):
                            jstart, jend = jchnk*max_chunk, (jchnk+1)*max_chunk
                            chnkmat = mat_to_use[istart:iend, jstart:jend]
                            tmpfiles.append(join_path(tmpdir, mat + "_{:03d}_{:03d}.png".format(ichnk, jchnk)))
                            print "        Running matshow for block {}, {}...".format(ichnk, jchnk)
                            tmp = pl.matshow(chnkmat, cmap=cmap, norm=norm)
                            tmp.write_png(tmpfiles[-1], noscale=True)
                        tmp_rows.append(join_path(tmpdir, mat + "_{:03d}_row.png".format(ichnk)))
                        print "      Compositing row chunks for row chunk {}...".format(ichnk)
                        os.system("{convert} {images} +append {rowimg}".format(
                            images=" ".join(tmpfiles), rowimg=tmp_rows[-1], convert=convert
                        ))
                        for f in tmpfiles:
                            os.unlink(f)
                    print "      Compositing all rows..."
                    os.system("{convert} {images} -append {outimg}".format(
                        images=" ".join(tmp_rows), outimg=png, convert=convert
                    ))
                    for f in tmp_rows:
                        os.unlink(f)
                    print "      Completed multi-pass matrix render!"
                else:
                    print "  Running matshow from pylab..."
                    tmp = pl.matshow(mat_to_use, cmap=cmap, norm=norm)
                    tmp.write_png(png, noscale=True)

            if not exists(rsz) or is_newer(png, rsz):
                image_size = image_width # for now
                print "  Converting image of {} to size {image_size}x{col_size}...".format(
                    mat, col_size = image_size if not do_scale else image_size + image_size // scale_denom,
                    **globals())
                os.system("{convert} {pngname} -scale {image_size}x{col_size} {resized_name}".format(
                    resized_name=rsz,
                    col_size = image_size if not do_scale else image_size + image_size // scale_denom,
                    pngname=png,
                    **globals()))

            if not exists(grd) or is_newer(rsz, grd):
                if isinstance(m, BasisFunctionMatrix):
                    print "  Drawing grid for {}...".format(mat)
                    os.system("{convert} {resized_name} -density 300 -fill none -stroke \"{grid_color}\" -strokewidth 1 -draw \"stroke-opacity 1 path '{drawstr}'\" {grid_version}".format(
                        resized_name=rsz,
                        grid_version=grd,
                        drawstr=m.draw_lines_string(image_width), **globals()
                    ))

            print "  Done with {}\n  {}".format(mat, '-'*50)

    if doing_printed:
        print "Done with {}\n{}".format(stuff, '='*60)

print "Processed {} different data groups of data ({} total matrix files),".format(len(sets), num_mats),
if num_processed == 0:
    print "nothing needed to be done for any of them."
else:
    print "produced or updated png files for {} matrices.".format(num_processed)
