//
// my_macros.h
//
// Copyright (C) 2013 David Hollman
//
// Author: David S. Hollman <dhollman@uga.edu>
// Maintainer: DSH
//
// This file contains some useful macros
// I realize that this sort of practice, in general, makes
// code harder to read.  Thus, I've tried to keep the custom
// macro usage to a minimum and documented the few that I do
// use thoroughly in the my_macros.h file.


//############################################################################//
// Macros to be used locally (i.e. only in this file).  Ignore these
//   if you've come here to figure out what the heck I'm doing with a macro
//   used elsewhere.
//----------------------------------------------------------------------------//
// The PP_NARG macro returns the number of arguments that have been passed in
#define PP_NARG(...) \
         PP_NARG_(__VA_ARGS__,PP_RSEQ_N())
#define PP_NARG_(...) \
         PP_ARG_N(__VA_ARGS__)
#define PP_ARG_N( \
          _1, _2, _3, _4, _5, _6, _7, _8, _9,_10, \
         _11,_12,_13,_14,_15,_16,_17,_18,_19,_20, \
         _21,_22,_23,_24,_25,_26,_27,_28,_29,_30, \
         _31,_32,_33,_34,_35,_36,_37,_38,_39,_40, \
         _41,_42,_43,_44,_45,_46,_47,_48,_49,_50, \
         _51,_52,_53,_54,_55,_56,_57,_58,_59,_60, \
         _61,_62,_63,N,...) N
#define PP_RSEQ_N() \
         63,62,61,60,                   \
         59,58,57,56,55,54,53,52,51,50, \
         49,48,47,46,45,44,43,42,41,40, \
         39,38,37,36,35,34,33,32,31,30, \
         29,28,27,26,25,24,23,22,21,20, \
         19,18,17,16,15,14,13,12,11,10, \
         9,8,7,6,5,4,3,2,1,0


#define _for2(x,y) \
    for(int x = 0; x < y; ++x)
#define _for4(x1,y1,x2,y2) \
    for(int x1 = 0; x1 < y1; ++x1) \
        for(int x2 = 0; x2 < y2; ++x2)
#define _for6(x1,y1,x2,y2,x3,y3) \
    for(int x1 = 0; x1 < y1; ++x1) \
        for(int x2 = 0; x2 < y2; ++x2) \
            for(int x3 = 0; x3 < y3; ++x3)
#define _for8(x1,y1,x2,y2,x3,y3,x4,y4) \
    for(int x1 = 0; x1 < y1; ++x1) \
        for(int x2 = 0; x2 < y2; ++x2) \
            for(int x3 = 0; x3 < y3; ++x3) \
                for(int x4 = 0; x4 < y4; ++x4)

#define macro_dispatcher(func, ...) \
            macro_dispatcher_(func, PP_NARG(__VA_ARGS__))
#define macro_dispatcher_(func, nargs) \
            macro_dispatcher__(func, nargs)
#define macro_dispatcher__(func, nargs) \
            func ## nargs
//############################################################################//


// Create nested for loops on one line.  For instance,
//     for_each(i,imax,j,jmax)
//   expands to:
//	   for(int i = 0; i < imax; ++i)
//	       for(int j = 0; j < imax; ++j)
#define for_each(...) macro_dispatcher(_for, __VA_ARGS__)(__VA_ARGS__)

// Iterate over something that implements the begin() and end() methods with
//   an existing iterator object (i.e. this macro does not declare the iterator,
//   only initializes it).
#define iterate(iter, over_what) \
	for(iter = over_what.begin(); iter != over_what.end(); iter++)

// For use with "assert(not_implemented)".  This is probably a bit too clever
//    and should be avoided in public code...
#define not_implemented false


#define PRINT_LIST(list, length, width, nperline) \
		for_each(ind, length){ \
			cout << setw(width) << list[ind]; \
			if((ind+1) % nperline == 0){ \
				cout << endl; \
			} \
		} \
		if (length % nperline != 0) \
			cout << endl;
