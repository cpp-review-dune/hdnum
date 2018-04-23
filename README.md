Heidelberg Educational Numerics Library (HDNUM)
===============================================

A simple to use and yet efficient C++ library for teaching numerical methods.

License
-------

HDNUM is licensed under version 2 of the GNU
General Public License, with the so-called "runtime exception", as
follows:

>   As a special exception, you may use the HDNUM source files as part
>   of a software library or application without restriction.
>   Specifically, if other files instantiate templates or use macros or
>   inline functions from one or more of the HDNUM source files, or you
>   compile one or more of the HDNUM source files and link them with
>   other files to produce an executable, this does not by itself cause
>   the resulting executable to be covered by the GNU General Public
>   License.  This exception does not however invalidate any other
>   reasons why the executable file might be covered by the GNU General
>   Public License.

This licence clones the one of the libstdc++ library. For further
implications of this library please see their [license page][0].

See the file [COPYING.md][] for full copying permissions.

Installation
------------

There is no installation. Just include the header file `hdnum.hh` and
that is it.

Building the documentation
--------------------------

In the `hdnum` top-level directory just run `doxygen`. This will build
the Doxygen documentation into the directory `doc/html` or `doc/latex`,
respectively.

History
-------

-    Version 0.11 Revision 1620
-    Version 0.12 from November, 5 2009.
-    Version 0.20 from April, 21 2011.
-    Version 0.22 from May, 11 2011. (Add more methods and documentation.)
-    Version 0.23 from June, 14 2011. (Do not pass by reference for element-wise operations!)
-    Version 0.24 from September, 9 2011. (import methods to solve odes/pdes)
-    Version 0.25 from October, 20 2013. (add exceptions in linear algebra)
-    Version 0.26 from October, 24 2013. delete countingptr and arrays

Links
-----

[0]: https://gcc.gnu.org/onlinedocs/libstdc++/faq.html#faq.license
[COPYING.md]: COPYING.md
