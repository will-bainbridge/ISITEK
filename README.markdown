ISITEK
======

* Author:	Will Bainbridge
* Date:		January, 2012
* GitHub:	<https://github.com/maninthemail/ISITEK>

LICENSE
-------

ISITEK is provided under the terms of the MIT License. See LICENSE or <http://www.opensource.org/licenses/MIT> for details.

DESCRIPTION
-----------

* ISITEK is a 2D unstructured discontinuous galerkin PDE solver.

INSTALLATION
------------

Use the install script in the third party directory to obtain and build the necessary linear algebra packages, GotoBLAS and UMFPACK

	cd thirdparty
	./install.sh

Alternatively, edit the makefile to reference intel MKL libraries

Back in the root directory, build the solver

	cd ..
	make

USAGE
-----

Move into one of the example directories

	cd cavity

Run the program with the input file as the sole argument

	../isitek cavity.INS.input

Paraview can then be used to view the resulting vtu files
