*****************************************************************************************************************************
*************************************Polinomial Root Finder******************************************************************
*****************************************************************************************************************************

README

Polinomial root finder is a command line tool to find the roots of real and complex polynomials.
The program is able to work with float, double, complex float and complex double data types, and it
implements the polynomial deflation, Müller and Jenkis-Traub algorithms to find out the roots given a polynomial. The program should have a flag called start, that could't be implemented.

This program has a complete implementation of the polynomial deflation method, and it can be used with every polynomial.

The Müller implementation works propertly, but it gets a NaN result on the case where there are no roots on the given polynomial or this one only have complex roots and you want to get real results.

The Jenkins-Traub algorithm has a complete implementation and it should work propertly on every case.

You can find a file with different cases for test, using all types of data and Deflation and Müller methods with some preloaded polynomials.

CONTACT

If you have problems or comments with this program you
can contact please contact us.

This project can also be found at GitHub in:
https://github.com/luisf1997cr/ANPI_Proyecto1

______________________________________________________________________________________________________________________________
_________________________________________________Run Instructions_____________________________________________________________
______________________________________________________________________________________________________________________________


Unzip the project, open a terminal and change your working directory to the unzipped folder

Create a directory build:

> mkdir build;

Go into that directory

> cd build;

You can choose to build a release version with:

> cmake ../code -DCMAKE_BUILD_TYPE=Release

or a debug version with

> cmake ../code -DCMAKE_BUILD_TYPE=Debug

And build everything with

> make

The executables will be stored at build/bin.

To execute the program go to the /build/bin directory

> cd build/bin

To execute the program, running muller method with default values you can use

> ./proyecto1

If you want to configure the use of the program

> ./proyecto1 -help

__________________________________________________________________________________________________________________________
______________________________________________________Dependencies________________________________________________________
__________________________________________________________________________________________________________________________

You need CMAKE and Boost to build the program

> sudo apt-get install libboost-all-dev
> sudo apt-get -y install cmake



