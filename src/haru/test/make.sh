#!/bin/sh

cd .. make -f makefile.cygwin
cd test
gcc -o Tests Tests.cc -DDEBUG -g -lharu -lz -lstdc++ -I.. -L..
./Tests
exit $?

