#!/bin/bash

javac -cp ..:../lib/core.jar:../lib/vecmath.jar:../lib/objimport.jar:../lib/bsim3.0.jar BSimLacOperon.java

java -cp ..:../lib/core.jar:../lib/vecmath.jar:../lib/objimport.jar:../lib/bsim3.0.jar BSimLacOperon.BSimLacOperon > /dev/null

rm *.class
