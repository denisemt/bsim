# Consortium oscillator example

The simulation shows the system presented in "Emergent genetic oscillations in a synthetic microbial consortium" (Ye Chen et al., Science, 2015: http://science.sciencemag.org/content/349/6251/986.long). 
We adapted the DDE model presented in the above article and implemented the model in BSim, in order to demonstrate new features available in the software (DDE simulation, growth and collision of capsular bacteria).

A video of the simulation can be seen [here](https://www.youtube.com/watch?v=FpG7EgIC5yI)


##### How to run it
# Running the example

Assuming the BSim repository root directory is `$BSIM_SRC`.

First, build everything:
```
cd $BSIM_SRC

git pull

ant -f bsim-build-tree.xml
```

Then go to the build dir:
```
cd out/production/$BSIM_SRC/

java  -cp .:../../../lib/core.jar:../../../lib/vecmath.jar:../../../lib/objimport.jar:../../../lib/jcommander-1.49-SNAPSHOT.jar:../../../lib/bsim-osp.jar BSimChenOscillator.BSimChenOscillatorExample
```