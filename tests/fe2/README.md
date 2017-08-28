This is a simple regression test for the FE2 execution of the *PGFem3D* code. It consists of a patch test with two macro-scale cohesive elements. The macro-scale force-displacement curve is evaluated over three loading steps and 'diffed' against a reference solution.

The repository contains the required scripts to regenerate the input files, however they are not executed in the regression test scripts, opting for the pre-built files instead. Thus, any differences can be attributed directly to the *PGFem3D* code.

See the [example script](scripts/local-run-test.sh) for hints on formating the Jenkins test.