
# Paths
SRC_MAIN=src/flexpart10.4
SRC_EXTRA=src/extras
BUILDDIR=build/flexpart
MAKEFILE=makefile.singularity.gfortran


build:
	python -m pip install loguru
	python runflex/compile.py --build ${BUILDDIR} --src ${SRC_MAIN} --extra ${SRC_EXTRA} --makefile ${MAKEFILE}

install:
	python -m pip install -e .[interactive]

clean:
	rm -Rf pyflex.egg-info
	rm -Rf build
	rm -f share/flexpart.x
	pip uninstall runflex

container:
	singularity build --fakeroot runflex.sif runflex.def

envcontainer:
	singularity build --fakeroot runflexenv.sif runflexenv.def
