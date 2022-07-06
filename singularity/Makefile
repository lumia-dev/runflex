BIN_DIR?=~/.local/bin
CONTAINER?=runflex.sif
ENVCONTAINER?=runflexenv.sif
METEO_DIR?=/data/FLEXPART/meteo
SCRATCH_DIR?=/data/FLEXPART/tmp

build: runflex.def
	[ ! -f ${CONTAINER}] || rm -f ${CONTAINER}
	singularity build ${CONTAINER} runflex.def

buildenv: runflexenv.def
	[ ! -f ${CONTAINER}] || rm -f ${CONTAINER}
	singularity build ${ENVCONTAINER} runflexenv.def

config:
	python genconfig.py --exec ${BIN_DIR} \
	--container ${CONTAINER} \
	--meteo ${METEO_DIR} \
	--scratch ${SCRATCH_DIR}

install: runflex
	cp runflex ${BIN_DIR}
