#!/bin/bash
EXEC=$SCRATCH/100k_sequential

echo ${0} ${1} ${2}
$EXEC ${1} >${2}

