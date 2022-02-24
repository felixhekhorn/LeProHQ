#!/bin/bash

MODE=2
FILES='hqcoef.f'

if [ "$MODE" -eq "1" ];
then
 f2py3 --overwrite-signature -m hqcoef -h hqcoef.template.pyf $FILES \
  only: c2log cllog c2nlog clnlog c2nlobarg clnlobarg c2nloq clnloq c2nlobarq clnlobarq d2nloq dlnloq   : ;
fi;

if [ "$MODE" -eq "2" ];
 then
  f2py3 -m hqcoef -c hqcoef.pyf $FILES --fcompiler=gfortran;
fi;
