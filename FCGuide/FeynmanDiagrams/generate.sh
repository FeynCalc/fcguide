#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

math -nopromt -script QCDOneLoopGluonSelfEnergy.m
sed -i 's|(432,504)|(432,150)|g' QCDOneLoopGluonSelfEnergy.tex
sed -i 's|u_k|q|g' QCDOneLoopGluonSelfEnergy.tex
sed -i 's|u_g|gh|g' QCDOneLoopGluonSelfEnergy.tex

math -nopromt -script QCDOneLoopQuarkSelfEnergy.m
sed -i 's|(432,504)|(432,150)|g' QCDOneLoopQuarkSelfEnergy.tex

math -nopromt -script QCDOneLoopGhostSelfEnergy.m
sed -i 's|(432,504)|(432,150)|g' QCDOneLoopGhostSelfEnergy.tex
sed -i 's|u_g|gh|g' QCDOneLoopGhostSelfEnergy.tex

math -nopromt -script EWOneLoopElectronElectronToHiggs.m

math -nopromt -script EWOneLoopWToFermionFermion.m
sed -i 's|{$W$}|{}|g' EWOneLoopWToFermionFermion.tex
sed -i 's|{$\\gamma$}|{}|g' EWOneLoopWToFermionFermion.tex
sed -i 's|{$Z$}|{}|g' EWOneLoopWToFermionFermion.tex
sed -i 's|{$e$}|{}|g' EWOneLoopWToFermionFermion.tex
sed -i 's|{$\\nu_e$}|{}|g' EWOneLoopWToFermionFermion.tex
sed -i 's|(432,504)|(432,150)|g' EWOneLoopWToFermionFermion.tex
