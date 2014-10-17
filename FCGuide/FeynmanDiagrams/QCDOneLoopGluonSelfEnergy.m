(* ::Package:: *)

$LoadFeynArts=True;
$FeynCalcStartupMessages=False;
<<HighEnergyPhysics`FeynCalc`
$FAVerbose=0;
Paint[DiagramDelete[InsertFields[
    CreateTopologies[1, 1 -> 1,
      ExcludeTopologies -> {Tadpoles}], {V[5]} -> {V[5]},
    InsertionLevel -> {Classes},
    Model -> "SMQCD",ExcludeParticles->{S[1],S[2],S[3],V[1],V[2],V[3]}],3], ColumnsXRows -> {4, 1}, SheetHeader -> False,
PaintLevel -> {Classes},Numbering -> None]//
Export["QCDOneLoopGluonSelfEnergy.tex",#,"TeX"]&;
