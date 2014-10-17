(* ::Package:: *)

$LoadFeynArts=True;
$FeynCalcStartupMessages=False;
<<HighEnergyPhysics`FeynCalc`
$FAVerbose=0;
Paint[InsertFields[
    CreateTopologies[1, 1 -> 1,
      ExcludeTopologies -> {Tadpoles}], {F[3,{1}]} -> {F[3,{1}]},
    InsertionLevel -> {Classes},
    Model -> "SMQCD",ExcludeParticles->{S[1],S[2],S[3],V[1],V[2],V[3]}], ColumnsXRows -> {2, 1}, SheetHeader -> False,
PaintLevel -> {Classes},Numbering -> None]//
Export["QCDOneLoopQuarkSelfEnergy.tex",#,"TeX"]&;
