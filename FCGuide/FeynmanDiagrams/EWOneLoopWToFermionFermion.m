(* ::Package:: *)

$LoadFeynArts=True;
$FeynCalcStartupMessages=False;
<<HighEnergyPhysics`FeynCalc`
$FAVerbose=0;
Paint[DiagramDelete[InsertFields[
    CreateTopologies[1, 1 -> 2,
      ExcludeTopologies -> {Tadpoles, WFCorrections}],{V[3]} -> {F[2,{1}], -F[1,{1}]},
    InsertionLevel -> {Classes},
    Model -> "SM",ExcludeParticles->{S[1],S[2],S[3]}],3,4], ColumnsXRows -> {2, 1}, SheetHeader -> False,
PaintLevel -> {Classes},Numbering -> None]//Export["EWOneLoopWToFermionFermion.tex",#,"TeX"]&;
