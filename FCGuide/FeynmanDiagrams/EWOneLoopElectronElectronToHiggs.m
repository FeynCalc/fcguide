(* ::Package:: *)

$LoadFeynArts=True;
$FeynCalcStartupMessages=False;
<<HighEnergyPhysics`FeynCalc`
$FAVerbose=0;
tops = CreateTopologies[1, 2 -> 2, Adjacencies -> {3},
ExcludeTopologies -> {SelfEnergies, WFCorrections, Tadpoles,
Boxes[3]}];
inserttops =
InsertFields[tops, {F[2, {1}], -F[2, {1}]} -> {V[2], S[1]}, Model -> "SM",
GenericModel -> "Lorentz", InsertionLevel -> Particles,
ExcludeFieldPoints -> {FieldPoint[F, -F, S]}];
graphs = Paint[inserttops, PaintLevel -> {Particles}, AutoEdit -> False,
SheetHeader -> False, Numbering -> False, ColumnsXRows -> {3, 2}]//
Export["EWOneLoopElectronElectronToHiggs.tex",#,"TeX"]&;
