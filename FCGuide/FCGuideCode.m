<< "HighEnergyPhysics`FeynCalc`"

ScalarProduct[p, q]

MetricTensor[\[Mu], \[Nu]]

FourVector[p, \[Mu]]

FourVector[p - 2q, \[Mu]]

FourVector[q, \[Mu], Dimension -> D]

PolarizationVector[k, \[Mu]]

Conjugate[PolarizationVector[k, \[Mu]]]

(DiracSlash[p + q] + m).DiracMatrix[\[Mu]]

DiracMatrix[\[Mu], 5, \[Nu], 6, \[Rho], 7]

DiracSlash[2b, a, 2(d - c), 6q - 3p]

DiracSlash[Polarization[k]]

DiracMatrix[\[Mu], Dimension -> D]

DiracMatrix[\[Mu], \[Mu], Dimension -> D - 4]

SpinorUBar[p, m].DiracSlash[p]

DiracSlash[p].SpinorU[p, m]

SpinorVBar[p, m].DiracSlash[q, p]

DiracSlash[p].SpinorV[p, m]

SpinorU[p, 0]

SpinorU[p]

SUNT[a] .SUNT[b].SUNT[c]

SUNT[a, b, c]

SUNT[a, b, c]SUNDelta[c, d]

SUNF[a, b, c]

SUND[a, b, c]

SUNF[a, b, c, Explicit -> True]

SUNF[a, b, c, Explicit -> True] + SUNF[a, c, b, Explicit -> True] // Expand

SUND[a, b, c, Explicit -> True] - SUND[a, c, b, Explicit -> True] // Expand

FeynAmpDenominator[PropagatorDenominator[q, m1], 
  PropagatorDenominator[q + p, m2]]

QuantumField[\[Psi]]

QuantumField[\[Psi], {i}]

QuantumField[g, {\[Mu]}, {i}]

QuantumField[PartialD[\[Nu]], g, {\[Mu]}, {i}]

GluonPropagator[p, \[Mu], a,\[NoBreak]\[Nu], b]

GP[q, 1, 2]

GP[q, 1, 2, OPE -> True]

QuarkPropagator[p]

Options[QP]



DataType[f, g, NonCommutative] = True;
t = f.g - g.(2a).f

DotSimplify[t]

DataType[f, g, NonCommutative] = False;
DataType[m, odd] = DataType[a, even] = True;
ptest1[x_] := x /. (-1)^ n_ :> -1 /; DataType[n, odd] ;
ptest2[x_] := x /. (-1)^ n_ :> 1 /; DataType[n, even] ;
t = (-1)^m + (-1)^a + (-1)^z

ptest1[t]

ptest2[%]

Clear[t]
SetOptions[$FrontEnd, 
  "CommonDefaultFormatTypes" -> {"Output" -> StandardForm}]
DiracMatrix[\[Mu]]

DiracMatrix[\[Mu], Dimension -> D]

DiracSlash[q].DiracSlash[q, Dimension -> D]

DiracSimplify[%]

MetricTensor[\[Mu], \[Nu], Dimension -> D]

% /. D -> 4

FourVector[p, \[Mu]]

FourVector[p, \[Mu]] // FeynCalcInternal

FourVector[p, \[Mu]]FourVector[p, \[Mu]] // Contract

FourVector[p, \[Mu], Dimension -> D - 4] // FeynCalcInternal

% /. D -> 4

PolarizationVector[k, \[Mu]]

SpinorV[-p, m] - SpinorU[p, m] // FeynCalcInternal

FeynCalcExternal[DiracMatrix[\[Mu]], FinalSubstitutions -> {\[Mu] -> mu}]

FeynCalcInternal[%, FinalSubstitutions -> {mu -> \[Sigma]}]

SetOptions[$FrontEnd, 
  "CommonDefaultFormatTypes" -> {"Output" -> TraditionalForm}]
?FeynRule


Contract[MetricTensor[\[Mu], \[Mu]]]

Contract[MetricTensor[\[Mu], \[Mu], Dimension -> D]]

Contract[MetricTensor[\[Alpha], \[Beta]]FourVector[p, \[Beta]]]

Contract[FourVector[q, \[Alpha]]FourVector[p - q, \[Alpha]]]

FourVector[2p, \[Mu]]FourVector[2p, \[Mu]] // Contract

Contract[MetricTensor[\[Alpha], \[Beta]]DiracMatrix[\[Alpha]]]

Contract[FourVector[q, \[Alpha]]DiracMatrix[\[Alpha]]]

Contract[LeviCivita[\[Alpha], \[Nu], \[Rho], \[Sigma]]LeviCivita[\[Beta], \
\[Nu], \[Rho], \[Sigma]]]

Contract[MetricTensor[\[Alpha], \[Sigma]]*
    FourVector[p, \[Alpha]] FourVector[
      p, \[Sigma]]*(FourVector[q, \[Beta]] + 
        FourVector[r, \[Beta]])*(FourVector[p, \[Beta]] - 
        FourVector[q, \[Beta]]), Expanding -> False]

Contract[FourVector[k, \[Mu]]PolarizationVector[k, \[Mu]]]

ExpandScalarProduct[ScalarProduct[a + b, c - 2 d]]

MomentumCombine[%]

Contract[FourVector[q, \[Alpha]]FourVector[p - q, \[Alpha]]]

% /. ScalarProduct[q, q] -> 0

ScalarProduct[q, q] = 0

Contract[FourVector[q, \[Alpha]]FourVector[p - q, \[Alpha]]]

DownValues[Pair] = Select[DownValues[Pair], FreeQ[#, q] &];
ExpandScalarProduct[FourVector[a + b, \[Mu]]]

(# == DiracEquation[#]) &[DiracSlash[p].Spinor[p, m]]

(# == Chisholm[#]) &[DiracMatrix[\[Mu], \[Nu], \[Rho]]]

DiracOrder[DiracMatrix[\[Beta], \[Alpha]]]

DiracOrder[%, {\[Beta], \[Alpha]}]

DiracOrder[DiracMatrix[\[Mu], \[Mu]], DiracSlash[p, p]]

DiracOrder[DiracMatrix[\[Alpha], \[Mu], \[Alpha], Dimension -> D]]

DiracOrder[DiracSlash[-p, q, p]]

DiracSlash[2b, a, 2(d - c), 6q - 3p] // DotSimplify

GA[\[Mu]].(a GS[p] - b GS[q]).GS[q].GA[\[Nu]]

DotSimplify[%, Expanding -> True]

DotSimplify[%, DotPower -> True]

DotSimplify[%, DotSimplifyRelations -> {GS[q]^2 -> 1}, DotPower -> True]

DeclareNonCommutative[a, b, c];
Commutator[a, c] = 1;
DotSimplify[a.(b - z c).a]

UndeclareNonCommutative[a, b, c];
Commutator[a, c] = 0;
DiracSimplify[DiracMatrix[\[Mu], \[Mu], Dimension -> D]]

DiracSimplify[DiracMatrix[\[Mu], \[Nu], \[Rho], \[Sigma], \[Mu]]]

DiracSimplify[
  1/2DiracMatrix[\[Mu], \[Alpha], \[Beta], \[Gamma], \[Delta], \[Mu]]]

DiracSimplify[DiracSlash[p], DiracSlash[-q] + m, DiracSlash[p]]

DiracSimplify[DiracMatrix[5], DiracMatrix[\[Mu]]]

DiracSimplify[DiracMatrix[6, \[Nu], 7, \[Mu]]]

DiracSimplify[(DiracSlash[p] - m).SpinorU[p, m]]

DiracSimplify[(DiracSlash[p] + m).SpinorV[p, m]]

DiracSimplify[SpinorUBar[p, m].(DiracSlash[p] - m)]

DiracSimplify[SpinorVBar[p, m].DiracSlash[q].(DiracSlash[p] - m)]

DiracSimplify[SpinorVBar[p, m1].DiracSlash[q, p].SpinorU[q, m2]]

DiracReduce[DiracMatrix[\[Mu], \[Nu]]]

DiracReduce[DiracMatrix[\[Mu], \[Nu], \[Rho]]]

DiracReduce[DiracMatrix[\[Mu], \[Nu], \[Rho], \[Sigma]]]

Contract[DiracSimplify[%.%]]

Calc[%%.%%]

$BreitMaison = True;
DiracMatrix[5]. DiracMatrix[\[Mu], Dimension -> D]

DiracSimplify[%]

DiracMatrix[6] .DiracMatrix[\[Mu], Dimension -> D]

% // DiracSimplify

$BreitMaison = False;
Tr[DiracMatrix[\[Alpha], \[Beta]]]

Tr[DiracSlash[a, b, c, d]]

Tr[DiracMatrix[\[Alpha], \[Beta], \[Gamma], \[Delta], 5]]

Tr[MetricTensor[\[Alpha], \[Beta]]/
        4 DiracMatrix[\[Mu]].DiracMatrix[\[Alpha]]FourVector[
        p, \[Mu]]] // Contract

pps = DiracSlash[p'];
ps = DiracSlash[p];
ks = DiracSlash[k];
a = DiracMatrix[\[Alpha]]; b = DiracMatrix[\[Beta]];
Tr[(pps + m).a.(ps + ks + m).b.(ps + m).b.(ps + ks + m).a/16] // Expand

Clear[pps, ps, ks, a, b];
a = DiracMatrix[\[Alpha], Dimension -> D];
b = DiracMatrix[\[Beta], Dimension -> D];
n = DiracMatrix[\[Nu], Dimension -> D];
{ps1, ps2, ps3} = Map[DiracSlash, {p1, p2, p3}];
Tr[b.a.ps1.ps2.n.b.ps2.ps3.a.ps1.n.ps3.ps1.ps2, PairCollect -> True]

Tr[b.a.ps1.ps2.n.b.ps2.ps3.a.ps1.n.ps3.ps1.ps2 /. D -> 4, 
  PairCollect -> True]

Clear[a, b, n, ps1, ps2, ps3]
DiracTrace[DiracMatrix[\[Alpha], \[Beta], \[Rho], \[Sigma]]]

Contract[%MetricTensor[\[Alpha], \[Beta]]]

% /. DiracTrace -> Tr

SetMandelstam[s, t, u, k1, k2, -p1, -p2, 0, 0, m, m]; {ks1, ps1, ps2} = 
  Map[DiracSlash, {k1, p1, p2}];
{si, ro} = Map[DiracMatrix, {\[Sigma], \[Rho]} ];
polsum1 = PolarizationSum[\[Alpha], \[Rho], k1, p1 - p2];
polsum2 = PolarizationSum[\[Beta], \[Sigma], k2, p1 - p2];
p1al = FourVector[p1, \[Alpha]];
p2be = FourVector[p2, \[Beta]];
Tr[(polsum1 polsum2 p1al p2be si). (ks1 - ps1 - m).ro.(ps1 + m).(ps2 - m), 
  Mandelstam -> {s, t, u, 2m^2}]

Tr[si.(ks1 - ps1 - m).ro.(ps1 + m).(ps2 - m), PairCollect -> True]

TrickMandelstam[
  ExpandScalarProduct[Contract[% polsum1 polsum2 p1al p2be]], {s, t, u, 
    2 m^2}]

Clear[ks1, ps1, ps2, si, ro, polsum1, polsum2, p1al, p2be]
t[n_] := t[n] = 
      Block[{gammas, calc}, 
        gammas = 
          Dot @@ Table[DiracMatrix[a[i], Dimension -> (d - 4)], {i, 1, n}]; 
        calc = Timing[Tr[gammas.gammas] // Expand]; Print["Time=", calc[[1]]];
         calc[[2]]];
t[8]


t[9]


t[10]


t[11]


Clear[t];
SUNF[a, b, c, Explicit -> True]

SUNTrace[SUNT[a].SUNT[b]]

SUNTrace[SUNT /@ (a.b.a.b)]

SUNTrace[SUNT /@ (a.b.c)]

SUNTrace[SUNF[a, b, c] SUNT /@ (a.b.c)] // SUNSimplify

SUNTrace[SUNF[a, r, s]SUNF[b, r, s]] // SUNSimplify

16SUNF[a, r, s]SUNF[b, r, s] // SUNSimplify

SUNSimplify[
  SUNF[a, b, r]SUNF[r, c, s] + SUNF[b, c, r]SUNF[r, a, s] + 
    SUNF[c, a, r]SUNF[r, b, s], SUNFJacobi -> True]

SUNT[a, b, a] // SUNSimplify

SUNF[c, a, b]SUNT[b, c] // SUNSimplify

e QuantumField[\[Gamma], {\[Mu]}] . QuantumField[\[Psi]\&_] . 
      DiracMatrix[\[Mu]] . QuantumField[\[Psi]]

FeynRule[%, {QuantumField[\[Psi]][
        p1], QuantumField[\[Psi]\&_][
        p2], QuantumField[\[Gamma], {\[Mu]3}][p3]}]

aa = QuantumField[PartialD[LorentzIndex[\[Beta]]], GaugeField, 
      LorentzIndex[\[Alpha]], SUNIndex[a]].QuantumField[
      PartialD[LorentzIndex[\[Beta]]], GaugeField, LorentzIndex[\[Alpha]], 
      SUNIndex[a]]

dd = FunctionalD[
    aa, {QuantumField[GaugeField, {\[Mu]1}, {i1}], 
      QuantumField[GaugeField, {\[Mu]2}, {i2}]}]

dd // StandardForm

dd.QuantumField[\[Psi], {\[Mu]1}, {i1}] // ExpandPartialD

% // StandardForm

FunctionalD[aa, {QuantumField[GaugeField, {\[Mu]1}, {i1}][p1], 
    QuantumField[GaugeField, {\[Mu]2}, {i2}][p2]}]

FeynRule[aa, {QuantumField[GaugeField, {\[Mu]1}, {i1}][p1], 
    QuantumField[GaugeField, {\[Mu]2}, {i2}][p2]}]

Clear[aa, dd];
Lagrangian["QCD"]

Lagrangian["QCD"] /. 
    FieldStrength[x__] :> FieldStrength[x, Explicit -> True] // DotExpand

FeynRule[Lagrangian["QCD"], {QuantumField[GaugeField, {\[Mu]1}, {i1}][p1], 
    QuantumField[GaugeField, {\[Mu]2}, {i2}][p2], 
    QuantumField[GaugeField, {\[Mu]3}, {i3}][p3], 
    QuantumField[GaugeField, {\[Mu]4}, {i4}][p4]}]


B0[s, mw2, mz2]

D[%, s]

SetOptions[A0, A0ToB0 -> False];
SetOptions[{B1, B00, B11}, BReduce -> False];
PaVe[0, {}, {m02}]

PaVe[0, {pp}, {m02, m12}]

PaVe[1, {pp}, {m12, m22}]

PaVe[0, 0, {pp}, {m02, m12}]

PaVe[1, 1, {p10}, {m12, m22}]

PaVe[0, {p10, p12, p20}, {m12, m22, m32}]

PaVe[0, {p10, p12, p23, p30, p20, p13}, {m12, m22, m32, m42}]

B1[pp, m12, m22]

B1[pp, SmallVariable[me2], m22]

B1[SmallVariable[me2], SmallVariable[me2], 0]

% // StandardForm

PaVeReduce[
  PaVe[2, {SmallVariable[me2], mw2, t}, {SmallVariable[me2], 0, mw2}]]

c12 = PaVeReduce[PaVe[1, 2, {s, m2, m2}, {m2, m2, w2}], IsolateNames -> c]

c[31]

SetOptions[B1, BReduce -> True];
Simplify /@ Collect[FixedPoint[ReleaseHold, c12], {_B0, _C0}]

d122 = PaVeReduce[
    PaVe[1, 2, 
      2, {SmallVariable[me2], mw2, mw2, SmallVariable[me2], s, t}, {0, 
        SmallVariable[me2], 0, SmallVariable[me2]}], 
    Mandelstam -> {s, t, u, 2mw2}, IsolateNames -> f]

Write2["d122.for", d122res = d122, FormatType -> FortranForm];
!! "d122.for"

DeleteFile["d122.for"];
a0 = -I/Pi^ 2FeynAmpDenominator[PropagatorDenominator[q, m]]

OneLoop[q, a0]

SetOptions[A0, A0ToB0 -> True];
OneLoop[q, a0]

SetOptions[OneLoop, Factoring -> True];
OneLoop[q, (I el^2)/(16Pi^4)/(1 - D)*
          FeynAmpDenominator[PropagatorDenominator[q, mf], 
            PropagatorDenominator[q - k, mf]]*
          DiracTrace[(mf + DiracSlash[q - k]).DiracMatrix[
                mu].(mf + DiracSlash[q]).DiracMatrix[mu]]] /. 
      ScalarProduct[k, k] -> k2 /. mf^2 -> mf2 // Simplify

gc[i_] := g[i, "-"]DiracMatrix[7] +
g[i, "+"]DiracMatrix[6];
ScalarProduct[p1, p1] = 0;
ScalarProduct[p2, p2] = 0;
ScalarProduct[p1, p2] = k2/2;
MakeBoxes[g[i_, j_], TraditionalForm] := SubsuperscriptBox[g, i, j];
wff1 = OneLoop[q, 
        I/(2Pi)^ 4FeynAmpDenominator[PropagatorDenominator[q, m], 
            PropagatorDenominator[q + p1], PropagatorDenominator[q - p2]]*
          Spinor[p1].DiracMatrix[nu].gc[1].DiracSlash[q + p1].DiracSlash[
              Polarization[k]].gc[3].DiracSlash[q - p2].DiracMatrix[nu].gc[
              2].Spinor[p2]] /. (m^ 2) -> m2 // Simplify

wff1a = PaVeReduce[wff1] // Simplify

var = Select[Variables[wff1a], (Head[#] === StandardMatrixElement) &]

Set @@ {var, {ma[1], ma[2]}};
Write2["wff1a.for", 
    vert = wff1a /. g[i_, "+"] -> gp[i] /. g[i_, "-"] -> gm[i], 
    FormatType -> FortranForm];
!! "wff1a.for"

DeleteFile["wff1a.for"];
Clear[gc, g, wff1, wff1a, vert];
r = DiracMatrix[6]; l = DiracMatrix[7];
ScalarProduct[q1, q1] = ScalarProduct[q2, q2] = 0;
ScalarProduct[q1, q2] = k2/2;
SetOptions[OneLoop, Factoring -> True, FormatType -> FortranForm, 
    ReduceToScalars -> True, WriteOut -> True, 
    FinalSubstitutions -> {g[i_, "+"] -> gp[i], g[i_, "-"] -> gm[i], 
        StandardMatrixElement -> ma}];
mt = MetricTensor;
fv = FourVector;
feynden[x : {_, _} ..] := FeynAmpDenominator @@
Map[Apply[PropagatorDenominator, #] &, {x}];
SetStandardMatrixElements[{(Spinor[q1].DiracSlash[Polarization[k]].r.Spinor[
              q2]) -> {1}, (Spinor[q1].DiracSlash[Polarization[k]].l.Spinor[
              q2]) -> {2}}];
OneLoop[wff2, q, 
      I/(2Pi)^ 4*feynden[{q, 0}, {q + q1, m1}, {q - q2, m2}]*
        Spinor[q1].DiracMatrix[
            nu].(g[1, "-"]l + g[1, "+"]r).DiracSlash[-q].DiracMatrix[
            ro].(g[2, "-"]l + g[2, "+"]r).Spinor[q2]*
        g3(mt[ro, mu]fv[q1 + 2q2 - q, nu] - mt[mu, nu]fv[2q1 + q2 + q, ro] + 
            mt[nu, ro]fv[2q + q1 - q2, mu])*
        PolarizationVector[k, mu]] /. (m1^2) -> m12 /. (m2^2) -> m22

!! "wff2.for"

DeleteFile["wff2.for"];
DeleteFile /@ FileNames["PaVe*"];
Clear[r, l, mt, fv, feynden, wff2];

Quit[];
<< FeynArts`;



tops = CreateTopologies[1, 2 -> 2, Adjacencies -> {3}, 
      ExcludeTopologies -> {SelfEnergies, WFCorrections, Tadpoles, 
          Boxes[3]}];


Paint[tops, ColumnsXRows -> {3, 1}];


inserttops = 
    InsertFields[tops, {F[2, {1}], -F[2, {1}]} -> {V[2], S[1]}, Model -> "SM",
       GenericModel -> "Lorentz", InsertionLevel -> Particles, 
      ExcludeFieldPoints -> {FieldPoint[F, -F, S]}];


graphs = Paint[inserttops, PaintLevel -> {Particles}, AutoEdit -> False, 
      SheetHeader -> False, Numbering -> False, ColumnsXRows -> {3, 2}];


eezhb = CreateFeynAmp[inserttops, AmplitudeLevel -> Particles];


PickLevel[Particles][eezhb] // ToFA1Conventions >> "eezhb.amp";

Quit[];
<< HighEnergyPhysics`FeynCalc`

  SetMandelstam[s, t, u, p1, p2, -k1, -k2, SmallVariable[ME], 
    SmallVariable[ME], MZ, MH];


SetOptions[PaVeOrder, PaVeOrderList -> {{s, t}, {s, u}, {t, u}}];


SetOptions[OneLoop, Mandelstam -> {s, t, u, MH^2 + MZ^2}, 
    Prefactor -> 1/ALPHA2, 
    InitialSubstitutions -> {k2 -> p1 + p2 - k1, CW -> MW/MZ, 
        EL -> Sqrt[4 Pi Sqrt[ALPHA2]]}, SmallVariables -> {ME, ME2}];


SetOptions[OneLoopSum, Prefactor -> 2 ALP4PI FLUFAC, 
    Mandelstam -> {s, t, u, MH^2 + MZ^2}, 
    FinalSubstitutions -> {SW -> Sqrt[SW2], ME -> Sqrt[ME2], MW -> Sqrt[MW2], 
        MZ -> Sqrt[MZ2], MH -> Sqrt[MH2], 
        ME2^n_ :> ME^(2 n) /; Head[n] =!= Integer, 
        MZ2^n_ :> MZ^(2 n) /; Head[n] =!= Integer, 
        MW2^n_ :> MW^(2 n) /; Head[n] =!= Integer, 
        MH2^n_ :> MH^(2 n) /; Head[n] =!= Integer, 
        SW2^n_ :> SW^(2 n) /; Head[n] =!= Integer, 
        StandardMatrixElement -> MBM}, WriteOutPaVe -> ""];


SetStandardMatrixElements[{Spinor[p1].DiracSlash[
            Conjugate[Polarization[k1]]].ChiralityProjector[+1].Spinor[
            p2] -> {0, 1}, 
      Spinor[p1].DiracSlash[
            Conjugate[Polarization[k1]]].ChiralityProjector[-1].Spinor[
            p2] -> {0, 2}, 
      ScalarProduct[Conjugate[Polarization[k1]], p1]*
          Spinor[p1].DiracSlash[k1].ChiralityProjector[+1].Spinor[p2] -> {1, 
          1}, 
      ScalarProduct[Conjugate[Polarization[k1]], p1]*
          Spinor[p1].DiracSlash[k1].ChiralityProjector[-1].Spinor[p2] -> {1, 
          2}, ScalarProduct[Conjugate[Polarization[k1]], p2]*
          Spinor[p1].DiracSlash[k1].ChiralityProjector[+1].Spinor[p2] -> {2, 
          1}, ScalarProduct[Conjugate[Polarization[k1]], p2]*
          Spinor[p1].DiracSlash[k1].ChiralityProjector[-1].Spinor[p2] -> {2, 
          2}}, {k2 -> (p1 + p2 - k1)}];


eezhamp = << eezhb.amp;


eezhboxes = OneLoopSum[eezhamp, CombineGraphs -> {1, 2, 3, 4, 5, 6}];


Write2["eezhb.m", EEZHBOXES = FRH[eezhboxes]];


Write2["eezhb.s", EEZHBOXES = eezhboxes];


Write2["eezhb.for", EEZHBOXES = eezhboxes, FormatType -> FortranForm];
<< "eezhb.m";
EEZHBOXES // LeafCount

eezhoxes = 
    Drop[EEZHBOXES, -1]((WriteString["stdout", "."]; 
              Collect[# // Expand, {_MBM, _D0, _C0, _B0,}, 
                If[FreeQ[#, _MBM | _D0 | _C0 | _B0 | _PaVe], 
                    FullSimplify[#], #] &]) & /@ EEZHBOXES[[-1]]);

eezhoxes = 
    Collect[EEZHBOXES // Expand, {ALP4PI, FLUFAC, _MBM, _D0, _C0, _B0}, 
      If[FreeQ[#, _MBM | _D0 | _C0 | _B0 | _PaVe], FullSimplify[#], #] &];
eezhoxes // LeafCount

DeleteFile /@ FileNames["eezhb*"];
DeleteFile /@ FileNames["PaVe*"];

Quit[];
<< HighEnergyPhysics`FeynArts`;





tops = CreateTopologies[0, 2 -> 2, Adjacencies -> {3}, 
      ExcludeTopologies -> {SelfEnergies, WFCorrections}];


Paint[tops, ColumnsXRows -> {3, 1}, AutoEdit -> False];


inserttops = 
    InsertFields[tops, {F[1, {1}], F[1, {1}]} -> {F[1, {1}], F[1, {1}]}, 
      Model -> "QED", GenericModel -> "QED", InsertionLevel -> Particles];


treegraphs = 
    Paint[inserttops, PaintLevel -> {Particles}, AutoEdit -> False, 
      SheetHeader -> False, Numbering -> False, ColumnsXRows -> {2, 1}];


amps = CreateFeynAmp[inserttops, AmplitudeLevel -> Particles];


PickLevel[Classes][amps] // ToFA1Conventions >> "moelleramps.m"

Quit[];
<< HighEnergyPhysics`FeynCalc`


  SetMandelstam[s, t, u, p1, p2, -k1, -k2, ME, ME, ME, ME];


amps = << moelleramps.m;


amp = (OneLoopSum[amps, CombineGraphs -> {1, 2}] // FRH) /. D :> Sequence[] //
     Contract


squaredamp = 
  FermionSpinSum[amp ComplexConjugate[amp /. li2 -> li1] // Expand] /. 
          DiracTrace -> Tr // DiracSimplify // 
      TrickMandelstam[#, {s, t, u, 4ME^2}] & // Simplify
squaredamp1 = 
  squaredamp // Contract // PropagatorDenominatorExplicit // Simplify


kinfac = 1/(64 \[Pi]^2s);


dcrosssection = 1/4*kinfac*squaredamp1 // Simplify

dc = dcrosssection /. u -> 4ME^2 - s - t /. {s -> 4 \[Omega]^2, 
          t -> -2 q2(1 - Sqrt[1 - sin[\[Theta]]^2])} /. 
      q2 -> \[Omega]^2 - ME^2 // Simplify


DeleteFile["moelleramps.m"];

Quit[];
$LoadPhi = True;
$LoadFeynArts = True;
$Configuration = "ChPTVirtualPhotons2";
$Lagrangians = {"ChPTVirtualPhotons2"[2], "ChPTVirtualPhotons2"[4]};
Get["HighEnergyPhysics`FeynCalc`"];

Lagrangian[ChPTVirtualPhotons2[2]]

FieldDerivative[QuantumField[p][x], x, LorentzIndex[\[Mu]]]

% // StandardForm

IsoDot[IsoVector[a], IsoVector[b]]

UVector[\[CapitalPsi]\&_] . 
    NM[IsoDot[IsoVector[a], IsoVector[UMatrix[b]]], Adjoint[UMatrix[c]]] . 
    UVector[\[CapitalPsi]]

IsoDot[IsoVector[a], IsoVector[b]] // WriteOutIsoVectors

UMatrix[a] // WriteOutUMatrices

UMatrix[UGenerator[1]] // WriteOutUMatrices

IsoVector[UMatrix[UGenerator[]]]

% // WriteOutIsoVectors

% // WriteOutUMatrices

IsoDot[IsoVector[a], IsoVector[b]] // IsoIndicesSupply

id = IsoDot;
iv = IsoVector;
ug = IsoVector[UMatrix[UGenerator[]]];
UTrace[NM[id[iv[a], ug] + id[iv[d], ug], id[iv[b], ug], id[iv[c], ug]] - 
    NM[id[iv[b], ug], id[iv[c], ug], id[iv[a], ug]]]

% // NMExpand // CycleUTraces

% // ExpandU

% // IsoIndicesSupply // CommutatorReduce

k(IsoDot[IsoCross[IsoVector[a], IsoVector[b]], 
        IsoCross[IsoVector[a], IsoVector[b]]] + 
      IsoDot[IsoCross[IsoVector[a], IsoVector[c]], 
        IsoCross[IsoVector[a], IsoVector[c]]])

% // IsoIndicesSupply

% // CommutatorReduce

% // IndicesCleanup

% // SUNReduce // Simplify

SUNReduce[%, FullReduce -> True] // Simplify

% // IndicesCleanup // Simplify

Clear[a, b, c, d, id, iv, ug];
UTrace[NM[UChiMatrix, Adjoint[MM]]]

UTrace[NM[UChiMatrix[x], Adjoint[MM[x, ExpansionOrder -> 2]]]]

CovariantFieldDerivative[QuantumField[Particle[Pion]][x]^2, x, 
    LorentzIndex[\[Mu]]] // CommutatorReduce

CovariantFieldDerivative[
    IsoDot[IsoVector[QuantumField[Particle[Pion]]][x], 
      IsoVector[UMatrix[UGenerator[]]]], x, 
    LorentzIndex[\[Mu]]] // CommutatorReduce

UVector[DiracBar[QuantumField[Particle[Nucleon]]]][x].DiracGamma[
    LorentzIndex[\[Mu]]].NM[UMatrix[SMM][x], 
    FieldDerivative[Adjoint[UMatrix[SMM][x]], x, LorentzIndex[\[Mu]]]].
UVector[QuantumField[Particle[Nucleon]]][x]

ArgumentsSupply[Lagrangian[ChPTVirtualPhotons2[2]], x, 
  RenormalizationState[0], ExpansionOrder -> 0, DropOrder -> 0]

l = ArgumentsSupply[Lagrangian[ChPTVirtualPhotons2[2]], x, 
      RenormalizationState[0], DiagonalToU -> True, ExpansionOrder -> 2, 
      DropOrder -> 2];
DiscardTerms[l, 
          Retain -> {ParticleField[Pion , RenormalizationState[0]] -> 2, 
              ParticleField[Photon, RenormalizationState[0]] -> 2}, 
          CommutatorReduce -> True, Method -> Expand] /. $Substitutions // 
      NMExpand // CycleUTraces // Expand

ExpandU[%, CommutatorReduce -> True] // Simplify

ll = % // IsoIndicesSupply // SUNReduce // IndicesCleanup // 
      CommutatorReduce // Simplify

fields = {QuantumField[Particle[PseudoScalar[2], RenormalizationState[0]], 
        SUNIndex[I1]][p1], 
    QuantumField[Particle[PseudoScalar[2], RenormalizationState[0]], 
        SUNIndex[I2]][p2], 
    QuantumField[Particle[Vector[1], RenormalizationState[0]], 
        LorentzIndex[\[Mu]3]][p3], 
    QuantumField[Particle[Vector[1], RenormalizationState[0]], 
        LorentzIndex[\[Mu]4]][p4]}

m = FeynRule[ll, fields] // SUNReduce[#, FullReduce -> True] & // 
        IndicesCleanup // CommutatorReduce // Simplify

mfa = MomentaCollect[m // Expand, PerturbationOrder -> 2]

gencoup = GenericCoupling[mfa]

classcoup = ClassesCoupling[mfa] // Together; classcoup // StandardForm

CheckF[gencoup, 
    XName[PhiModel -> ChPTVirtualPhotons2, 
        VertexFields -> {Pion[0], Pion[0], Photon[0], Photon[0]}, 
        PerturbationOrder -> 2] <> ".Gen"];
CheckF[classcoup, 
    XName[PhiModel -> ChPTVirtualPhotons2, 
        VertexFields -> {Pion[0], Pion[0], Photon[0], Photon[0]}, 
        PerturbationOrder -> 2] <> ".Mod"];
Clear[l, ll, m, mfa];
Amplitude["ChPTVirtualPhotons2P20P20V10V10o2"]



mesonstop = CreateTopologies[1, 1 -> 1, Adjacencies -> {3, 4}];
mesontreeinsert = 
    InsertFields[mesonstop, {Pion[0, {i1}]} -> {Pion[0, {i2}]}, 
      Model -> "Automatic", GenericModel -> "Automatic", 
      InsertionLevel -> Classes];
mesontreegraph = 
    Paint[mesontreeinsert, PaintLevel -> {Classes}, AutoEdit -> False, 
      SheetHeader -> False, Numbering -> False, ColumnsXRows -> {2, 1}];
$ConstantIsoIndices = Union[$ConstantIsoIndices, {i1, i2, I1}];
amplFC = CreateFCAmp[mesontreeinsert, AmplitudeLevel -> Classes, 
          EqualMasses -> False] /. i2 -> i1 // SUNReduce // Simplify
ampreduced = OneLoop[q1, #, Sum -> Explicit] & /@ amplFC;
ampsimple = ampreduced // IsoToChargedMasses // Simplify

CheckF // Options // StandardForm


Collect2[(b - 2) d f[x] + f[x] + c p[y] + p[y], {f, y}]

Combine[(a - b) (c - d)/e + g]

test = (a - b) x + (b - a) y

Map[Factor, test]

Map[Factor2, test]

Factor[test]

Factor2[test]

Isolate[a + b]

test = Isolate[(a + b) f + (c + d) f + e, f]

FullForm[test]

{KK[1], test, ReleaseHold[test]}

Isolate[a[z](b + c (y + z)) + d[z] (y + z), {a, d}, IsolateNames -> g]

und[x__] := Isolate[Plus[x], IsolateNames -> h] /;
FreeQ2[{x}, {a, d}];
(a[z] (b + c (y + z)) + d[z] (y + z)) /.
Plus -> und /. und -> Plus

ReleaseHold[%]

Isolate[a - b - c - d - e, IsolateNames -> l, IsolateSplit -> 15]

{l[2], l[1]}

Clear[test];
FreeQ2[w^2 + m^2 B0[pp, m1, m2], {w, B0}]

NumericalFactor[-137x]

PartitHead[f[m] (s - u), f]

PartitHead[s^2 + m^2 - f[m], f]

tpol = z + Isolate[Isolate[2 Pi I + f[x] (a - b), f]]

Write2["test.m", test = tpol];
!! test.m

Write2["test.for", test = tpol, FormatType -> FortranForm];
!! test.for

DeleteFile /@ {"test.m", "test.for"};
Contract[LeviCivita[\[Mu], \[Nu], \[Rho], \[Sigma]] FourVector[
      p + q, \[Sigma]], EpsContract -> False]

EpsEvaluate[%]

EpsChisholm[
  LeviCivita[\[Alpha], \[Beta], \[Gamma], \[Delta]] DiracMatrix[\[Mu], \
\[Alpha]] // FCI]

LeviCivita[\[Mu], \[Nu], \[Rho], \[Sigma]] // StandardForm

Contract[% FourVector[p, \[Sigma]]] // StandardForm

TrickMandelstam[(s + t - u) (2 mw2 - t - u),
  {s, t, u, 2 mw2}]

TrickMandelstam[m^ 2 s - s^2 + m^2 t - s t + m^2 u - s u,
  {s, t, u, 2 m^2}]

PolarizationSum[\[Mu], \[Nu], k]

PolarizationSum[\[Mu], \[Nu], k, p1 - p2]

PaVeOrder[D0[me2, me2, mw2, mw2, t, s, me2, 0, me2, 0], 
  PaVeOrderList -> {me2, me2, 0, 0}]

PaVeOrder[D0[a, b, c, d, e, f, m12, m22, m32, m42], PaVeOrderList -> {f, e}]

PaVeOrder[
  D0[a, b, c, d, e, f, m12, m22, m32, m42] +  
    D0[me2, me2, mw2, mw2, t, s, me2, 0, me2, 0], 
  PaVeOrderList -> { {me2, me2, 0, 0}, {f, e}}]



?FortranForm

?ChisholmSpinor

DiracSigma[GA[\[Alpha]], GA[\[Beta]]] // FCI

DiracSigmaExplicit[%]

% // StandardForm


?EpsEvaluate

?*Fierz*

?*Schouten*

SUNReduce[
  SUNF[a, b, r]SUNF[r, c, s] + SUNF[b, c, r]SUNF[r, a, s] + 
      SUNF[c, a, r]SUNF[r, b, s] // FCI]

