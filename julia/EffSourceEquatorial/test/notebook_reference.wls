nb = Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "notebooks", "ScalarSingularField.nb"}]];
Do[Quiet[ToExpression[b, StandardForm]], {b, Cases[nb, Cell[b_BoxData, "Input", ___] :> b, Infinity]}];
digits[x_] := ToString[Round[N[x, 90] 10^60]] <> "e-60";
rat[x_] := Rationalize[N[x, 100], 10^-95];
particles = {{"circular", 1/2, 10, 0, 10, 0}, {"A_peri", 1/2, 10, 3/10, 100/13, 0}, {"A_out9", 1/2, 10, 3/10, 9, 1},
   {"A_in9", 1/2, 10, 3/10, 9, -1}, {"A_out12", 1/2, 10, 3/10, 12, 1}, {"A_apo", 1/2, 10, 3/10, 100/7, 0},
   {"B_peri", 9/10, 7, 1/2, 14/3, 0}, {"B_out7", 9/10, 7, 1/2, 7, 1}, {"B_in10", 9/10, 7, 1/2, 10, -1}};
offsets = {{1/10, 1/100, 1/50}, {1/2, 1/20, 3/10}, {-3/10, 1/10, -1/5}, {2, 3/10, 1/2}, {1/1000, 1/1000, 1/1000}, {1/10^6, 1/10^6, 1/10^6}};
modeParticles = {"A_peri", "A_out9", "A_in9", "B_in10"};
out = OpenWrite[FileNameJoin[{DirectoryName[$InputFileName], "notebook_reference.txt"}]];
Do[
  {label, aa, p, e, r, sgn} = part;
  g = KerrEquatorialConstants[1, aa, p, e, r, 100];
  tw = TwistedPuncture[{a -> aa, M -> 1, r1 -> r, e1 -> rat[e1 /. g], l1 -> rat[l1 /. g], ur1 -> If[sgn == 0, 0, sgn rat[ur1 /. g]]}];
  WriteString[out, "particle ", label, " ", StringRiffle[ToString[#, InputForm] & /@ {aa, p, e, r, sgn}, " "], "\n"];
  f = TwistedPhiS[tw, dr, dth, dph];
  Do[WriteString[out, "real ", label, " ", StringRiffle[Join[ToString[#, InputForm] & /@ o,
       digits /@ ({f, D[f, dr], D[f, dth], D[f, dph]} /. Thread[{dr, dth, dph} -> o])], " "], "\n"], {o, offsets}];
  If[MemberQ[modeParticles, label],
   Do[With[{gg = TwistedPhiS[tw, o[[1]], o[[2]], ph], pk = tw["c"] o[[1]]},
      val = NIntegrate[gg Exp[-I m ph], {ph, pk - Pi, pk, pk + Pi}, WorkingPrecision -> 70, PrecisionGoal -> 50, MaxRecursion -> 40];
      WriteString[out, "mode ", label, " ", m, " ", StringRiffle[ToString[#, InputForm] & /@ o, " "], " ", digits[Re[val]], " ", digits[Im[val]], "\n"]],
     {o, {{1/2, 1/20}, {2, 3/10}, {-3/10, 1/10}, {1/100, 1/100}, {1/1000, 1/1000}}}, {m, {0, 2, 10, 20}}]],
  {part, particles}];
Close[out];
