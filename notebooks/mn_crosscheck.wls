dir = DirectoryName[$InputFileName];
nbfile = FileNameJoin[{dir, "ScalarSingularField.nb"}];
mode = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "run"];
out = Which[mode === "selftest", FileNameJoin[{$TemporaryDirectory, "mn_crosscheck_selftest.txt"}], Length[$ScriptCommandLine] >= 3, $ScriptCommandLine[[3]], True, FileNameJoin[{dir, "mn_crosscheck.txt"}]];
loadDefs[] := With[{nb = Get[nbfile]}, Do[Quiet[ToExpression[b, StandardForm]], {b, Cases[nb, Cell[b_BoxData, "Input", ___] :> b, Infinity]}]];
rules[] := (RuleParameters["xcheck16"] = {16, 1/16, 16}; RuleParameters["xcheck20"] = {20, 1/24, 20});
loadDefs[];
rules[];
ms = Range[0, 20]; ns = Range[-20, 20];
interior = {{"11", "1/100"}, {"11", "1/1000"}, {"11", "1/10000"}, {"9", "1/10"}, {"13", "1/1000"}, {"13/2", "1/100"}, {"33/2", "1/5"}};
faces = {{"34099538200381/10000000000000", "1/10"}, {"9274383314342119/500000000000000", "1/100"}, {"11", "ArcSin[1/2]"}};
label[pt_] := "r=" <> pt[[1]] <> ",dth=" <> pt[[2]];
orbit[prec_] := orbit[prec] = KerrEquatorialOrbit[1/2, 10, 3/10, prec];
spec[pt_] := Which[MemberQ[interior, pt], {110, {"box"}}, pt[[2]] === "ArcSin[1/2]", {110, {"Phi", "r", "th"}}, True, {50, {"Phi", "r", "th"}}];
fmtp[x_, scale_] := Module[{k = Min[20, Floor[Precision[x]]], d, e},
  If[! TrueQ[Accuracy[x] + Log10[scale] >= 15], Throw[{"precision", x}, "fmtp"]];
  If[k < 1 || x == 0, "0", {d, e} = RealDigits[x, 10, k]; If[x < 0, "-", ""] <> ToString[First[d]] <> "." <> StringJoin[ToString /@ Rest[d]] <> "E" <> ToString[e - 1]]];
pointDigits[res_, qs_] := Min[Table[Module[{v = res["values"][q], d = res["digits"][q], scale, diff},
    scale = Max[Abs[Flatten[v]]];
    diff = Max[Table[If[d[[mi]] === Infinity, 0, 10^-d[[mi]] Max[Abs[v[[mi]]]]], {mi, Length[v]}]];
    If[diff == 0, Infinity, -Log10[diff/scale]]], {q, qs}]];
pointLines[pt_, res_, qs_] := Module[{lab = label[pt], col, scale, pd = pointDigits[res, qs]},
  If[! MatchQ[pd, _Real | _Integer | Infinity], Throw[{"digits", pd}, "fmtp"]];
  scale = Association[Table[q -> Max[Abs[Flatten[res["values"][q]]]], {q, qs}]];
  col[q_, mi_, ni_] := If[MemberQ[qs, q], fmtp[res["values"][q][[mi, ni]], scale[q]], "-"];
  Join[{"point " <> lab <> " " <> pt[[1]] <> " " <> pt[[2]] <> " -"},
    Flatten[Table["mn " <> lab <> " " <> ToString[ms[[mi]]] <> " " <> ToString[ns[[ni]]] <> " " <> StringRiffle[{col["Phi", mi, ni], col["r", mi, ni], col["th", mi, ni], col["box", mi, ni], "-"}, " "], {mi, Length[ms]}, {ni, Length[ns]}]],
    Table["digits " <> lab <> " " <> ToString[ms[[mi]]] <> " " <> StringRiffle[Table[If[MemberQ[qs, q], ToString[N[res["digits"][q][[mi]], 4]], "-"], {q, {"Phi", "r", "th", "box"}}], " "] <> " -", {mi, Length[ms]}],
    {"digits " <> lab <> " all " <> ToString[N[pd, 4]], "done " <> lab}]];
writeLines[lines_] := With[{st = OpenAppend[out]}, WriteString[st, StringRiffle[lines, "\n"] <> "\n"]; Close[st]];
header[] := With[{orb = orbit[110]}, {"orbit 1/2 10 3/10", "constants " <> fmtp[orb["E"], orb["E"]] <> " " <> fmtp[orb["L"], orb["L"]], "freq " <> StringRiffle[fmtp[#, #] & /@ {orb["OmegaR"], orb["OmegaPhi"], orb["TR"], orb["LambdaR"]}, " "]}];
prepare[] := Module[{lines, last, keep},
  If[! FileExistsQ[out], writeLines[header[]]; Return[{}, Module]];
  lines = ReadList[out, String];
  last = Max[0, Flatten[Position[lines, s_String /; StringStartsQ[s, "done "]]]];
  keep = If[last == 0, header[], Take[lines, last]];
  With[{st = OpenWrite[out]}, WriteString[st, StringRiffle[keep, "\n"] <> "\n"]; Close[st]];
  StringDrop[#, 5] & /@ Select[keep, StringStartsQ[#, "done "] &]];
safeMap[f_, list_] := Module[{res = ParallelMap[f, list, Method -> "FinestGrained"], bad},
  bad = Flatten[Position[res, x_ /; ! AssociationQ[x], {1}, Heads -> False]];
  If[bad =!= {}, res[[bad]] = f /@ list[[bad]]]; res];
computePoint[pt_] := Module[{prec, qs, orb, r, th, res, t},
  {prec, qs} = spec[pt]; orb = orbit[prec];
  r = N[ToExpression[pt[[1]]], prec]; th = Pi/2 + N[ToExpression[pt[[2]]], prec];
  {t, res} = AbsoluteTiming[MNModes[orb, r, th, ms, ns, qs, "Rules" -> {"xcheck16", "xcheck20"}, "Map" -> safeMap]];
  <|"res" -> res, "qs" -> qs, "digits" -> pointDigits[res, qs], "seconds" -> t|>];
stepDigits[] := Module[{orb = orbit[110], r = N[11, 110], dth = N[1/10000, 110], lc, pa, w, src, scale},
  lc = CrossingTime[orb, r]; pa = PunctureAt[orb, lc];
  w = Sqrt[pa["alpha02"]/pa["alpha20"]] dth/Abs[orb["dr"][lc]];
  src[l_, h_] := Module[{pt = ParticleAt[orb, l]}, EffectiveSource4D[orb, l, N[h, 110]][r, Pi/2 + dth, pt["phip"] + PunctureAt[orb, l]["c"] (r - pt["rp"])]];
  scale = Max[Abs[src[#, 10^-12] & /@ Join[N[Range[1, 63]/64 orb["LambdaR"]/2, 110], {lc, lc + w}]]];
  Table[N[-Log10[Abs[src[l, 10^-12] - src[l, 5 10^-13]]/scale], 4], {l, {lc, lc + w}}]];
report[name_, ok_, val_] := (Print["check ", name, " ", If[TrueQ[ok], "PASS", "FAIL"], " ", val]; TrueQ[ok]);
runPoints[pts_] := Module[{done = prepare[], ok = True, c, lines},
  Do[If[! MemberQ[done, label[pt]],
      c = computePoint[pt];
      lines = Catch[pointLines[pt, c["res"], c["qs"]], "fmtp"];
      If[MatchQ[lines, {"precision" | "digits", _}],
        ok = report[lines[[1]] <> ":" <> label[pt], False, lines[[2]]],
        writeLines[lines];
        ok = report["digits:" <> label[pt], c["digits"] >= 14, {N[c["digits"], 4], Round[c["seconds"]]}] && ok]], {pt, pts}];
  ok];
If[mode === "selftest",
  Module[{ok = True, syn, lines},
    ok = report["rules", RuleParameters["xcheck16"] === {16, 1/16, 16} && RuleParameters["xcheck20"] === {20, 1/24, 20} && RuleParameters["primary"] === {28, 1/32, 28} && RuleParameters["check"] === {32, 1/48, 32}, "-"] && ok;
    ok = report["format", StringLength[StringDelete[First[StringSplit[fmtp[N[Pi/7, 32], 1], "E"]], {"-", "."}]] == 20 && Abs[ToExpression[StringReplace[fmtp[N[Pi/7, 32], 1], "E" -> "*^"]] - Pi/7] < 10^-19 && StringMatchQ[fmtp[N[1/3, 17], 1], "3.3333333333333333" ~~ "E-1"] && fmtp[SetAccuracy[77/10^14, 25], 1] === "7.70000000000E-13" && fmtp[0``32, 1] === "0" && MatchQ[Catch[fmtp[N[1, 10], 1], "fmtp"], {"precision", _}] && MatchQ[Catch[fmtp[SetAccuracy[77/10^14, 25], 10^-11], "fmtp"], {"precision", _}], fmtp[SetAccuracy[77/10^14, 25], 1]] && ok;
    ok = report["label", label[{"11", "ArcSin[1/2]"}] === "r=11,dth=ArcSin[1/2]" && Length[Join[interior, faces]] == 10 && spec[{"11", "1/10000"}] === {110, {"box"}} && spec[{"11", "ArcSin[1/2]"}] === {110, {"Phi", "r", "th"}} && spec[faces[[1]]] === {50, {"Phi", "r", "th"}}, "-"] && ok;
    syn = <|"values" -> <|"box" -> {{1.`30, 2.`30}, {0.5`30, 0.1`30}}|>, "digits" -> <|"box" -> {3, Infinity}|>|>;
    ok = report["point-digits", Abs[pointDigits[syn, {"box"}] - 3] < 10^-10, pointDigits[syn, {"box"}]] && ok;
    syn = <|"values" -> <|"box" -> Table[N[1/(1 + m + n^2), 30], {m, ms}, {n, ns}]|>, "digits" -> <|"box" -> Table[16, {m, ms}]|>|>;
    lines = pointLines[{"11", "1/100"}, syn, {"box"}];
    ok = report["lines", Length[lines] == 1 + 21 41 + 21 + 2 && lines[[1]] === "point r=11,dth=1/100 11 1/100 -" && StringMatchQ[lines[[2]], "mn r=11,dth=1/100 0 -20 - - - " ~~ __ ~~ " -"] && Last[lines] === "done r=11,dth=1/100", Length[lines]] && ok;
    syn = <|"values" -> <|"box" -> Table[If[m == 3 && n == 0, N[1/3, 10], N[1/(1 + m + n^2), 30]], {m, ms}, {n, ns}]|>, "digits" -> <|"box" -> Table[16, {m, ms}]|>|>;
    ok = report["precision-guard", MatchQ[Catch[pointLines[{"11", "1/100"}, syn, {"box"}], "fmtp"], {"precision", _}], "-"] && ok;
    syn = <|"values" -> <|"box" -> Table[N[1/(1 + m + n^2), 30], {m, ms}, {n, ns}]|>, "digits" -> <|"box" -> Table[If[m == 5, Indeterminate, 16], {m, ms}]|>|>;
    lines = Catch[pointLines[{"11", "1/100"}, syn, {"box"}], "fmtp"];
    syn = <|"values" -> <|"box" -> Table[N[1/(1 + m + n^2), 30], {m, ms}, {n, ns}]|>, "digits" -> <|"box" -> Table[Infinity, {m, ms}]|>|>;
    ok = report["digits-guard", MatchQ[lines, {"digits", _}] && Last[pointLines[{"11", "1/100"}, syn, {"box"}]] === "done r=11,dth=1/100", lines[[1]]] && ok;
    If[FileExistsQ[out], DeleteFile[out]];
    writeLines[Join[header[], {"point a 1 1 -", "done a", "point b 1 1 -", "mn b 0 0 - - - 1.0E0 -"}]];
    ok = report["resume", prepare[] === {"a"} && Last[ReadList[out, String]] === "done a", "-"] && ok;
    DeleteFile[out];
    Exit[If[ok, 0, 1]]]];
LaunchKernels[8];
DistributeDefinitions[nbfile, loadDefs, rules];
ParallelEvaluate[loadDefs[]; rules[]];
step = stepDigits[];
stepOk = report["step", Min[step] >= 16, step];
If[mode === "time",
  Module[{t0, t1, ok},
    t0 = AbsoluteTime[]; ok = runPoints[{{"11", "1/10000"}}]; t1 = AbsoluteTime[] - t0;
    t0 = AbsoluteTime[]; ok = runPoints[{faces[[1]]}] && ok; t0 = AbsoluteTime[] - t0;
    Print["time interior ", Round[t1], " s, face ", Round[t0], " s, projection ", N[(7 t1 + 3 t0)/3600, 3], " h"];
    Exit[If[ok && stepOk, 0, 1]]]];
Exit[If[runPoints[Join[interior, faces]] && stepOk, 0, 1]];
