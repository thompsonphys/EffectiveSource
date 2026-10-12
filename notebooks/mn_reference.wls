dir = DirectoryName[$InputFileName];
nbfile = FileNameJoin[{dir, "ScalarSingularField.nb"}];
mode = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "run"];
defaultOut[m_] := If[MemberQ[{"selftest", "test"}, m], FileNameJoin[{$TemporaryDirectory, "mn_reference_" <> m <> ".txt"}], FileNameJoin[{dir, "mn_reference.txt"}]];
out = If[Length[$ScriptCommandLine] >= 3, $ScriptCommandLine[[3]], defaultOut[mode]];
loadDefs[] := With[{nb = Get[nbfile]}, Do[Quiet[ToExpression[b, StandardForm]], {b, Cases[nb, Cell[b_BoxData, "Input", ___] :> b, Infinity]}]];
loadDefs[];
LaunchKernels[8];
DistributeDefinitions[nbfile, loadDefs];
ParallelEvaluate[loadDefs[]];
prec = 110;
orb = KerrEquatorialOrbit[1/2, 10, 3/10, prec];
ms = {0, 2, 10}; ns = Range[-20, 20];
radii = {87/13, 9, 11, 13, 107/7}; angles = {2/25, 1/50, 1/200, ArcSin[2/5]};
checkPoints = {{107/7, 2/25}, {11, 1/50}, {13, 1/200}};
points = Join[Tuples[{radii, angles}], {{9, 0}, {12, 0}, {107/7, 0}}];
spots = {{{107/7, 2/25}, 2, 0, "Phi", 25}, {{11, 1/50}, 2, 3, "r", 25}, {{12, 0}, 2, -1, "r", 20}};
inRange[r_] := 100/13 < r < 100/7;
quantities[pt_] := Which[pt[[2]] === 0 && inRange[pt[[1]]], {"Phi", "r"}, MemberQ[checkPoints, pt], {"Phi", "r", "th", "rr", "thth", "box"}, True, {"Phi", "r", "th", "box"}];
If[mode === "test", points = {{107/7, 0}, {12, 0}}; ms = {2}; ns = {-1, 0, 1};
  quantities[pt_] := If[pt[[2]] === 0 && inRange[pt[[1]]], {"Phi", "r"}, {"Phi", "r", "th"}]];
label[pt_] := "r=" <> ToString[pt[[1]], InputForm] <> ",dth=" <> ToString[pt[[2]], InputForm];
theta[d_] := If[d === 0, Pi/2, Pi/2 + N[d, prec]];
fmt[x_, k_] := If[x == 0, "0", Module[{d, e}, {d, e} = RealDigits[x, 10, k]; If[x < 0, "-", ""] <> ToString[First[d]] <> "." <> StringJoin[ToString /@ Rest[d]] <> "E" <> ToString[e - 1]]];
num[s_] := If[s === "0", 0, ToExpression[StringReplace[s, "E" -> "*^"]]];
safeMap[f_, list_] := Module[{res = ParallelMap[f, list, Method -> "FinestGrained"], bad},
  bad = Flatten[Position[res, x_ /; ! AssociationQ[x], {1}, Heads -> False]];
  If[bad =!= {}, res[[bad]] = f /@ list[[bad]]]; res];
writeLines[lines_] := With[{st = OpenAppend[out]}, WriteString[st, StringRiffle[lines, "\n"] <> "\n"]; Close[st]];
passes[x_, t_] := TrueQ[x >= t];
finish[lines_, ok_] := writeLines[If[TrueQ[ok], Append[lines, "end"], lines]];
adaptiveValue[sp_, limit_] := TimeConstrained[MNModeAdaptive[orb, N[sp[[1, 1]], prec], theta[sp[[1, 2]]], sp[[2]], sp[[3]], sp[[4]], sp[[5]]], limit, $Aborted];
adaptiveLines[sp_, val_, ref_] := Module[{dig = If[NumericQ[val], N[-Log10[Abs[SetPrecision[val, Infinity] - SetPrecision[ref[sp[[3]]], Infinity]]/Max[Abs[SetPrecision[Values[ref], Infinity]]]], 4], "-"]},
  {{"adaptive " <> label[sp[[1]]] <> " " <> ToString[sp[[2]]] <> " " <> ToString[sp[[3]]] <> " " <> sp[[4]] <> " " <> If[NumericQ[val], fmt[val, 35], "aborted"] <> " " <> ToString[dig],
    "check adaptive:" <> label[sp[[1]]] <> " " <> If[passes[dig, sp[[5]]], "PASS", "FAIL"] <> " " <> ToString[dig]}, passes[dig, sp[[5]]]}];
header = {"orbit 1/2 10 3/10", "constants " <> fmt[orb["E"], 40] <> " " <> fmt[orb["L"], 40], "freq " <> StringRiffle[fmt[#, 40] & /@ {orb["OmegaR"], orb["OmegaPhi"], orb["TR"], orb["LambdaR"]}, " "]};
prepare[] := Module[{lines, last, keep},
  If[! FileExistsQ[out], writeLines[header]; Return[{}, Module]];
  lines = ReadList[out, String];
  If[MemberQ[lines, "end"], Print["complete: ", out]; Exit[0]];
  last = Max[0, Flatten[Position[lines, s_String /; StringStartsQ[s, "done "]]]];
  keep = If[last == 0, header, Take[lines, last]];
  With[{st = OpenWrite[out]}, WriteString[st, StringRiffle[keep, "\n"] <> "\n"]; Close[st]];
  StringDrop[#, 5] & /@ Select[keep, StringStartsQ[#, "done "] &]];
stored[lab_, m_, q_] := Module[{col = Position[{"Phi", "r", "th", "box", "boxmode"}, q][[1, 1]] + 4},
  Association[(ToExpression[#[[4]]] -> num[#[[col]]]) & /@ Select[StringSplit /@ ReadList[out, String], #[[1]] === "mn" && #[[2]] === lab && #[[3]] === ToString[m] &]]];
If[mode === "selftest",
  Module[{ok = True, chk, blockA, v = N[-Pi 10^-7, 50]},
    chk[name_, cond_] := (Print["check ", name, " ", If[TrueQ[cond], "PASS", "FAIL"]]; ok = ok && TrueQ[cond]);
    chk["safemap", safeMap[If[# == 2 && $KernelID > 0, $Failed, <|"x" -> #|>] &, Range[4]] === Table[<|"x" -> k|>, {k, 4}]];
    chk["format", Abs[num[fmt[v, 35]] - v] < 10^-40 && fmt[0, 35] === "0" && StringMatchQ[fmt[v, 35], "-3." ~~ DigitCharacter .. ~~ "E-7"]];
    chk["label", label[{87/13, ArcSin[2/5]}] === "r=87/13,dth=ArcSin[2/5]" && label[{12, 0}] === "r=12,dth=0"];
    blockA = {"point a 1 0 -", "mn a 2 0 " <> fmt[v, 35] <> " - - - -", "done a"};
    If[FileExistsQ[out], DeleteFile[out]];
    writeLines[Join[header, blockA, {"point b 1 0 -", "mn b 2 0 1.0E0 - - - -"}]];
    chk["resume", prepare[] === {"a"} && ReadList[out, String] === Join[header, blockA]];
    chk["stored", Abs[stored["a", 2, "Phi"][0] - v] < 10^-40];
    chk["default-out", defaultOut["selftest"] =!= FileNameJoin[{dir, "mn_reference.txt"}] && defaultOut["test"] =!= FileNameJoin[{dir, "mn_reference.txt"}] && defaultOut["run"] === FileNameJoin[{dir, "mn_reference.txt"}] && defaultOut["measure"] === FileNameJoin[{dir, "mn_reference.txt"}]];
    chk["passes", passes[Indeterminate, 30] === False && passes[$Failed, 28] === False && passes[Infinity, 30] && passes[31, 30] && ! passes[29, 30]];
    finish[{"check reality FAIL 1"}, False];
    finish[{"check adaptive:a FAIL 1"}, Indeterminate >= 25];
    chk["finish-fail", ReadList[out, String] === Join[header, blockA, {"check reality FAIL 1", "check adaptive:a FAIL 1"}] && prepare[] === {"a"} && ReadList[out, String] === Join[header, blockA]];
    finish[{"check reality PASS 0"}, True];
    chk["finish-pass", Last[ReadList[out, String]] === "end"];
    chk["adaptive-digits", With[{w = adaptiveLines[spots[[1]], N[1/7, 50], <|0 -> num[fmt[N[1/7, 50], 35]], 1 -> N[1/3, 35]|>]}, w[[2]] && ! StringContainsQ[StringJoin[w[[1]]], "Indeterminate"]]];
    chk["spot-goals", spots[[All, 5]] === {25, 25, 20}];
    chk["spot-timeout", adaptiveValue[spots[[1]], 1/10] === $Aborted];
    chk["adaptive-lines", With[{x = adaptiveLines[spots[[3]], $Aborted, <|-1 -> 1|>], y = adaptiveLines[spots[[3]], N[1 + 10^-22, 40], <|-1 -> 1|>], z = adaptiveLines[spots[[1]], N[1 + 10^-22, 40], <|0 -> 1|>]},
      ! x[[2]] && StringContainsQ[x[[1, 1]], " aborted "] && StringContainsQ[x[[1, 2]], " FAIL "] && y[[2]] && StringContainsQ[y[[1, 2]], " PASS "] && ! z[[2]] && StringContainsQ[z[[1, 2]], " FAIL "]]];
    DeleteFile[out];
    Exit[If[ok, 0, 1]]]];
runPoint[pt_] := Module[{q = quantities[pt], lab = label[pt], res, lines, col, fails = {}, routes = {}},
  res = MNModes[orb, N[pt[[1]], prec], theta[pt[[2]]], ms, ns, q, "Map" -> safeMap];
  col[name_, mi_, ni_] := Which[! KeyExistsQ[res["values"], name], "-", Max[Abs[res["values"][name][[mi]]]] == 0, "0", True, fmt[res["values"][name][[mi, ni]], 35]];
  lines = {"point " <> lab <> " " <> ToString[pt[[1]], InputForm] <> " " <> ToString[pt[[2]], InputForm] <> " " <> If[inRange[pt[[1]]], fmt[CrossingTime[orb, N[pt[[1]], prec]], 40], "-"]};
  lines = Join[lines, Flatten[Table[StringRiffle[{"mn", lab, ToString[ms[[mi]]], ToString[ns[[ni]]], col["Phi", mi, ni], col["r", mi, ni], col["th", mi, ni], col["box", mi, ni], col["boxmode", mi, ni]}, " "], {mi, Length[ms]}, {ni, Length[ns]}]]];
  lines = Join[lines, Table[StringRiffle[Join[{"digits", lab, ToString[ms[[mi]]]}, Table[If[KeyExistsQ[res["digits"], name] && res["digits"][name][[mi]] =!= Infinity, ToString[N[Round[res["digits"][name][[mi]], 1/10]]], "-"], {name, {"Phi", "r", "th", "box", "boxmode"}}]], " "], {mi, Length[ms]}]];
  Do[If[! passes[res["digits"][name][[mi]], 30], AppendTo[fails, "check refinement:" <> lab <> ":m=" <> ToString[ms[[mi]]] <> ":" <> name <> " FAIL " <> ToString[N[res["digits"][name][[mi]], 4]]]], {name, Keys[res["digits"]]}, {mi, Length[ms]}];
  If[KeyExistsQ[res["values"], "boxmode"],
    routes = Table[N[-Log10[Max[Abs[res["values"]["box"][[mi]] - res["values"]["boxmode"][[mi]]]]/Max[Abs[res["values"]["box"][[mi]]]]], 4], {mi, Length[ms]}];
    If[! passes[Min[routes], 28], AppendTo[fails, "check routes:" <> lab <> " FAIL " <> ToString[Min[routes]]]]];
  If[MemberQ[checkPoints, pt] && ! KeyExistsQ[res["values"], "boxmode"], AppendTo[fails, "check routes:" <> lab <> " FAIL missing"]];
  If[fails =!= {}, writeLines[fails]; Print /@ fails; Exit[1]];
  writeLines[Join[lines, If[routes === {}, {}, {"check routes:" <> lab <> " PASS " <> ToString[Min[routes]]}], {"done " <> lab}]];
  Print[lab, " done, digits ", N[res["digits"], 3]]];
finals[] := Module[{full, ratios, lines = {}, adapt, ok, ref, dig},
  full = MNModes[orb, N[11, prec], theta[1/50], ms, ns, {"Phi"}, "FullPeriod" -> True, "Rules" -> {"primary"}, "Map" -> safeMap];
  ratios = Table[Max[Abs[Im[full["values"]["Phi"][[mi]]]]]/Max[Abs[full["values"]["Phi"][[mi]]]], {mi, Length[ms]}];
  lines = Table["imag r=11,dth=1/50 " <> ToString[ms[[mi]]] <> " " <> fmt[ratios[[mi]], 5], {mi, Length[ms]}];
  ok = TrueQ[Max[ratios] <= 10^-30];
  AppendTo[lines, "check reality " <> If[ok, "PASS", "FAIL"] <> " " <> fmt[Max[ratios], 5]];
  adapt = ParallelMap[adaptiveValue[#, 7200] &, spots];
  Do[With[{al = adaptiveLines[spots[[k]], adapt[[k]], stored[label[spots[[k, 1]]], spots[[k, 2]], spots[[k, 4]]]]}, lines = Join[lines, al[[1]]]; ok = ok && al[[2]]], {k, Length[spots]}];
  finish[lines, ok];
  ok];
If[mode === "measure",
  Module[{classes = {{11, 1/200}, {107/7, 2/25}, {9, 0}}, per = <||>, total = 0, nodes, sample, t},
    Do[Do[nodes = MNNodes[orb, N[c[[1]], prec], theta[c[[2]]], rule][[1]];
        sample = nodes[[Round[Range[1, Length[nodes], (Length[nodes] - 1)/31]]]];
        t = First[AbsoluteTiming[safeMap[MModeSamples[orb, #, N[c[[1]], prec], theta[c[[2]]], ms, quantities[c], RuleParameters[rule][[3]]] &, sample]]];
        per[{c, rule}] = 8 t/Length[sample];
        Print["class ", c, " rule ", rule, ": ", N[per[{c, rule}], 3], " kernel-seconds per node"], {rule, {"primary", "check"}}], {c, classes}];
    Do[With[{cls = Which[pt[[2]] === 0 && inRange[pt[[1]]], {9, 0}, inRange[pt[[1]]], {11, 1/200}, True, {107/7, 2/25}]},
        total += Sum[Length[MNNodes[orb, N[pt[[1]], prec], theta[pt[[2]]], rule][[1]]] per[{cls, rule}], {rule, {"primary", "check"}}]], {pt, points}];
    Print["projected hours on 8 kernels for the reference points: ", N[total/8/3600, 3]];
    Print["projected hours including 1 h for the final checks: ", N[total/8/3600 + 1, 3]]];
  Exit[0]];
done = prepare[];
Do[If[! MemberQ[done, label[pt]], runPoint[pt]], {pt, points}];
If[mode === "test", writeLines[{"end"}], If[TrueQ[finals[]], Print["all checks PASS"], Print["final checks FAIL"]; Exit[1]]];
Print["wrote ", out];
