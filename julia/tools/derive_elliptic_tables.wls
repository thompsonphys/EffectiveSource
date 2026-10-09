nb = Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "notebooks", "ScalarSingularField.nb"}]];
Do[Quiet[ToExpression[b, StandardForm]], {b, Cases[nb, Cell[b_BoxData, "Input", ___] :> b, Infinity]}];
out = OpenWrite[$ScriptCommandLine[[2]]];
Do[WriteString[out, row[[1]], " ", StringRiffle[ToString[#, InputForm] & /@ Rest[row], " "], "\n"], {row, EllipticTableEntries[20]}];
Close[out];
