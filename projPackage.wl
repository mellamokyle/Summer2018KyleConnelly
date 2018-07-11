(* ::Package:: *)

BeginPackage["TimeSeriesDescription`"]

semanticInfo::usage = "Get semantic information of graph, in numerical form:
						{{{length, {relative slope}}, ..., {length, {relative slope}}}, global trend ratio}"
plotPiecewiseFit::usage = "Calculate fit of piecewise polynomials and plot it."
dataToWords::usage = "Get words to describe data."
fixTimeSeries::usage = "Put time series data into pure numerical form for processing."

Begin["`Private`"]
fixTimeSeries[ts_] := Apply[{QuantityMagnitude[#1-Normal[ts][[1,1]]], QuantityMagnitude[#2]}&, Normal[ts], 1]

base[center_, support_, x_] := N[BSplineBasis[3, (1./support)*(x-center) + 1/2]];

piecewiseCoeffs[data_, nPts_]:= 
	Block[{range, length, basisMat},
		range = data[[-1,1]] - data[[1,1]];
		length = Length@data;
		basisMat = SparseArray[Table[base[range*(c-1)/(nPts-1)+data[[1,1]], 4 range/nPts, data[[r, 1]]], {r, length}, {c, nPts}]];
		LeastSquares[basisMat, data[[All, 2]]];
	]

plotPiecewiseFit[data_, nPts_]:=
	Block[{range, length, basisMat, consts},
		range = data[[-1,1]] - data[[1,1]];
		length = Length@data;
		basisMat = SparseArray[Table[base[range*(c-1)/(nPts-1)+data[[1,1]], 4 range/nPts, data[[r, 1]]], {r, length}, {c, nPts}]];
		consts = LeastSquares[basisMat, data[[All, 2]]];
		Plot[consts.Table[base[range (c-1)/(nPts-1) + data[[1,1]], 4range/nPts, x], {c, nPts}], {x, data[[1,1]], data[[-1,1]]}, PlotStyle->Red]
	]
	
(*it does not like to have the c.Table stuff returned as a function*)

piecewiseCoeffsWLinear[data_, nPts_]:= 
	Block[{dom, length, basisMat, linFit, fixedData},
		linFit = LinearModelFit[data, xp, xp];
		fixedData = {#1, #2 - linFit[#1]}&@@@data;
		dom = data[[-1,1]] - data[[1,1]];
		length = Length@data;
		basisMat = SparseArray[Table[base[dom*(c-1)/(nPts-1)+data[[1,1]], 4 dom/nPts, data[[r, 1]]], {r, length}, {c, nPts}]];
		{LeastSquares[basisMat, fixedData[[All, 2]]], Abs[linFit["BestFitParameters"][[2]]*dom]}
	]

membership[x_, width_] := Ramp[x/width + 1] - UnitStep[x]Ramp[x/width] - UnitStep[x-width]Ramp[x/width - 1] + UnitStep[x-2width]Ramp[x/width - 2]
classify[x_, n_, min_, max_]:= 
	With[{range = max - min},
		membership[x - (#-1.)/(n) range - min, range/(n)]&/@Range[n]
	]

delta[x_, y_] := 
With [{len = Length@x},
	Map[Max, 
		Map[x[[#]]&, Table[i, {i, Max[1, 1-#], Min[len, len-#]}]&/@(Range[2len-1]-len)]*
		Map[y[[#]]&, Table[i + #, {i, Max[1, 1-#], Min[len, len-#]}]&/@(Range[2len-1]-len)]
	]
]

heightWords = <|{-2}-> "fast decline", {-1} -> "slow decline", {0}->"flat", {1}-> "growth", {2}-> "fast growth"|>;
lengthWordsTemplate[n_] := AssociationThread[Range[n],
	Join[Table["Brief period", {i, Floor[n/10.]}], 
			Table["Period", {i, Floor[n/10.]+1, Floor[n/5.]}], 
			Table["Long period", {i, Floor[n/5.]+1, Floor[n/2.]}], 
			Table["Mostly", {i, Floor[n/2.]+1, n}]]];

semanticInfo[data_, detail_]:=
	Block[{info, coeffs},
		coeffs = piecewiseCoeffsWLinear[data, detail];
		info = Map[classify[#, 7, Min@coeffs[[1]], Max@coeffs[[1]]]&, coeffs[[1]]];
		info = Map[(Range[13]-7).#&, Map[delta[info[[#]], info[[#+1]]]&, Range[Length@info-1]]];
		info = Map[Round, Map[{Length[#], Mean[#]}&, Split[info, (Abs[#1- #2] <=  2.1)&]]];
		info = Apply[{#1, classify[#2, 5, Min[info[[All,2]]], Max[info[[All,2]]]]}&, info, {1}];
		info = Apply[{#1, Position[#2,Max@#2][[1]]-3}&, info, {1}];
		lengthWords = lengthWordsTemplate[detail];
		{info, coeffs[[2]]/Max[Abs[coeffs[[1]]]]}
	]

wordify[semantics_]:= 
	Block[{firstpass},
		firstPass = SequenceReplace[test[[1]], {Repeated[{{_, {2}}, {_, {-2}}}, {2, Infinity}]-> "large oscillation", 
												Repeated[{{_, {-2}}, {_, {2}}}, {2, Infinity}]-> "large oscillation",
												Repeated[{{_, {1}}, {_, {-1}}}, {2, Infinity}]-> "small oscillation",
												Repeated[{{_, {-1}}, {_, {1}}}, {2, Infinity}]-> "small oscillation",   
												{{_, {2}}, RepeatedNull[{1|2, {-1|0|1}}, 1], {_, {-2}}}-> "peak", 
												{{_, {-2}}, RepeatedNull[{1|2, {-1|0|1}}, 1], 
												{_, {2}}} -> "valley"}];
												
		secondPass = {Replace[Replace[firstPass, heightWords, {2}], lengthWords, All], semantics[[2]]}
	]

ratioDesc[x_] := 
	Piecewise[{{"Strong upward trend", x >= 4}, 
			    {"Upward trend", 2 <= x < 4}, 
				{"Flat", 12 <= x < 2}, 
				{"Downward trend", -4 <= x < -2}, 
				{"Strong downward trend", x < -4}}
	]

dataToWords[data_, detail_]:=
	Block[{info, lengthWords, coeffs, secondPass, firstPass},
		coeffs = piecewiseCoeffsWLinear[data, detail];
		info = Map[classify[#, 7, Min@coeffs[[1]], Max@coeffs[[1]]]&, coeffs[[1]]];
		info = Map[(Range[13]-7).#&, Map[delta[info[[#]], info[[#+1]]]&, Range[Length@info-1]]];
		info = Map[Round, Map[{Length[#], Mean[#]}&, Split[info, (Abs[#1- #2] <=  2.1)&]]];
		info = Apply[{#1, classify[#2, 5, Min[info[[All,2]]], Max[info[[All,2]]]]}&, info, {1}];
		info = Apply[{#1, Position[#2,Max@#2][[1]]-3}&, info, {1}];
		lengthWords = lengthWordsTemplate[detail];
		firstPass = SequenceReplace[info, {Repeated[{{_, {2}}, {_, {-2}}}, {2, Infinity}]-> "Large oscillation", 
												Repeated[{{_, {-2}}, {_, {2}}}, {2, Infinity}]-> "Large oscillation",
												Repeated[{{_, {1}}, {_, {-1}}}, {2, Infinity}]-> "Small oscillation",
												Repeated[{{_, {-1}}, {_, {1}}}, {2, Infinity}]-> "Small oscillation",   
												{{_, {2}}, RepeatedNull[{1|2, {-1|0|1}}, 1], {_, {-2}}}-> "Peak", 
												{{_, {-2}}, RepeatedNull[{1|2, {-1|0|1}}, 1], 
												{_, {2}}} -> "Valley"}];
		secondPass = Replace[Replace[firstPass, heightWords, {2}], lengthWords, All];
		Join[Quiet[secondPass /. {a_String,b_String}->a<> " of "<>b], {ratioDesc[coeffs[[2]]/Max[Abs[coeffs[[1]]]]]}]
	]

End[]
EndPackage[]









