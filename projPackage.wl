(* ::Package:: *)

BeginPackage["TimeSeriesDescription`"]

piecewiseCoeffs::usage = "Get coefficients of piecewise polynomial fit of data: nPts = number of polynomials in fit; data = time series data"
plotPiecewiseFit::usage = "Calculate fit and plot it."

Begin["`Private`"]
base[center_, support_, x_] := N[BSplineBasis[3, (1./support)*(x-center) + 1/2]];

piecewiseCoeffs[nPts_, data_]:= 
	Block[{range, length, basisMat},
		range = data[[-1,1]] - data[[1,1]];
		length = Length@data;
		basisMat = SparseArray[Table[base[range*(c-1)/(nPts-1)+data[[1,1]], 4 range/nPts, data[[r, 1]]], {r, length}, {c, nPts}]];
		LeastSquares[basisMat, data[[All, 2]]];
	]

plotPiecewiseFit[nPts_, data_]:=
	Block[{range, length, basisMat, consts},
		range = data[[-1,1]] - data[[1,1]];
		length = Length@data;
		basisMat = SparseArray[Table[base[range*(c-1)/(nPts-1)+data[[1,1]], 4 range/nPts, data[[r, 1]]], {r, length}, {c, nPts}]];
		consts = LeastSquares[basisMat, data[[All, 2]]];
		Plot[consts.Table[base[range (c-1)/(nPts-1) + data[[1,1]], 4range/nPts, x], {c, nPts}], {x, data[[1,1]], data[[-1,1]]}, PlotStyle->Red]
	]
	
(*it does not like to have the c.Table stuff returned as a function*)

piecewiseCoeffsWLinear[nPts_, data_]:= 
	Block[{dom, length, basisMat, linFit, fixedData},
		linFit = LinearModelFit[data, xp, xp];
		fixedData = {#1, #2 - linFit[#1]}&@@@data;
		dom = data[[-1,1]] - data[[1,1]];
		length = Length@data;
		basisMat = SparseArray[Table[base[dom*(c-1)/(nPts-1)+data[[1,1]], 4 dom/nPts, data[[r, 1]]], {r, length}, {c, nPts}]];
		{LeastSquares[basisMat, fixedData[[All, 2]]], Abs[linFit["BestFitParameters"][[2]]*dom]}
	]

membership[x_, width_] := Piecewise[{{x/width + 1, -width <= x <= 0}, {1, 0 <= x <= width}, {-x/width + 2, width <= x <= 2 width}}]
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

semanticInfo[data_, detail_]:=
	Block[{info, coeffs},
		coeffs = piecewiseCoeffsWLinear[detail, data];
		info = Map[classify[#, 7, Min@coeffs[[1]], Max@coeffs[[1]]]&, coeffs[[1]]];
		info = Map[(Range[13]-7).#&, Map[delta[info[[#]], info[[#+1]]]&, Range[Length@info-1]]];
		info = Map[Round, Map[{Length[#], Mean[#]}&, Split[info, (Abs[#1- #2] <=  2.1)&]]];
		info = Apply[{#1, classify[#2, 5, Min[info[[All,2]]], Max[info[[All,2]]]]}&, info, {1}];
		info = Apply[{#1, Position[#2,Max@#2][[1]]-3}&, info, {1}];
		{info, coeffs[[2]]}
	]

End[]
EndPackage[]



