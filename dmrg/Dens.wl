(* ::Package:: *)

ClearAll["Global`*"]
SetDirectory[NotebookDirectory[]];
Needs["ErrorBarPlots`"]


names=FileNames["*.out","N_*",Infinity];


names


ff1=Compile[{L,{mat,_Complex,4}},
4\[Pi]^2/L^4*
Sum[
N[
p*k*
Exp[-2*\[Pi]*I*p*(m-l)/L]*
Exp[-2*\[Pi]*I*k*(n-j)/L]*
mat[[m,l,n,j]]
]
,{p,-L/2+1,L/2,1}
,{k,-L/2+1,L/2,1}
,{m,1,L}
,{l,1,L}
,{n,1,L}
,{j,1,L}],
RuntimeAttributes->{Listable},Parallelization->True,CompilationTarget->"C"
];

ff2=Compile[{L,{mat,_Complex,2}},
-4\[Pi]^2/L^4*
(
Sum[
N[k*
Exp[-2*\[Pi]*I*k*(m-l)/L]*
mat[[m,l]]
]
,{k,-L/2+1,L/2,1}
,{m,1,L}
,{l,1,L}]
)^2,
RuntimeAttributes->{Listable},Parallelization->True,CompilationTarget->"C"
];


NN=DeleteDuplicates[StringSplit[names,{"_","/"}][[All,2]]]


mumb={};
Do[
L=10;

data0=ReadList["N_"<>nn<>"/uga.prp",        {Number,Number,Number,Number}][[All,{3,4}]];
data1=ReadList["N_"<>nn<>"/cdc.prp",        {Number,Number,Number,Number}][[All,{3,4}]];
cur  =ReadList["N_"<>nn<>"/corrente.prp",   {Number,Number,Number,Number}][[All,{1,3}]];

dt=cur[[2,1]]-cur[[1,1]];
Nstep=Dimensions[cur][[1]];

mat0 = ArrayReshape[data0[[All,1]],{Nstep,L,L,L,L}]+I*ArrayReshape[data0[[All,2]],{Nstep,L,L,L,L}];
mat1 = ArrayReshape[data1[[All,1]],{Nstep,L,L}]+I*ArrayReshape[data1[[All,2]],{Nstep,L,L}];

fish=Re[Table[
{(t-1)dt,
ff1[L,mat0[[t]]]+
ff2[L,mat1[[t]]]}
,{t,1,Nstep,20}]];

Export["N_"<>nn<>"/fisher.dat",fish];

,{nn,NN}];
