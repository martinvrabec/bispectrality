z = {u, v};
X = {x, y};

(* Vectors *)
b1 = {Sqrt[2.], 0.};
b2 = {Sqrt[2.]/2., Sqrt[6.]/2.};
b3 = {-Sqrt[2.]/2., Sqrt[6.]/2.};
a3 = {3. Sqrt[2.]/2., Sqrt[6.]/2.};
a2 = {-3. Sqrt[2.]/2., Sqrt[6.]/2.};
a1 = {0., Sqrt[6.]};

R   = {0.5 b1, -0.5 b1, 0.5 b2, -0.5 b2, 0.5 b3, -0.5 b3};
Rplus = {0.5 b1, 0.5 b2, 0.5 b3}; 
m = Association[0.5 b1 -> 2, -0.5 b1 -> 2, 0.5 b2 -> 2, -0.5 b2 -> 2, 0.5 b3 -> 2, -0.5 b3 -> 2];

(* Coroots *)
check[\[Gamma]_] := 2.0 \[Gamma]/(\[Gamma].\[Gamma]); 
(* Reflection *)
s[\[Gamma]_, x_] := x - (check[\[Gamma]].x) \[Gamma];

(* A shift operator *)
T[F_, \[Gamma]_] := 
  F /. {u -> u + \[Gamma][[1]], v -> v + \[Gamma][[2]] };


(* The Macdonald operator for \pi = 2/3 a[1] *)
Wpi = { 2./3. a1, -2./3. a2,  -2./3. a3}; 
Dpi[F_] := Sum[Product[If[\[Gamma].\[Tau] == 1, (1.0 - m[\[Gamma]]/(\[Gamma].z)), 1], {\[Gamma], R}] T[F, \[Tau]], {\[Tau] , Wpi}];


\[Mu] = Sum[Exp[\[Tau].X], {\[Tau] , Wpi}];
M = Sum[m[\[Gamma]], {\[Gamma], Rplus}];
Q = Product[
	    Product[((\[Gamma].z)^2 - s^2), {s, 1, m[\[Gamma]]}], {\[Gamma], 
								   Rplus}];
c = M! Product[
	       Sum[(\[Gamma].\[Tau]) Exp[\[Tau].X], {\[Tau] , Wpi}]^ 
	       m[\[Gamma]], {\[Gamma], Rplus}];

DiffOp[F_] := Dpi[F] - \[Mu] F;

PsiA2 = c^-1  Nest[DiffOp, Q Exp[X.z], M]; 

Print["The size of PsiA2 is: " <> ToString[ByteCount[PsiA2]]]

Print[N[PsiA2 /. {u -> 0.5, v -> 0.1, x -> 0.2, y -> -1}]]
Print[N[PsiA2 /. {u -> 0.5, v -> 0.1, x -> 0.22, y -> -0.11}]]
zero = Nest[DiffOp, Q Exp[X.z], M + 1] // AbsoluteTiming; 
Print[N[zero /. {u -> 0.5, v -> 0.1, x -> 0.2, y -> -1}]]
Print[N[zero /. {u -> 0.5, v -> 0.1, x -> 0.22, y -> -0.11}]]
zero = Nest[DiffOp, Q Exp[X.z], M + 2] // AbsoluteTiming;
Print[N[zero /. {u -> 0.5, v -> 0.1, x -> 0.2, y -> -1}]]
Print[N[zero /. {u -> 0.5, v -> 0.1, x -> 0.22, y -> -0.11}]]

(* H[F_] := -Laplacian[F, {x, y}]  + Sum[m[\[Gamma]] (m[\[Gamma]] + 1) check[\[Gamma]].check[\[Gamma]]/4/Sinh[1/2 check[\[Gamma]].X]^2, {\[Gamma], Rplus }] F;
Print[N[H[PsiA2] + (u^2 + v^2) PsiA2 /. {u -> 0.5, v -> 0.1, x -> 0.22, y -> -0.11}] // AbsoluteTiming]
Print[N[(u^2 + v^2) PsiA2 /. {u -> 0.5, v -> 0.1, x -> 0.22, y -> -0.11}]] *)


R   = {0.5 b1, -0.5 b1, 0.5 b2, -0.5 b2, 0.5 b3, -0.5 b3, 1./6. a1, -1./6. a1, 1./6. a2, -1./6. a2, 1./6. a3, -1./6. a3};
m = Association[0.5 b1 -> 3, -0.5 b1 -> 3, 0.5 b2 -> 3, -0.5 b2 -> 3, 
		      0.5 b3 -> 3, -0.5 b3 -> 3,  1./6. a1 -> 1, -1./6. a1 -> 1, 
		      1./6. a2 -> 1, -1./6. a2 -> 1, 1./6. a3 -> 1, -1./6. a3 -> 1];
Rplus = {0.5 b1, 0.5 b2, 0.5 b3, 1./6. a1, 1./6. a2, 1./6. a3};

Wpi = { 2 b1, 2 b2, 2 b3, -2 b1, -2 b2,  -2 b3}; 
					    
Dpi[F_] := 12 + Sum[Product[If[\[Gamma].\[Tau] == 2, (1.0 - m[\[Gamma]]/(\[Gamma].z)) (1.0 - m[\[Gamma]]/(\[Gamma].z + 1)), 1], {\[Gamma], R}] Product[If[\[Gamma].\[Tau] == 1, (1.0 - m[\[Gamma]]/(\[Gamma].z)), 1], {\[Gamma], R}] T[F, \[Tau]], {\[Tau] , Wpi}]; 

\[Mu] = Sum[Exp[\[Tau].X], {\[Tau] , Wpi}];
M = Sum[m[\[Gamma]], {\[Gamma], Rplus}];
Print["M is equal to: " <> ToString[M]]
Q = Product[
	    Product[((\[Gamma].z)^2 - s^2), {s, 1, m[\[Gamma]]}], {\[Gamma], 
								   Rplus}];
c = M! Product[
	       Sum[(\[Gamma].\[Tau]) Exp[\[Tau].X], {\[Tau] , Wpi}]^ 
	       m[\[Gamma]], {\[Gamma], Rplus}];

DiffOp[F_] := Dpi[F] - \[Mu] F;
Qalt = Product[((\[Gamma].z)^2 - m[\[Gamma]]^2), {\[Gamma], Rplus}];
PsiG2alt = c^-1 Nest[DiffOp, Qalt PsiA2, 5];

Print["The size of PsiG2alt is: " <> ToString[ByteCount[PsiG2alt]]]
(* Print[N[PsiG2alt /. {u -> 0.5, v -> 0.01, x -> 0.2, y -> -0.02}]] *)

Print["MaxMemoryUsed was: " <> ToString[MaxMemoryUsed[]]]

Exit[] 
