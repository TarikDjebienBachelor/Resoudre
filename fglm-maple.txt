with (Groebner):
with (LinearAlgebra):

regle_reec := proc (P, ordre)
    local paire, lc, lt;
    paire := leadmon (P, ordre);
    lc, lt := paire;
    lt = (lc * lt - expand (P)) / lc
end:

p1 := 2*x*y^2 - 5*y^2 - y + 2;
p2 := x^2*y + x - y - 4;

ordre1 := plex (y, x);
base1 := gbasis ([p1, p2], ordre1);
map (regle_reec, base1, ordre1);

ordre2 := tdeg (y, x);
base2 := [];


T := [1, x, x^2, x^3];
AEnumerer := [1, x, y, x^2, x*y];

for i from 1 to nops (AEnumerer) do
    normalf (AEnumerer [i], base1, ordre1)
od;

m := nops (AEnumerer);
n := nops (T);
M := Matrix (m, n):
for lig from 1 to m do
    nf := normalf (AEnumerer [lig], base1, ordre1):
    koeffs := [ coeffs (nf, indets (nf), 'termes') ]:
    termes := [ termes ]:
    for j from 1 to nops (termes) do
        member (termes [j], T, 'col'):
	M [lig, col] := koeffs [j]:
    od:
od:
M;

GaussianElimination (M);

M1 := < M | IdentityMatrix (m) >;
M2 := GaussianElimination (M1);
lambdas := Row (SubMatrix (M2, 1..m, n+1..n+m), -1);

q1 := lambdas . Vector (m, AEnumerer [1..m]);
base2 := [ op (base2), q1 ];
map (regle_reec, base2, ordre2);

AEnumerer := [1, x, y, x^2];
AEnumerer := [1, x, y, x^2, y^2];

m := nops (AEnumerer);
n := nops (T);
M := Matrix (m, n):
for lig from 1 to m do
    nf := normalf (AEnumerer [lig], base1, ordre1):
    koeffs := [ coeffs (nf, indets (nf), 'termes') ]:
    termes := [ termes ]:
    for j from 1 to nops (termes) do
        member (termes [j], T, 'col'):
        M [lig, col] := koeffs [j]:
    od:
od:
M;

GaussianElimination (M);

M1 := < M | IdentityMatrix (m) >;
M2 := GaussianElimination (M1);
lambdas := Row (SubMatrix (M2, 1..m, n+1..n+m), -1);

q2 := lambdas . Vector (m, AEnumerer [1..m]);
base2 := [ op (base2), q2 ];
map (regle_reec, base2, ordre2);


AEnumerer := [1, x, y, x^2];
AEnumerer := [1, x, y, x^2, x^3];

m := nops (AEnumerer);
n := nops (T);
M := Matrix (m, n):
for lig from 1 to m do
    nf := normalf (AEnumerer [lig], base1, ordre1):
    koeffs := [ coeffs (nf, indets (nf), 'termes') ]:
    termes := [ termes ]:
    for j from 1 to nops (termes) do
        member (termes [j], T, 'col'):
        M [lig, col] := koeffs [j]:
    od:
od:
M;

GaussianElimination (M);

M1 := < M | IdentityMatrix (m) >;
M2 := GaussianElimination (M1);
lambdas := Row (SubMatrix (M2, 1..m, n+1..n+m), -1);

q3 := lambdas . Vector (m, AEnumerer [1..m]);
base2 := [ op (base2), q3 ];
map (regle_reec, base2, ordre2);

basis := gbasis ([p1,p2], ordre2):
map (regle_reec, basis, ordre2);

