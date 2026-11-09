# Supplementary exact certificates for the finite detector hierarchy

This supplement belongs to *A Positive Reciprocal-Product Criterion for the
Riemann Hypothesis*. All entries below are exact rational numbers. They are
reconstructed independently by verify_detector_certificates.py; no finite
grid or floating-point sign decision is used.

For a detector $V_d$ of degree $d$, write $y=N+r$, where
$N=\lfloor y\rfloor\ge1$ and $0\le r<1$. The certificates use

$$
P_{V_d}(N+r)=\frac{N(N+1)}{(N+r)^d}Q_N(r)
$$

and expand $Q_N$ in the Bernstein basis
$\beta_{k,d-2}(r)=\binom{d-2}{k}r^k(1-r)^{d-2-k}$.
Nonnegative Bernstein coefficients imply $P_{V_d}(y)\ge0$ on the complete
half-line.

## Table S1. Degree-twelve detector

Let $m=N-1$. For

$$
V_{12}(u)=u^2(1-u)^7(91u^3-3u^2-7u-1),
$$

write

$$
Q_N(r)=\sum_{k=0}^{10}b_k(m)\beta_{k,10}(r),\qquad
b_k(m)=\frac1{L_k}\sum_{j=0}^{3}c_{k,j}m^{\underline j}.
$$

The exact coefficients are

$$
\begin{array}{c|r|rrrr}
k&L_k&c_{k,0}&c_{k,1}&c_{k,2}&c_{k,3}\\ \hline
0&30&0&245&1170&234\\
1&30&0&395&1245&234\\
2&135&0&2810&5840&1033\\
3&45&0&1445&2015&331\\
4&45&0&2160&2120&316\\
5&135&0&9280&6960&928\\
6&45&0&4160&2620&316\\
7&45&15&5185&2950&331\\
8&135&330&18070&9655&1033\\
9&30&285&4475&2265&234\\
10&30&735&4925&2340&234
\end{array}
$$

Every entry is nonnegative.

## Table S2. Degree-fifteen detector

Let $m=N-1$. For

$$
V_{15}(u)=u^2(1-u)^9
(1024u^4+207u^3-57u^2-27u-3),
$$

write

$$
Q_N(r)=\sum_{k=0}^{13}b_k(m)\beta_{k,13}(r),\qquad
b_k(m)=\frac1{L_k}\sum_{j=0}^{4}c_{k,j}m^{\underline j}.
$$

The exact coefficients are

$$
\begin{array}{c|r|rrrrr}
k&L_k&c_{k,0}&c_{k,1}&c_{k,2}&c_{k,3}&c_{k,4}\\ \hline
0&6&0&946&37133&21996&2444\\
1&6&0&1567&43973&23364&2444\\
2&15&0&6425&128005&62127&6143\\
3&15&0&10420&146420&66015&6209\\
4&15&0&16684&164492&69840&6296\\
5&3&0&5264&36336&14660&1276\\
6&5&0&13600&65960&25376&2144\\
7&5&0&20640&71240&26080&2144\\
8&3&0&18304&46116&15964&1276\\
9&15&12&130932&250172&81264&6296\\
10&15&120&180360&273815&83001&6209\\
11&15&645&238635&301840&85305&6143\\
12&6&975&121497&133433&35292&2444\\
13&6&2838&149478&147113&36660&2444
\end{array}
$$

Every entry is nonnegative.

## Table S3. Degree-eighteen detector

For

$$
V_{18}(u)=\frac{u^2(1-u)^{11}}{8151}
(8151u^5+3645u^4+198u^3-262u^2-77u-7),
$$

the $N=1$ Bernstein coefficients satisfy

$$
b_0=\cdots=b_{10}=0,
$$

$$
(b_{11},\ldots,b_{16})=
\left(
\frac4{24453},\frac{28}{13585},\frac{189}{13585},
\frac{8122}{122265},\frac{1357}{5434},\frac{4227}{5434}
\right).
$$

For $N\ge2$, put $t=N-2$. For $k=1,\ldots,16$, $b_k(t+2)$
is $D_k^{-1}$ times a polynomial whose coefficients, ordered as
$(t^5,t^4,t^3,t^2,t,1)$, are

$$
\begin{array}{c|r|rrrrrr}
k&D_k&c_5&c_4&c_3&c_2&c_1&c_0\\ \hline
1&570570&1746848&16398812&45847764&33396073&2448668&248395\\
2&12837825&39240380&383636520&1146319999&1012925500&220326411&9324770\\
3&4279275&13037660&132666490&423688589&441930505&143021591&5150845\\
4&4279275&12978580&137486620&469237708&564266320&228009592&8472940\\
5&2567565&7747276&85497408&311556284&422928132&197417400&8295420\\
6&855855&2570152&29568736&114798368&172629344&89301720&4472160\\
7&1426425&4269120&51211920&211087376&345846000&193638384&11936960\\
8&12837825&38375040&479688000&2089268992&3677317440&2199529728&170168320\\
9&1426425&4269120&55516080&254128976&475076640&301771584&29576960\\
10&855855&2570152&34685064&165961648&326725776&220242880&27363840\\
11&2567565&7747276&108184492&538427124&1109899928&796390880&124481280\\
12&4279275&12978580&186977880&964150308&2074531300&1592926992&308555520\\
13&4279275&13037660&193275010&1029773789&2310499830&1909803216&448899360\\
14&12837825&39240380&597372980&3283684599&7686871485&6872998836&1912532160\\
15&570570&1746848&27272388&154583524&378058587&366564238&117706120\\
16&570570&1746848&27949568&163386864&418062432&439374073&158939620
\end{array}
$$

For $k=0$,

$$
b_0(N)=
\frac{(N-1)(1746848N^4-14791632N^2+31364905)}{570570}>0
\qquad(N\ge2).
$$

Writing the quartic factor as
$F(N)=1746848N^4-14791632N^2+31364905$, its positivity is closed by
the exact identity

$$
54589F(N)=2(218356N^2-924477)^2+2863351987>0.
$$

## Minimality matrices

At level $\ell=1,\ldots,5$, the endpoint argument forces

$$
V(u)=u^2(1-u)^{2\ell+1}\sum_{j=0}^{\ell-1}a_ju^j
$$

whenever $\deg V\le3\ell+2$. The $\ell-1$ conditions
$[u^{2j+1}]V=0$, $1\le j\le\ell-1$, together with
$\int_0^1V=0$, form an $\ell\times\ell$ rational system.
The verifier recomputes these matrices and confirms rank $\ell$ for each
$\ell=1,\ldots,5$. Thus no nonzero detector in the stated endpoint-flat
polynomial class exists below degrees $6,9,12,15,18$, respectively.

## Reproduction

Run in an environment containing SymPy

    python3 verify_detector_certificates.py

or use `uv run --with sympy python verify_detector_certificates.py` when
`uv` is available.

The expected output is five detector PASS lines followed by
minimality ranks m=1,...,5: PASS.
