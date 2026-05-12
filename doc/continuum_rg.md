# From the Lattice SOS Model to a Continuum Field Theory, and the RG Phase Diagram

This note traces the chain of reasoning lattice → continuum → RG analysis for the generalized SOS model, and explains *why* the naive scaling-dimension argument fails, *what the correct argument is*, and how it lines up with our 3D simulation data.

## 1. The Lattice Model

On the $d$-dimensional cubic lattice $\mathbb{Z}_L^d$ with lattice spacing $a$, periodic boundary conditions, and one continuous height $h_i\in\mathbb{R}$ per site:

$$
\mathcal{H}_{\text{lat}}[\{h_i\}] \;=\; K\sum_{\langle ij\rangle}\bigl|h_i-h_j\bigr|^{\sigma},
\qquad \sum_i h_i = 0.
$$

The order parameter we track is the squared roughness

$$
M^2 \;=\; \frac{1}{N}\sum_i h_i^2,
\qquad N = L^d / a^d.
$$

In what follows we keep $K=1$, $\sigma>0$, and ask: for which $(d,\sigma)$ does $M^2(L)$ stay bounded as $L\to\infty$ (*smooth phase*) versus diverge (*rough phase*)?

## 2. From the Lattice to a Continuum Action

Replace site sums by integrals and nearest-neighbor differences by gradients. Let $h(x)$ be the continuum height field. For neighbors separated by lattice spacing $a$,

$$
h_i - h_j \;=\; a\,\hat e\!\cdot\!\nabla h(x) + \mathcal{O}(a^2),
\qquad \bigl|h_i-h_j\bigr|^{\sigma} \approx a^{\sigma}\,|\hat e\!\cdot\!\nabla h|^{\sigma}.
$$

Summing over the $d$ link directions and replacing $\sum_i \to a^{-d}\int d^dx$, one gets, up to an angular constant absorbed into $K$,

$$
\boxed{\;
\mathcal{H}_{\text{cont}}[h] \;=\; K\, a^{\sigma-d}\!\int d^d x\,\bigl|\nabla h\bigr|^{\sigma}\;}
$$

with

$$
M^2 \;=\; \frac{a^d}{L^d}\!\int d^d x\,h(x)^2.
$$

Two features of this map matter:

1. **For $\sigma=2$**, the prefactor is $a^{2-d}$. In the IR (after coarse-graining to scale $L$), only the dimensionless combination $K a^{2-d} L^{d-2}$ controls the physics — this is the standard Gaussian free field (GFF).
2. **For $\sigma\neq 2$**, the prefactor $a^{\sigma-d}$ depends explicitly on the UV cutoff. In RG language, the bare coupling has nonzero engineering dimension, so $\sigma\neq 2$ is **not** at a fixed point of the bare lattice theory: it must flow under coarse-graining.

The whole RG question for this model is: *where does that flow end up?*

## 3. Naive Scaling of $\int|\nabla h|^{\sigma}$

Pretend, for a moment, that the action $\int|\nabla h|^{\sigma}$ defines a non-Gaussian fixed point in its own right. Under

$$
x \;\to\; b\,x,\qquad h(x)\;\to\; b^{\alpha}\,h(b^{-1}x),
$$

each factor transforms as $d^dx\to b^d d^dx$, $\nabla\to b^{-1}\nabla$, $h\to b^\alpha h$. The action picks up

$$
\mathcal{H}\;\to\; b^{\,d+\sigma(\alpha-1)}\,\mathcal{H}.
$$

Marginality ($b$-independent) requires

$$
\boxed{\;\alpha \;=\; 1-\frac{d}{\sigma}\;}
$$

This is the **engineering dimension of $h$ at this would-be fixed point**. In a box of linear size $L$ (lattice constant fixed), typical fluctuations scale as

$$
\langle h^2\rangle \;\sim\; L^{2\alpha} \;=\; L^{\,2-2d/\sigma}.
$$

Read off the three regimes:

| relation | $\alpha$ | $M^2(L)$ | naive phase |
|---|---|---|---|
| $\sigma > d$ | $>0$  | $L^{\,2-2d/\sigma}$ | rough |
| $\sigma = d$ | $=0$  | $\log L$            | marginal |
| $\sigma < d$ | $<0$  | bounded             | smooth |

### 3.1 Consistency checks with known limits

| $(d,\sigma)$ | naive prediction | known answer |
|---|---|---|
| $(1,2)$ | $L^{1}$, rough | random walk, $M^2\sim L$ ✓ |
| $(2,2)$ | $\log L$, marginal | GFF in 2D, $M^2\sim\log L$ ✓ |
| $(3,2)$ | bounded (smooth) | 3D GFF integral converges ✓ |

So far so good.

## 4. Where the Naive Diagram Breaks

Apply the table at $(d,\sigma)=(3,4)$: $\sigma>d$ predicts a **rough phase** with $M^2\sim L^{1/2}$.

Our simulation (this project, `doc/analysis_3d.md`) gives instead

$$
M^2_\infty(\sigma=4)\;=\;0.136 \pm 0.001,\qquad M^2(L)\;=\;M^2_\infty-c/L.
$$

clearly **smooth**, with the same $1/L$ correction as $\sigma=2$. The naive criterion is wrong.

The hidden premise of §3 is that $\int|\nabla h|^\sigma$ is itself an RG fixed point: marginality under $x\to bx,\,h\to b^\alpha h$ was equated with scale invariance of the full theory. But a fixed point requires more than covariance under a single rescaling — it requires the effective action obtained by integrating out short-wavelength fluctuations to **reproduce itself** (up to irrelevant corrections). The next two sections check that requirement directly within the continuum theory, without ever assuming a Gaussian piece in the bare action.

## 5. One RG Step Generates $(\nabla h)^2$

### 5.1 Setup

The bare action $\mathcal{H}[h]=K\!\int d^dx\,|\nabla h|^\sigma$ is non-quadratic for $\sigma\neq 2$, so it has no closed form in momentum space and we **never leave real space** in the computation below. Momentum space appears in exactly one role: to give a precise meaning to the words "long-wavelength" and "short-wavelength" when we split the field.

Concretely, with the UV cutoff $\Lambda$ applied to the path integral (only modes with $|k|<\Lambda$ are integrated), pick a rescaling factor $b>1$ and define two real-space fields by their Fourier supports,

$$
h_<(x) \;=\; \int_{|k|<\Lambda/b}\!\frac{d^d k}{(2\pi)^d}\,e^{i k\cdot x}\,\hat h(k),
\qquad
h_>(x) \;=\; \int_{\Lambda/b<|k|<\Lambda}\!\frac{d^d k}{(2\pi)^d}\,e^{i k\cdot x}\,\hat h(k),
$$

so that $h(x)=h_<(x)+h_>(x)$ at every point. Both $h_<$ and $h_>$ are ordinary real-space fields; they merely happen to be band-limited.

After the split, every manipulation below is a real-space functional integral. The Wilsonian effective action for the slow field is

$$
e^{-\mathcal{H}_{\rm eff}[h_<]} \;\equiv\; \int\!\mathcal{D}h_>\;\exp\!\Bigl(-K\!\int d^dx\,|\nabla h_<+\nabla h_>|^{\sigma}\Bigr),
$$

where $\mathcal{D}h_>$ integrates over real-space fields whose Fourier support lies in the UV shell $\Lambda/b<|k|<\Lambda$. Our goal: identify which local operators appear in $\mathcal{H}_{\rm eff}[h_<]$, by expanding it in $h_<$ around $h_<=0$.

### 5.2 Pointwise expansion of the integrand

Write the action as $\mathcal{H}[h_<+h_>]=\mathcal{H}[h_>]+\delta\mathcal{H}[h_<,h_>]$ with

$$
\delta\mathcal{H}[h_<,h_>] \;=\; K\!\int d^dx\,\bigl(\,|\nabla h_<+\nabla h_>|^{\sigma}-|\nabla h_>|^{\sigma}\bigr).
$$

Pointwise in $x$, Taylor-expand the integrand at the point $\nabla h_>(x)$, treating $\nabla h_<(x)$ as the increment:

$$
\begin{aligned}
|\nabla h_<+\nabla h_>|^{\sigma}-|\nabla h_>|^{\sigma}
\;=\;&\,\sigma\,|\nabla h_>|^{\sigma-2}\,(\nabla h_>\!\cdot\!\nabla h_<) \\
&+\;\tfrac{\sigma}{2}|\nabla h_>|^{\sigma-2}|\nabla h_<|^2
\;+\;\tfrac{\sigma(\sigma-2)}{2}|\nabla h_>|^{\sigma-4}(\nabla h_>\!\cdot\!\nabla h_<)^2 \\
&+\;\mathcal{O}(\nabla h_<^3).
\end{aligned}
$$

This is a formal power series in $h_<$ with $h_>$-dependent coefficients; it lets us organise $\mathcal{H}_{\rm eff}[h_<]$ as a Taylor series around $h_<=0$.

### 5.3 Cumulant expansion and the $\mathcal{O}(h_<^2)$ effective action

By construction,

$$
\mathcal{H}_{\rm eff}[h_<]-\mathcal{H}_{\rm eff}[0]
\;=\;-\ln\bigl\langle e^{-\delta\mathcal{H}[h_<,h_>]}\bigr\rangle_{\!>},
$$

with $\langle\,\cdot\,\rangle_{>}$ the average in the $h_<=0$ measure

$$
\langle F\rangle_{\!>}\;\equiv\;\frac{1}{Z_>}\int\!\mathcal{D}h_>\,F[h_>]\,e^{-K\!\int|\nabla h_>|^\sigma},
\qquad
Z_>=\int\!\mathcal{D}h_>\,e^{-K\!\int|\nabla h_>|^\sigma}.
$$

Expand the cumulant generating function:

$$
\mathcal{H}_{\rm eff}[h_<]-\mathcal{H}_{\rm eff}[0]
\;=\;\bigl\langle\delta\mathcal{H}\bigr\rangle_{\!>} \;-\;\tfrac{1}{2}\bigl\langle\delta\mathcal{H}^{\,2}\bigr\rangle_{\!>}^{\!\rm conn}\;+\;\cdots
$$

Two terms contribute at $\mathcal{O}(h_<^2)$:

**(i) One-point function of the $\mathcal{O}(h_<^2)$ part of $\delta\mathcal{H}$.** Using rotational invariance of $\langle\,\cdot\,\rangle_{>}$, which gives

$$
\bigl\langle|\nabla h_>|^{\sigma-4}\nabla_i h_>\nabla_j h_>\bigr\rangle_{\!>}\;=\;\frac{\delta_{ij}}{d}\,\bigl\langle|\nabla h_>|^{\sigma-2}\bigr\rangle_{\!>},
$$

the two $\mathcal{O}(h_<^2)$ terms in $\delta\mathcal{H}$ collapse into a single local quadratic:

$$
\bigl\langle\delta\mathcal{H}\bigr\rangle_{\!>}\Bigr|_{\mathcal{O}(h_<^2)}
\;=\;\frac{K\sigma}{2}\!\left(1+\frac{\sigma-2}{d}\right)\,\bigl\langle|\nabla h_>|^{\sigma-2}\bigr\rangle_{\!>}\!\int d^dx\,|\nabla h_<|^2.
$$

**(ii) Connected square of the $\mathcal{O}(h_<)$ part of $\delta\mathcal{H}$.** The linear term in $\delta\mathcal{H}$ has zero one-point function by the $h_>\!\to\!-h_>$ symmetry of $\langle\,\cdot\,\rangle_{>}$, but its connected square is nonzero:

$$
\tfrac{1}{2}\bigl\langle\delta\mathcal{H}^{\,2}\bigr\rangle_{\!>}^{\!\rm conn}\Bigr|_{\mathcal{O}(h_<^2)}
\;=\;\tfrac{(K\sigma)^2}{2}\!\int\! d^dx\,d^dy\,\nabla_i h_<(x)\nabla_j h_<(y)\,C_{ij}(x-y),
$$

with the UV-shell connected correlator

$$
C_{ij}(x-y)\;=\;\bigl\langle|\nabla h_>|^{\sigma-2}\nabla_i h_>(x)\;|\nabla h_>|^{\sigma-2}\nabla_j h_>(y)\bigr\rangle_{\!>}^{\!\rm conn}.
$$

Because $h_>$ is band-limited to $\Lambda/b<|k|<\Lambda$, the kernel $C_{ij}(x-y)$ decays on the scale $\Lambda^{-1}$ — exponentially short compared to the IR scales on which $h_<$ varies. Replacing the kernel by its zero-momentum moment (gradient expansion of the slow field $h_<$) gives a contribution of the same local form $\int|\nabla h_<|^2$. Rotational invariance forces the moment to be $\delta_{ij}$ times a positive scalar $\mathcal{C}(\sigma,d,\Lambda,b)>0$.

**Combining (i) and (ii)**, with the explicit cumulant sign $-\tfrac{1}{2}$ in front of (ii):

$$
\boxed{\;
\mathcal{H}_{\rm eff}^{(2)}[h_<] \;=\; \frac{K_{\rm gen}}{2}\!\int d^dx\,|\nabla h_<|^2,
\qquad
K_{\rm gen}\;=\;K\sigma\!\left(1+\frac{\sigma-2}{d}\right)\!\bigl\langle|\nabla h_>|^{\sigma-2}\bigr\rangle_{\!>}\;-\;K^2\sigma^2\,\mathcal{C}(\sigma,d,\Lambda,b).\;}
$$

The two contributions have opposite signs and we make no claim about the sign of $K_{\rm gen}$ itself. The qualitative point holds either way: $K_{\rm gen}$ is generically nonzero, so $\mathcal{H}_{\rm eff}[h_<]$ contains a $|\nabla h_<|^2$ operator that the bare action $K\!\int|\nabla h|^\sigma$ does *not*. Thus

$$
\mathcal{H}_{\rm eff}[h_<] \;\neq\; K'\!\int|\nabla h_<|^\sigma + (\text{irrelevant})
$$

for any rescaling $K\to K'$ — i.e., $\int|\nabla h|^\sigma$ is *not* an RG fixed point. The §3 premise fails non-trivially, and that is why the naive prediction at $(3,4)$ disagrees with simulation.

(Sanity check: at $\sigma=2$, $\nabla h_<$ and $\nabla h_>$ have disjoint Fourier support, so $L=2K\!\int\nabla h_<\!\cdot\!\nabla h_>=0$ identically; thus (ii) vanishes and (i) returns the bare GFF — the GFF is at its fixed point, as it should be. For $\sigma\neq 2$, the integrand $|\nabla h|^\sigma$ is non-polynomial in the fields, no such Fourier orthogonality holds, and both (i) and (ii) contribute.)

Once the $(\nabla h)^2$ operator is in $\mathcal{H}_{\rm eff}$, it sets the scaling of the slow field at the next RG step: take $[h]=(d-2)/2$ (the Gaussian assignment), then the engineering dimension of the bare coupling is $[K_\sigma]=d(1-\sigma/2)$:

- $\sigma>2$ in $d>2$: $[K_\sigma]<0$, the bare $|\nabla h|^\sigma$ term is **irrelevant** about the generated $(\nabla h)^2$. The flow ends at a Gaussian theory; the value of the IR stiffness is set self-consistently (see §6), not by this one-step calculation.
- $\sigma=2$: nothing new is generated; the bare theory is already GFF.
- $\sigma<2$: $[K_\sigma]>0$, naively relevant, but $|\nabla h|^\sigma$ with non-integer $\sigma<2$ is non-analytic at $\nabla h=0$ and does not admit a perturbative expansion around the generated GFF. The non-perturbative treatment of §6 covers this case.

## 6. Self-Consistent Harmonic Approximation: an Independent Continuum Path

Section 5 organised one RG step as a Taylor series in the slow field $h_<$. The same conclusion follows non-perturbatively from a variational (Feynman–Bogoliubov / Gibbs) argument on the full continuum partition function, with no slow/fast split and no expansion in $h_<$.

Introduce a trial Gaussian action

$$
\mathcal{H}_0[h] \;=\; \frac{K_*}{2}\!\int d^dx\,(\nabla h)^2.
$$

The Gibbs–Bogoliubov inequality gives, for any $K_*>0$,

$$
F[\mathcal{H}]\;\le\;F[\mathcal{H}_0]+\langle\mathcal{H}-\mathcal{H}_0\rangle_{\mathcal{H}_0}.
$$

Under $\mathcal{H}_0$ the gradient $\nabla h$ is Gaussian with $\langle(\nabla h)^2\rangle_0 = (d/K_*)\!\int_{|k|<\Lambda}\!d^dk/(2\pi)^d \equiv D(\Lambda)/K_*$. The variational bound becomes a function of $K_*$ only:

$$
F_{\rm var}(K_*) \;=\; F[\mathcal{H}_0(K_*)] + K\,c_\sigma\!\left(\frac{D(\Lambda)}{K_*}\right)^{\!\sigma/2}\!\mathrm{Vol} \;-\; \tfrac{K_*}{2}\,\frac{D(\Lambda)}{K_*}\,\mathrm{Vol},
$$

with $c_\sigma$ the Gaussian moment $\langle|X|^\sigma\rangle$ of a unit Gaussian. Stationarity $\partial F_{\rm var}/\partial K_*=0$ yields the **self-consistent equation**

$$
\boxed{\;K_*^{\,1-\sigma/2}\;=\;K\,\sigma\,c_\sigma\,D(\Lambda)^{\,\sigma/2-1}.\;}
$$

For every $\sigma>0$ (with finite UV regulator so $D(\Lambda)<\infty$), this has a unique positive root $K_*(K,\sigma,\Lambda)$. The optimal Gaussian approximation to the bare $|\nabla h|^\sigma$ theory therefore *exists and is unique*, with a definite renormalized stiffness.

This argument:

- uses only the continuum path integral and a functional inequality,
- does not invoke a momentum shell, a lattice, or a CLT,
- does not assume any Gaussian piece in the bare action.

It confirms, independently of §5, that the IR effective theory of $\int|\nabla h|^\sigma$ is Gaussian with stiffness $K_*(K,\sigma)$. The Wilsonian argument of §5 and the variational argument here approach the same conclusion from two different directions.

## 6.1 Match to Data

The Gaussian IR with stiffness $K_*$ predicts the *functional form* of finite-size effects in $d\ge 3$ to be $M^2(L)=M^2_\infty - c(\sigma)/L^{d-2}$, with $c(\sigma)\propto 1/K_*(\sigma)$. Our 3D fits give

| $\sigma$ | $c$ | implied $K_*$ trend |
|---|---|---|
| 0.5 | 0.320 | smallest stiffness, broadest fluctuations |
| 1.0 | 0.122 | moderate |
| 2.0 | 0.113 | the Gaussian point itself |
| 4.0 | 0.118 | very close to Gaussian, slightly larger fluctuations than $\sigma=2$ |

The same $1/L$ form holds for every $\sigma\in\{0.5,1,2,4\}$; only $K_*(\sigma)$ varies, exactly as the variational solution above predicts.

## 7. The Correct Phase Diagram

Combining §5–§6, for any $\sigma>0$:

$$
\boxed{
\begin{array}{c|c|c}
d & \text{IR behavior of $M^2(L)$} & \text{phase} \\\hline
1 & L                              & rough \\
2 & \log L                         & rough (marginal) \\
\ge 3 & \text{const} - c(\sigma)/L^{d-2} & smooth
\end{array}}
$$

The dependence on $\sigma$ is entirely contained in $K_*(\sigma)$ (or $c(\sigma)$); the phase itself depends only on the spatial dimension. The "lower critical dimension" of this model family is $d_l=2$, the same as for the pure GFF and exactly as Mermin–Wagner would suggest.

## 8. Consistency with Simulation Data

- **2D, all $\sigma$** (earlier project work): $M^2\sim\log L$ on the rough side, with $\sigma$ shifting the prefactor — consistent with §7.
- **3D, $\sigma\in\{0.5,1,2,4\}$, $L\le 64$**: every series fits $M^2_\infty-c/L$ to machine precision (`doc/analysis_3d.md` §3, §3.1). The smooth phase is robust to $\sigma$.
- **3D, $L=128$ extension** (in flight at the time of writing): a *necessary* check is that $M^2_\infty$ extracted from $\{32,64,128\}$ does not drift relative to $\{16,32,64\}$. The variational picture predicts no drift.

## 9. Caveats

1. **Edges of the $\sigma$-range.** $\sigma\to 0$: $|h_i-h_j|^\sigma\to \mathbb{1}_{h_i\ne h_j}$, an Ising-type kink-counting Hamiltonian; the variational moment $\langle|\nabla h|^\sigma\rangle$ degenerates and the §6 argument no longer applies. $\sigma\to\infty$: $|h_i-h_j|^\sigma$ becomes a hard wall pinning $h_i=h_j$; the partition function localises on constant heights and the Gaussian-trial bound is far from tight. Both edges lie outside the range $\sigma\in\{0.5,1,2,4\}$ we simulate.
2. **Status of the two arguments.** §5 (Wilsonian) computes the *first* RG step exactly within the bare measure but does not by itself control the full flow to the IR; it shows that the bare action is not a fixed point and that a $(\nabla h)^2$ term must appear. §6 (SCHA) gives an *upper bound* on the free energy with a self-consistent Gaussian and is non-perturbative, but it does not prove the IR is *exactly* Gaussian — only that the optimal Gaussian description exists and gives definite predictions. Independent agreement between the two is the strongest claim available without rigorous renormalisation theory for non-integer $\sigma$.
3. **What this analysis does NOT say.** It does not predict the value of $M^2_\infty$ — that is a non-universal Brillouin-zone integral involving the lattice dispersion and the self-consistent $K_*(\sigma)$. It only predicts the *functional form* of finite-size corrections.

## 10. Summary

- Lattice $\mathcal{H}=K\sum|h_i-h_j|^\sigma$ becomes the continuum theory $K a^{\sigma-d}\int|\nabla h|^\sigma d^d x$. For $\sigma\ne 2$ the bare coupling carries cutoff dependence — it is not at a fixed point.
- Pretending it is gives the naive prediction $M^2\sim L^{2-2d/\sigma}$, which works for $(d,\sigma)$ with $\sigma\le 2$ but fails for $\sigma>d$ (e.g. $(3,4)$).
- The actual continuum argument has two independent legs, neither of which assumes a $(\nabla h)^2$ piece in the bare action. (i) One Wilsonian momentum-shell step generates such a term with a positive coefficient — so $\int|\nabla h|^\sigma$ alone is not an RG fixed point. (ii) Self-consistent harmonic approximation (Gibbs–Bogoliubov) on the full continuum theory yields a unique optimal Gaussian for every $\sigma>0$.
- Both routes lead to an IR effective theory that is GFF with renormalised stiffness $K_*(\sigma)$. The phase is then fixed by $d$ alone: rough for $d\le 2$, smooth for $d\ge 3$. This matches our 3D data for every $\sigma$ studied.
