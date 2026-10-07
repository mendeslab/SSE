# Origin conditioning in QuaSSE

[QuaSSE_origin.xml](QuaSSE_origin.xml) starts one lineage at an estimated origin age and includes
an ancestral fossil before the first observed bifurcation. It demonstrates age proposals and
sampling parameters in a short MCMC run, not a converged inference or a topology-search example.

The initial tree is `(S:0,(A:1,B:1):1);`: S is a fossil ancestor at age 2, the ordinary split is at
age 1, and A and B are living samples. The origin starts at age 3. Its Uniform(0,4) prior is proper;
proposals younger than the sampled root have zero QuaSSE support. The tree ScaleOperator changes
ordinary internal heights while preserving the fossil date and its zero-length observation branch.
The root here is the dated fossil event, so its age stays fixed. This operator does not explore
sampled-ancestor topology changes. Such analyses need compatible sampled-ancestor operators.

Run from an output directory, substituting the absolute checkout path:

```sh
export BEAST2_HOME=/path/to/beast
cd /path/to/SSE
ant native
export LD_LIBRARY_PATH="$PWD/build/native${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
mkdir -p build/origin-run
cd build/origin-run
/path/to/SSE/beast-sse -seed 42 /path/to/SSE/examples/QuaSSE_origin.xml
```

For Java-only integration, copy the XML and set `fftBackend="sst"` and
`integrationBackend="java"`. See the [documentation index](../docs/README.md) for building the
[guide](../docs/quasse/guide.tex) and [mathematical reference](../docs/quasse/theory.tex).

## Starting and sampling options

| Input | Meaning and default |
| --- | --- |
| `conditionOnRoot` | Defaults to `true`: start at an ordinary bifurcation. Set `false` for an origin. |
| `origin` | Scalar age parameter, required only in origin mode; must be at least the sampled root age. |
| `conditionOnSampling` | Require fossil or living samples: on each daughter in root mode, anywhere in origin mode. |
| `conditionOnRhoSampling` | Require living samples: on each daughter in root mode, anywhere in origin mode. |

If both sampling flags are omitted, `conditionOnSampling=true`. If only one flag is supplied,
the omitted flag is false. Both true is an error; both false removes sampling-success conditioning.
Origin and living-sample conditioning require `dynDt=true`. This preserves the historical
root/any-sample default; BEAST SA uses different defaults, and its v2.1.1 implementation does not
expose the root/no-sampling combination. Names alone do not establish matching likelihood settings.

Root mode rejects a sampled-ancestor root and an explicit `origin`. Root/living conditioning also
rejects an all-fossil daughter. Origin mode allows sampled-ancestor roots and multiple stem fossils,
but still requires at least one living sample at age zero. Fossil-only trees, an offset to the
present, and single-leaf trees are not supported. A zero-length origin stem is allowed; its
likelihood need not equal the root-conditioned likelihood because the starting conditions differ.

## Trait weighting and numerical scope

The example uses `Flat` over a finite useful grid interval. Root mode retains SSE's convention
that the effective trait density at a starting split is proportional to the weights times λ.
Origin weights are normalized on the useful bins at the origin age. `Observed` remains available
for comparison and uses the endpoint likelihood to construct its weights; it is not a generative
Bayesian prior. Automatically chosen grid support can itself depend on observations.

Changing the origin does not automatically widen the grid. Check likelihood stability when
refining spatial spacing and timestep and when widening the useful interval separately. With Flat weighting, widening the useful interval changes the trait prior itself, so the resulting
likelihood difference is not solely numerical error. Padding changes can also change useful support;
keeping `nX` fixed alone does not hold the prior interval fixed. A long stem can require a wider domain. These checks establish numerical adequacy for the
analysis being run, not a universal accuracy guarantee.

## Age priors and the posterior

The XML puts `origin` in `State`, assigns a proper prior, and supplies a random-walk proposal.
Conditioning on an age in the likelihood does not mean fixing it in the analysis.

For an ordinary root in root-conditioned mode, use the tree's root height rather than introducing
an unrelated scalar `rootAge`. An all-taxa MRCA prior can provide an age density; for example, with
samples A, B, F, S in the fixed-tree fossil example:

```xml
<distribution spec="beast.base.evolution.tree.MRCAPrior" tree="@tree" monophyletic="true">
  <taxonset spec="beast.base.evolution.alignment.TaxonSet">
    <taxon spec="beast.base.evolution.alignment.Taxon" id="A"/>
    <taxon spec="beast.base.evolution.alignment.Taxon" id="B"/>
    <taxon spec="beast.base.evolution.alignment.Taxon" id="F"/>
    <taxon spec="beast.base.evolution.alignment.Taxon" id="S"/>
  </taxonset>
  <distr spec="Exponential" mean="3.0"/>
</distribution>
```

If those taxon objects are already defined elsewhere, use `idref` instead of defining their IDs
again. Supply an explicit tree taxon set when estimating trees: this avoids a BEAST 2.7.8 parser
ordering issue that can associate incorrect tip labels with restored ages after rejection. To estimate ages, put the tree in `State` and choose compatible
operators. The original fixed-tree example deliberately does neither. Applying an age density to
a fixed tree only evaluates that density; it does not cause the age to move.

QuaSSE supplies the joint diversification/trait factor. Do not add an ordinary FBD tree prior on
top of it: that would count diversification twice. Parameter and age priors and additional data
likelihoods belong in the usual BEAST posterior. Age-prior factors specify an analysis and can
interact with the diversification factor; they do not redefine the local conditional likelihood.
