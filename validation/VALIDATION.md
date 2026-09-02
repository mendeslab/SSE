# SSE validation

This will guide you through the validation procedures we conducted (as well as generate all figures).
You are going to need the following R packages:

* ape
* diversitree
* HDInterval
* sjPlot
* gridExtra
* ggplot2

The calibrated-validation commands in section 1 are historical. Several input
templates, Python scripts, and generated data directories referenced there are
not present in this checkout, so that section is not currently reproducible
without obtaining those files separately.

## (1) Calibrated validation of BiSSE and ClaSSE    

We will simulate 2000 trees ("-n 2000"; some will not be considered for being too small or too large) with a simulation stop time of 50 ("-st 50") in step (1.1).
Simulation parameters will come from the specified priors and their respective parameters (e.g., "-pt 'exp,exp,exp,exp,exp,exp'" and "-pp '20;20;80;80;20;20'").
In step (1.2) we then parse the results.

Simulations require BiSSE and ClaSSE templates, which we provide and pass in with "-xt bisse_beast_template.xml".
We also specify a working directory with "-pd /path/to/working/directory/", an R script directory with "-rd /path/to/rscripts/directory", and an output directory with "-od /path/to/output/directory".

### (1.1) BiSSE (note that the value after "-n" will affect the results, even if setting the seed in R)

```python /path/to/SSE/validation/scripts/simulate_prep4beast.py -rd /path/to/SSE/validation/r_scripts/ -od /path/to/SSE/validation/calibrated/ -pt 'exp,exp,exp,exp,exp,exp' -pp '20;20;80;80;20;20' -pn l0,l1,m0,m1,q01,q10 -p bisse -xd /path/to/SSE/validation/calibrated/bisse_xmls/ -xt /path/to/SSE/validation/bisse_beast_template.xml -n 2000 -pd /path/to/SSE/validation/ -st 50 -b```

It should be necessary to ignore 10 too-large simulations, and have a n-tip median of 10; if you don't get this, something went wrong with the seeding and your beast_outputs won't match.

### (1.1.1) Add "_bisse" to all files inside "/path/to/SSE/validation/calibrated/" to diferentiate them from the ClaSSE files produced below.

### (1.2) ClaSSE

```
python /path/to/SSE/validation/scripts/generate_priors_table.py -n 3 -od /path/to/SSE/validation/calibrated/ -s SS,V,S -pt exp,exp,exp,exp,exp -pp '80;20;20;20;20' -p classe
python /path/to/SSE/validation/scripts/simulate_prep4beast.py -rd /path/to/SSE/validation/r_scripts/ -od /path/to/SSE/validation/calibrated/ -p classe -xd /path/to/SSE/validation/calibrated/classe_xmls/ -xt /path/to/SSE/validation/classe_beast_template.xml -n 2000 -pd /path/to/SSE/validation/ -ppf /path/to/SSE/validation/calibrated/classe_prior_params.csv -e2t /path/to/SSE/validation/calibrated/spec_event_to_triplet.csv -st 50
```

### (1.2.1) Add "_classe" to all files inside "/path/to/SSE/validation/calibrated/" to diferentiate them from the BiSSE files produced above.

### (1.3) After running all .xml files on the cluster, we need to parse the .log files

We do BiSSE first:

```
python /path/to/SSE/validation/scripts/parse_beast_logs.py -bd /path/to/SSE/validation/bisse_beast_outputs/ -rd /path/to/SSE/validation/r_scripts/ -cd /path/to/SSE/validation/calibrated/ -b 500000 -n 100 -n1 l0,l1,m0,m1,q01,q10 -n2 Lambda1,Lambda2,Mu1,Mu2,FlatQMatrix1,FlatQMatrix2
```

Then ClaSSE:

```
python /path/to/SSE/validation/scripts/parse_beast_logs.py -bd /path/to/SSE/validation/classe_beast_outputs/ -rd /path/to/SSE/validation/r_scripts/ -cd /path/to/SSE/validation/calibrated/ -b 500000 -n 100 -n1 l_111,l_313,l_312,m1,m2,m3,q01,q02,q10,q12,q20,q21 -n2 SympatricRate,SubsympatricRate,VicariantRate,Mu1,Mu2,Mu3,FlatQMatrix1,FlatQMatrix2,FlatQMatrix3,FlatQMatrix4,FlatQMatrix5,FlatQMatrix6
```

## (2) Plotting calibrated validation graphs

``Rscript /path/to/SSE/validation/r_scripts/calibrated_validation_plots.R /path/to/SSE/validation/ /path/to/SSE/validation/calibrated/``

## (3) Ancestral state reconstruction via stochastic character mapping    
In the SSE package, we implement the same stochastic character mapping approach proposed in Freyman and H&ouml;hna (2017).
We compare our implementation to diversitree's like they do (and produce the same figure), and also compare the two mapping methods (sampling at internal nodes only vs. sampling along the entire tree in fixed-size chunks).

### (3.1) Preparing input for stochastic character mapping (tree, tip data, internal node data, diversitree reconstructions) in R

The following command will simulate a tree and tip states given parameters
λ<sub>0</sub>=0.2, λ<sub>1</sub>=0.4, µ<sub>0</sub>=0.001,
µ<sub>1</sub>=0.1, q<sub>01</sub>=0.1, and q<sub>10</sub>=0.4. We use a fixed
seed to obtain the same graph in Freyman and H&ouml;hna (2017). The R script also
performs ancestral state reconstruction using diversitree.

Below, "/path/to/SSE/validation" should be replaced with your local path to the validation folder that comes with this git repository.

```
Rscript /path/to/SSE/validation/r_scripts/prepare_BiSSE_ASR_input.R /path/to/SSE/validation/ asr
```

The following files will be produced:

* asr/asr_piechart.pdf
* asr/simulated_tree.txt
* asr/simulated_tip_states.txt
* asr/true_ancestral_states.txt
* asr/diversitree_ancestral_states.txt

### (3.2) Running stochastic character mapping (both only on internal nodes and along branches as well) in Java

Compile the test and validation classes, then run the restored source directly.
The output argument is relative to the repository root in this example.

```
cd /path/to/SSE
ant compile-test
java -cp 'build:build-test:../beast2/build/dist/launcher.jar:'\
'../beast2/build/dist/BEAST.base.jar:'\
'../beast2/lib/junit/junit-platform-console-standalone-1.8.2.jar:lib/*' \
  validation.SDSEPValidationBiSSEASR validation/asr/
```

The following files will be produced:

* asr/SSE_asr_joint.csv
* asr/SSE_asr_stoc.csv

### (3.3) Plotting diversitree's ancestral state reconstruction vs. SSE packages's
Now we show that our ancestral state reconstructions matches diversitree's for BiSSE.

```
Rscript /path/to/SSE/validation/r_scripts/plot_diversitree_vs_SSEjoint.R /path/to/SSE/validation/asr/diversitree_ancestral_states.txt /path/to/SSE/validation/asr/SSE_asr_joint.csv /path/to/SSE/validation/plots/ bisse_joint
Rscript /path/to/SSE/validation/r_scripts/plot_diversitree_vs_SSEjoint.R /path/to/SSE/validation/asr/diversitree_ancestral_states.txt /path/to/SSE/validation/asr/SSE_asr_stoc.csv /path/to/SSE/validation/plots/ bisse_stoc
```

The R script above will produce two graph in the "plots/" folder:

* plots/bisse_joint_sse_vs_diversitree_asr.pdf
* plots/bisse_stoc_sse_vs_diversitree_asr.pdf
