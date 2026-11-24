# Central Dogma Interaction Equation: implemented equations

These equations summarize the dimensionless illustrative models in [compartments.py](compartments.py) and [model.py](model.py). See [Methods](METHODS.md) for assumptions and biological references, [Data sources](DATA_SOURCES.md) for dataset provenance, and [Reproducibility](REPRODUCIBILITY.md) for the validation protocol.

## Replication and inheritance

- Nuclear fork arrival: $t(x)=\min_i[t_i+|x-x_i|/v]$, using scheduled origin positions $x_i$ and initiation times $t_i$.
- Mitochondrial replication: $t_L=d_{\rm circ}(O_H,O_L)/v_H+\tau$ and $t_{\rm end}=\max(L/v_H, t_L+L/v_L)$. The compartment example scales fork speed as $v_{\max}f(\mathrm{POLG})f(\mathrm{TWINKLE})f(\mathrm{SSBP1})$ and uses $\tau=1/f(\mathrm{POLRMT}_{\rm prime})$, with $f(x)=x/(K+x)$ and default $K=1$.
- Offspring mtDNA fraction: $K_v\sim\mathrm{Binomial}(n,h_m)$ and $h_{\rm child}=K_v/n$. This requires a maternal germline fraction $h_m$ and an effective bottleneck size $n$.
- At an XX locus the maternal and paternal activity weights are $(p+e(1-p),1-p+ep)$, with maternal X activation fraction $p$ and escape fraction $e$. The code enumerates phased offspring genotypes for a single autosomal or non-pseudoautosomal locus.
- Tissue or mosaic mixture: $Y_{\rm tissue}=\sum_c w_c Y_c$ with $\sum_c w_c=1$. Each group can have its own intracellular mtDNA fraction $h_c$.

## Nuclear transcription, RNA processing and protein

- Precursor flux: $J=d A f(\mathrm{PolII})\exp[b_{\rm CCAAT}f(\mathrm{NFY})+b_{\rm TATA}f(\mathrm{TBP})+E C-R]$. Here $d$ is dosage, $A$ accessibility, $b$ indicates an available promoter motif, $E$ enhancer activity, $C$ promoter contact and $R$ repressor drive.
- Exportable flux: $J_{\rm exp}=J f_{\rm cap}f_{\rm poly(A)}$. The two factors describe effective 5′ capping and 3′ cleavage/polyadenylation.
- For splice isoform $i$, $M_i=J_{\rm exp}p_i/(\delta_i+I_i\,\delta_{N,i})$ and $\sum_i p_i=1$. Frame changes follow the coding-length difference modulo three. Here the indicator $I_i$ denotes a premature termination codon and $\delta_{N,i}$ is the user-supplied extra decay rate. A separate heuristic checks premature-stop distance from the last junction.
- Translation flux is $J_i=M_i\alpha_i/(1+\alpha_i L_i/v_i)$, where effective initiation $\alpha_i$ includes a Kozak factor or a ribosome availability factor depending on the entry point. Delivered folded protein is $P_i=J_i f_{\rm fold,i}f_{\rm delivery,i}/\delta_{\rm protein,i}$, with functional output $P_i a_i$.

## Mitochondrial transcription, translation and assembly

- Primary transcript fluxes are $J_s=N f(\mathrm{POLRMT})f(\mathrm{TFAM})f(\mathrm{TFB2M})b_s$, with mtDNA copies $N$ and distinct LSP/HSP strengths $b_s$.
- For the illustrative mature RNA pool, $M=(J_L+J_H)f(\mathrm{TEFM})[w f(\mathrm{RNaseP})f(\mathrm{ELAC2})+(1-w)e_{\rm other}]$. The stop-completion term is $q=1-u+u f(\mathrm{MTPAP})$, with incomplete-stop fraction $u$.
- Usable tRNA per allele is abundance times 3′-CCA completion, conditional charging, and average decoding efficiency across modification states. The MT-TL1 example mixes those efficiencies as $(1-h)t_{\rm WT}+h t_{\rm var}$, then scales by heavy-strand output, elongation and cleavage to obtain $T$.
- Mitochondrial subunit flux is $M q r_{28}r_{39}T/(K_T+T)$. The nuclear contribution multiplies precursor supply by TOM, TIM23 and HSPA9 efficiencies, membrane potential and HSP60/HSP10 folding factors. Assembly capacity is the smaller of that contribution and the insertion-adjusted mitochondrial subunit pool.
- The interaction contrast is $\Gamma_Y=Y_{AB}-Y_A-Y_B+Y_0$, applied to the same endpoint across all four arms.

The checks and example outputs cover synthetic scenarios. These formulas have no patient-level calibration, and the complex-assembly proxy does not measure respiration or ATP.
