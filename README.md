

R scripts that accompanies the following paper: 

## Challenging Additivity: Comparing Predicted and Observed AhR Activity of Polycyclic Aromatic Compound (PAC) Mixtures Containing Active and Inactive Constituents

Class-based cumulative risk assessment approaches have been applied to high-priority environmental contaminants such as polycyclic aromatic compounds (PACs), yet uncertainties remain in their application. In this study, we evaluated the influence of inactive chemicals on mixture modeling outcomes and explored strategies for predicting the aryl hydrocarbon receptor (AhR)-mediated toxicity of PAC mixtures. Using an in vitro AhR reporter gene assay, we tested seven defined mixtures composed of six active and seven inactive PACs. Observed concentration-response curves were compared to predictions from three established mixture models: concentration addition (CA), independent action (IA), and generalized concentration addition (GCA), using both effective concentration eliciting 10% response (EC10) and benchmark concentration (BMC10) approaches. Including inactive chemicals without scaling led to consistent overestimation of potency, especially in models assuming equal efficacy. Predictive accuracy improved across all models when mixtures were limited to active chemicals and contributions were scaled to 100%, excluding inactives. Among approaches, GCA consistently produced the best agreement with measured responses, particularly when paired with BMC modeling. BMC10 values better accommodate partial agonists. Our findings support a pragmatic, mechanism-based framework for modeling environmental mixtures, one that prioritizes active components, scales their contributions, and adopts BMC-based methods to estimate potency. 

-----------------------------------------------
This repository provides code to:
- fit concentration–response models for individual PACs
- analyze in vitro responses of defined PAC mixtures
- generate mixture model predictions
- compare predicted and observed concentration–response relationships
- summarize potency metrics such as EC10 and benchmark concentrations

## Repository contents
These scripts must be run in sequential order

| File | Description |
|-----|------------|
| `01-individual_chemical_dr.R` | Fit concentration–response models for individual PACs |
| `02-mixtures_modeling.R` | Implement mixture models (e.g., CA, IA) |
| `03-mixtures_dr.R` | concentration–response analysis for PAC mixtures |
| `04-invitro_EC10_mixture.R` | Estimate EC10 values for mixtures |
| `05-compare_dr.R` | Compare observed vs. predicted responses |
| `06-dr_actives_fig.R` | Plot concentration–response curves for active chemicals |
| `07-pie_plot_fig.R` | Visualize mixture composition |
| `08-plot_model_params.R` | Plot mixture model parameters |
| `09-summarize_bmc.R` | Summarize benchmark concentration metrics |
| `LICENSE` | MIT license |

