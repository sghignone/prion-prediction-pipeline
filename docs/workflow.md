## Step 1: Data Collection
* Download fungal proteomes from UniProt Trembl fungi dataset (uniprot_trembl_fungi.dat.gz) to have a broad and high-quality starting point.
* (optional) Extract proteomes of target species from metadata:
    AMF species (Rhizophagus, Gigaspora, Diversispora, etc.)
    Other guild representatives (ECM fungi, pathogens, saprotrophs) using taxonomic and ecological metadata.

## Step 2: Identification of Prion-Like Domains (Q/N-rich)
* Use **PrionScan**, an online database and accompanying software based on Espinosa Angarica's work, which applies probabilistic models (HMMs) trained on experimentally validated prion domains to scan full proteomes for candidate prion domains with Q/N enrichment and disorder characteristics (score threshold recommended by original papers ~0.2 for candidate confidence). 
[in short: Use PrionScan, which applies probabilistic models of Q/N-rich prion domains for high-confidence candidates.]  

    Reference papers:  
        - "PrionScan: an online database of predicted prion domains in complete proteomes"  
        - "Discovering putative prion sequences in complete proteomes using probabilistic representations of Q/N-rich domains"

* Include **PAPA** (Prion Aggregation Prediction Algorithm) to predict prion-like domains based on amino acid composition and aggregation propensity from Toombs et al. 2012. PAPA outperforms traditional predictors in some cases by capturing sequence features related to synthetic and natural prions. Default thresholds from the original work (aggregation scores above 0.05) guide candidate selection.
* Alternatively or complementarily, run **PLAAC** (Prion-Like Amino Acid Composition) to identify prion-like domains using a log-likelihood scoring function for Q/N-rich sequences relative to a background proteome. PLAAC offers a score cutoff of 20 for high-confidence hits, but tuning may be needed. [in short: PLAAC complementary scans for prion-like domains based on log-likelihood of Q/N composition.]
* Supplement with **fLPS** to identify low-complexity Q/N enriched regions in proteins, using default thresholds but focusing on those overlapping predicted disordered regions. [in short: Use fLPS to annotate low-complexity Q/N-rich regions overlapping with predicted disordered sections.]

## Step 3: Identification of Non-Q/N Prion-Like or Amyloid-Forming Domains
* Run **WALTZ** (a sequence-based amyloid motif predictor) on the proteomes or selected prion-like candidates to capture short amyloidogenic peptides that may not be Q/N-rich but still form amyloids/prions. Thresholds around 0.5-0.6 stringency usually balance sensitivity and specificity.
* Use **ArchCandy** or **AmyloGram** to further confirm amyloid propensity of predicted domains.
* Analyze predicted intrinsic disorder using **IUPred2A** to correlate prion predictions with disorder, a strong hallmark of prion domains. Prion domains typically occur in disordered segments. [in short: Use IUPred2A for disorder prediction, key for validating prion domain candidates from both Q/N-rich and amyloid predictors.]

## Additional Integration
* Combine predictions from PrionScan, PAPA, and PLAAC into a union set of prion-like proteins to maximize sensitivity while providing confidence levels.
* Score overlap and consensus among these tools to prioritize the most likely prion-like candidates for further analysis.

## Step 4: Domain Annotation and Functional Context
* Integrate **InterProScan** results to annotate predicted prion-like proteins with domain functions, signal peptides, transmembrane regions, and GO terms.
* Look for fusion of prion-like domains with NLR signaling domains or secretion signals as these may indicate immune or ecological functions.

## Step 5: Quantification and Comparative Analysis
* Calculate per-proteome metrics:
    - Number and fraction (%) of proteins with prion-like domains (Q/N-rich and amyloid).
    - Average prion domain length and log-likelihood scores.
    - Distribution of prion-like proteins by functional category (secreted, cell-wall, signaling).
* Normalize by proteome size for comparison across species.
* Compare AMF to other guilds using statistical tests:
    - Non-parametric tests (Mann-Whitney U) for distribution differences in prion protein counts.
    - Chi-square or Fisher’s exact tests for enrichment of prion-like domains in functional groups or guilds.

## Step 6: Ecological Guild Mapping and Interpretation
* Get or curate fungal ecological guild annotations based on FUNGuild or literature.
* Visualize prion domain expansions/clusters in AMF relative to saprotrophs, pathogens, and ECM fungi.
* Highlight prion domain architectures uniquely expanded in AMF or symbiotic fungi.

## Updated Quantification and Comparative Analysis
* Track which proteins are predicted by which tools (e.g., Venn distribution of PrionScan, PAPA, PLAAC hits).
* Quantify prion-like protein counts, lengths, and functional categories per tool and combined.
* Comparative statistical testing as before, with detailed breakdown by prediction method, to robustly demonstrate expansion or enrichment patterns in AMF.

## Summary Table of Prediction Tools
| Tool      | Purpose                       | Key Features                                                         | Threshold Example                 |
| --------- | ----------------------------- | -------------------------------------------------------------------- | --------------------------------- |
| PrionScan | Probabilistic Q/N-rich PrLDs  | HMM-based on experimentally validated prions                         | Score > 0.2 (default)             |
| PAPA      | Aggregation-prone PrLDs       | Composition and aggregation propensity-based prediction              | Aggregation score > 0.05          |
| PLAAC     | Log-likelihood Q/N PrLDs      | Statistical Q/N enrichment relative to background proteome           | Score > 20 (high-confidence hits) |
| fLPS      | Low-complexity Q/N regions    | Identifies simple sequence repeats and Q/N-rich low complexity       | Default settings                  |
| WALTZ     | Amyloid motif predictor       | Detects short amyloidogenic peptides independently of Q/N enrichment | Score ~0.5-0.6                    |
| ArchCandy | Amyloid propensity predictor  | Complementary sequence-based amyloid prediction                      | Default parameters                |
| IUPred2A  | Intrinsic disorder prediction | Validates prion candidates by disorder context                       | Default threshold (IUPred > 0.5)  |