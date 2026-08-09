# 🦠 Discerning Ecological Behaviour and Interactions of Minimal Reactomes

![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-blue.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

> **Project 2026/WSAI/018** | Wadhwani School of Data Science and AI (WSAI), IIT Madras

This repository contains the computational pipeline, MATLAB scripts, and empirical datasets for the systematic evaluation of metabolic redundancy and ecosystem dynamics in **52 human gut bacterial species**. 

## 🔬 Project Overview
The human gut microbiome relies on complex metabolic interactions to maintain homeostasis. For organisms to survive sudden nutritional perturbations, they utilize redundant metabolic pathways (backup safety nets). 

This project computationally isolates the bare essential pathways (the **Minimal Reactome**) for 52 gut microbes using Mixed-Integer Linear Programming (MILP) and parsimonious Flux Balance Analysis (pFBA). By simulating 1,326 pairwise microbial communities, the pipeline maps how ecological behaviors shift when microbes are stripped of these redundant pathways. The findings mathematically demonstrate that metabolic reduction systematically forces cooperative microbial guilds into competitive states, highlighting the critical role of redundancy in sustaining mutualism.

---

## 🧪 Methodology & Simulation Parameters

### 1. Dietary Environments
The models were evaluated across six distinct dietary profiles to test adaptation to nutrient availability:
* **No Diet:** Unrestricted nutrient availability (Baseline)
* **High Fiber:** Rich in complex carbohydrates
* **Mediterranean:** Balanced, plant-rich, and healthy fats
* **Vegetarian:** Strict plant-based nutrient bounds
* **Western:** High fat and high sugar composition
* **Unhealthy:** Highly processed, nutrient-poor bounds

### 2. Interaction Classification (α Thresholds)
Community interaction types were determined by calculating the relative growth change (α) for each organism against its monoculture baseline:

α_i = (X_i(community) - X_i(monoculture)) / X_i(monoculture)

Using a tolerance band of ± 0.1 for an "unaffected" state, the 1,326 pairs were classified into six ecological categories:

| Condition on α_1 | Condition on α_2 | Interaction Type |
| :---: | :---: | :---: |
| -0.1 ≤ α_1 ≤ 0.1 <br> α_1 ≤ -0.1 | α_2 ≤ -0.1 <br> -0.1 ≤ α_2 ≤ 0.1 | **Amensalism** |
| -0.1 ≤ α_1 ≤ 0.1 <br> α_1 ≥ 0.1 | α_2 ≥ 0.1 <br> -0.1 ≤ α_2 ≤ 0.1 | **Commensalism** |
| α_1 ≤ -0.1 <br> α_1 ≥ 0.1 | α_2 ≥ 0.1 <br> α_2 ≤ -0.1 | **Exploitation** |
| α_1 ≤ -0.1 | α_2 ≤ -0.1 | **Competition** |
| α_1 ≥ 0.1 | α_2 ≥ 0.1 | **Mutualism** |
| -0.1 ≤ α_1 ≤ 0.1 | -0.1 ≤ α_2 ≤ 0.1 | **Neutralism** |

---

## 📂 Repository Structure

    Gut-Microbiome-Minimal-Reactomes/
    ├── data/
    │   ├── agora_2_models/                  # 52 Wild-Type .mat models (NCBI mapped)
    │   ├── diets/                           # 6 dietary constraint profiles (.txt)
    │   ├── minReactModels_WideSummary.csv   # Monoculture baseline datasets
    │   └── Community_Results_Wide_52Models.csv # Dual-objective community results
    ├── scripts/
    │   ├── Minimal_Reactomes.m              # Core MILP minimization engine
    │   ├── Monoculture52.m                  # FBA baseline generation
    │   ├── Pairwise_Community.m             # Combinatorial ecosystem simulator
    │   ├── Plot1_InteractionCounts.m        # Visualization suite (Plots 1-7)
    │   ├── Plot2_GrowthDynamics.m
    │   ├── Plot3_MacroComposition.m
    │   ├── Plot4_Heatmaps.m
    │   ├── Plot5_SankeyTransitions.m
    │   ├── Plot6_PhylumGrids.m
    │   └── Plot7_JaccardOverlap.m
    ├── figures/                             # 32 High-resolution .png ecological plots
    └── docs/                                
        ├── Final_Report.pdf                 # Compiled project report
        └── Final_Report.tex                 # LaTeX source code

---

## 💾 Data Availability
Due to GitHub's repository file size limitations, the complete database of algorithmically generated minimal reactome models (>1 GB of .mat files across 6 diets) is hosted externally on Google Drive. 

* **Minimal Reactome Models:** [Download Complete .zip Dataset Here](https://drive.google.com/file/d/1Pxx1Lv3hPc6DUgRE9lAYgn_EvOzOl9Zt/view?usp=sharing)
* The foundational Wild-Type models were sourced directly from the [AGORA2 database](https://www.vmh.life/).

---

## ⚙️ Computational Requirements
* **MATLAB** (R2019b or newer recommended)
* **The COBRA Toolbox v2.0**
* **Gurobi Optimizer** (Solver configured with 10^-6 feasibility tolerance margin)

---

## 🚀 Execution Pipeline
To reproduce the empirical findings, execute the scripts sequentially:

1. **Minimal_Reactomes.m**: Generates the minimal reactomes for all 52 models across the 6 diets using MILP/pFBA. *(Note: Requires parallel pool `parpool` initialization; runtime varies heavily by CPU architecture).*
2. **Monoculture52.m**: Applies dietary constraints and computes baseline monoculture growth rates.
3. **Pairwise_Community.m**: Assembles the 1,326 combinatorial pairs in a shared microenvironment and executes dual-objective maximization.
4. **Visualization Suite**: Execute `Plot1` through `Plot7` sequentially to classify the α interaction vectors and render the 32 ecological visualizations.

---



## 📄 License
This project is licensed under the MIT License - see the `LICENSE` file for details.
