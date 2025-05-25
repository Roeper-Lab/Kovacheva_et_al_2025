# Kovacheva_et_al_2025

This repository provides analysis scripts and example datasets accompanying the manuscript:

**"Recovery of the full in vivo firing range in post-lesion surviving DA SN neurons associated with Kv4.3-mediated pacemaker plasticity"**  
Kovacheva et al., 2025

The dataset includes MATLAB and Python tools for electrophysiological, behavioral, and immunohistochemical analysis of dopamine (DA) neurons in the substantia nigra (SN) in a mouse model of Parkinson's disease.

---

## 📁 Repository Structure

### 🧪 In Vitro Electrophysiology (MATLAB)

| Script | Description | Example Data (provided in this repository) |
|--------|-------------|--------------|
| `in_vitro_Kv4_activation.m` | Fits activation curves of Kv4.3 K+ currents | `in_vitro_Kv4_activation_example.mat` |
| `in_vitro_Kv4_inactivation.m` | Quantifies Kv4.3 current inactivation | `in_vitro_Kv4_inactivation_example.mat` |
| `in_vitro_HCN.m` | Analyzes HCN currents from hyperpolarizing steps | `in_vitro_HCN_example.mat` |
| `in_vitro_contCC.m` | Compares pacemaker variability pre/post whole-cell access | `in_vitro_contCC_example.mat` |
| `in_vitro_CChyp.m` | Evaluates firing behavior under Kv4.3-blockade (AmmTx3) | `in_vitro_CChyp_example.mat` |
| `in_vitro_onCell.m` | Analyzes spontaneous firing in on-cell configuration | `in_vito_onCell_example.mat` |

---

### 🧠 In Vivo Spike Analysis (MATLAB)

| Script | Description | Input |
|--------|-------------|-------|
| `in_vivo_TimeStamps_1-Import.m` | Imports spike timestamps and computes ISIs | `in_vivo_late_6OHDA_ISIs.txt` |
| `in_vivo_TimeStamps_2-Stat.m` | Calculates ISI distribution, CV, bursts, autocorrelation | Uses output from above |

---

### 🧍 Behavioral Analysis (MATLAB)

| Script | Description | Example Data (provided in this repository) |
|--------|-------------|--------------|
| `behavior_rotSequence.m` | Analyzes directional turning behavior | `behavior_rotSequence.mat` |
| `behavior_PercentRightRotations.m` | Computes % rightward rotations across time | `behavior_GroupsRotAnalysis.mat` |

---

### 🔬 Immunohistochemistry Analysis (Python)

| Script | Description | Example Data (provided in this repository) |
|--------|-------------|--------------|
| `Kv4_immuno.py` | Quantifies Kv4.3 intensity in SNc from confocal images | `Kv4_immuno_example_c1.jpeg`, `Kv4_immuno_example_c2.jpeg` |

---

## ▶️ Usage Instructions

All scripts are designed to run with the provided example data.  

### MATLAB

1. Open MATLAB and navigate to the script directory.
2. Open the desired `.m` file.
3. Ensure the matching `*_example.mat` or `.txt` file is in the same directory.
4. Run the script to generate plots, fits, and statistics.

### Python (for immuno)

```bash
python Kv4_immuno.py
```

Make sure the example `.jpeg` images are in the same folder.

---

## 🛠 Requirements

### MATLAB
- R2019a or newer
- Curve Fitting Toolbox
- Signal Processing Toolbox

### Python
- Python 3.8+
- numpy, scipy, pandas, matplotlib, scikit-image, opencv-python

Install via:

```bash
pip install numpy scipy pandas matplotlib scikit-image opencv-python
```

---

## 📜 Citation

If you use this repository, please cite:

> Kovacheva L., Shin J., Zaldivar-Diez J., et al. (2025). *Recovery of the full in vivo firing range in post-lesion surviving DA SN neurons associated with Kv4.3-mediated pacemaker plasticity*. In Review.

---

## 📝 License

This repository is distributed for academic use under the **MIT License**.  
See `LICENSE` file for details.

---

## ✉️ Contact

For questions or contributions, please contact:  
**Dr. Lora Kovacheva** – lorask@gmail.com  
**Prof. Dr. Jochen Roeper** – roeper@em.uni-frankfurt.de
