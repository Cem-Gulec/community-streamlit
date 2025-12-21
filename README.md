# 🧬 Community Web

This is a Streamlit-based web interface for the [Community](https://github.com/SoloveyMaria/community) R package, which enables interactive exploration of cell-to-cell communication differences using single-cell RNA sequencing (scRNA-seq) data.

## 🧠 About

**Community** detects and compares ligand-receptor interactions across case and control samples. This web frontend enhances accessibility through visual dashboards, data upload interfaces, dynamic tables, and visualizations including:

- Volcano plots
- Heatmaps
- Network plots
- Forest plots

All analysis runs locally, integrating Python and R.

## 🚀 Features

- 📂 Upload or use demo datasets for analysis.
- 🔍 Filter and visualize interactions using thresholds and statistical tests.
- 🖼️ Generate dynamic visualizations (plots, heatmaps, networks).
- 📊 Explore differential interactions and component-wise effects.
- 📁 Customize ligand-receptor database interactively.
- 🎨 Modify visualization parameters and styles.

## 🛠️ Tech Stack

- **Frontend**: Streamlit
- **Backend**: Python & R (via `subprocess`)
- **Visualization**: Plotly, Streamlit-AgGrid
- **Analysis Engine**: R package [`community`](https://github.com/SoloveyMaria/community)

## 📦 Requirements

Install dependencies from `requirements.txt`:

```bash
pip install -r requirements.txt
````

You also need R installed with the following R packages:

* `community`
* `data.table`
* `tidyverse`
* `gridExtra`

## ▶️ Run the App

```bash
streamlit run app.py
```

Ensure R is properly installed and available in your system path, as it is called from Python using `subprocess`.

## 📁 Directory Structure

```
├── app.py                  # Main Streamlit app
├── backend.R               # Main R backend logic
├── heatmap.R               # Heatmap generation script
├── visualization.R         # Network/volcano plots
├── visualizations/         # Data and helper scripts for visualizations
├── plots/                  # Generated plot images
├── each_component_values/  # Interaction-level CSV data
├── input_data/             # Toy demo datasets (loaded dynamically)
└── requirements.txt
```

## 📄 License

This project is intended for research purposes. Refer to the source R package [`community`](https://github.com/SoloveyMaria/community) for licensing.

## 🙌 Acknowledgements

* Developed on top of the [`community`](https://github.com/SoloveyMaria/community) R package.
* Web interface built by [Cem Güleç](https://github.com/Cem-Gulec) with the guidance of Muhammet A. Celik and Dr. Maria Solovey.
