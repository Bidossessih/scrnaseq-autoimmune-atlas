# Single-cell autoimmune atlas

Analysis code for building and interpreting a single-cell RNA-seq atlas of autoimmune disease samples. The repository combines Python/Jupyter workflows for preprocessing, integration, label transfer, and Milo analysis with R and Quarto workflows for statistical analysis, visualization, trajectory analysis, and gene-target interpretation.

## Workflow

1. **Preprocessing** (`notebooks/1-preprocessing/`)
   - Convert and standardize study-specific inputs.
   - Perform quality control, filtering, normalization, concatenation, and doublet or ambient-RNA correction where required.
   - The `remapped/` and `ummapped/` folders contain separate workflows for remapped and unmapped studies.
2. **Integration** (`notebooks/2-Integration/`)
   - Integrate PBMC and tissue datasets.
   - Apply CellTypist and scANVI-based annotation or label transfer.
   - Reintegrate B-cell, T-cell, and myeloid subsets for downstream analyses.
3. **Milo** (`notebooks/3-Milo/`)
   - Identify differential abundance across conditions and cell populations.
   - Use the exported cell-level tables in `notebooks/data/` for downstream summaries and plots.
4. **R and Quarto analyses** (`notebooks/4-R/`)
   - Produce figures, differential-expression summaries, trajectory analyses, supplementary tables, and drug or gene-target analyses.

## Repository layout

```text
notebooks/
├── 1-preprocessing/   Study preprocessing and QC notebooks
├── 2-Integration/     Integration, annotation, and label transfer
├── 3-Milo/            Differential-abundance analysis
├── 4-R/               R and Quarto reporting workflows
└── data/              Exported tables used by downstream analyses
```

The notebooks are organized by analysis stage rather than as a single executable pipeline. Some notebooks are study-specific, and filenames containing `Copy`, numbered variants, or dates represent alternate analysis iterations.

## Python environment

The main Python environment is defined in [`notebooks/2-Integration/environment.yml`](notebooks/2-Integration/environment.yml). It includes the Scanpy, AnnData, scVI-tools, CellTypist, GPU, and plotting dependencies used by the notebooks.

Create the environment with Conda or Mamba:

```bash
conda env create -f notebooks/2-Integration/environment.yml
conda activate autoimmune
```

The integration environment includes CUDA packages and is intended for a Linux system with a compatible NVIDIA driver when GPU notebooks are used. CPU execution may require adapting the relevant notebook settings.

## Data and reproducibility

The input data and exported analysis tables are available from [Zenodo (10.5281/zenodo.21514797)](https://doi.org/10.5281/zenodo.21514797). Download the required files from Zenodo and place them where the relevant notebook expects them. The CSV files under `notebooks/data/` and the archive `./notebooks/script_autoimmune.zip` are ignored by Git because they are generated or distributed data files.

Because this is a notebook-based research workflow, results depend on the selected study inputs, package versions, random seeds, GPU availability, and the order in which notebooks are run. Record any changed paths, parameters, environment details, and notebook outputs when creating a reproducible analysis run.

## Citation and contact

Add the associated manuscript, dataset accession numbers, and maintainer contact here when they are finalized.
