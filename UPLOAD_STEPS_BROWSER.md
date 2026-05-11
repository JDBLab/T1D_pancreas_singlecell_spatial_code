# Browser upload instructions for GitHub and GitHub Pages

These steps avoid command line tools.

## 1. Create the repository

1. Open GitHub in your browser.
2. Click `+` in the top-right corner.
3. Click `New repository`.
4. Suggested name: `t1d-pancreas-singlecell-spatial-code`.
5. Add a description: `Code and notebooks for single-cell and spatial transcriptomic analysis of human diabetic pancreas.`
6. Choose `Public` only after approval from your PI/collaborators.
7. Do not initialize with a README if you are uploading this package.
8. Click `Create repository`.

## 2. Upload files

If your browser supports folder drag-and-drop:

1. Open the unzipped repository folder on your computer.
2. Select all contents inside the folder, not the folder itself.
3. Drag them into GitHub using `Add file` > `Upload files`.
4. Commit with message: `Initial upload of manuscript analysis code`.

If GitHub does not let you upload folders:

1. Upload root files first: `README.md`, `_config.yml`, `.gitignore`.
2. Create folders manually using `Add file` > `Create new file`.
3. To create a folder, type a path such as `docs/index.md`, `notebooks/README.md`, or `R/README.md` in the file-name box.
4. Paste the content of the corresponding file and commit.
5. Then open each created folder and use `Add file` > `Upload files` for files in that folder.

## 3. Add your real notebooks

1. Open the `notebooks` folder in GitHub.
2. Click `Add file` > `Upload files`.
3. Upload your revised `.ipynb` notebooks.
4. Commit with message: `Add stLearn and SIRV scVelo PAGA notebooks`.

Recommended notebook names:

- `stLearn_spatial_trajectory_CCI.ipynb`
- `SIRV_scVelo_PAGA_spatial_velocity.ipynb`

## 4. Enable GitHub Pages

1. Go to repository `Settings`.
2. Click `Pages` in the left menu.
3. Under `Build and deployment`, select `Deploy from a branch`.
4. Select branch `main`.
5. Select folder `/docs`.
6. Click `Save`.

After a few minutes, GitHub will show a website link like:

`https://YOUR_USERNAME.github.io/t1d-pancreas-singlecell-spatial-code/`

## 5. Update the manuscript

Use the repository URL for code availability and the GitHub Pages URL as a clean landing page.
