# Final GitHub and Zenodo Checklist

## Before Making the Repository Public

1. Confirm that the working tree is clean:

   ```bash
   git status
   ```

2. Confirm that no private data are tracked:

   ```bash
   git ls-files
   ```

3. Confirm that no large RDS/raw-data files are staged:

   ```bash
   git status --short
   ```

4. Confirm the final license choice.

5. Confirm whether analysis-ready Uruguay inputs can be redistributed. If not,
   keep `data/analysis_ready/` empty except for its README and document the
   required inputs in `data/metadata/`.

6. Replace placeholders in `CITATION.cff`:

   - affiliation;
   - GitHub URL;
   - Zenodo DOI after release.

## GitHub Release

1. Create a public GitHub repository.

2. Push the final branch:

   ```bash
   git remote add origin https://github.com/OWNER/REPOSITORY.git
   git push -u origin main
   ```

3. Create a release:

   - tag: `v1.0-submission`
   - title: `Replication repository for Biostatistics submission`
   - description: include the manuscript title and state that this is the
     submission replication archive.

## Zenodo

1. Log in to Zenodo.

2. Connect GitHub in Zenodo.

3. Enable archiving for the GitHub repository.

4. Create or re-create the GitHub release `v1.0-submission`.

5. Wait for Zenodo to archive the release and assign a DOI.

6. Update:

   - `CITATION.cff`;
   - manuscript reproducibility statement;
   - supplementary material;
   - cover letter.

7. Commit the DOI update and optionally create a follow-up release:

   ```bash
   git add CITATION.cff README.md
   git commit -m "Add Zenodo DOI"
   git push
   ```
